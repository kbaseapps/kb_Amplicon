
import errno
import json
import logging
import os
import shutil
import uuid
import zipfile
import re
import subprocess

import pandas as pd
import plotly.express as px
from plotly.offline import plot
import plotly.graph_objs as go

from installed_clients.DataFileUtilClient import DataFileUtil
from installed_clients.KBaseReportClient import KBaseReport


class MDSUtils:

    R_BIN = '/kb/deployment/bin'
    MDS_OUT_DIR = 'mds_output'
    PARAM_IN_WS = 'workspace_name'
    PARAM_IN_MATRIX = 'input_obj_ref'
    PARAM_OUT_MATRIX = 'mds_matrix_name'

    def _mkdir_p(self, path):
        """
        _mkdir_p: make directory for given path
        """
        if not path:
            return
        try:
            os.makedirs(path)
        except OSError as exc:
            if exc.errno == errno.EEXIST and os.path.isdir(path):
                pass
            else:
                raise

    def _validate_run_mds_params(self, params):
        """
        _validate_run_mds_params:
            validates params passed to run_mds method
        """

        logging.info('start validating run_mds params')

        # check for required parameters
        for p in [self.PARAM_IN_MATRIX, self.PARAM_IN_WS, self.PARAM_OUT_MATRIX]:
            if p not in params:
                raise ValueError('"{}" parameter is required, but missing'.format(p))

    def _build_rMDS_script(self, params):
        """
        _build_rMDS_script: build a sequence of R command calls according to params
        Note: To run the NMDS, we will use the function metaMDS from the vegan package.
        # The metaMDS function requires only a community-by-species matrix.
        """
        data_file_path = params.get('datafile')
        if not data_file_path:
            return ''

        exists = os.path.isfile(os.path.join(self.output_dir, os.path.basename(data_file_path)))
        if not exists:
            shutil.copyfile(data_file_path,
                            os.path.join(self.output_dir, os.path.basename(data_file_path)))

        associated_matrix_file = params.get('associated_matrix_file')

        n_components = params.get('n_components', 2)
        max_iter = params.get('max_iter', 300)
        run_metric = True if params.get('metric', 0) else False
        dist_metric = params.get('distance_metric', 'bray')

        mds_cfg = 'distance="' + dist_metric + '",try=20,trymax=' + str(max_iter) + \
                  ',autotransform=TRUE,noshare=0.1,expand=TRUE,trace=1,' + \
                  'plot=FALSE,engine=c("monoMDS","isoMDS"),k=' + str(n_components)
        if run_metric:
            mds_cfg += 'metric=True'

        mds_scrpt = 'library(vegan)\n'
        mds_scrpt += 'library(jsonlite)\n'
        mds_scrpt += 'vg_data <- read.table("' + data_file_path + \
                     '",header=TRUE,row.names=1,sep="")\n'
        # remove the last (taxonomy) column
        # mds_scrpt += 'vg_data<-vg_data[,1:dim(vg_data)[2]-1]\n'
        # Function metaMDS returns an object of class metaMDS.
        mds_scrpt += 'vg_data.mds <- metaMDS(vg_data,' + mds_cfg + ')\n'
        mds_scrpt += 'vg_data.mds\n'

        # save the results in the memory
        # 1) store species ordination
        mds_scrpt += 'variableScores <- vg_data.mds$species\n'
        # 2) store site ordination
        mds_scrpt += 'sampleScores <- vg_data.mds$points\n'
        # 3) store other ordination results
        mds_scrpt += 'stress <- vg_data.mds$stress\n'
        mds_scrpt += 'dist_metric <- vg_data.mds$distance\n'
        mds_scrpt += 'dist_matrix <- vg_data.mds$diss\n'
        mds_scrpt += 'dist_call <- vg_data.mds$distcall\n'
        mds_scrpt += 'converged <- vg_data.mds$converged\n'
        mds_scrpt += 'dims <- vg_data.mds$ndim\n'
        mds_scrpt += 'tries <- vg_data.mds$tries\n'
        mds_scrpt += 'maxits <- vg_data.mds$maxits\n'

        # save the results to the current dir
        # Write CSV in R
        mds_scrpt += 'write.csv(dist_matrix,file="dist_matrix.csv",row.names=TRUE,na="")\n'
        mds_scrpt += 'write.csv(variableScores,file="species_ordination.csv",' + \
                     'row.names=TRUE,na="")\n'
        mds_scrpt += 'write.csv(sampleScores,file="site_ordination.csv",row.names=TRUE,na="")\n'

        if associated_matrix_file:
            mds_scrpt += 'chem_data <- read.table("' + associated_matrix_file + \
                 '",header=TRUE,row.names=1,sep="")\n'
            mds_scrpt += '(fit <- envfit(vg_data.mds,chem_data,perm=999))\n'
            mds_scrpt += 'vectors <- scores(fit, "vectors")\n'
            mds_scrpt += 'write.csv(vectors,file="vectors.csv",row.names=TRUE,na="")\n'

        # Write JSON in R
        mds_scrpt += 'item_name=c("stress","distance_metric","dist_call","converged",' + \
                     '"dimesions","trials","maxits")\n'
        mds_scrpt += 'item_value=c(stress,dist_metric,dist_call,converged,dims,tries,maxits)\n'
        mds_scrpt += 'df <- data.frame(item_name,item_value,stringsAsFactors=FALSE)\n'
        mds_scrpt += 'write_json(toJSON(df),path="others.json",pretty=TRUE,auto_unbox=FALSE)\n'

        # If there is user input plotting script:
        plt_scrpt = params.get('plot_script', '').lower()
        if plt_scrpt and re.match("^plot\(\s*[a-zA-Z]+.*\)$", plt_scrpt):
            arr_plt = plt_scrpt.split(',')
            arr_plt[0] = 'plot(vg_data.mds'  # make sure to pass the correct data
            plt_scrpt = (',').join(arr_plt)
            if len(arr_plt) == 1:
                plt_scrpt += ')'
            plt_type = params.get('plot_type', 'pdf').lower()
            if not plt_type:
                plt_type = 'pdf'

            plt_name = params.get('plot_name', 'usr_plt_name').lower()
            if not plt_name:
                plt_name = 'usr_plt_name'
            plt_name += '.' + plt_type

            if plt_type == 'jpg':
                plt_type = 'jpeg'
            if plt_type == 'ps':
                plt_type = 'postscript'
                mds_scrpt += plt_type
                mds_scrpt += '(file="' + plt_name + '")\n'
            if plt_type == 'tiff':
                mds_scrpt += plt_type
                mds_scrpt += '(file="' + plt_name + '",width=4,height=4,units="in",' + \
                             'compression="lzw",res=300)\n'
            if plt_type in ['jpg', 'jpeg', 'bmp', 'png']:
                mds_scrpt += plt_type
                mds_scrpt += '(file="' + plt_name + '",width=580,height=580,units="px",' + \
                             'res=100, pointsize=12)\n'

            mds_scrpt += plt_scrpt + '\n'

            if associated_matrix_file:
                mds_scrpt += 'plot(fit)\n'
            mds_scrpt += 'dev.off()\n'

        logging.info('R script: {}'.format(mds_scrpt))

        mds_rscript = 'mds_script.R'
        rscrpt_file_path = os.path.join(self.output_dir, mds_rscript)

        with open(rscrpt_file_path, 'w') as r_file:
            r_file.write(mds_scrpt)
        return rscrpt_file_path

    def _execute_r_script(self, rfile_name):
        """
        _execute_r_script: Calling the Rscript executable to run the R script in rfile_name
        """
        logging.info('Calling R......')

        result_dir = os.path.dirname(rfile_name)
        if not result_dir:
            result_dir = self.working_dir

        rcmd = [os.path.join(self.R_BIN, 'Rscript')]
        rcmd.append(rfile_name)

        logging.info('Running metaMDS script in current working directory: {}'.format(result_dir))
        exitCode = 0
        try:
            complete_proc = subprocess.run(rcmd, cwd=result_dir, stdin=subprocess.PIPE,
                                           stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                                           close_fds=True)
            exitCode = complete_proc.returncode
            if (exitCode == 0):
                logging.info('\n{}'.format(complete_proc.stdout))
                logging.info('\n{} was executed successfully, exit code was: {}'.format(
                    ' '.join(rcmd), str(exitCode)))
                logging.info("Finished calling R.")
            else:
                logging.info('Error running command: {} Exit Code: {}'.format(' '.join(rcmd),
                                                                              str(exitCode)))
                logging.info('\n{}'.format(complete_proc.stderr))
                logging.info('\n{}'.format(complete_proc.stdout))
        except subprocess.CalledProcessError as sub_e:
            exitCode = -99
            logging.info('Caught subprocess.CalledProcessError {}'.format(sub_e))

        logging.info('created files in {}:\n{}'.format(result_dir, os.listdir(result_dir)))

        return exitCode

    def _df_to_list(self, df):
        """
        _df_to_list: convert Dataframe to FloatMatrix2D matrix data
        """

        df.index = df.index.astype('str')
        df.columns = df.columns.astype('str')
        df.fillna(0, inplace=True)
        matrix_data = {'row_ids': df.index.tolist(),
                       'col_ids': df.columns.tolist(),
                       'values': df.values.tolist()}

        return matrix_data

    def _save_mds_matrix(self, workspace_name, input_obj_ref, mds_matrix_name,
                         distance_df, mds_params_df, site_ordin_df, species_ordin_df):

        logging.info('Saving MDSMatrix...')

        if not isinstance(workspace_name, int):
            ws_name_id = self.dfu.ws_name_to_id(workspace_name)
        else:
            ws_name_id = workspace_name

        mds_data = {}

        mds_data.update({'distance_matrix': self._df_to_list(distance_df)})
        mds_data.update({'site_ordination': self._df_to_list(site_ordin_df)})
        mds_data.update({'species_ordination': self._df_to_list(species_ordin_df)})
        mds_data.update({'mds_parameters': self._df_to_list(mds_params_df)})
        mds_data.update({'original_matrix_ref': input_obj_ref})
        mds_data.update({'rotation_matrix': self._df_to_list(distance_df)})

        obj_type = 'KBaseExperiments.PCAMatrix'
        info = self.dfu.save_objects({
            "id": ws_name_id,
            "objects": [{
                "type": obj_type,
                "data": mds_data,
                "name": mds_matrix_name
            }]
        })[0]

        return "%s/%s/%s" % (info[6], info[0], info[4])

    def _zip_folder(self, folder_path, output_path):
        """
        _zip_folder: Zip the contents of an entire folder (with that folder included in the
         archive). Empty subfolders could be included in the archive as well if the 'Included
         all subfolders, including empty ones' portion.
         portion is used.
        """
        with zipfile.ZipFile(output_path, 'w',
                             zipfile.ZIP_DEFLATED,
                             allowZip64=True) as ziph:
            for root, folders, files in os.walk(folder_path):
                # Include all subfolders, including empty ones.
                for folder_name in folders:
                    absolute_fpath = os.path.join(root, folder_name)
                    relative_fpath = os.path.join(os.path.basename(root), folder_name)
                    logging.info("Adding folder {} to archive.".format(absolute_fpath))
                    ziph.write(absolute_fpath, relative_fpath)
                for f in files:
                    absolute_path = os.path.join(root, f)
                    relative_path = os.path.join(os.path.basename(root), f)
                    logging.info("Adding file {} to archive.".format(absolute_path))
                    ziph.write(absolute_path, relative_path)

        logging.info("{} created successfully.".format(output_path))

    def _generate_output_file_list(self, out_dir):
        """
        _generate_output_file_list: zip result files and generate file_links for report
        """

        logging.info('Start packing result files from MDS...')

        output_files = list()

        output_dir = os.path.join(self.working_dir, str(uuid.uuid4()))
        self._mkdir_p(output_dir)
        mds_output = os.path.join(output_dir, 'metaMDS_output.zip')
        self._zip_folder(out_dir, mds_output)

        output_files.append({'path': mds_output,
                             'name': os.path.basename(mds_output),
                             'label': os.path.basename(mds_output),
                             'description': 'Output file(s) generated by metaMDS'})
        return output_files

    def _generate_mds_html_report(self, mds_outdir, n_components):

        logging.info('Start generating html report for MDS results...')
        html_report = list()

        mds_plots = list()
        for root, folders, files in os.walk(mds_outdir):
            # Find the image files by their extensions.
            for f in files:
                if re.match('^[a-zA-Z]+.*.(html)$', f):  # jpeg|jpg|bmp|png|tiff|pdf|ps|
                    absolute_path = os.path.join(root, f)
                    logging.info("Adding file {} to plot archive.".format(absolute_path))
                    mds_plots.append(absolute_path)

        result_dir = os.path.join(self.working_dir, str(uuid.uuid4()))
        self._mkdir_p(result_dir)
        result_file_path = os.path.join(result_dir, 'mds_result.html')

        visualization_content = ''

        for mds_plot in mds_plots:
            shutil.copy2(mds_plot,
                         os.path.join(result_dir, os.path.basename(mds_plot)))
            visualization_content += '<iframe height="900px" width="100%" '
            visualization_content += 'src="{}" '.format(os.path.basename(mds_plot))
            visualization_content += 'style="border:none;"></iframe>\n<p></p>\n'

        with open(result_file_path, 'w') as result_file:
            with open(os.path.join(os.path.dirname(__file__), 'templates', 'mds_template.html'),
                      'r') as report_template_file:
                report_template = report_template_file.read()
                report_template = report_template.replace('<p>Visualization_Content</p>',
                                                          visualization_content)
                report_template = report_template.replace('n_components',
                                                          '{} Components'.format(n_components))
                result_file.write(report_template)

        report_shock_id = self.dfu.file_to_shock({'file_path': result_dir,
                                                  'pack': 'zip'})['shock_id']

        html_report.append({'shock_id': report_shock_id,
                            'name': os.path.basename(result_file_path),
                            'label': os.path.basename(result_file_path),
                            'description': 'HTML summary report for MDS Matrix App'
                            })
        return html_report

    def _generate_mds_report(self, mds_ref, output_dir, workspace_name, n_components):
        logging.info('Creating MDS report...')

        output_files = self._generate_output_file_list(output_dir)
        output_html_files = self._generate_mds_html_report(output_dir, n_components)

        objects_created = list()
        objects_created.append({'ref': mds_ref,
                                'description': 'MDS Matrix'})

        report_params = {'message': '',
                         'workspace_name': workspace_name,
                         'file_links': output_files,
                         'objects_created': objects_created,
                         'html_links': output_html_files,
                         'direct_html_link_index': 0,
                         'html_window_height': 666,
                         'report_object_name': 'kb_mds_report_' + str(uuid.uuid4())}

        kbase_report_client = KBaseReport(self.callback_url)
        output = kbase_report_client.create_extended_report(report_params)

        report_output = {'report_name': output['name'], 'report_ref': output['ref']}

        return report_output

    def _get_metadata_from_obj(self, dimension,
                               associated_matrix_obj_ref, attribute_mapping_obj_ref,
                               color_marker_by, scale_size_by):
        logging.info('Retrieving metadata..')

        # build color_marker_by only
        if color_marker_by is not None:

            attr_obj = self.dfu.get_objects({'object_refs': [attribute_mapping_obj_ref]})
            attr_l = attr_obj['data'][0]['data']['attributes']

            color_index = None
            for i in range(len(attr_l)):
                if attr_l[i]['attribute'] == color_marker_by:
                    color_index = i
                    break

            color_data = []
            mdf_indx = attr_obj['data'][0]['data']['instances'].keys()
            for sample in mdf_indx:
                color_data.append(attr_obj['data'][0]['data']['instances'][sample][color_index])

            mdf = pd.DataFrame(index=mdf_indx, columns=[color_marker_by, scale_size_by])
            if color_index is not None:
                mdf[color_marker_by] = color_data

            if scale_size_by is not None:
                if associated_matrix_obj_ref is not None:
                    matrix_obj = self.dfu.get_objects({
                        'object_refs': [associated_matrix_obj_ref]})['data'][0]['data']
                    matrix_data = matrix_obj['data']

                    size_data = list()
                    if dimension == 'col':
                        size_index = matrix_data['row_ids'].index(scale_size_by)
                        for sample in mdf_indx:
                            if sample in matrix_data['col_ids']:
                                idx = matrix_data['col_ids'].index(sample)
                                size_data.append(matrix_data['values'][size_index][idx])
                            else:
                                size_data.append(None)
                    else:
                        size_index = matrix_data['col_ids'].index(scale_size_by)
                        for sample in mdf_indx:
                            if sample in matrix_data['row_ids']:
                                idx = matrix_data['row_ids'].index(sample)
                                size_data.append(matrix_data['values'][idx][size_index])
                            else:
                                size_data.append(None)

                    mdf[scale_size_by] = size_data
                else:
                    size_index = None
                    for i in range(len(attr_l)):
                        if attr_l[i]['attribute'] == scale_size_by:
                            size_index = i
                            break

                    size_data = []
                    for sample in mdf_indx:
                        try:
                            size_data.append(
                                float(
                                    attr_obj['data'][0]['data']['instances'][sample][size_index]))
                        except Exception:
                            logging.info(
                                'ERROR: scaling is not int or float. scaling has been dropped')
                            scale_size_by = None
                            size_index = None

                    if size_index is not None:
                        mdf[scale_size_by] = size_data
        # build scale_size_by only
        else:
            if associated_matrix_obj_ref is not None:
                matrix_data = self.dfu.get_objects({
                    'object_refs': [associated_matrix_obj_ref]})['data'][0]['data']['data']
                if dimension == 'col':
                    size_index = matrix_data['row_ids'].index(scale_size_by)
                    size_data = matrix_data['values'][size_index]

                    mdf = pd.DataFrame(index=matrix_data['col_ids'],
                                       columns=[color_marker_by, scale_size_by])
                    mdf[scale_size_by] = size_data

                else:
                    size_index = matrix_data['col_ids'].index(scale_size_by)
                    size_data = list()
                    for value in matrix_data['values']:
                        size_data.append(value[size_index])

                    mdf = pd.DataFrame(index=matrix_data['col_ids'],
                                       columns=[color_marker_by, scale_size_by])
                    mdf[scale_size_by] = size_data
            else:
                attr_obj = self.dfu.get_objects({'object_refs': [attribute_mapping_obj_ref]})
                attr_l = attr_obj['data'][0]['data']['attributes']

                size_index = None
                for i in range(len(attr_l)):
                    if attr_l[i]['attribute'] == scale_size_by:
                        size_index = i
                        break

                size_data = []
                mdf_indx = attr_obj['data'][0]['data']['instances'].keys()
                for sample in mdf_indx:
                    try:
                        size_data.append(float(
                            attr_obj['data'][0]['data']['instances'][sample][size_index]))
                    except Exception:
                        err_msg = 'ERROR: scaling is not int or float. scaling has been dropped'
                        logging.info(err_msg)
                        scale_size_by = None
                        size_index = None

                mdf = pd.DataFrame(index=mdf_indx,
                                   columns=[color_marker_by, scale_size_by])
                if size_index is not None:
                    mdf[scale_size_by] = size_data

        logging.info('created metadata df:\n{}'.format(mdf))

        return mdf

    def _get_metadata_from_file(self, metadata_file, color_marker_by, scale_size_by):
        """
        Get metadata from file and return simplified pd.DataFrame
        :return:
        """

        logging.info('Retrieving metadata..')

        mdf = pd.read_csv(metadata_file, sep='\t', index_col=0)

        logging.info('MDF: {}'.format(mdf))

        mdf = mdf[[color_marker_by, scale_size_by]]

        return mdf

    def _plot_without_grouping(self, dimension):

        # Get site data from previously saved file
        site_ordin_df = pd.read_csv(os.path.join(self.output_dir, "site_ordination.csv"),
                                    index_col=0)
        logging.info('SITE_ORDIN_DF:\n {}'.format(site_ordin_df))
        site_ordin_df.fillna('na', inplace=True)

        fig = px.scatter(site_ordin_df, x="MDS1", y="MDS2", hover_name=site_ordin_df.index)

        # Save plotly_fig.html and return path
        plotly_html_file_path = os.path.join(self.output_dir, "plotly_fig.html")
        plot(fig, filename=plotly_html_file_path)
        return plotly_html_file_path

    def _plot_with_grouping(self, dimension, associated_matrix_obj_ref, attribute_mapping_obj_ref,
                            metadata_file, color_marker_by, scale_size_by, highlight,
                            only_highlight):
        logging.info('Plotting with grouping: "{}", and "{}"'.format(color_marker_by,
                                                                     scale_size_by))

        # Both can not be the same right now.. mdf is now new pd would lead to problems
        if color_marker_by == scale_size_by:
            logging.info('ERROR: both color and scale are same field. scale set to None')
            scale_size_by = None

        if (attribute_mapping_obj_ref is not None or
                associated_matrix_obj_ref is not None):
            mdf = self._get_metadata_from_obj(dimension,
                                              associated_matrix_obj_ref,
                                              attribute_mapping_obj_ref,
                                              color_marker_by,
                                              scale_size_by)
        elif metadata_file is not None:
            mdf = self._get_metadata_from_file(metadata_file, color_marker_by, scale_size_by)
        else:
            raise ValueError('No metadata file was specified')

        grouping_meta_file = os.path.join(self.output_dir, 'grouping_meta.csv')
        with open(grouping_meta_file, 'w') as m_file:
            mdf.to_csv(m_file, sep='\t')

        # Get site data from previously saved file
        site_ordin_file = os.path.join(self.output_dir, "site_ordination.csv")

        if not os.path.exists(site_ordin_file):
            raise ValueError('failed to generate metaMDS points')
        site_ordin_df = pd.read_csv(site_ordin_file, index_col=0)
        logging.info('SITE_ORDIN_DF:\n {}'.format(site_ordin_df))

        # Check if metadata file is valid for this method
        for sample in site_ordin_df.index:
            try:
                mdf.loc[sample]
            except KeyError:
                raise KeyError('One or more samples in site_ordination is not found in chosen '
                               'metadata obj. If you ran this using files, you might need to '
                               'transpose the data in your files so samples are rows and OTU '
                               'are columns.')

        # Fill site_ordin_df with metadata from mdf
        site_ordin_df['color'] = None
        site_ordin_df['size'] = None
        for ID in site_ordin_df.index:
            site_ordin_df['color'].loc[ID] = mdf[color_marker_by].loc[ID]
            site_ordin_df['size'].loc[ID] = mdf[scale_size_by].loc[ID]

        site_ordin_df.fillna('na', inplace=True)

        # Plot
        if color_marker_by is not None and scale_size_by is not None and all(
                isinstance(x, (int, float)) for x in list(site_ordin_df['size'])):
            fig = px.scatter(site_ordin_df, x="MDS1", y="MDS2", color="color", size="size",
                             hover_name=site_ordin_df.index)
        elif color_marker_by is not None:
            fig = px.scatter(site_ordin_df, x="MDS1", y="MDS2", color="color",
                             hover_name=site_ordin_df.index)
        elif scale_size_by is not None:
            fig = px.scatter(site_ordin_df, x="MDS1", y="MDS2", size="size",
                             hover_name=site_ordin_df.index)

        # add vectors
        vector_file = os.path.join(self.output_dir, "vectors.csv")

        if os.path.exists(vector_file):
            vector_df = pd.read_csv(vector_file, index_col=0)
            logging.info('VECTOR_DF:\n {}'.format(vector_df))
            loading_x, loading_y, loading_text = list(), list(), list()
            highlight_x, highlight_y, highlight_text = list(), list(), list()
            for idx, row in vector_df.iterrows():
                x, y, name = row[0], row[1], idx

                if name in highlight:
                    highlight_x.extend([0, x])
                    highlight_y.extend([0, y])
                    highlight_text.extend(['0', name])

                    fig.add_annotation(x=x, y=y, ax=0, ay=0, xanchor="center", yanchor="bottom",
                                       text=name, font=dict(color="mediumvioletred"))
                else:
                    loading_x.extend([0, x])
                    loading_y.extend([0, y])
                    loading_text.extend(['0', name])

            if not (highlight and only_highlight):
                fig.add_trace(go.Scatter(
                                x=loading_x,
                                y=loading_y,
                                mode="lines+markers",
                                name="environmental vectors",
                                text=loading_text,
                                textposition="bottom center",
                                line=dict(color="RoyalBlue", width=0.5)
                            ))

            fig.add_trace(go.Scatter(
                            x=highlight_x,
                            y=highlight_y,
                            mode="lines+markers",
                            name="selected environmental vectors",
                            text=highlight_text,
                            textposition="bottom center",
                            line=dict(color="mediumvioletred", width=1.5)
                        ))

        # Save plotly_fig.html and return path
        plotly_html_file_path = os.path.join(self.output_dir, "plotly_fig.html")
        plot(fig, filename=plotly_html_file_path)
        return plotly_html_file_path

    def __init__(self, config):

        self.ws_url = config["workspace-url"]
        self.callback_url = config['SDK_CALLBACK_URL']
        self.token = config['KB_AUTH_TOKEN']
        self.scratch = config['scratch']
        self.working_dir = self.scratch
        self.dfu = DataFileUtil(self.callback_url)
        self.output_dir = os.path.join(self.working_dir, self.MDS_OUT_DIR)
        self._mkdir_p(self.output_dir)

    def run_metaMDS(self, params):
        """
        run_metaMDS: perform metaMDS analysis on matrix
        :param input_obj_ref: object reference of a matrix
        :param workspace_name: the name of the workspace
        :param mds_matrix_name: name of MDS (KBaseExperiments.MDSMatrix) object
        :param n_components - dimentionality of the reduced space (default 2)
        :param max_iter: maximum iterations allowed
        :param metric: indication of running metric or non-metric MDS
        :param distance_metric: distance the ordination will be performed on, default to "bray"
        """

        logging.info('--->\nrunning metaMDS with input\n' +
                     'params:\n{}'.format(json.dumps(params, indent=1)))

        self._validate_run_mds_params(params)

        input_obj_ref = params.get(self.PARAM_IN_MATRIX)
        workspace_name = params.get(self.PARAM_IN_WS)
        mds_matrix_name = params.get(self.PARAM_OUT_MATRIX)
        n_components = int(params.get('n_components', 2))
        dimension = params.get('dimension', 'col')

        res = self.dfu.get_objects({'object_refs': [input_obj_ref]})['data'][0]
        obj_data = res['data']
        obj_name = res['info'][1]
        obj_type = res['info'][2]

        max_size = len(obj_data['data']['col_ids'])
        if n_components > max_size:
            raise ValueError('Number of components should be less than number of samples')

        exitCode = -99
        if "KBaseMatrices" in obj_type:
            # create the input file from obj_data
            matrix_tab = obj_data['data']['values']
            row_ids = obj_data['data']['row_ids']
            col_ids = obj_data['data']['col_ids']
            matrix_df = pd.DataFrame(matrix_tab, index=row_ids, columns=col_ids)
            # Transpose DataFrame
            if dimension == 'col':
                matrix_df = matrix_df.T
            elif dimension != 'row':
                err_msg = 'Input dimension [{}] is not available.\n'.format(dimension)
                err_msg += 'Please choose either "col" or "row"'
                raise ValueError(err_msg)
            logging.info('input matrix:\n {}'.format(matrix_df))

            all_zero_rows = matrix_df.index[(matrix_df == 0).all(1)].tolist()

            if all_zero_rows:
                err_msg = 'Please alter the input matrix so that no observation has all 0s\n'
                err_msg += 'Observations with all zero values:\n{}'.format(all_zero_rows)
                raise ValueError(err_msg)

            matrix_data_file = os.path.join(self.output_dir, obj_name + '.csv')
            with open(matrix_data_file, 'w') as m_file:
                matrix_df.to_csv(m_file, sep='\t')

            params['datafile'] = matrix_data_file
            exitCode = self.run_metaMDS_with_file(params, matrix_df)
        else:
            err_msg = 'Ooops! [{}] is not supported.\n'.format(obj_type)
            err_msg += 'Please provide a KBaseMatrices object'
            raise ValueError(err_msg)

        if exitCode == -99:
            raise ValueError('Caught subprocess.CalledProcessError while calling R.')

        # saving the mds_matrix object
        # read metaMDS results from files into data frames
        dist_matrix_df = pd.read_csv(os.path.join(self.output_dir, "dist_matrix.csv"))
        mds_params_df = pd.read_json(os.path.join(self.output_dir, "others.json"))
        site_ordin_df = pd.read_csv(os.path.join(self.output_dir, "site_ordination.csv"))
        species_ordin_df = pd.read_csv(os.path.join(self.output_dir, "species_ordination.csv"))

        mds_ref = self._save_mds_matrix(workspace_name, input_obj_ref, mds_matrix_name,
                                        dist_matrix_df, mds_params_df, site_ordin_df,
                                        species_ordin_df)
        return_val = {'mds_ref': mds_ref}

        # generating report
        report_output = self._generate_mds_report(mds_ref, self.output_dir,
                                                  workspace_name, n_components)

        return_val.update(report_output)
        return return_val

    def run_metaMDS_with_file(self, params, matrix_df):
        """
        run_metaMDS_with_file: perform metaMDS analysis on matrix
        :param datafile: a file that contains the matrix data
        :param workspace_name: the name of the workspace
        :param mds_matrix_name: name of MDS (KBaseExperiments.MDSMatrix) object
        :param n_components - dimentionality of the reduced space (default 2)
        :param max_iter: maximum iterations allowed
        :param metric: indication of running metric or non-metric MDS
        :param distance_metric: distance the ordination will be performed on, default to "bray"
        """

        # Variables for Grouping Features

        logging.info('--->\nrunning run_metaMDS_with_file with input\n' +
                     'params:\n{}'.format(json.dumps(params, indent=1)))

        attribute_mapping_obj_ref = params.get('attribute_mapping_obj_ref')
        metadata_file = params.get('metadata_file')
        color_marker_by = params.get('color_marker_by')
        associated_matrix_obj_ref = params.get('associated_matrix_obj_ref')

        dimension = params.get('dimension', 'col')

        if color_marker_by is not None:
            try:
                color_marker_by = color_marker_by['attribute_color'][0]
            except KeyError:
                raise KeyError('Expected dictionary with key "attribute_color" containing a list '
                               'of one element. Instead found: {}'.format(color_marker_by))
        scale_size_by = params.get('scale_size_by')
        highlight = list()
        if scale_size_by is not None:
            if scale_size_by.get('attribute_size'):
                scale_size_by = scale_size_by['attribute_size'][0]
            elif scale_size_by.get('row_size'):
                highlight = scale_size_by.get('highlight_row', list())
                scale_size_by = scale_size_by['row_size'][0]
                if dimension != 'col':
                    err_msg = 'Please choose Column dimension in order for the plot size to be '
                    err_msg += 'associated with Matrix row'
                    raise ValueError(err_msg)

            elif scale_size_by.get('col_size'):
                highlight = scale_size_by.get('highlight_col', list())
                scale_size_by = scale_size_by['col_size'][0]

                if dimension != 'row':
                    err_msg = 'Please choose Row dimension in order for the plot size to be '
                    err_msg += 'associated with Matrix column'
                    raise ValueError(err_msg)
            else:
                raise KeyError('Expected dictionary with key "attribute_size" containing a list '
                               'of one element. Instead found: {}'.format(scale_size_by))

        if associated_matrix_obj_ref:

            associated_matrix_obj = self.dfu.get_objects({
                            'object_refs': [associated_matrix_obj_ref]})['data'][0]

            associated_matrix_data = associated_matrix_obj['data']
            associated_matrix_name = associated_matrix_obj['info'][1]

            values = associated_matrix_data['data']['values']
            row_ids = associated_matrix_data['data']['row_ids']
            col_ids = associated_matrix_data['data']['col_ids']
            associated_matrix_df = pd.DataFrame(values, index=row_ids, columns=col_ids)
            # Transpose DataFrame
            if dimension == 'col':
                associated_matrix_df = associated_matrix_df.T
            logging.info('input associated matrix:\n {}'.format(associated_matrix_df))

            common_idx = associated_matrix_df.index.intersection(matrix_df.index)

            if len(common_idx):
                associated_matrix_df = associated_matrix_df.loc[common_idx]
                associated_matrix_data_file = os.path.join(self.output_dir,
                                                           associated_matrix_name + '.csv')
                with open(associated_matrix_data_file, 'w') as m_file:
                    associated_matrix_df.to_csv(m_file, sep='\t')

                params['associated_matrix_file'] = associated_matrix_data_file
            else:
                logging.info('associated matrix and input matrix share no common observations')

        rscrpt_file = self._build_rMDS_script(params)
        logging.info('--->\nR script file has been written to {}'.format(rscrpt_file))

        result = self._execute_r_script(rscrpt_file)

        # Make and save plotly fig
        if color_marker_by is not None or scale_size_by is not None:
            self._plot_with_grouping(dimension,
                                     associated_matrix_obj_ref,
                                     attribute_mapping_obj_ref,
                                     metadata_file,
                                     color_marker_by, scale_size_by,
                                     highlight, params.get('only_highlight', True))
        else:
            self._plot_without_grouping(dimension)

        return result
