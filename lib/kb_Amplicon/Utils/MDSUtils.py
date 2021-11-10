
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
import plotly.offline as pyo
import plotly.graph_objs as go
import plotly.express as px
from plotly.offline import plot

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
        data_file_path = params.get('datafile', None)
        if not data_file_path:
            return ''

        exists = os.path.isfile(os.path.join(self.output_dir, os.path.basename(data_file_path)))
        if not exists:
            shutil.copyfile(data_file_path,
                            os.path.join(self.output_dir, os.path.basename(data_file_path)))

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
        mds_scrpt += 'func_call <- vg_data.mds$call\n'
        mds_scrpt += 'mds_data <- vg_data.mds$data\n'

        # save the results to the current dir
        # Write CSV in R
        mds_scrpt += 'write.csv(dist_matrix,file="dist_matrix.csv",row.names=TRUE,na="")\n'
        mds_scrpt += 'write.csv(variableScores,file="species_ordination.csv",' + \
                     'row.names=TRUE,na="")\n'
        mds_scrpt += 'write.csv(sampleScores,file="site_ordination.csv",row.names=TRUE,na="")\n'

        # Write JSON in R
        mds_scrpt += 'write_json(toJSON(dist_matrix),path="dist_matrix.json",pretty=TRUE,' + \
                     'auto_unbox=FALSE)\n'
        mds_scrpt += 'write_json(toJSON(variableScores),path="species_ordination.json",' + \
                     'pretty=TRUE,auto_unbox=FALSE)\n'
        mds_scrpt += 'write_json(toJSON(sampleScores),path="site_ordination.json",' + \
                     'pretty=TRUE,auto_unbox=FALSE)\n'
        mds_scrpt += 'item_name=c("stress","distance_metric","dist_call","converged",' + \
                     '"dimesions","trials","maxits")\n'
        mds_scrpt += 'item_value=c(stress,dist_metric,dist_call,converged,dims,tries,maxits)\n'
        mds_scrpt += 'df <- data.frame(item_name,item_value,stringsAsFactors=FALSE)\n'
        mds_scrpt += 'write_json(toJSON(df),path="others.json",pretty=TRUE,auto_unbox=FALSE)\n'

        # save mds plots
        '''
        mds_scrpt += 'bmp(file="saving_mds_plot.bmp",width=580,height=580,units="px",' + \
                     'res=100, pointsize=12)\n'
        mds_scrpt += 'plot(vg_data.mds,type="n",display="sites")\n'
        mds_scrpt += 'points(vg_data.mds)\n'
        mds_scrpt += 'dev.off()\n'
        mds_scrpt += 'pdf(file="saving_mds_plot.pdf",width=6,height=6)\n'
        mds_scrpt += 'plot(vg_data.mds,type="n",display="sites")\n'
        mds_scrpt += 'points(vg_data.mds)\n'
        mds_scrpt += 'dev.off()\n'
        mds_scrpt += 'pdf(file="mds_plot_withlabel.pdf",width=6,height=6)\n'
        mds_scrpt += 'plot(vg_data.mds,type="n",display="sites")\n'
        mds_scrpt += 'ordilabel(vg_data.mds,dis="sites",cex=1.2,font=3,fill="hotpink",col="blue")\n'
        mds_scrpt += 'dev.off()\n'
        mds_scrpt += 'pdf(file="mds_plot_withcolor.pdf",width=6,height=6)\n'
        mds_scrpt += 'fig <- ordiplot(vg_data.mds,type="none")\n'
        mds_scrpt += 'points(fig,"sites",pch=21,col="red",bg="yellow")\n'
        mds_scrpt += 'points(fig,"species",pch=21,col="green",bg="blue")\n'
        # mds_scrpt += 'text(fig, "species", col="blue", cex=0.9)\n'
        mds_scrpt += 'dev.off()\n'
        '''
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
                logging.info('Error running command: {} Exit Code: '.format(
                    ' '.join(rcmd), str(exitCode)))
                logging.info('\n{}'.format(complete_proc.stderr))
        except subprocess.CalledProcessError as sub_e:
            exitCode = -99
            logging.info('Caught subprocess.CalledProcessError {}'.format(sub_e))

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

    def _mds_df_to_excel(self, mds_df, distance_df, result_dir, mds_matrix_ref):
        """
        write MDS matrix df into excel
        """
        logging.info('writting mds data frame to excel file')
        mds_matrix_obj = self.dfu.get_objects({'object_refs': [mds_matrix_ref]})['data'][0]
        mds_matrix_info = mds_matrix_obj['info']
        mds_matrix_name = mds_matrix_info[1]

        file_path = os.path.join(result_dir, mds_matrix_name + ".xlsx")
        writer = pd.ExcelWriter(file_path)

        mds_df.to_excel(writer, "mds_matrix", index=True)
        if distance_df:
            distance_df.to_excel(writer, "mds_distance_matrix", index=True)

        writer.close()

    def _Matrix2D_to_df(self, Matrix2D):
        """
        _Matrix2D_to_df: transform a FloatMatrix2D to data frame
        """

        index = Matrix2D.get('row_ids')
        columns = Matrix2D.get('col_ids')
        values = Matrix2D.get('values')

        df = pd.DataFrame(values, index=index, columns=columns)

        return df

    def _mds_to_df(self, mds_matrix_ref):
        """
        retrieve MDS matrix ws object to mds_df
        """
        logging.info('converting mds matrix to data frame')
        mds_data = self.dfu.get_objects({'object_refs': [mds_matrix_ref]})['data'][0]['data']

        rotation_matrix_data = mds_data.get('rotation_matrix')
        distance_matrix_data = mds_data.get('distance_matrix')
        original_matrix_ref = mds_data.get('original_matrix_ref')
        dimension = mds_data.get('mds_parameters').get('n_components')

        mds_df = self._Matrix2D_to_df(rotation_matrix_data)
        distance_df = None
        if distance_matrix_data:
            distance_df = self._Matrix2D_to_df(distance_matrix_data)

        if original_matrix_ref:
            logging.info('appending instance group information to mds data frame')
            obj_data = self.dfu.get_objects(
                {'object_refs': [original_matrix_ref]})['data'][0]['data']

            attributemapping_ref = obj_data.get('{}_attributemapping_ref'.format(dimension))

            am_data = self.dfu.get_objects(
                {'object_refs': [attributemapping_ref]})['data'][0]['data']

            attributes = am_data.get('attributes')
            instances = am_data.get('instances')
            am_df = pd.DataFrame(data=list(instances.values()),
                                 columns=list(map(lambda x: x.get('attribute'), attributes)),
                                 index=instances.keys())

            mds_df = mds_df.merge(am_df, left_index=True, right_index=True, how='left',
                                  validate='one_to_one')

        return mds_df, distance_df

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

    def _get_metadata_from_obj(self):
        """
        Get metadata from obj and return simplified pd.DataFrame
        :return:
        """

        logging.info('Retrieving metadata..')

        # KBase obj data
        mdf = self.dfu.get_objects({'object_refs': [self.attribute_mapping_obj_ref]})
        attr_l = mdf['data'][0]['data']['attributes']

        # Get index location in mdf(metadata) of chosen color and scale
        color_index = None
        size_index = None
        for i in range(len(attr_l)):
            if attr_l[i]['attribute'] == self.color_marker_by:
                color_index = i
            if attr_l[i]['attribute'] == self.scale_size_by:
                size_index = i

        # Make list of color and scale data
        color_data = []
        size_data = []
        mdf_indx = mdf['data'][0]['data']['instances'].keys()
        for sample in mdf_indx:
            if color_index is not None:
                color_data.append(mdf['data'][0]['data']['instances'][sample][color_index])
            if size_index is not None:
                try:
                    size_data.append(float(mdf['data'][0]['data']['instances'][sample][size_index]))
                except:
                    logging.info('ERROR: scaling is not int or float. scaling has been dropped')
                    self.scale_size_by = None
                    size_index = None

        # mdf is now new pd.DataFrame that only includes needed data
        mdf = pd.DataFrame(index=mdf_indx, columns=[self.color_marker_by, self.scale_size_by])
        if color_index is not None:
            mdf[self.color_marker_by] = color_data
        if size_index is not None:
            mdf[self.scale_size_by] = size_data

        return mdf

    def _get_metadata_from_file(self):
        """
        Get metadata from file and return simplified pd.DataFrame
        :return:
        """

        logging.info('Retrieving metadata..')

        mdf = pd.read_csv(self.metadata_file, sep='\t', index_col=0)

        logging.info('MDF: {}'.format(mdf))

        mdf = mdf[[self.color_marker_by, self.scale_size_by]]

        return mdf

    def _plot_with_grouping(self):
        logging.info('Plotting with grouping: "{}", and "{}"'.format(self.color_marker_by, self.scale_size_by))

        # Both can not be the same right now.. mdf is now new pd would lead to problems
        if self.color_marker_by == self.scale_size_by:
            logging.info('ERROR: both color and scale are same field. scale set to None')
            self.scale_size_by = None

        if self.attribute_mapping_obj_ref is not None:
            mdf = self._get_metadata_from_obj()
        elif self.metadata_file is not None:
            mdf = self._get_metadata_from_file()
        else:
            FileNotFoundError('No metadata file was specified')

        # Get site data from previously saved file
        site_ordin_df = pd.read_csv(os.path.join(self.output_dir, "site_ordination.csv"), index_col=0)
        logging.info('SITE_ORDIN_DF:\n {}'.format(site_ordin_df))

        # Check if metadata file is valid for this method
        for sample in site_ordin_df.index:
            try:
                mdf.loc[sample]
            except KeyError:
                raise KeyError('One or more samples in site_ordination is not found in chosen metadata obj. If you ran '
                               'this using files, you might need to transpose the data in your files so samples are '
                               'rows and OTU are columns.')

        # Fill site_ordin_df with metadata from mdf
        site_ordin_df['color'] = None
        site_ordin_df['size'] = None
        for ID in site_ordin_df.index:
            site_ordin_df['color'].loc[ID] = mdf[self.color_marker_by].loc[ID]
            site_ordin_df['size'].loc[ID] = mdf[self.scale_size_by].loc[ID]

        site_ordin_df.fillna('na', inplace=True)

        # Plot
        if self.color_marker_by is not None and self.scale_size_by is not None and all(
                isinstance(x, (int, float)) for x in list(site_ordin_df['size'])):
            fig = px.scatter(site_ordin_df, x="MDS1", y="MDS2", color="color", size="size",
                             hover_name=site_ordin_df.index)
        elif self.color_marker_by is not None:
            fig = px.scatter(site_ordin_df, x="MDS1", y="MDS2", color="color", hover_name=site_ordin_df.index)
        elif self.scale_size_by is not None:
            fig = px.scatter(site_ordin_df, x="MDS1", y="MDS2", size="size", hover_name=site_ordin_df.index)

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

        # If input is from files, then pd.DataFrame needs to be transposed in run_metaMDS_with_file method
        self.need_to_transpose = True


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
            matrix_df = matrix_df.T
            self.need_to_transpose = False

            matrix_data_file = os.path.join(self.output_dir, obj_name + '.csv')
            with open(matrix_data_file, 'w') as m_file:
                matrix_df.to_csv(m_file, sep='\t')

            params['datafile'] = matrix_data_file
            exitCode = self.run_metaMDS_with_file(params)
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
        returnVal = {'mds_ref': mds_ref}

        # generating report
        report_output = self._generate_mds_report(mds_ref, self.output_dir,
                                                  workspace_name, n_components)

        returnVal.update(report_output)
        return returnVal

    def run_metaMDS_with_file(self, params):
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
        self.attribute_mapping_obj_ref = params.get('attribute_mapping_obj_ref')
        self.metadata_file = params.get('metadata_file')
        self.color_marker_by = params.get('color_marker_by')
        if self.color_marker_by is not None:
            try:
                self.color_marker_by = self.color_marker_by['attribute_color'][0]
            except KeyError:
                raise KeyError('Expected dictionary with key "attribute_color" containing a list of one element. '
                               'Instead found: {}'.format(self.color_marker_by))
        self.scale_size_by = params.get('scale_size_by')
        if self.scale_size_by is not None:
            try:
                self.scale_size_by = self.scale_size_by['attribute_size'][0]
            except KeyError:
                raise KeyError('Expected dictionary with key "attribute_size" containing a list of one element. '
                               'Instead found: {}'.format(self.scale_size_by))

        logging.info('--->\nrunning metaMDS with input \n' +
                     'params:\n{}'.format(json.dumps(params, indent=1)))

        rscrpt_file = self._build_rMDS_script(params)
        logging.info('--->\nR script file has been written to {}'.format(rscrpt_file))

        result = self._execute_r_script(rscrpt_file)

        # Make and save plotly fig
        if self.color_marker_by is not None or self.scale_size_by is not None:
            self._plot_with_grouping()

        return result

    def export_mds_matrix_excel(self, params):
        """
        export MDSMatrix as Excel
        """
        logging.info('start exporting mds matrix')
        mds_matrix_ref = params.get('input_ref')

        mds_df, components_df = self._mds_to_df(mds_matrix_ref)

        result_dir = os.path.join(self.scratch, str(uuid.uuid4()))
        self._mkdir_p(result_dir)

        self._mds_df_to_excel(mds_df, components_df, result_dir, mds_matrix_ref)

        package_details = self.dfu.package_for_download({
            'file_path': result_dir,
            'ws_refs': [mds_matrix_ref]
        })

        return {'shock_id': package_details['shock_id']}
