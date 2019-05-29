import errno
import itertools
import json
import logging
import os
import shutil
import sys
import uuid
import zipfile
import re
import subprocess

import pandas as pd
import numpy as np
import plotly.graph_objs as go
from matplotlib import pyplot as plt
from plotly.offline import plot
from sklearn.manifold import MDS
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
# from skbio.stats.distance import DistanceMatrix

from kb_Vegan.Utils.DataUtil import DataUtil
from installed_clients.DataFileUtilClient import DataFileUtil
from installed_clients.KBaseReportClient import KBaseReport


class MDSUtils:

    R_BIN = '/kb/deployment/bin'
    VEGAN_OUT_DIR = 'Vegan_output'
    PARAM_IN_WS = 'output_workspace'
    PARAM_IN_OTU_FILE = 'otu_file'
    PARAM_IN_MATRIX = 'matrix_ref'

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
        for p in ['input_obj_ref', 'workspace_name', 'mds_matrix_name']:
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
        mds_scrpt += 'vg_data <- read.table("' + data_file_path + '",header=TRUE,row.names=1,sep="")\n'
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
        mds_scrpt += 'write_json(dist_matrix,path="dist_matrix.json",pretty=TRUE,auto_unbox=FALSE)\n'
        mds_scrpt += 'write.csv(variableScores,file="species_ordination.csv",row.names=TRUE,na="")\n'
        mds_scrpt += 'write_json(variableScores,path="species_ordination.json",pretty=TRUE,auto_unbox=FALSE)\n'
        mds_scrpt += 'write.csv(sampleScores,file="site_ordination.csv",row.names=TRUE,na="")\n'
        mds_scrpt += 'write_json(sampleScores,path="site_ordination.json",pretty=TRUE,auto_unbox=FALSE)\n'
        mds_scrpt += 'item_name=c("stress","distance_metric","dist_call","converged","dimesions","trials","maxits")\n'
        mds_scrpt += 'item_value=c(stress,dist_metric,dist_call,converged,dims,tries,maxits)\n' 
        mds_scrpt += 'df <- data.frame(item_name,item_value,stringsAsFactors=FALSE)\n'
        mds_scrpt += 'write_json(toJSON(df),path="others.json",pretty=TRUE,auto_unbox=FALSE)\n'

        # save mds plot
        mds_scrpt += 'bmp(file="saving_mds_plot.bmp",width=6,height=4,units="in",res=100)\n'
        mds_scrpt += 'plot(vg_data.mds,type="n",display="sites")\n'
        mds_scrpt += 'points(vg_data.mds)\n'
        mds_scrpt += 'dev.off()\n'
        mds_scrpt += 'pdf(file="saving_mds_plot.pdf",width=6,height=4)\n'
        mds_scrpt += 'plot(vg_data.mds,type="n",display="sites")\n'
        mds_scrpt += 'points(vg_data.mds)\n'
        mds_scrpt += 'dev.off()\n'

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
                logging.info('\n' + ' '.join(rcmd) + ' was executed successfully, exit code was: ' +
                    str(exitCode))
                logging.info("Finished calling R.")
            else:
                logging.info('Error running command: ' + ' '.join(rcmd) + 'Exit Code: ' +
                    str(exitCode))
                logging.info('\n{}'.format(complete_proc.stderr))
        except subprocess.CalledProcessError as subE:
            pass
            exitCode = -99
            logging.info('Caught subprocess.CalledProcessError {}'.format(subE))

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

        mds_df = self._Matrix2D_to_df(rotation_matrix_data)
        distance_df = None
        if distance_matrix_data:
            distance_df = self._Matrix2D_to_df(distance_matrix_data)

        if original_matrix_ref:
            logging.info('appending instance group information to mds data frame')
            obj_data = self.dfu.get_objects({'object_refs': [original_matrix_ref]})['data'][0]['data']

            attributemapping_ref = obj_data.get('{}_attributemapping_ref'.format(dimension))

            am_data = self.dfu.get_objects({'object_refs': [attributemapping_ref]})['data'][0]['data']

            attributes = am_data.get('attributes')
            instances = am_data.get('instances')
            am_df = pd.DataFrame(data=list(instances.values()),
                                 columns=list(map(lambda x: x.get('attribute'), attributes)),
                                 index=instances.keys())

            mds_df = mds_df.merge(am_df, left_index=True, right_index=True, how='left',
                                  validate='one_to_one')

        return mds_df, distance_df

    def _save_mds_matrix(self, workspace_name, input_obj_ref, mds_matrix_name, rotation_matrix_df,
                         distance_df, n_components, distance_metric):

        logging.info('saving MDSMatrix')

        if not isinstance(workspace_name, int):
            ws_name_id = self.dfu.ws_name_to_id(workspace_name)
        else:
            ws_name_id = workspace_name

        mds_data = {}

        mds_data.update({'rotation_matrix': self._df_to_list(rotation_matrix_df)})
        mds_data.update({'distance_matrix': self._df_to_list(distance_df)})
        mds_data.update({'mds_parameters': {'n_components': str(n_components),
                                            'distance_metric': str(distance_metric)}})
        mds_data.update({'original_matrix_ref': input_obj_ref})

        obj_type = 'MDSMatrix'
        info = self.dfu.save_objects({
            "id": ws_name_id,
            "objects": [{
                "type": obj_type,
                "data": mds_data,
                "name": mds_matrix_name
            }]
        })[0]

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
                    logging.info("Adding {} to archive.".format(absolute_fpath))
                    ziph.write(absolute_fpath, relative_fpath)
                for f in files:
                    absolute_path = os.path.join(root, f)
                    relative_path = os.path.join(os.path.basename(root), f)
                    logging.info("Adding {} to archive.".format(absolute_path))
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

        result_dir = os.path.join(self.working_dir, str(uuid.uuid4()))
        self._mkdir_p(result_dir)
        result_file_path = os.path.join(result_dir, 'mds_result.html')

        mds_plots = list()
        for root, folders, files in os.walk(mds_outdir):
            # Find the image files by their extensions.
            for f in files:
                if re.match("^[a-zA-Z]+.*\.(jpeg|jpg|bmp|tiff|pdf|ps)$", f):
                # for p in ['.jpeg', '.jpg', '.bmp', '.tiff','.pdf', '.ps']:
                #    if p not in f:
                    absolute_path = os.path.join(root, f)
                    logging.info("Adding {} to plot archive.".format(absolute_path))
                    mds_plots.append(absolute_path)

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

    def _generate_mds_report(self, mds_ref, mds_plots, workspace_name, n_components):
        logging.info('creating report')

        output_html_files = self._generate_mds_html_report(mds_plots, n_components)

        objects_created = list()
        objects_created.append({'ref': mds_ref,
                                'description': 'MDS Matrix'})

        report_params = {'message': '',
                         'workspace_name': workspace_name,
                         'objects_created': objects_created,
                         'html_links': output_html_files,
                         'direct_html_link_index': 0,
                         'html_window_height': 666,
                         'report_object_name': 'kb_mds_report_' + str(uuid.uuid4())}

        kbase_report_client = KBaseReport(self.callback_url)
        output = kbase_report_client.create_extended_report(report_params)

        report_output = {'report_name': output['name'], 'report_ref': output['ref']}

        return report_output

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

    def _append_instance_group(self, plot_mds_matrix, obj_data, dimension):
        plot_mds_matrix = plot_mds_matrix.copy()

        if dimension == 'row':
            attribute_mapping = obj_data.get('row_mapping')
        elif dimension == 'col':
            attribute_mapping = obj_data.get('col_mapping')
        else:
            raise ValueError('Unexpected dimension')

        if not attribute_mapping:
            logging.warning('Matrix object does not have {}_mapping attribute'.format(dimension))
            # build matrix with unify color and shape
            return plot_mds_matrix
        else:
            # append instance col mapping from row/col_mapping
            plot_mds_matrix['instance'] = plot_mds_matrix.index.map(attribute_mapping)

        return plot_mds_matrix

    def _build_size_mds_matrix(self, plot_mds_matrix, obj_data, dimension, attribute_name):
        """
        _build_size_mds_matrix: append attribute value to rotation_matrix
        """
        logging.info('appending attribute value for sizing to rotation matrix')

        plot_mds_matrix = plot_mds_matrix.copy()

        if dimension == 'row':
            attribute_mapping = obj_data.get('row_mapping')
            attribute_mapping_ref = obj_data.get('row_attributemapping_ref')
        elif dimension == 'col':
            attribute_mapping = obj_data.get('col_mapping')
            attribute_mapping_ref = obj_data.get('col_attributemapping_ref')
        else:
            raise ValueError('Unexpected dimension')

        if not attribute_mapping:
            logging.warning('Matrix object does not have {}_mapping attribute'.format(dimension))
            # build matrix with unify color and shape
            return plot_mds_matrix
        else:
            # append instance col mapping from row/col_mapping
            plot_mds_matrix['instance'] = plot_mds_matrix.index.map(attribute_mapping)

        res = self.dfu.get_objects({'object_refs': [attribute_mapping_ref]})['data'][0]
        attri_data = res['data']
        attri_name = res['info'][1]

        attributes = attri_data.get('attributes')

        attr_pos = None
        for idx, attribute in enumerate(attributes):
            if attribute.get('attribute') == attribute_name:
                attr_pos = idx
                break

        if attr_pos is None:
            raise ValueError('Cannot find attribute [{}] in [{}]'.format(attribute_name,
                                                                         attri_name))

        instances = attri_data.get('instances')

        plot_mds_matrix['attribute_value_size'] = None
        for instance_name, attri_values in instances.items():
            plot_mds_matrix.loc[plot_mds_matrix.instance == instance_name,
                                ['attribute_value_size']] = attri_values[attr_pos]

        return plot_mds_matrix

    def _build_color_mds_matrix(self, plot_mds_matrix, obj_data, dimension, attribute_name):
        """
        _build_color_mds_matrix: append attribute value to rotation_matrix
        """
        logging.info('appending attribute value for grouping color to rotation matrix')

        plot_mds_matrix = plot_mds_matrix.copy()

        if dimension == 'row':
            attribute_mapping = obj_data.get('row_mapping')
            attribute_mapping_ref = obj_data.get('row_attributemapping_ref')
        elif dimension == 'col':
            attribute_mapping = obj_data.get('col_mapping')
            attribute_mapping_ref = obj_data.get('col_attributemapping_ref')
        else:
            raise ValueError('Unexpected dimension')

        if not attribute_mapping:
            logging.warning('Matrix object does not have {}_mapping attribute'.format(dimension))
            # build matrix with unify color and shape
            return plot_mds_matrix
        else:
            # append instance col mapping from row/col_mapping
            plot_mds_matrix['instance'] = plot_mds_matrix.index.map(attribute_mapping)

        res = self.dfu.get_objects({'object_refs': [attribute_mapping_ref]})['data'][0]
        attri_data = res['data']
        attri_name = res['info'][1]

        attributes = attri_data.get('attributes')

        attr_pos = None
        for idx, attribute in enumerate(attributes):
            if attribute.get('attribute') == attribute_name:
                attr_pos = idx
                break

        if attr_pos is None:
            raise ValueError('Cannot find attribute [{}] in [{}]'.format(attribute_name,
                                                                         attri_name))

        instances = attri_data.get('instances')

        plot_mds_matrix['attribute_value_color'] = None
        for instance_name, attri_values in instances.items():
            plot_mds_matrix.loc[plot_mds_matrix.instance == instance_name,
                                ['attribute_value_color']] = attri_values[attr_pos]

        return plot_mds_matrix

    def _build_2_comp_trace(self, plot_mds_matrix, components_x, components_y):

        traces = []

        if 'attribute_value_color' in plot_mds_matrix.columns and 'attribute_value_size' in plot_mds_matrix.columns:

            maximum_marker_size = 10
            sizeref = 2.*float(max(plot_mds_matrix['attribute_value_size']))/(maximum_marker_size**2)

            for name in set(plot_mds_matrix.attribute_value_color):
                attribute_value_size = plot_mds_matrix.loc[plot_mds_matrix['attribute_value_color'].eq(name)].attribute_value_size
                size_list = list(map(abs, list(map(float, attribute_value_size))))
                for idx, val in enumerate(size_list):
                    if val == 0:
                        size_list[idx] = sys.float_info.min
                trace = go.Scatter(
                    x=list(plot_mds_matrix.loc[plot_mds_matrix['attribute_value_color'].eq(name)][components_x]),
                    y=list(plot_mds_matrix.loc[plot_mds_matrix['attribute_value_color'].eq(name)][components_y]),
                    mode='markers',
                    name=name,
                    text=list(plot_mds_matrix.loc[plot_mds_matrix['attribute_value_color'].eq(name)].index),
                    textposition='bottom center',
                    marker=go.Marker(symbol='circle', sizemode='area', sizeref=sizeref,
                                     size=size_list, sizemin=2,
                                     line=go.Line(color='rgba(217, 217, 217, 0.14)', width=0.5),
                                     opacity=0.8))
                traces.append(trace)
        elif 'attribute_value_color' in plot_mds_matrix.columns:
            for name in set(plot_mds_matrix.attribute_value_color):
                trace = go.Scatter(
                    x=list(plot_mds_matrix.loc[plot_mds_matrix['attribute_value_color'].eq(name)][components_x]),
                    y=list(plot_mds_matrix.loc[plot_mds_matrix['attribute_value_color'].eq(name)][components_y]),
                    mode='markers',
                    name=name,
                    text=list(plot_mds_matrix.loc[plot_mds_matrix['attribute_value_color'].eq(name)].index),
                    textposition='bottom center',
                    marker=go.Marker(size=10, opacity=0.8, line=go.Line(color='rgba(217, 217, 217, 0.14)',
                                                                        width=0.5)))
                traces.append(trace)
        elif 'attribute_value_size' in plot_mds_matrix.columns:

            maximum_marker_size = 10
            sizeref = 2.*float(max(plot_mds_matrix['attribute_value_size']))/(maximum_marker_size**2)

            for name in set(plot_mds_matrix.instance):
                attribute_value_size = plot_mds_matrix.loc[plot_mds_matrix['instance'].eq(name)].attribute_value_size
                size_list = list(map(abs, list(map(float, attribute_value_size))))
                for idx, val in enumerate(size_list):
                    if val == 0:
                        size_list[idx] = sys.float_info.min
                trace = go.Scatter(
                    x=list(plot_mds_matrix.loc[plot_mds_matrix['instance'].eq(name)][components_x]),
                    y=list(plot_mds_matrix.loc[plot_mds_matrix['instance'].eq(name)][components_y]),
                    mode='markers',
                    name=name,
                    text=list(plot_mds_matrix.loc[plot_mds_matrix['instance'].eq(name)].index),
                    textposition='bottom center',
                    marker=go.Marker(symbol='circle', sizemode='area', sizeref=sizeref,
                                     size=size_list, sizemin=2,
                                     line=go.Line(color='rgba(217, 217, 217, 0.14)', width=0.5),
                                     opacity=0.8))
                traces.append(trace)
        else:
            trace = go.Scatter(
                x=list(plot_mds_matrix[components_x]),
                y=list(plot_mds_matrix[components_y]),
                mode='markers',
                text=list(plot_mds_matrix.index),
                textposition='bottom center',
                marker=go.Marker(size=10, opacity=0.8,
                                 line=go.Line(color='rgba(217, 217, 217, 0.14)', width=0.5)))
            traces.append(trace)

        return traces

    def _plot_mds_matrix(self, plot_mds_matrix, n_components):

        output_directory = os.path.join(self.scratch, str(uuid.uuid4()))
        self._mkdir_p(output_directory)
        result_file_paths = []

        all_pairs = list(itertools.combinations(range(1, n_components+1), 2))

        for pair in all_pairs:
            first_component = pair[0]
            second_component = pair[1]
            result_file_path = os.path.join(output_directory, 'mds_plot_{}_{}.html'.format(
                                                                                first_component,
                                                                                second_component))

            traces = self._build_2_comp_trace(plot_mds_matrix,
                                              'principal_component_{}'.format(first_component),
                                              'principal_component_{}'.format(second_component))

            data = go.Data(traces)
            layout = go.Layout(xaxis=go.XAxis(title='PC{}'.format(first_component), showline=False),
                               yaxis=go.YAxis(title='PC{}'.format(second_component), showline=False))
            fig = go.Figure(data=data, layout=layout)

            plot(fig, filename=result_file_path)

            result_file_paths.append(result_file_path)

        return result_file_paths

    def _validate_mds_matrix(self, obj_data, dimension,
                             color_marker_by, scale_size_by):

        if dimension == 'row':
            attribute_mapping = obj_data.get('row_mapping')
        elif dimension == 'col':
            attribute_mapping = obj_data.get('col_mapping')
        else:
            raise ValueError('Unexpected dimension')

        if not attribute_mapping:
            if (color_marker_by and color_marker_by.get('attribute_color')[0]) or \
               (scale_size_by and scale_size_by.get('attribute_size')[0]):
                raise ValueError('Matrix object is not associated with any {} attribute mapping'.format(dimension))

    def __init__(self, config):

        self.ws_url = config["workspace-url"]
        self.callback_url = config['SDK_CALLBACK_URL']
        self.token = config['KB_AUTH_TOKEN']
        self.scratch = config['scratch']
        self.dfu = DataFileUtil(self.callback_url, service_ver='release')
        self.working_dir = self.scratch

        self.data_util = DataUtil(config)
        self.dfu = DataFileUtil(self.callback_url)
        self.output_dir = os.path.join(self.working_dir, self.VEGAN_OUT_DIR)
        self._mkdir_p(self.output_dir)


    def run_mds(self, params):
        """
        run_mds: perform MDS analysis on matrix
        :param input_obj_ref: object reference of a matrix
        :param workspace_name: the name of the workspace
        :param mds_matrix_name: name of MDS (KBaseExperiments.MDSMatrix) object
        :param n_components - dimentionality of the reduced space (default 2)
        :param max_iter: maximum iterations allowed
        :param metric: indication of running metric or non-metric MDS
        :param distance_metric: distance the ordination will be performed on, default to "bray"

        The metaMDS routine in R vegan library has the useful default behavior of following
        the ordination with a rotation via principal components analysis such that MDS
        axis 1 reflects the principal source of variation, and so on, as is characteristic of
        eigenvalue methods.
        Procrustes rotation--The minimum sum of squares is called the root mean square error (rmse).
        The smaller the rmse, the more similar the two configurations are.
        Two final configurations are considered to have converged (arrived at essentially the same
        solution) when the rmse is less than 0.01, and no single residual value exceeds 0.005.
        Procrustes analysis thereby provides a mechanism for determining when to stop repeatedly
        re-running the analysis - stop when there is convergence as measured by procrustes rmse.
        # for ecological data, samples should be standardized by sample size to avoid ordinations
        # that reflect primarily sample size unless the input matrix has already been standardized
        # The Wisconsin double standarization is done automatically following the square transformation. 
        # Generate a distance (dissimiliarity) matrix from the multivariate data
        # Bray-Curtis dissimilarity matrix
        # if distance_metric == 'bray_curtis', call metaMDS with the default Bray distance setting
        # else: # euclidean
        """

        logging.info('--->\nrunning mds with input\n' +
                     'params:\n{}'.format(json.dumps(params, indent=1)))

        self._validate_run_mds_params(params)

        input_obj_ref = params.get('input_obj_ref')
        workspace_name = params.get('workspace_name')
        mds_matrix_name = params.get('mds_matrix_name')
        n_components = int(params.get('n_components', 2))

        res = self.dfu.get_objects({'object_refs': [input_obj_ref]})['data'][0]
        obj_data = res['data']
        obj_name = res['info'][1]
        obj_type = res['info'][2]

        max_size = len(obj_data['data']['col_ids'])
        if n_components > max_size:
            raise ValueError('Number of components should be less than number of samples')

        if "KBaseMatrices" in obj_type:
            # create the input file from obj_data
            matrix_tab = obj_data['data']['values']
            row_ids = obj_data['data']['row_ids']
            col_ids = obj_data['data']['col_ids']
            matrix_df = pd.DataFrame(matrix_tab, index=row_ids, columns=col_ids)

            matrix_data_file = os.path.join(self.output_dir, obj_name + '.csv')
            with open(matrix_data_file, 'w') as m_file:
                matrix_df.to_csv(m_file, sep='\t')

            params['datafile'] = matrix_data_file
            self.run_mds_with_file(params)
        else:
            err_msg = 'Ooops! [{}] is not supported.\n'.format(obj_type)
            err_msg += 'Please provide a KBaseMatrices object'
            raise ValueError("err_msg")
   
        # saving the mds_matrix object
        """
        mds_ref = self._save_mds_matrix(workspace_name, input_obj_ref, mds_matrix_name,
                                        rotation_matrix_df, components_df, explained_variance,
                                        explained_variance_ratio, singular_values,
                                        n_components, dimension)
        """
        mds_ref = None
        returnVal = {'mds_ref': mds_ref}

        # generating report
        report_output = self._generate_mds_report(mds_ref, self.output_dir,
                                                  workspace_name, n_components)

        returnVal.update(report_output)
        return returnVal

    def run_mds_with_file(self, params):
        """
        run_mds_with_file: perform MDS analysis on matrix
        :param datafile: a file that contains the matrix data
        :param workspace_name: the name of the workspace
        :param mds_matrix_name: name of MDS (KBaseExperiments.MDSMatrix) object
        :param n_components - dimentionality of the reduced space (default 2)
        :param max_iter: maximum iterations allowed
        :param metric: indication of running metric or non-metric MDS
        :param distance_metric: distance the ordination will be performed on, default to "bray"
        """

        logging.info('--->\nrunning mds with input \n' +
                     'params:\n{}'.format(json.dumps(params, indent=1)))

        rscrpt_file = self._build_rMDS_script(params)
        logging.info('--->\nR script file has been written to {}'.format(rscrpt_file))

        exitCode = self._execute_r_script(rscrpt_file)

        returnVal = {'mds_ref': None,
		     'report_name': None,
                     'report_ref': None}

        return returnVal

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

