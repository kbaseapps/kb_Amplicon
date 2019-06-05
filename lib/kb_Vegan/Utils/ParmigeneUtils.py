
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

from kb_Vegan.Utils.DataUtil import DataUtil
from installed_clients.DataFileUtilClient import DataFileUtil
from installed_clients.KBaseReportClient import KBaseReport


class ParmigeneUtils:

    R_BIN = '/kb/deployment/bin'
    PARMI_OUT_DIR = 'parmigene_output'
    PARAM_IN_WS = 'workspace_name'
    PARAM_IN_MATRIX = 'input_obj_ref'
    PARAM_OUT_MATRIX = 'parmigene_matrix_name'
    OMP_NUM_THREADS = 'num_threads'

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

    def _validate_run_mi_params(self, params):
        """
        _validate_run_mi_params:
            validates params passed to run_mi method
        """

        logging.info('start validating run_mi params')

        # check for required parameters
        for p in [self.PARAM_IN_MATRIX, self.PARAM_IN_WS, self.PARAM_OUT_MATRIX]:
            if p not in params:
                raise ValueError('"{}" parameter is required, but missing'.format(p))

    def _build_rParmigene_script(self, mi, algthm, num_threads, eps):
        """
        _build_rParmigene_script: build a sequence of R command calls according to params
        Note: To run the Parmigene functions, we will call different functions from the parmigene
        package, that requires a mutual information matrix (mi), the algorithm name (algthm),
        actual number of threads used (num_threads) and, sometimes, a positive numeric criteria
        (eps) to remove the weakest edge of each triple of nodes.
        """
        parmi_scrpt = 'library(parmigene)\n'
        parmi_scrpt += algthm + '(' + mi + ',' + eps + ')\n'
        # save the results in the memory
        # 1) store species ordination
        parmi_scrpt += 'variableScores <- vg_data.parmi$species\n'
        # 2) store site ordination
        parmi_scrpt += 'sampleScores <- vg_data.parmi$points\n'

        # save the results to the current dir
        # Write CSV in R
        parmi_scrpt += 'write.csv(dist_matrix,file="dist_matrix.csv",row.names=TRUE,na="")\n'
        parmi_scrpt += 'write.csv(variableScores,file="species_ordination.csv",' + \
                       'row.names=TRUE,na="")\n'

        # Write JSON in R
        parmi_scrpt += 'write_json(toJSON(dist_matrix),path="dist_matrix.json",pretty=TRUE,' + \
                       'auto_unbox=FALSE)\n'

        # save Parmigene plot
        parmi_scrpt += 'bmp(file="saving_mi_plot.bmp",width=6,height=4,units="in",res=100)\n'
        parmi_scrpt += 'plot(vg_data.parmi,type="n",display="sites")\n'
        parmi_scrpt += 'points(vg_data.parmi)\n'
        parmi_scrpt += 'dev.off()\n'

        parmi_rscript = 'parmi_script.R'
        rscrpt_file_path = os.path.join(self.output_dir, parmi_rscript)

        with open(rscrpt_file_path, 'w') as r_file:
            r_file.write(parmi_scrpt)
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

        logging.info('Running Parmigene script in current working directory: {}'.format(result_dir))
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

    def _mi_df_to_excel(self, mi_df, distance_df, result_dir, mi_matrix_ref):
        """
        write mutual information matrix df into excel
        """
        logging.info('writting mi data frame to excel file')
        mi_matrix_obj = self.dfu.get_objects({'object_refs': [mi_matrix_ref]})['data'][0]
        mi_matrix_info = mi_matrix_obj['info']
        mi_matrix_name = mi_matrix_info[1]

        file_path = os.path.join(result_dir, mi_matrix_name + ".xlsx")
        writer = pd.ExcelWriter(file_path)

        mi_df.to_excel(writer, "mi_matrix", index=True)
        if distance_df:
            distance_df.to_excel(writer, "mi_distance_matrix", index=True)

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

    def _mi_to_df(self, mi_matrix_ref):
        """
        retrieve mutual information matrix ws object to mi_df
        """
        logging.info('converting mutual information matrix to data frame')
        mi_data = self.dfu.get_objects({'object_refs': [mi_matrix_ref]})['data'][0]['data']

        rotation_matrix_data = mi_data.get('rotation_matrix')
        distance_matrix_data = mi_data.get('distance_matrix')
        original_matrix_ref = mi_data.get('original_matrix_ref')
        dimension = mi_data.get('mi_parameters').get('n_components')

        mi_df = self._Matrix2D_to_df(rotation_matrix_data)
        distance_df = None
        if distance_matrix_data:
            distance_df = self._Matrix2D_to_df(distance_matrix_data)

        if original_matrix_ref:
            logging.info('appending instance group information to mutual information data frame')
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

            mi_df = mi_df.merge(am_df, left_index=True, right_index=True, how='left',
                                validate='one_to_one')

        return mi_df, distance_df

    def _save_mi_matrix(self, workspace_name, input_obj_ref, mi_matrix_name,
                        distance_df, mi_params_df, site_ordin_df, species_ordin_df):

        logging.info('Saving MIMatrix...')

        if not isinstance(workspace_name, int):
            ws_name_id = self.dfu.ws_name_to_id(workspace_name)
        else:
            ws_name_id = workspace_name

        mi_data = {}

        mi_data.update({'distance_matrix': self._df_to_list(distance_df)})
        mi_data.update({'site_ordination': self._df_to_list(site_ordin_df)})
        mi_data.update({'species_ordination': self._df_to_list(species_ordin_df)})
        mi_data.update({'mi_parameters': self._df_to_list(mi_params_df)})
        mi_data.update({'original_matrix_ref': input_obj_ref})
        mi_data.update({'rotation_matrix': self._df_to_list(distance_df)})

        obj_type = 'KBaseExperiments.PCAMatrix'
        info = self.dfu.save_objects({
            "id": ws_name_id,
            "objects": [{
                "type": obj_type,
                "data": mi_data,
                "name": mi_matrix_name
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

        logging.info('Start packing result files from Parmigene...')

        output_files = list()

        output_dir = os.path.join(self.working_dir, str(uuid.uuid4()))
        self._mkdir_p(output_dir)
        mi_output = os.path.join(output_dir, 'metami_output.zip')
        self._zip_folder(out_dir, mi_output)

        output_files.append({'path': mi_output,
                             'name': os.path.basename(mi_output),
                             'label': os.path.basename(mi_output),
                             'description': 'Output file(s) generated by Parmigene'})
        return output_files

    def _generate_mi_html_report(self, mi_outdir, n_components):

        logging.info('Start generating html report for Parmigene results...')
        html_report = list()

        result_dir = os.path.join(self.working_dir, str(uuid.uuid4()))
        self._mkdir_p(result_dir)
        result_file_path = os.path.join(result_dir, 'mi_result.html')

        mi_plots = list()
        for root, folders, files in os.walk(mi_outdir):
            # Find the image files by their extensions.
            for f in files:
                if re.match('^[a-zA-Z]+.*.(jpeg|jpg|bmp|tiff|pdf|ps)$', f):
                    absolute_path = os.path.join(root, f)
                    logging.info("Adding {} to plot archive.".format(absolute_path))
                    mi_plots.append(absolute_path)

        visualization_content = ''

        for mi_plot in mi_plots:
            shutil.copy2(mi_plot,
                         os.path.join(result_dir, os.path.basename(mi_plot)))
            visualization_content += '<iframe height="900px" width="100%" '
            visualization_content += 'src="{}" '.format(os.path.basename(mi_plot))
            visualization_content += 'style="border:none;"></iframe>\n<p></p>\n'

        with open(result_file_path, 'w') as result_file:
            with open(os.path.join(os.path.dirname(__file__), 'templates', 'mi_template.html'),
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
                            'description': 'HTML summary report for Parmigene Matrix App'
                            })
        return html_report

    def _generate_mi_report(self, mi_ref, output_dir, workspace_name, n_components):
        logging.info('Creating Parmigene report...')

        output_files = self._generate_output_file_list(output_dir)
        output_html_files = self._generate_mi_html_report(output_dir, n_components)

        objects_created = list()
        objects_created.append({'ref': mi_ref,
                                'description': 'Mutual Information Matrix'})

        report_params = {'message': '',
                         'workspace_name': workspace_name,
                         'file_links': output_files,
                         'objects_created': objects_created,
                         'html_links': output_html_files,
                         'direct_html_link_index': 0,
                         'html_window_height': 666,
                         'report_object_name': 'kb_mi_report_' + str(uuid.uuid4())}

        kbase_report_client = KBaseReport(self.callback_url)
        output = kbase_report_client.create_extended_report(report_params)

        report_output = {'report_name': output['name'], 'report_ref': output['ref']}

        return report_output

    def __init__(self, config):

        self.ws_url = config["workspace-url"]
        self.callback_url = config['SDK_CALLBACK_URL']
        self.token = config['KB_AUTH_TOKEN']
        self.scratch = config['scratch']
        self.dfu = DataFileUtil(self.callback_url, service_ver='release')
        self.working_dir = self.scratch

        self.data_util = DataUtil(config)
        self.dfu = DataFileUtil(self.callback_url)
        self.output_dir = os.path.join(self.working_dir, self.PARMI_OUT_DIR)
        self._mkdir_p(self.output_dir)

    def run_mi(self, params):
        """
        run_mi: perform Parmigene analysis on matrix
        :param input_obj_ref: object reference of a matrix
        :param workspace_name: the name of the workspace
        :param mi_matrix_name: name of Parmigene (KBaseExperiments.MIMatrix) object
        :param n_components - dimentionality of the reduced space (default 2)
        :param max_iter: maximum iterations allowed
        :param distance_metric: distance the ordination will be performed on, default to "bray"
        """

        logging.info('--->\nrunning Parmigene with input\n' +
                     'params:\n{}'.format(json.dumps(params, indent=1)))

        self._validate_run_mi_params(params)

        input_obj_ref = params.get(self.PARAM_IN_MATRIX)
        workspace_name = params.get(self.PARAM_IN_WS)
        mi_matrix_name = params.get(self.PARAM_OUT_MATRIX)
        n_threads = int(params.get(self.OMP_NUM_THREADS, 2))

        res = self.dfu.get_objects({'object_refs': [input_obj_ref]})['data'][0]
        obj_data = res['data']
        obj_name = res['info'][1]
        obj_type = res['info'][2]

        exitCode = -99
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
            exitCode = self.run_mi_with_file(params)
        else:
            err_msg = 'Ooops! [{}] is not supported.\n'.format(obj_type)
            err_msg += 'Please provide a KBaseMatrices object'
            raise ValueError("err_msg")

        if exitCode == -99:
            raise ValueError('Caught subprocess.CalledProcessError while calling R.')

        # saving the mi_matrix object
        # read Parmigene results from files into data frames
        dist_matrix_df = pd.read_csv(os.path.join(self.output_dir, "dist_matrix.csv"))
        mi_params_df = pd.read_json(os.path.join(self.output_dir, "others.json"))
        site_ordin_df = pd.read_csv(os.path.join(self.output_dir, "site_ordination.csv"))
        species_ordin_df = pd.read_csv(os.path.join(self.output_dir, "species_ordination.csv"))

        mi_ref = self._save_mi_matrix(workspace_name, input_obj_ref, mi_matrix_name,
                                      dist_matrix_df, mi_params_df, site_ordin_df,
                                      species_ordin_df)
        returnVal = {'mi_ref': mi_ref}

        # generating report
        report_output = self._generate_mi_report(mi_ref, self.output_dir,
                                                 workspace_name, n_threads)

        returnVal.update(report_output)
        return returnVal

    def run_mi_with_file(self, params):
        """
        run_mi_with_file: perform Parmigene analysis on matrix
        :param datafile: a file that contains the matrix data
        :param workspace_name: the name of the workspace
        :param mi_matrix_name: name of Parmigene (KBaseExperiments.MIMatrix) object
        :param n_components - dimentionality of the reduced space (default 2)
        :param max_iter: maximum iterations allowed
        :param distance_metric: distance the ordination will be performed on, default to "bray"
        """

        logging.info('--->\nrunning Parmigene with input \n' +
                     'params:\n{}'.format(json.dumps(params, indent=1)))

        rscrpt_file = self._build_rmi_script(params)
        logging.info('--->\nR script file has been written to {}'.format(rscrpt_file))

        return self._execute_r_script(rscrpt_file)

    def export_mi_matrix_excel(self, params):
        """
        export MIMatrix as Excel
        """
        logging.info('start exporting Parmigene matrix')
        mi_matrix_ref = params.get('input_ref')

        mi_df, components_df = self._mi_to_df(mi_matrix_ref)

        result_dir = os.path.join(self.scratch, str(uuid.uuid4()))
        self._mkdir_p(result_dir)

        self._mi_df_to_excel(mi_df, components_df, result_dir, mi_matrix_ref)

        package_details = self.dfu.package_for_download({
            'file_path': result_dir,
            'ws_refs': [mi_matrix_ref]
        })

        return {'shock_id': package_details['shock_id']}
