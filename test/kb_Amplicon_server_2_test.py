# -*- coding: utf-8 -*-
import os
import time
import unittest
import inspect
from mock import patch
import requests
from configparser import ConfigParser

from kb_Amplicon.kb_AmpliconImpl import kb_Amplicon
from kb_Amplicon.kb_AmpliconServer import MethodContext
from installed_clients.authclient import KBaseAuth as _KBaseAuth
from kb_Amplicon.Utils.MDSUtils import MDSUtils
from installed_clients.DataFileUtilClient import DataFileUtil
from installed_clients.WorkspaceClient import Workspace
from installed_clients.AbstractHandleClient import AbstractHandle as HandleService


class kb_AmpliconTest(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.token = os.environ.get('KB_AUTH_TOKEN', None)
        config_file = os.environ.get('KB_DEPLOYMENT_CONFIG', None)
        cls.cfg = {}
        config = ConfigParser()
        config.read(config_file)
        for nameval in config.items('kb_Amplicon'):
            cls.cfg[nameval[0]] = nameval[1]
        # Getting username from Auth profile for token
        auth_service_url = cls.cfg['auth-service-url']
        auth_client = _KBaseAuth(auth_service_url)
        user_id = auth_client.get_user(cls.token)
        # WARNING: don't call any logging methods on the context object,
        # it'll result in a NoneType error
        cls.ctx = MethodContext(None)
        cls.ctx.update({'token': cls.token,
                        'user_id': user_id,
                        'provenance': [
                            {'service': 'kb_Amplicon',
                             'method': 'please_never_use_it_in_production',
                             'method_params': []
                             }],
                        'authenticated': 1})
        cls.wsURL = cls.cfg['workspace-url']
        cls.wsClient = Workspace(cls.wsURL)
        cls.shockURL = cls.cfg['shock-url']
        cls.serviceImpl = kb_Amplicon(cls.cfg)
        cls.scratch = cls.cfg['scratch']
        cls.callback_url = os.environ['SDK_CALLBACK_URL']

        cls.dfu = DataFileUtil(cls.callback_url)
        cls.mds_util = MDSUtils(cls.cfg)
        cls.hs = HandleService(url=cls.cfg['handle-service-url'],
                               token=cls.token)

        suffix = int(time.time() * 1000)
        cls.wsName = "test_kb_Amplicon_" + str(suffix)
        ret = cls.wsClient.create_workspace({'workspace': cls.wsName})
        cls.wsId = ret[0]

        small_file = os.path.join(cls.scratch, 'test.txt')
        with open(small_file, "w") as f:
            f.write("empty content")
        cls.test_shock = cls.dfu.file_to_shock({'file_path': small_file, 'make_handle': True})
        cls.handles_to_delete = []
        cls.nodes_to_delete = []
        cls.handles_to_delete.append(cls.test_shock['handle']['hid'])
        cls.nodes_to_delete.append(cls.test_shock['shock_id'])

    @classmethod
    def tearDownClass(cls):
        if hasattr(cls, 'wsName'):
            cls.wsClient.delete_workspace({'workspace': cls.wsName})
            print('Test workspace was deleted')
        if hasattr(cls, 'nodes_to_delete'):
            for node in cls.nodes_to_delete:
                cls.delete_shock_node(node)
        if hasattr(cls, 'handles_to_delete'):
            cls.hs.delete_handles(cls.hs.hids_to_handles(cls.handles_to_delete))
            print('Deleted handles ' + str(cls.handles_to_delete))

    @classmethod
    def delete_shock_node(cls, node_id):
        header = {'Authorization': 'Oauth {0}'.format(cls.token)}
        requests.delete(cls.shockURL + '/node/' + node_id, headers=header,
                        allow_redirects=True)
        print('Deleted shock node ' + node_id)

    def getMDSUtil(self):
        return self.__class__.mds_util

    def mock_file_to_shock(params):
        print('Mocking DataFileUtilClient.file_to_shock')
        print(params)

        return kb_AmpliconTest().test_shock

    def loadExpressionMatrix(self):
        if hasattr(self.__class__, 'expr_matrix_ref'):
            return self.__class__.expr_matrix_ref

        # matrix_file_name = 'test_import.xlsx'
        col_attribute = {'attributes': [{'attribute': 'test_attribute_1',
                                         'attribute_ont_id': 'OBI_0500020',
                                         'source': 'upload',
                                         'unit': 'Hour',
                                         'unit_ont_id': 'UO_0000032'},
                                        {'attribute': 'test_attribute_2',
                                         'attribute_ont_id': 'CHEBI:9168',
                                         'source': 'upload',
                                         'unit': 'nanogram per milliliter',
                                         'unit_ont_id': 'UO_0000275'},
                                        {'attribute': 'test_attribute_3',
                                         'attribute_ont_id': 'CHEBI:9168',
                                         'source': 'upload',
                                         'unit': 'nanogram per milliliter',
                                         'unit_ont_id': 'UO_0000275'}],
                         'instances': {'instance_1': ['1', '5', '9'],
                                       'instance_2': ['2', '6', '10'],
                                       'instance_3': ['3', '7', '11'],
                                       'instance_4': ['4', '8', '12']},
                         'ontology_mapping_method': 'User Curation'}

        info = self.dfu.save_objects({
            'id': self.wsId,
            'objects': [{'type': 'KBaseExperiments.AttributeMapping',
                         'data': col_attribute,
                         'name': 'test_ExpressionMatrix_col_attribute_mapping'}]})[0]

        col_attributemapping_ref = "%s/%s/%s" % (info[6], info[0], info[4])

        self.__class__.col_attributemapping_ref = col_attributemapping_ref
        print('Loaded Col AttributeMapping: ' + col_attributemapping_ref)

        row_attribute = {'attributes': [{'attribute': 'test_attribute_1',
                                         'attribute_ont_id': 'OBI_0500020',
                                         'source': 'upload',
                                         'unit': 'Hour',
                                         'unit_ont_id': 'UO_0000032'},
                                        {'attribute': 'test_attribute_2',
                                         'attribute_ont_id': 'CHEBI:9168',
                                         'source': 'upload',
                                         'unit': 'nanogram per milliliter',
                                         'unit_ont_id': 'UO_0000275'},
                                        {'attribute': 'test_attribute_3',
                                         'attribute_ont_id': 'CHEBI:9168',
                                         'source': 'upload',
                                         'unit': 'nanogram per milliliter',
                                         'unit_ont_id': 'UO_0000275'}],
                         'instances': {'WRI_RS00010_CDS_1': ['1', '4', '7'],
                                       'WRI_RS00015_CDS_1': ['3', '4', '8'],
                                       'WRI_RS00025_CDS_1': ['3', '6', '7'],
                                       'WRI_RS00030_CDS_1': ['3', '6', '7'],
                                       'WRI_RS00035_CDS_1': ['3', '6', '7']},
                         'ontology_mapping_method': 'User Curation'}

        info = self.dfu.save_objects({
            'id': self.wsId,
            'objects': [{'type': 'KBaseExperiments.AttributeMapping',
                         'data': row_attribute,
                         'name': 'test_ExpressionMatrix_row_attribute_mapping'}]})[0]

        row_attributemapping_ref = "%s/%s/%s" % (info[6], info[0], info[4])

        self.__class__.row_attributemapping_ref = row_attributemapping_ref
        print('Loaded Row AttributeMapping: ' + row_attributemapping_ref)

        matrix_data = {'attributes': {'Instrument': 'Old Faithful',
                                      'Scientist': 'Marie Currie'},
                       'col_attributemapping_ref': col_attributemapping_ref,
                       'col_mapping': {'instance_1': 'instance_1',
                                       'instance_2': 'instance_2',
                                       'instance_3': 'instance_3',
                                       'instance_4': 'instance_4'},
                       'col_normalization': 'test_col_normalization',
                       'data': {'col_ids': ['instance_1', 'instance_2', 'instance_3',
                                            'instance_4'],
                                'row_ids': ['WRI_RS00010_CDS_1', 'WRI_RS00015_CDS_1',
                                            'WRI_RS00025_CDS_1', 'WRI_RS00030_CDS_1',
                                            'WRI_RS00035_CDS_1'],
                                'values': [[1, 2, 3, 4],
                                           [50, 60, 70, 80],
                                           [9, 10, 11, 12],
                                           [9, 10, 11, 12],
                                           [9, 10, 11, 12]]},
                       'description': 'test_desc',
                       'row_attributemapping_ref': row_attributemapping_ref,
                       'row_mapping': {'WRI_RS00010_CDS_1': 'WRI_RS00010_CDS_1',
                                       'WRI_RS00015_CDS_1': 'WRI_RS00015_CDS_1',
                                       'WRI_RS00025_CDS_1': 'WRI_RS00025_CDS_1',
                                       'WRI_RS00030_CDS_1': 'WRI_RS00030_CDS_1',
                                       'WRI_RS00035_CDS_1': 'WRI_RS00035_CDS_1'},
                       'row_normalization': 'test_row_normalization',
                       'scale': 'log2',
                       'search_attributes': ['Scientist | Marie Currie',
                                             'Instrument | Old Faithful']}

        info = self.dfu.save_objects({'id': self.wsId,
                                      'objects': [{'type': 'KBaseMatrices.ExpressionMatrix',
                                                   'data': matrix_data,
                                                   'name': 'test_ExpressionMatrix'}]})[0]

        expr_matrix_ref = "%s/%s/%s" % (info[6], info[0], info[4])

        self.__class__.expr_matrix_ref = expr_matrix_ref
        print('Loaded ExpressionMatrix: ' + expr_matrix_ref)

        # load associated matrix
        matrix_data = {'attributes': {'Instrument': 'Old Faithful',
                                      'Scientist': 'Marie Currie'},
                       'col_attributemapping_ref': col_attributemapping_ref,
                       'col_mapping': {'instance_1': 'instance_1',
                                       'instance_2': 'instance_2',
                                       'instance_3': 'instance_3',
                                       'instance_4': 'instance_4'},
                       'col_normalization': 'test_col_normalization',
                       'data': {'col_ids': ['instance_1', 'instance_2', 'instance_3',
                                            'instance_4'],
                                'row_ids': ['WRI_RS00010_CDS_1', 'WRI_RS00015_CDS_1',
                                            'WRI_RS00025_CDS_1', 'WRI_RS00030_CDS_1',
                                            'WRI_RS00035_CDS_1'],
                                'values': [[0.1, 0.2, 0.3, 0.4],
                                           [0.5, 0.6, 0.7, 0.8],
                                           [0.9, 1, 1.1, 1.2],
                                           [0.9, 1, 1.1, 1.2],
                                           [0.9, 1, 1.1, 1.2]]},
                       'description': 'test_desc',
                       'row_attributemapping_ref': row_attributemapping_ref,
                       'row_mapping': {'WRI_RS00010_CDS_1': 'WRI_RS00010_CDS_1',
                                       'WRI_RS00015_CDS_1': 'WRI_RS00015_CDS_1',
                                       'WRI_RS00025_CDS_1': 'WRI_RS00025_CDS_1',
                                       'WRI_RS00030_CDS_1': 'WRI_RS00030_CDS_1',
                                       'WRI_RS00035_CDS_1': 'WRI_RS00035_CDS_1'},
                       'row_normalization': 'test_row_normalization',
                       'scale': 'log2',
                       'search_attributes': ['Scientist | Marie Currie',
                                             'Instrument | Old Faithful']}

        info = self.dfu.save_objects({
            'id': self.wsId,
            'objects': [{'type': 'KBaseMatrices.ExpressionMatrix',
                         'data': matrix_data,
                         'name': 'test_associated_ExpressionMatrix'}]})[0]

        asso_matrix_ref = "%s/%s/%s" % (info[6], info[0], info[4])

        self.__class__.asso_matrix_ref = asso_matrix_ref
        print('Loaded Associated ExpressionMatrix: ' + asso_matrix_ref)

    def start_test(self):
        testname = inspect.stack()[1][3]
        print('\n*** starting test: ' + testname + ' **')

    def test_init_ok(self):
        self.start_test()
        class_attri = ['ws_url', 'callback_url', 'token', 'scratch', 'dfu', 'working_dir',
                       'output_dir']

        mds_util = self.getMDSUtil()
        self.assertTrue(set(class_attri) <= set(mds_util.__dict__.keys()))
        self.assertEqual(mds_util.scratch, self.cfg.get('scratch'))

    @patch.object(DataFileUtil, "file_to_shock", side_effect=mock_file_to_shock)
    def test_run_metaMDS_scale_by_attri_plot_associated_matrix_without_color(self, file_to_shock):
        self.start_test()
        self.loadExpressionMatrix()

        # testing col dimension with linked matrix
        params = {'workspace_name': self.wsName,
                  'input_obj_ref': self.expr_matrix_ref,
                  'n_components': 3,
                  'max_iter': 20,
                  'plot_script': 'plot(my_data.mds,type="t",display="sites")',
                  'plot_type': 'ps',
                  'plot_name': '',
                  'attribute_mapping_obj_ref': self.col_attributemapping_ref,
                  'associated_matrix_obj_ref': self.asso_matrix_ref,
                  'scale_size_by': {'attribute_size': ["test_attribute_1"]},
                  # 'color_marker_by': {'attribute_color': ['test_attribute_2']},
                  'mds_matrix_name': 'output_mds_from_obj',
                  'dimension': 'col'}

        ret = self.serviceImpl.run_metaMDS(self.ctx, params)[0]

        self.assertTrue('report_name' in ret)
        self.assertTrue('report_ref' in ret)
        self.assertTrue('mds_ref' in ret)

        pca_matrix_ref = ret.get('mds_ref')

        pca_data = self.dfu.get_objects({"object_refs": [pca_matrix_ref]})['data'][0]['data']

        expected_values = ['distance_matrix', 'mds_parameters', 'original_matrix_ref',
                           'rotation_matrix', 'site_ordination', 'species_ordination']
        self.assertTrue(set(expected_values) <= set(pca_data.keys()))

        expected_row_ids = ['WRI_RS00010_CDS_1', 'WRI_RS00015_CDS_1', 'WRI_RS00025_CDS_1',
                            'WRI_RS00030_CDS_1', 'WRI_RS00035_CDS_1']
        expected_col_ids = ['instance_1', 'instance_2', 'instance_3', 'instance_4']

        result_row_ids = [value[0] for value in pca_data.get('species_ordination').get('values')]
        result_col_ids = [value[0] for value in pca_data.get('site_ordination').get('values')]
        self.assertCountEqual(result_row_ids, expected_row_ids)
        self.assertCountEqual(result_col_ids, expected_col_ids)

        mds_dir = '/kb/module/work/tmp/mds_output'
        expected_files = ['dist_matrix.csv', 'mds_script.R', 'others.json',
                          'plotly_fig.html', 'site_ordination.csv', 'species_ordination.csv',
                          'test_ExpressionMatrix.csv',
                          'usr_plt_name.ps']
        self.assertTrue(set(expected_files) <= set(os.listdir(mds_dir)))

    @patch.object(DataFileUtil, "file_to_shock", side_effect=mock_file_to_shock)
    def test_run_metaMDS_scale_by_attri_plot_associated_matrix(self, file_to_shock):
        self.start_test()
        self.loadExpressionMatrix()

        # testing col dimension with linked matrix
        params = {'workspace_name': self.wsName,
                  'input_obj_ref': self.expr_matrix_ref,
                  'n_components': 3,
                  'max_iter': 20,
                  'plot_script': 'plot(my_data.mds,type="t",display="sites")',
                  'plot_type': 'ps',
                  'plot_name': '',
                  'attribute_mapping_obj_ref': self.col_attributemapping_ref,
                  'associated_matrix_obj_ref': self.asso_matrix_ref,
                  'scale_size_by': {'attribute_size': ["test_attribute_1"]},
                  'color_marker_by': {'attribute_color': ['test_attribute_2']},
                  'mds_matrix_name': 'output_mds_from_obj',
                  'dimension': 'col'}

        ret = self.serviceImpl.run_metaMDS(self.ctx, params)[0]

        self.assertTrue('report_name' in ret)
        self.assertTrue('report_ref' in ret)
        self.assertTrue('mds_ref' in ret)

        pca_matrix_ref = ret.get('mds_ref')

        pca_data = self.dfu.get_objects({"object_refs": [pca_matrix_ref]})['data'][0]['data']

        expected_values = ['distance_matrix', 'mds_parameters', 'original_matrix_ref',
                           'rotation_matrix', 'site_ordination', 'species_ordination']
        self.assertTrue(set(expected_values) <= set(pca_data.keys()))

        expected_row_ids = ['WRI_RS00010_CDS_1', 'WRI_RS00015_CDS_1', 'WRI_RS00025_CDS_1',
                            'WRI_RS00030_CDS_1', 'WRI_RS00035_CDS_1']
        expected_col_ids = ['instance_1', 'instance_2', 'instance_3', 'instance_4']

        result_row_ids = [value[0] for value in pca_data.get('species_ordination').get('values')]
        result_col_ids = [value[0] for value in pca_data.get('site_ordination').get('values')]
        self.assertCountEqual(result_row_ids, expected_row_ids)
        self.assertCountEqual(result_col_ids, expected_col_ids)

        mds_dir = '/kb/module/work/tmp/mds_output'
        expected_files = ['dist_matrix.csv', 'mds_script.R', 'others.json',
                          'plotly_fig.html', 'site_ordination.csv', 'species_ordination.csv',
                          'test_ExpressionMatrix.csv',
                          'usr_plt_name.ps']
        self.assertTrue(set(expected_files) <= set(os.listdir(mds_dir)))

    @patch.object(DataFileUtil, "file_to_shock", side_effect=mock_file_to_shock)
    def test_run_metaMDS_with_linked_matrix_ok(self, file_to_shock):
        self.start_test()
        self.loadExpressionMatrix()

        # testing col dimension with linked matrix
        params = {'workspace_name': self.wsName,
                  'input_obj_ref': self.expr_matrix_ref,
                  'n_components': 3,
                  'max_iter': 20,
                  'plot_script': 'plot(my_data.mds,type="t",display="sites")',
                  'plot_type': 'ps',
                  'plot_name': '',
                  'attribute_mapping_obj_ref': self.col_attributemapping_ref,
                  'associated_matrix_obj_ref': self.asso_matrix_ref,
                  'scale_size_by': {'row_size': ['WRI_RS00010_CDS_1']},
                  'color_marker_by': {'attribute_color': ['test_attribute_2']},
                  'mds_matrix_name': 'output_mds_from_obj',
                  'dimension': 'col'}

        ret = self.serviceImpl.run_metaMDS(self.ctx, params)[0]

        self.assertTrue('report_name' in ret)
        self.assertTrue('report_ref' in ret)
        self.assertTrue('mds_ref' in ret)

        pca_matrix_ref = ret.get('mds_ref')

        pca_data = self.dfu.get_objects({"object_refs": [pca_matrix_ref]})['data'][0]['data']

        expected_values = ['distance_matrix', 'mds_parameters', 'original_matrix_ref',
                           'rotation_matrix', 'site_ordination', 'species_ordination']
        self.assertTrue(set(expected_values) <= set(pca_data.keys()))

        expected_row_ids = ['WRI_RS00010_CDS_1', 'WRI_RS00015_CDS_1', 'WRI_RS00025_CDS_1',
                            'WRI_RS00030_CDS_1', 'WRI_RS00035_CDS_1']
        expected_col_ids = ['instance_1', 'instance_2', 'instance_3', 'instance_4']

        result_row_ids = [value[0] for value in pca_data.get('species_ordination').get('values')]
        result_col_ids = [value[0] for value in pca_data.get('site_ordination').get('values')]
        self.assertCountEqual(result_row_ids, expected_row_ids)
        self.assertCountEqual(result_col_ids, expected_col_ids)

        mds_dir = '/kb/module/work/tmp/mds_output'
        expected_files = ['dist_matrix.csv', 'mds_script.R', 'others.json',
                          'plotly_fig.html', 'site_ordination.csv', 'species_ordination.csv',
                          'test_ExpressionMatrix.csv',
                          'usr_plt_name.ps']
        self.assertTrue(set(expected_files) <= set(os.listdir(mds_dir)))

    @patch.object(DataFileUtil, "file_to_shock", side_effect=mock_file_to_shock)
    def test_run_metaMDS_with_row_linked_matrix_ok(self, file_to_shock):
        self.start_test()
        self.loadExpressionMatrix()

        # testing row dimension with linked matrix
        params = {'workspace_name': self.wsName,
                  'input_obj_ref': self.expr_matrix_ref,
                  'n_components': 3,
                  'max_iter': 20,
                  'plot_script': 'plot(my_data.mds,type="t",display="sites")',
                  'plot_type': 'ps',
                  'plot_name': '',
                  'attribute_mapping_obj_ref': self.row_attributemapping_ref,
                  'associated_matrix_obj_ref': self.asso_matrix_ref,
                  'scale_size_by': {'col_size': ['instance_2']},
                  'color_marker_by': {'attribute_color': ['test_attribute_2']},
                  'mds_matrix_name': 'output_mds_from_obj',
                  'dimension': 'row'}

        ret = self.serviceImpl.run_metaMDS(self.ctx, params)[0]

        self.assertTrue('report_name' in ret)
        self.assertTrue('report_ref' in ret)
        self.assertTrue('mds_ref' in ret)

        pca_matrix_ref = ret.get('mds_ref')

        pca_data = self.dfu.get_objects({"object_refs": [pca_matrix_ref]})['data'][0]['data']

        expected_values = ['distance_matrix', 'mds_parameters', 'original_matrix_ref',
                           'rotation_matrix', 'site_ordination', 'species_ordination']
        self.assertTrue(set(expected_values) <= set(pca_data.keys()))

        expected_row_ids = ['instance_1', 'instance_2', 'instance_3', 'instance_4']
        expected_col_ids = ['WRI_RS00010_CDS_1', 'WRI_RS00015_CDS_1', 'WRI_RS00025_CDS_1',
                            'WRI_RS00030_CDS_1', 'WRI_RS00035_CDS_1']

        result_row_ids = [value[0] for value in pca_data.get('species_ordination').get('values')]
        result_col_ids = [value[0] for value in pca_data.get('site_ordination').get('values')]
        self.assertCountEqual(result_row_ids, expected_row_ids)
        self.assertCountEqual(result_col_ids, expected_col_ids)

        mds_dir = '/kb/module/work/tmp/mds_output'
        expected_files = ['dist_matrix.csv', 'mds_script.R', 'others.json',
                          'plotly_fig.html', 'site_ordination.csv', 'species_ordination.csv',
                          'test_ExpressionMatrix.csv',
                          'usr_plt_name.ps']
        self.assertTrue(set(expected_files) <= set(os.listdir(mds_dir)))

    @patch.object(DataFileUtil, "file_to_shock", side_effect=mock_file_to_shock)
    def test_run_metaMDS_with_linked_matrix_ok_only_scale_size(self, file_to_shock):
        self.start_test()
        self.loadExpressionMatrix()

        # testing only scale_size_by with linked matrix
        params = {'workspace_name': self.wsName,
                  'input_obj_ref': self.expr_matrix_ref,
                  'n_components': 3,
                  'max_iter': 20,
                  'plot_script': 'plot(my_data.mds,type="t",display="sites")',
                  'plot_type': 'ps',
                  'plot_name': '',
                  'associated_matrix_obj_ref': self.asso_matrix_ref,
                  'scale_size_by': {'row_size': ['WRI_RS00010_CDS_1']},
                  'mds_matrix_name': 'output_mds_from_obj',
                  'dimension': 'col'}

        ret = self.serviceImpl.run_metaMDS(self.ctx, params)[0]

        self.assertTrue('report_name' in ret)
        self.assertTrue('report_ref' in ret)
        self.assertTrue('mds_ref' in ret)

        pca_matrix_ref = ret.get('mds_ref')

        pca_data = self.dfu.get_objects({"object_refs": [pca_matrix_ref]})['data'][0]['data']

        expected_values = ['distance_matrix', 'mds_parameters', 'original_matrix_ref',
                           'rotation_matrix', 'site_ordination', 'species_ordination']
        self.assertTrue(set(expected_values) <= set(pca_data.keys()))

        expected_row_ids = ['WRI_RS00010_CDS_1', 'WRI_RS00015_CDS_1', 'WRI_RS00025_CDS_1',
                            'WRI_RS00030_CDS_1', 'WRI_RS00035_CDS_1']
        expected_col_ids = ['instance_1', 'instance_2', 'instance_3', 'instance_4']

        result_row_ids = [value[0] for value in pca_data.get('species_ordination').get('values')]
        result_col_ids = [value[0] for value in pca_data.get('site_ordination').get('values')]
        self.assertCountEqual(result_row_ids, expected_row_ids)
        self.assertCountEqual(result_col_ids, expected_col_ids)

        mds_dir = '/kb/module/work/tmp/mds_output'
        expected_files = ['dist_matrix.csv', 'mds_script.R', 'others.json',
                          'plotly_fig.html', 'site_ordination.csv', 'species_ordination.csv',
                          'test_ExpressionMatrix.csv',
                          'usr_plt_name.ps']
        self.assertTrue(set(expected_files) <= set(os.listdir(mds_dir)))

    @patch.object(DataFileUtil, "file_to_shock", side_effect=mock_file_to_shock)
    def test_run_metaMDS_ok_col_dimension(self, file_to_shock):
        self.start_test()
        self.loadExpressionMatrix()

        # testing col dimension
        params = {'workspace_name': self.wsName,
                  'input_obj_ref': self.expr_matrix_ref,
                  'n_components': 3,
                  'max_iter': 20,
                  'plot_script': 'plot(my_data.mds,type="t",display="sites")',
                  'plot_type': 'ps',
                  'plot_name': '',
                  'attribute_mapping_obj_ref': self.col_attributemapping_ref,
                  'scale_size_by': {'attribute_size': ["test_attribute_1"]},
                  'color_marker_by': {'attribute_color': ['test_attribute_2']},
                  'mds_matrix_name': 'output_mds_from_obj',
                  'dimension': 'col'}

        ret = self.serviceImpl.run_metaMDS(self.ctx, params)[0]

        self.assertTrue('report_name' in ret)
        self.assertTrue('report_ref' in ret)
        self.assertTrue('mds_ref' in ret)

        pca_matrix_ref = ret.get('mds_ref')

        pca_data = self.dfu.get_objects({"object_refs": [pca_matrix_ref]})['data'][0]['data']

        expected_values = ['distance_matrix', 'mds_parameters', 'original_matrix_ref',
                           'rotation_matrix', 'site_ordination', 'species_ordination']
        self.assertTrue(set(expected_values) <= set(pca_data.keys()))

        expected_row_ids = ['WRI_RS00010_CDS_1', 'WRI_RS00015_CDS_1', 'WRI_RS00025_CDS_1',
                            'WRI_RS00030_CDS_1', 'WRI_RS00035_CDS_1']
        expected_col_ids = ['instance_1', 'instance_2', 'instance_3', 'instance_4']

        result_row_ids = [value[0] for value in pca_data.get('species_ordination').get('values')]
        result_col_ids = [value[0] for value in pca_data.get('site_ordination').get('values')]
        self.assertCountEqual(result_row_ids, expected_row_ids)
        self.assertCountEqual(result_col_ids, expected_col_ids)

        mds_dir = '/kb/module/work/tmp/mds_output'
        expected_files = ['dist_matrix.csv', 'mds_script.R', 'others.json',
                          'plotly_fig.html', 'site_ordination.csv', 'species_ordination.csv',
                          'test_ExpressionMatrix.csv',
                          'usr_plt_name.ps']
        self.assertTrue(set(expected_files) <= set(os.listdir(mds_dir)))

    @patch.object(DataFileUtil, "file_to_shock", side_effect=mock_file_to_shock)
    def test_run_metaMDS_ok_row_dimension(self, file_to_shock):
        self.start_test()
        self.loadExpressionMatrix()

        # testing row dimension
        params = {'workspace_name': self.wsName,
                  'input_obj_ref': self.expr_matrix_ref,
                  'n_components': 3,
                  'max_iter': 20,
                  'plot_script': 'plot(my_data.mds,type="t",display="sites")',
                  'plot_type': 'ps',
                  'plot_name': '',
                  'attribute_mapping_obj_ref': self.row_attributemapping_ref,
                  'scale_size_by': {'attribute_size': ["test_attribute_1"]},
                  'color_marker_by': {'attribute_color': ['test_attribute_2']},
                  'mds_matrix_name': 'output_mds_from_obj',
                  'dimension': 'row'}

        ret = self.serviceImpl.run_metaMDS(self.ctx, params)[0]

        self.assertTrue('report_name' in ret)
        self.assertTrue('report_ref' in ret)
        self.assertTrue('mds_ref' in ret)

        pca_matrix_ref = ret.get('mds_ref')

        pca_data = self.dfu.get_objects({"object_refs": [pca_matrix_ref]})['data'][0]['data']

        expected_values = ['distance_matrix', 'mds_parameters', 'original_matrix_ref',
                           'rotation_matrix', 'site_ordination', 'species_ordination']
        self.assertTrue(set(expected_values) <= set(pca_data.keys()))

        expected_row_ids = ['instance_1', 'instance_2', 'instance_3', 'instance_4']
        expected_col_ids = ['WRI_RS00010_CDS_1', 'WRI_RS00015_CDS_1', 'WRI_RS00025_CDS_1',
                            'WRI_RS00030_CDS_1', 'WRI_RS00035_CDS_1']

        result_row_ids = [value[0] for value in pca_data.get('species_ordination').get('values')]
        result_col_ids = [value[0] for value in pca_data.get('site_ordination').get('values')]
        self.assertCountEqual(result_row_ids, expected_row_ids)
        self.assertCountEqual(result_col_ids, expected_col_ids)

        mds_dir = '/kb/module/work/tmp/mds_output'
        expected_files = ['dist_matrix.csv', 'mds_script.R', 'others.json',
                          'plotly_fig.html', 'site_ordination.csv', 'species_ordination.csv',
                          'test_ExpressionMatrix.csv',
                          'usr_plt_name.ps']
        self.assertTrue(set(expected_files) <= set(os.listdir(mds_dir)))

    @patch.object(DataFileUtil, "file_to_shock", side_effect=mock_file_to_shock)
    def test_run_metaMDS_with_linked_matrix_ok_only_color_by(self, file_to_shock):
        self.start_test()
        self.loadExpressionMatrix()

        # testing only color_marker_by
        params = {'workspace_name': self.wsName,
                  'input_obj_ref': self.expr_matrix_ref,
                  'n_components': 3,
                  'max_iter': 20,
                  'plot_script': 'plot(my_data.mds,type="t",display="sites")',
                  'plot_type': 'ps',
                  'plot_name': '',
                  'attribute_mapping_obj_ref': self.col_attributemapping_ref,
                  'scale_size_by': None,
                  'color_marker_by': {'attribute_color': ['test_attribute_2']},
                  'mds_matrix_name': 'output_mds_from_obj',
                  'dimension': 'col'}

        ret = self.serviceImpl.run_metaMDS(self.ctx, params)[0]

        self.assertTrue('report_name' in ret)
        self.assertTrue('report_ref' in ret)
        self.assertTrue('mds_ref' in ret)

        pca_matrix_ref = ret.get('mds_ref')

        pca_data = self.dfu.get_objects({"object_refs": [pca_matrix_ref]})['data'][0]['data']

        expected_values = ['distance_matrix', 'mds_parameters', 'original_matrix_ref',
                           'rotation_matrix', 'site_ordination', 'species_ordination']
        self.assertTrue(set(expected_values) <= set(pca_data.keys()))

        expected_row_ids = ['WRI_RS00010_CDS_1', 'WRI_RS00015_CDS_1', 'WRI_RS00025_CDS_1',
                            'WRI_RS00030_CDS_1', 'WRI_RS00035_CDS_1']
        expected_col_ids = ['instance_1', 'instance_2', 'instance_3', 'instance_4']

        result_row_ids = [value[0] for value in pca_data.get('species_ordination').get('values')]
        result_col_ids = [value[0] for value in pca_data.get('site_ordination').get('values')]
        self.assertCountEqual(result_row_ids, expected_row_ids)
        self.assertCountEqual(result_col_ids, expected_col_ids)

        mds_dir = '/kb/module/work/tmp/mds_output'
        expected_files = ['dist_matrix.csv', 'mds_script.R', 'others.json',
                          'plotly_fig.html', 'site_ordination.csv', 'species_ordination.csv',
                          'test_ExpressionMatrix.csv',
                          'usr_plt_name.ps']
        self.assertTrue(set(expected_files) <= set(os.listdir(mds_dir)))

    @patch.object(DataFileUtil, "file_to_shock", side_effect=mock_file_to_shock)
    def test_run_metaMDS_ok_without_grouping(self, file_to_shock):
        self.start_test()
        self.loadExpressionMatrix()

        # testing without grouping
        params = {'workspace_name': self.wsName,
                  'input_obj_ref': self.expr_matrix_ref,
                  'n_components': 3,
                  'max_iter': 20,
                  'plot_script': 'plot(my_data.mds,type="t",display="sites")',
                  'plot_type': 'ps',
                  'plot_name': '',
                  'mds_matrix_name': 'output_mds_from_obj'}

        ret = self.serviceImpl.run_metaMDS(self.ctx, params)[0]

        self.assertTrue('report_name' in ret)
        self.assertTrue('report_ref' in ret)
        self.assertTrue('mds_ref' in ret)

        pca_matrix_ref = ret.get('mds_ref')

        pca_data = self.dfu.get_objects({"object_refs": [pca_matrix_ref]})['data'][0]['data']

        expected_values = ['distance_matrix', 'mds_parameters', 'original_matrix_ref',
                           'rotation_matrix', 'site_ordination', 'species_ordination']
        self.assertTrue(set(expected_values) <= set(pca_data.keys()))

        expected_row_ids = ['WRI_RS00010_CDS_1', 'WRI_RS00015_CDS_1', 'WRI_RS00025_CDS_1',
                            'WRI_RS00030_CDS_1', 'WRI_RS00035_CDS_1']
        expected_col_ids = ['instance_1', 'instance_2', 'instance_3', 'instance_4']

        result_row_ids = [value[0] for value in pca_data.get('species_ordination').get('values')]
        result_col_ids = [value[0] for value in pca_data.get('site_ordination').get('values')]
        self.assertCountEqual(result_row_ids, expected_row_ids)
        self.assertCountEqual(result_col_ids, expected_col_ids)

        mds_dir = '/kb/module/work/tmp/mds_output'
        expected_files = ['dist_matrix.csv', 'mds_script.R', 'others.json',
                          'plotly_fig.html', 'site_ordination.csv', 'species_ordination.csv',
                          'test_ExpressionMatrix.csv',
                          'usr_plt_name.ps']
        self.assertTrue(set(expected_files) <= set(os.listdir(mds_dir)))
