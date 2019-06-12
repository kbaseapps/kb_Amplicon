# -*- coding: utf-8 -*-
import os
import time
import unittest
from configparser import ConfigParser

from kb_Amplicon.kb_AmpliconImpl import kb_Amplicon
from kb_Amplicon.kb_AmpliconServer import MethodContext
from kb_Amplicon.authclient import KBaseAuth as _KBaseAuth
from kb_Amplicon.Utils.MDSUtils import MDSUtils

from installed_clients.WorkspaceClient import Workspace


class kb_AmpliconTest(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        token = os.environ.get('KB_AUTH_TOKEN', None)
        config_file = os.environ.get('KB_DEPLOYMENT_CONFIG', None)
        cls.cfg = {}
        config = ConfigParser()
        config.read(config_file)
        for nameval in config.items('kb_Amplicon'):
            cls.cfg[nameval[0]] = nameval[1]
        # Getting username from Auth profile for token
        authServiceUrl = cls.cfg['auth-service-url']
        auth_client = _KBaseAuth(authServiceUrl)
        user_id = auth_client.get_user(token)
        # WARNING: don't call any logging methods on the context object,
        # it'll result in a NoneType error
        cls.ctx = MethodContext(None)
        cls.ctx.update({'token': token,
                        'user_id': user_id,
                        'provenance': [
                            {'service': 'kb_Amplicon',
                             'method': 'please_never_use_it_in_production',
                             'method_params': []
                             }],
                        'authenticated': 1})
        cls.wsURL = cls.cfg['workspace-url']
        cls.wsClient = Workspace(cls.wsURL)
        cls.serviceImpl = kb_Amplicon(cls.cfg)
        cls.scratch = cls.cfg['scratch']
        cls.callback_url = os.environ['SDK_CALLBACK_URL']
        suffix = int(time.time() * 1000)
        cls.wsName = "test_kb_Amplicon_" + str(suffix)
        ret = cls.wsClient.create_workspace({'workspace': cls.wsName})  # noqa

    @classmethod
    def tearDownClass(cls):
        if hasattr(cls, 'wsName'):
            cls.wsClient.delete_workspace({'workspace': cls.wsName})
            print('Test workspace was deleted')


    # NOTE: According to Python unittest naming rules test method names should start from 'test'. # noqa
    # Uncomment to skip this test
    # @unittest.skip("skipped test_run_mds_with_objref")
    def test_run_mds_with_objref(self):
        # Prepare test objects in workspace if needed using
        # self.getWsClient().save_objects({'workspace': self.getWsName(),
        #                                  'objects': []})
        #
        # Run your method by
        # ret = self.getImpl().your_method(self.getContext(), parameters...)
        #
        # Check returned data with
        # self.assertEqual(ret[...], ...) or other unittest methods
        ret = self.serviceImpl.run_mds(self.ctx, {'workspace_name': self.wsName,
                                                  'input_obj_ref': '40925/Incubation-16S',
                                                  'n_components': 3,
                                                  'max_iter': 20,
                                                  'mds_matrix_name': 'output_mds_from_obj'})
        self.assertTrue(ret[0]['mds_ref'])
        self.assertTrue(ret[0]['report_name'])
        self.assertTrue(ret[0]['report_ref'])
        mds_dir = '/kb/module/work/tmp/mds_output'
        self.assertTrue(os.path.isfile(os.path.join(mds_dir, 'dist_matrix.csv')))
        self.assertTrue(os.path.isfile(os.path.join(mds_dir, 'others.json')))
        self.assertTrue(os.path.isfile(os.path.join(mds_dir, 'site_ordination.csv')))
        self.assertTrue(os.path.isfile(os.path.join(mds_dir, 'species_ordination.csv')))
        self.assertTrue(os.path.isfile(os.path.join(mds_dir, 'Incubation-16S.csv')))
        self.assertTrue(os.path.isfile(os.path.join(mds_dir, 'mds_script.R')))
        self.assertTrue(os.path.isfile(os.path.join(mds_dir, 'saving_mds_plot.bmp')))
        self.assertTrue(os.path.isfile(os.path.join(mds_dir, 'saving_mds_plot.pdf')))

    # Uncomment to skip this test
    # @unittest.skip("skipped test_MDSUtilsrun_mds_with_file")
    def test_MDSUtils_run_mds_with_file(self):
        # Prepare test objects in workspace if needed using
        # self.getWsClient().save_objects({'workspace': self.getWsName(),
        #                                  'objects': []})
        #
        # Run your method by
        # ret = self.getImpl().your_method(self.getContext(), parameters...)
        #
        # Check returned data with
        # self.assertEqual(ret[...], ...) or other unittest methods
        mds_util = MDSUtils(self.cfg)

        ret = mds_util.run_mds_with_file({'workspace_name': self.wsName,
                                          'datafile': 'smpl_16s.csv',
                                          'n_components': 3,
                                          'max_iter': 20,
                                          'mds_matrix_name': 'output_mds_from_file'})
        self.assertEqual(ret, 0)
        mds_dir = '/kb/module/work/tmp/mds_output'
        self.assertTrue(os.path.isfile(os.path.join(mds_dir, 'dist_matrix.csv')))
        self.assertTrue(os.path.isfile(os.path.join(mds_dir, 'others.json')))
        self.assertTrue(os.path.isfile(os.path.join(mds_dir, 'site_ordination.csv')))
        self.assertTrue(os.path.isfile(os.path.join(mds_dir, 'species_ordination.csv')))
        self.assertTrue(os.path.isfile(os.path.join(mds_dir, 'smpl_16s.csv')))
        self.assertTrue(os.path.isfile(os.path.join(mds_dir, 'mds_script.R')))
        self.assertTrue(os.path.isfile(os.path.join(mds_dir, 'saving_mds_plot.bmp')))
        self.assertTrue(os.path.isfile(os.path.join(mds_dir, 'saving_mds_plot.pdf')))

