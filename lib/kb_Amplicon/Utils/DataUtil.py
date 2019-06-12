import logging
import re
from collections import defaultdict

from dotmap import DotMap

from installed_clients.DataFileUtilClient import DataFileUtil
from installed_clients.GenericsServiceClient import GenericsService
from installed_clients.WorkspaceClient import Workspace as workspaceService

GENERICS_TYPE = ['FloatMatrix2D']  # add case in _convert_data for each additional type
GENERICS_MODULES = ['KBaseMatrices']


class DataUtil:

    @staticmethod
    def _find_between(s, start, end):
        """
        _find_between: find string in between start and end
        """

        return re.search('{}(.*){}'.format(start, end), s).group(1)

    def _find_constraints(self, obj_type):
        """
        _find_constraints: retrieve constraints (@contains, rowsum, unique, conditionally_required)
        """

        type_info = self.wsClient.get_type_info(obj_type)
        type_desc = type_info.get('description')
        constraints = {}

        for tag in ('contains', 'rowsum', 'unique', 'conditionally_required'):
            constraints[tag] = [line.strip().split()[1:] for line in type_desc.split("\n")
                                if line.startswith(f'@{tag}')]

        return constraints

    def _filter_constraints(self, constraints, data):
        """filters out constraints with missing keys"""
        contains_constraints = constraints.get('contains')

        filtered_constraints = []
        for contains_constraint in contains_constraints:
            in_values = contains_constraint[1:]
            missing_key = True
            for in_value in in_values:
                if in_value.startswith('values'):
                    search_value = re.search('{}(.*){}'.format('\(', '\)'), in_value).group(1)
                    unique_list = search_value.split('.')
                    key = unique_list[0]
                elif ':' in in_value:
                    key = in_value.split(':')[0]
                else:
                    unique_list = in_value.split('.')
                    key = unique_list[0]

                if key in data:
                    missing_key = False
                    break

            if missing_key:
                filtered_constraints.append(contains_constraint)

        for x in filtered_constraints:
            contains_constraints.remove(x)

        return constraints

    def _retrieve_value(self, data, value):
        """Parse the provided 'data' object to retrieve the item in 'value'."""
        logging.info('Getting value for {}'.format(value))
        retrieve_data = []
        m_data = DotMap(data)
        if value.startswith('set('):
            retrieve_data = value[4:-1].split(",")
        elif value.startswith('values('):  # TODO: nested values e.g. values(values(ids))
            search_value = re.search('{}(.*){}'.format('\(', '\)'), value).group(1)
            unique_list = search_value.split('.')
            m_data_cp = m_data.copy()
            for attr in unique_list:
                m_data_cp = getattr(m_data_cp, attr)
            retrieve_data = list(m_data_cp.values())
        elif ':' in value:
            obj_ref = getattr(m_data, value.split(':')[0])
            if obj_ref:
                included = value.split(':')[1]
                included = '/' + included.replace('.', '/')
                ref_data = self.wsClient.get_objects2({'objects': [{'ref': obj_ref,
                                                       'included': [included]}]})['data'][0]['data']
                m_ref_data = DotMap(ref_data)
                if ref_data:
                    if '*' not in included:
                        for key in included.split('/')[1:]:
                            m_ref_data = getattr(m_ref_data, key)
                    else:
                        keys = included.split('/')[1:]
                        m_ref_data = [x.get(keys[2]) for x in ref_data.get(keys[0])]  # TODO: only works for 2 level nested data like '/features/[*]/id'

                retrieve_data = list(m_ref_data)
        else:
            unique_list = value.split('.')
            m_data_cp = m_data.copy()
            for attr in unique_list:
                m_data_cp = getattr(m_data_cp, attr)
            retrieve_data = list(m_data_cp)

        logging.info('Retrieved value (first 20):\n{}\n'.format(retrieve_data[:20]))

        return retrieve_data

    def _validate(self, constraints, data):
        """
        _validate: validate data
        """

        validated = True
        failed_constraints = defaultdict(list)

        unique_constraints = constraints.get('unique')
        for unique_constraint in unique_constraints:
            retrieved_value = self._retrieve_value(data, unique_constraint[0])
            if len(set(retrieved_value)) != len(retrieved_value):
                validated = False
                failed_constraints['unique'].append(unique_constraint[0])

        contains_constraints = constraints.get('contains')
        for contains_constraint in contains_constraints:
            value = contains_constraint[0]
            in_values = contains_constraint[1:]
            retrieved_in_values = []
            for in_value in in_values:
                retrieved_in_values += self._retrieve_value(data, in_value)
            if not (set(self._retrieve_value(data, value)) <= set(retrieved_in_values)):
                validated = False
                failed_constraints['contains'].append(" ".join(contains_constraint))

        conditional_constraints = constraints.get('conditionally_required')
        for conditional_constraint in conditional_constraints:
            trigger = conditional_constraint[0]
            required_keys = conditional_constraint[1:]
            if trigger in data:
                missing_keys = [key for key in required_keys if key not in data]
                if missing_keys:
                    validated = False
                    failed_constraints['conditionally_required'].append(
                        (trigger, required_keys, missing_keys))

        return validated, failed_constraints

    @staticmethod
    def _raise_validation_error(params, validate):
        """Raise a meaningful error message for failed validation"""
        logging.error('Data failed type checking')
        failed_constraints = validate.get('failed_constraints')
        error_msg = ['Object {} failed type checking:'.format(params.get('obj_name'))]
        if failed_constraints.get('unique'):
            unique_values = failed_constraints.get('unique')
            error_msg.append('Object should have unique field: {}'.format(unique_values))
        if failed_constraints.get('contains'):
            contained_values = failed_constraints.get('contains')
            for contained_value in contained_values:
                subset_value = contained_value.split(' ')[0]
                super_value = ' '.join(contained_value.split(' ')[1:])
                if 'col_mapping' in super_value:
                    error_msg.append('Column attribute mapping instances should contain all '
                                     'column index from original data')

                if 'row_mapping' in super_value:
                    error_msg.append('Row attribute mapping instances should contain all row '
                                     'index from original data')

                error_msg.append('Object field [{}] should contain field [{}]'.format(
                    super_value,
                    subset_value))
        for failure in failed_constraints.get('conditionally_required', []):
            error_msg.append('If object field "{}" is present than object field(s) {} should '
                             'also be present. Object is missing {}'.format(*failure))
        raise ValueError('\n'.join(error_msg))

    def __init__(self, config):
        self.ws_url = config["workspace-url"]
        self.callback_url = config['SDK_CALLBACK_URL']
        self.token = config['KB_AUTH_TOKEN']
        self.scratch = config['scratch']
        self.serviceWizardURL = config['srv-wiz-url']
        self.wsClient = workspaceService(self.ws_url, token=self.token)
        self.dfu = DataFileUtil(self.callback_url)
        self.generics_service = GenericsService(self.serviceWizardURL)

    def list_generic_types(self, params=None):
        """
        *Not yet exposed in spec*
        list_generic_types: lists the current valid generics types

        arguments:
            none

        return:
            A list of generic types in the current environment
        """
        returnVal = [x['type_def'] for module in GENERICS_MODULES
                     for x in self.wsClient.get_all_type_info(module)]
        return returnVal

    def fetch_data(self, params):
        """
        fetch_data: fetch generics data as pandas dataframe for a generics data object

        arguments:
        obj_ref: generics object reference

        optional arguments:
        generics_module: the generics data module to be retrieved from
                        e.g. for an given data type like below:
                        typedef structure {
                          FloatMatrix2D data;
                          condition_set_ref condition_set_ref;
                        } SomeGenericsMatrix;
                        generics_module should be
                        {'data': 'FloatMatrix2D',
                         'condition_set_ref': 'condition_set_ref'}

        return:
        data_matrix: a pandas dataframe in json format
        """
        for p in ['obj_ref']:
            if p not in params:
                raise ValueError('"{}" parameter is required, but missing'.format(p))

        return self.generics_service.fetch_data(params)

    def validate_data(self, params):
        """
        validate_data: validate data

        arguments:
        obj_type: obj type e.g.: 'KBaseMatrices.ExpressionMatrix-1.1'
        data: obj data to be validated

        return:
        validated: True or False
        """

        constraints = self._find_constraints(params.get('obj_type'))
        data = params.get('data')

        constraints = self._filter_constraints(constraints, data)

        validated, failed_constraints = self._validate(constraints, data)

        return {'validated': validated,
                'failed_constraints': failed_constraints}

    def save_object(self, params):
        """
        save_object: validate data constraints and save matrix object

        arguments:
        obj_type: saving object data type
        obj_name: saving object name
        data: data to be saved
        workspace_name: workspace name matrix object to be saved to

        return:
        obj_ref: object reference
        """
        logging.info('Starting saving object')

        obj_type = params.get('obj_type')

        module_name = obj_type.split('.')[0]
        type_name = obj_type.split('.')[1]

        types = self.wsClient.get_module_info({'mod': module_name}).get('types')

        for module_type in types:
            if self._find_between(module_type, '\.', '\-') == type_name:
                obj_type = module_type
                break

        data = dict((k, v) for k, v in params.get('data').items() if v)
        validate = self.validate_data({'obj_type': obj_type,
                                       'data': data})

        if not validate.get('validated'):
            self._raise_validation_error(params, validate)

        workspace_name = params.get('workspace_name')
        if not isinstance(workspace_name, int):
            ws_name_id = self.dfu.ws_name_to_id(workspace_name)
        else:
            ws_name_id = workspace_name

        info = self.dfu.save_objects({
            "id": ws_name_id,
            "objects": [{
                "type": obj_type,
                "data": data,
                "name": params.get('obj_name')
            }]
        })[0]

        return {"obj_ref": "%s/%s/%s" % (info[6], info[0], info[4])}

