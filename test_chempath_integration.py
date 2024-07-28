import unittest
from unittest.mock import patch
from chempath_api import ChemPathAPI

class TestChemPathIntegration(unittest.TestCase):

    def setUp(self):
        self.api = ChemPathAPI('test_chempath.db')

    def tearDown(self):
        self.api.close_connection()

    @patch('requests.get')
    def test_external_api_call(self, mock_get):
        mock_get.return_value.json.return_value = {'compound': 'Quercetin', 'additional_info': 'Test data'}
        result = self.api.get_external_compound_info('Quercetin')
        self.assertEqual(result['compound'], 'Quercetin')

    @patch('requests.get')
    def test_external_api_call_error(self, mock_get):
        mock_get.return_value.json.side_effect = ValueError('Mocked error')
        with self.assertRaises(ValueError):
            self.api.get_external_compound_info('Quercetin')

    def test_data_exchange(self):
        external_data = {
            'name': 'New Compound',
            'smiles': 'C1=CC=CC=C1',
            'molecular_weight': 78.11,
            'logp': 2.0,
            'plant_source': 'Test Plant',
            'biological_activities': 'Test Activity',
            'traditional_use': 'Test Use',
            'h_bond_donors': 1,
            'h_bond_acceptors': 1,
            'polar_surface_area': 0.0,
            'rotatable_bonds': 0
        }
        compound_id = self.api.insert_compound(external_data)
        self.assertIsNotNone(compound_id)
        retrieved_data = self.api.get_compound(compound_id)
        self.assertEqual(retrieved_data['name'], 'New Compound')

    def test_data_exchange_invalid_data(self):
        external_data = {'name': 'Invalid Compound', 'smiles': 'Invalid SMILES'}
        with self.assertRaises(ValueError):
            self.api.insert_compound(external_data)

    @unittest.skipIf(not hasattr(ChemPathAPI, 'generate_report'), "Report generation not implemented")
    def test_report_generation(self):
        report = self.api.generate_report('Quercetin')
        self.assertIn('Quercetin', report)
        self.assertIn('Molecular Weight', report)

    @unittest.skipIf(not hasattr(ChemPathAPI, 'generate_report'), "Report generation not implemented")
    def test_report_generation_invalid_compound(self):
        report = self.api.generate_report('Invalid Compound')
        self.assertNotIn('Molecular Weight', report)

if __name__ == '__main__':
    unittest.main()
