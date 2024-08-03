import unittest
from unittest.mock import patch
from sqlalchemy.orm import Session
from database import engine
from models import Compound, TherapeuticArea
from chempath_api import ChemPathAPI

class TestChemPathIntegration(unittest.TestCase):

    def setUp(self):
        self.api = ChemPathAPI()

    def tearDown(self):
        with Session(engine) as session:
            session.query(Compound).delete()
            session.query(TherapeuticArea).delete()
            session.commit()

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
            'h_bond_donors': 1,
            'h_bond_acceptors': 1,
            'ph': 7.0,
            'temperature': 25.0
        }
        compound_id = self.api.insert_compound(external_data)
        self.assertIsNotNone(compound_id)
        retrieved_data = self.api.get_compound(compound_id)
        self.assertEqual(retrieved_data['name'], 'New Compound')

    def test_data_exchange_invalid_data(self):
        external_data = {'name': 'Invalid Compound', 'smiles': 'Invalid SMILES'}
        with self.assertRaises(ValueError):
            self.api.insert_compound(external_data)

    def test_predict_therapeutic_areas(self):
        compound_data = {
            'smiles': 'CC(=O)Oc1ccccc1C(=O)O',
            'name': 'Aspirin'
        }
        result = self.api.predict_therapeutic_areas(compound_data)
        self.assertIn('predicted_therapeutic_areas', result)
        self.assertIsInstance(result['predicted_therapeutic_areas'], list)

    def test_retrosynthesis(self):
        compound_data = {
            'smiles': 'CC(=O)Oc1ccccc1C(=O)O',
            'name': 'Aspirin'
        }
        result = self.api.perform_retrosynthesis(compound_data)
        self.assertIn('steps', result)
        self.assertIsInstance(result['steps'], list)

    def test_advanced_retrosynthesis(self):
        compound_data = {
            'smiles': 'CC(=O)Oc1ccccc1C(=O)O',
            'name': 'Aspirin'
        }
        result = self.api.perform_advanced_retrosynthesis(compound_data)
        self.assertIn('tree_image', result)
        self.assertIsInstance(result['tree_image'], str)

    def test_optimize_compounds(self):
        compound_data = [
            {'smiles': 'CC(=O)Oc1ccccc1C(=O)O', 'name': 'Aspirin'},
            {'smiles': 'CC12CCC3C(C1CCC2O)CCC4=CC(=O)CCC34C', 'name': 'Testosterone'}
        ]
        result = self.api.optimize_compounds(compound_data)
        self.assertIsInstance(result, list)
        self.assertEqual(len(result), 2)
        self.assertIn('optimization_priority', result[0])

    def test_analyze_reagent(self):
        result = self.api.analyze_reagent('Sodium borohydride')
        self.assertIn('availability_score', result)
        self.assertIn('estimated_cost', result)
        self.assertIn('handling_requirements', result)
        self.assertIn('environmental_impact_score', result)

if __name__ == '__main__':
    unittest.main()
