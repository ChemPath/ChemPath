import unittest
import os
from chempath_api import ChemPathAPI
from database import engine
from sqlalchemy.orm import Session
from models import Compound, TherapeuticArea

class TestChemPathAPI(unittest.TestCase):
    def setUp(self):
        self.api = ChemPathAPI()
        self.populate_sample_data()

    def populate_sample_data(self):
        with Session(engine) as session:
            compound_data = {
                'name': 'Test Compound',
                'smiles': 'CC(=O)Nc1ccc(cc1)S(=O)(=O)N',
                'molecular_weight': 300.0,
                'logp': 2.5,
                'h_bond_donors': 2,
                'h_bond_acceptors': 5,
                'ph': 7.0,
                'temperature': 25.0
            }
            compound = Compound(**compound_data)
            session.add(compound)
            session.commit()

    def test_chempath_api(self):
        compound_data = {
            'name': 'New Compound',
            'smiles': 'CC(=O)Oc1ccccc1C(=O)O',
            'molecular_weight': 180.16,
            'logp': 1.19,
            'h_bond_donors': 1,
            'h_bond_acceptors': 4,
            'ph': 7.0,
            'temperature': 25.0
        }
        result = self.api.insert_compound(compound_data)
        self.assertIsNotNone(result)
        retrieved_compound = self.api.get_compound_by_smiles(compound_data['smiles'])
        self.assertIsNotNone(retrieved_compound)

    def test_predict_therapeutic_areas(self):
        smiles = 'CC(=O)Oc1ccccc1C(=O)O'
        result = self.api.predict_therapeutic_areas(smiles)
        self.assertIsInstance(result, list)

    def test_perform_retrosynthesis(self):
        smiles = 'CC(=O)Oc1ccccc1C(=O)O'
        result = self.api.perform_retrosynthesis(smiles)
        self.assertIn('steps', result)

    def test_optimize_structure(self):
        smiles = 'CC(=O)Oc1ccccc1C(=O)O'
        result = self.api.optimize_structure(smiles)
        self.assertIsInstance(result, list)

    def tearDown(self):
        with Session(engine) as session:
            session.query(Compound).delete()
            session.query(TherapeuticArea).delete()
            session.commit()

if __name__ == '__main__':
    unittest.main()
