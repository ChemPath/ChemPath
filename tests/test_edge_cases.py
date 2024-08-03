import unittest
from chempath_api import ChemPathAPI
from database import engine
from sqlalchemy.orm import Session
from models import Compound, TherapeuticArea

class TestChemPathEdgeCases(unittest.TestCase):
    def setUp(self):
        self.api = ChemPathAPI()
        self.all_areas = [
            "Area1", "Area2", "Area3", "Antimicrobial", "Antifungal", 
            "Enzyme Inhibitor", "Protein Stabilizer", "Cryoprotectant", 
            "Thermostable Enzyme"
        ]

    def tearDown(self):
        with Session(engine) as session:
            session.query(Compound).delete()
            session.query(TherapeuticArea).delete()
            session.commit()

    def test_01_large_molecule(self):
        fullerene_data = {
            "name": "Fullerene C60",
            "smiles": "C1=C2C=C3C=C4C=C5C=C6C=C7C=C8C=C9C=C%10C=C%11C=C%12C=C%13C=C%14C=C%15C=C1C2=C3C4=C5C6=C7C8=C9C%10=C%11C%12=C%13C%14=C%15",
            "molecular_weight": 720.66,
            "logp": 6.8,
            "h_bond_donors": 0,
            "h_bond_acceptors": 0,
            "ph": 7.0,
            "temperature": 25.0
        }
        compound_id = self.api.insert_compound(fullerene_data)
        self.assertIsNotNone(compound_id)

    def test_02_small_molecule(self):
        hydrogen_data = {
            "name": "Hydrogen",
            "smiles": "[H][H]",
            "molecular_weight": 2.02,
            "logp": 0.0,
            "h_bond_donors": 0,
            "h_bond_acceptors": 0,
            "ph": 7.0,
            "temperature": 25.0
        }
        compound_id = self.api.insert_compound(hydrogen_data)
        self.assertIsNotNone(compound_id)

    def test_03_invalid_structure(self):
        invalid_compound_data = {
            "name": "Invalid Compound",
            "smiles": "C1=CC=C",  # This is an invalid SMILES string
            "molecular_weight": 100.0,
            "logp": 1.5,
            "h_bond_donors": 0,
            "h_bond_acceptors": 0,
            "ph": 7.0,
            "temperature": 25.0
        }
        with self.assertRaises(ValueError):
            self.api.insert_compound(invalid_compound_data)

    def test_04_strained_molecule(self):
        cubane_data = {
            "name": "Cubane",
            "smiles": "C12C3C4C1C5C2C3C45",
            "molecular_weight": 104.15,
            "logp": 1.6,
            "h_bond_donors": 0,
            "h_bond_acceptors": 0,
            "ph": 7.0,
            "temperature": 25.0
        }
        compound_id = self.api.insert_compound(cubane_data)
        self.assertIsNotNone(compound_id)

    def test_05_extreme_ph(self):
        predictions_low_ph = self.api.predict_therapeutic_areas("C1=CC=CC=C1", ph=-1)
        self.assertIsInstance(predictions_low_ph, list)
        self.assertTrue(len(predictions_low_ph) > 0)

        predictions_high_ph = self.api.predict_therapeutic_areas("C1=CC=CC=C1", ph=14)
        self.assertIsInstance(predictions_high_ph, list)
        self.assertTrue(len(predictions_high_ph) > 0)

    def test_06_extreme_temperature(self):
        predictions_low_temp = self.api.predict_therapeutic_areas("C1=CC=CC=C1", temperature=-273.16)
        self.assertIsInstance(predictions_low_temp, list)
        self.assertTrue(len(predictions_low_temp) > 0)

        predictions_high_temp = self.api.predict_therapeutic_areas("C1=CC=CC=C1", temperature=1000)
        self.assertIsInstance(predictions_high_temp, list)
        self.assertTrue(len(predictions_high_temp) > 0)

    def test_07_extreme_molecular_weight(self):
        tiny_compound_data = {
            "name": "Tiny Compound",
            "smiles": "[H][H]",
            "molecular_weight": 0.01,
            "logp": 0.0,
            "h_bond_donors": 0,
            "h_bond_acceptors": 0,
            "ph": 7.0,
            "temperature": 25.0
        }
        tiny_id = self.api.insert_compound(tiny_compound_data)
        self.assertIsNotNone(tiny_id)

        huge_compound_data = {
            "name": "Huge Compound",
            "smiles": "CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC",
            "molecular_weight": 12008.0,
            "logp": 1000.0,
            "h_bond_donors": 0,
            "h_bond_acceptors": 0,
            "ph": 7.0,
            "temperature": 25.0
        }
        huge_id = self.api.insert_compound(huge_compound_data)
        self.assertIsNotNone(huge_id)

if __name__ == '__main__':
    unittest.main()
