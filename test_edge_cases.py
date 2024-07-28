
import unittest
from chempath_api import ChemPathAPI

class TestChemPathEdgeCases(unittest.TestCase):
    def setUp(self):
        self.api = ChemPathAPI(db_path=':memory:')
        self.all_areas = [
            "Area1", "Area2", "Area3", "Antimicrobial", "Antifungal", 
            "Enzyme Inhibitor", "Protein Stabilizer", "Cryoprotectant", 
            "Thermostable Enzyme"
        ]

    def tearDown(self):
        self.api.close_connection()

    def create_full_compound_data(self, name, smiles, molecular_weight, logp, 
                                 plant_source, biological_activities, traditional_use):
        return {
            'name': name,
            'smiles': smiles,
            'molecular_weight': molecular_weight,
            'logp': logp,
            'plant_source': plant_source,
            'biological_activities': biological_activities,
            'traditional_use': traditional_use,
            'h_bond_donors': 0,
            'h_bond_acceptors': 0,
            'polar_surface_area': 0.0,
            'rotatable_bonds': 0
        }

    def test_01_large_molecule(self):
        fullerene_data = self.create_full_compound_data(
            "Fullerene C60",
            "C1=C2C=C3C=C4C=C5C=C6C=C7C=C8C=C9C=C%10C=C%11C=C%12C=C%13C=C%14C=C%15C=C1C2=C3C4=C5C6=C7C8=C9C%10=C%11C%12=C%13C%14=C%15",
            molecular_weight=720.66,
            logp=6.8,
            plant_source="Synthetic",
            biological_activities="Antioxidant",
            traditional_use="Nanotechnology"
        )
        compound_id = self.api.insert_compound(fullerene_data)
        self.assertIsNotNone(compound_id)

    def test_02_small_molecule(self):
        hydrogen_data = self.create_full_compound_data(
            "Hydrogen",
            "[H][H]",
            molecular_weight=2.02,
            logp=0.0,
            plant_source="Various sources",
            biological_activities="Reducing agent",
            traditional_use="Fuel cells"
        )
        compound_id = self.api.insert_compound(hydrogen_data)
        self.assertIsNotNone(compound_id)

    def test_03_invalid_structure(self):
        invalid_compound_data = self.create_full_compound_data(
            "Invalid Compound",
            "C1=CC=C",  # This is an invalid SMILES string
            molecular_weight=100.0,
            logp=1.5,
            plant_source="Unknown",
            biological_activities="Unknown",
            traditional_use="Unknown"
        )
        with self.assertRaises(ValueError):
            self.api.insert_compound(invalid_compound_data)

    def test_04_strained_molecule(self):
        cubane_data = self.create_full_compound_data(
            "Cubane",
            "C12C3C4C1C5C2C3C45",
            molecular_weight=104.15,
            logp=1.6,
            plant_source="Synthetic",
            biological_activities="Potential fuel",
            traditional_use="Research"
        )
        compound_id = self.api.insert_compound(cubane_data)
        self.assertIsNotNone(compound_id)

    def test_05_extreme_ph(self):
        predictions_low_ph = self.api.predict_therapeutic_areas(
            "C1=CC=CC=C1", 
            all_therapeutic_areas=self.all_areas, 
            ph=-1
        )
        self.assertIsInstance(predictions_low_ph, list)
        self.assertTrue(len(predictions_low_ph) > 0)
        for area in predictions_low_ph:
            self.assertIn(area, self.all_areas)

        predictions_high_ph = self.api.predict_therapeutic_areas(
            "C1=CC=CC=C1", 
            all_therapeutic_areas=self.all_areas, 
            ph=14
        )
        self.assertIsInstance(predictions_high_ph, list)
        self.assertTrue(len(predictions_high_ph) > 0)
        for area in predictions_high_ph:
            self.assertIn(area, self.all_areas)

    def test_06_extreme_temperature(self):
        predictions_low_temp = self.api.predict_therapeutic_areas(
            "C1=CC=CC=C1", 
            all_therapeutic_areas=self.all_areas, 
            temperature=-273.16
        )
        self.assertIsInstance(predictions_low_temp, list)
        self.assertTrue(len(predictions_low_temp) > 0)
        for area in predictions_low_temp:
            self.assertIn(area, self.all_areas)

        predictions_high_temp = self.api.predict_therapeutic_areas
        predictions_high_temp = self.api.predict_therapeutic_areas("C1=CC=CC=C1", all_therapeutic_areas=self.all_areas, temperature=1000)
        self.assertIsInstance(predictions_high_temp, list)
        self.assertTrue(len(predictions_high_temp) > 0)
        for area in predictions_high_temp:
            self.assertIn(area, self.all_areas)

    def test_07_extreme_molecular_weight(self):
        tiny_compound_data = self.create_full_compound_data(
        "Tiny Compound",
        "[H][H]",  # SMILES string for a single hydrogen atom
        0.01,  # extremely low molecular weight
        1.5,
        "Unknown",
        "Unknown",
        "Unknown"
    )
        tiny_id = self.api.insert_compound(tiny_compound_data)
        self.assertIsNotNone(tiny_id)
        huge_compound_data = self.create_full_compound_data(
        "Huge Compound",
        "CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC",  # a valid SMILES string for a large compound
        12008.0, 1000.0, "Test", "Test huge MW", "MW testing"
    )
        huge_id = self.api.insert_compound(huge_compound_data)
        self.assertIsNotNone(huge_id)

if __name__ == '__main__':
    unittest.main()