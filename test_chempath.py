import unittest
import os
from chempath_api import ChemPathAPI
from chempath_database import populate_sample_data


class TestChemPathAPI(unittest.TestCase):
    def setUp(self):
        self.api = ChemPathAPI("test_database.db")
        populate_sample_data(self.api.conn)

    def test_chempath_api(self):
        self.api.create_tables()
        compound_data = {
            'name': 'Test Compound',
            'smiles': 'CC(=O)Nc1ccc(cc1)S(=O)(=O)N',
            'molecular_weight': 300.0,
            'logp': 2.5,
            'plant_source': 'Test Plant',
            'biological_activity': 'Test Activity',
            'traditional_use': 'Test Use',
            'h_bond_donors': 2,
            'h_bond_acceptors': 5,
            'polar_surface_area': 100.0,
            'rotatable_bonds': 5
        }
        self.api.insert_compound(compound_data)
        result = self.api.get_compound_by_smiles(compound_data['smiles'])
        self.assertIsNotNone(result)


class TestChemPathCore(unittest.TestCase):
    def setUp(self):
        self.api = ChemPathAPI("test_database.db")
        populate_sample_data(self.api.conn)

    def test_optimize_structure(self):
        compound = 'CC(=O)Nc1ccc(cc1)S(=O)(=O)N'
        goal = 'functional_group'
        params = {"target": "[OH]", "replacement": "[NH2]"}
        optimized_compounds = self.api.optimize_structure(compound, goal, params)
        self.assertIsInstance(optimized_compounds, str)

         # Test with invalid input compound
        invalid_compound = 'invalid_smiles_string'
        with self.assertRaises(ValueError):  # Change Exception to ValueError
            self.api.optimize_structure(invalid_compound, goal, params)



class TestChemPathDatabase(unittest.TestCase):
    def setUp(self):
        self.db_path = "test_database.db"
        if os.path.exists(self.db_path):
            try:
                os.remove(self.db_path)
            except PermissionError:
                pass  # File is in use, continue with existing database
        self.api = ChemPathAPI(self.db_path)
        populate_sample_data(self.api.conn)



    def test_create_table(self):
        self.api.create_tables()
        result = self.api.search_compounds("CC(=O)Nc1ccc(cc1)S(=O)(=O)N")
        self.assertIsNotNone(result)

    def test_insert_compound(self):
        compound_data = {
            'name': 'Test Compound',
            'smiles': 'CC(=O)Nc1ccc(cc1)S(=O)(=O)N',
            'molecular_weight': 300.0,
            'logp': 2.5,
            'plant_source': 'Test Plant',
            'biological_activity': 'Test Activity',
            'traditional_use': 'Test Use',
            'h_bond_donors': 2,
            'h_bond_acceptors': 5,
            'polar_surface_area': 100.0,
            'rotatable_bonds': 5
        }
        self.api.insert_compound(compound_data)
        result = self.api.get_compound_by_smiles(compound_data['smiles'])
        self.assertIsNotNone(result)

    def tearDown(self):
        self.api.conn.close()
        try:
            os.remove(self.api.db_name)
        except OSError:
            print(f"Unable to delete {self.api.db_name}. It may be in use by another process.")


if __name__ == '__main__':
    unittest.main()
