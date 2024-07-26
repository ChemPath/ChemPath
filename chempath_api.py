from structural_optimization import generate_analogs, chemical_space_exploration
from rdkit import Chem
from chempath_core import optimize_structure
from chempath_core import create_connection, search_compounds, insert_compound, get_therapeutic_areas, predict_therapeutic_areas, optimize_structure

class ChemPathAPI:
    def __init__(self, db_path):
        self.conn = create_connection(db_path)

    def search(self, query=None, filters=None, page=1, per_page=10):
        return search_compounds(self.conn, query, filters, page, per_page)

    def add_compound(self, compound_data):
        return insert_compound(self.conn, compound_data)

    def get_therapeutic_areas(self):
        return get_therapeutic_areas(self.conn)

    def predict_therapeutic_areas(self, smiles):
        all_therapeutic_areas = self.get_therapeutic_areas()
        return predict_therapeutic_areas(smiles, all_therapeutic_areas)

    def optimize_structure(self, smiles, optimization_type, params):
        return optimize_structure(smiles, optimization_type, params)

    def explore_chemical_space(self, smiles, num_iterations=10):
        return chemical_space_exploration(smiles, num_iterations)
    
    def retrosynthesis_informed_optimization(self, smiles):
        compound = get_compound_by_smiles(self.conn, smiles)
    def retrosynthesis_informed_optimization(self, smiles):
        compound = get_compound_by_smiles(self.conn, smiles)
        if compound:
            retrosynthesis_data = get_retrosynthesis_data(self.conn, compound[0])
            optimized_smiles = retrosynthesis_informed_optimization(smiles, retrosynthesis_data)
            store_retrosynthesis_informed_optimization(self.conn, compound[0], optimized_smiles)
            return optimized_smiles
        else:
            return None

    def close_connection(self):
        self.conn.close()

