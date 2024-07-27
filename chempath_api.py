from structural_optimization import generate_analogs, chemical_space_exploration
from rdkit import Chem
from chempath_core import optimize_structure
from chempath_core import create_connection, search_compounds, insert_compound, get_therapeutic_areas, predict_therapeutic_areas, optimize_structure
from chempath_core import train_ml_model
from chempath_ml_models import (train_retrosynthesis_model, train_reaction_prediction_model,
                                predict_retrosynthesis, predict_reaction, smiles_to_fingerprint)
import numpy as np

class ChemPathAPI:
    def table_exists(self, table_name):
        cursor = self.conn.cursor()
        cursor.execute(f"SELECT name FROM sqlite_master WHERE type='table' AND name='{table_name}'")
        return cursor.fetchone() is not None

    def __init__(self, db_path):
        self.conn = create_connection(db_path)
        print("About to create tables")
        create_tables(self.conn)
        print("Tables created successfully")
        create_indexes(self.conn)
        if not self.table_exists('plant_compounds'):
            print("plant_compounds table does not exist. Creating now...")
            create_tables(self.conn)
        self.model, self.scaler = train_ml_model(self.conn)



    def search(self, query=None, filters=None, page=1, per_page=10):
        return search_compounds(self.conn, query, filters, page, per_page)

    def add_compound(self, compound_data):
        from chempath_database import insert_compound
        return insert_compound(self.conn, compound_data)

    def get_therapeutic_areas(self):
        return get_therapeutic_areas(self.conn)

    def predict_therapeutic_areas(self, smiles):
        return predict_therapeutic_areas(smiles, self.model, self.scaler)

    def optimize_structure(self, smiles, optimization_type, params):
        return optimize_structure(smiles, optimization_type, params)

    def explore_chemical_space(self, smiles, num_iterations=10):
        return chemical_space_exploration(smiles, num_iterations)
    
    def retrosynthesis_informed_optimization(self, smiles):
        from chempath_core import get_compound_by_smiles
        compound = get_compound_by_smiles(self.conn, smiles)
    def retrosynthesis_informed_optimization(self, smiles):
        from chempath_core import get_compound_by_smiles, get_retrosynthesis_data, retrosynthesis_informed_optimization, store_retrosynthesis_informed_optimization
        compound = get_compound_by_smiles(self.conn, smiles)
        if compound:
            retrosynthesis_data = get_retrosynthesis_data(self.conn, compound[0])
            optimized_smiles = retrosynthesis_informed_optimization(smiles, retrosynthesis_data)
            store_retrosynthesis_informed_optimization(self.conn, compound[0], optimized_smiles)
            return optimized_smiles
        else:
            return None
        
    def train_ml_model(self, conn):
        from chempath_core import train_ml_model
        return train_ml_model(conn)
    
    def train_ml_models(self):
        # Fetch data from the database
        cursor = self.conn.cursor()
        cursor.execute("SELECT smiles, retrosynthesis_feasibility, reaction_class FROM compounds")
        data = cursor.fetchall()

        X = np.array([smiles_to_fingerprint(row[0]) for row in data])
        y_retro = np.array([row[1] for row in data])
        y_reaction = np.array([row[2] for row in data])

        self.retro_model, self.retro_scaler = train_retrosynthesis_model(X, y_retro)
        self.reaction_model, self.reaction_scaler = train_reaction_prediction_model(X, y_reaction)

    def predict_retrosynthesis_feasibility(self, smiles):
        return predict_retrosynthesis(self.retro_model, self.retro_scaler, smiles)

    def predict_reaction_class(self, smiles):
        return predict_reaction(self.reaction_model, self.reaction_scaler, smiles)




    
    def search_compounds(self, query):
        from chempath_core import search_compounds
        return search_compounds(self.conn, query)

    
    def __init__(self, db_path):
        self.conn = create_connection(db_path)
        self.model, self.scaler = train_ml_model(self.conn)

    def close_connection(self):
        self.conn.close()

