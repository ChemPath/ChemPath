from structural_optimization import generate_analogs, chemical_space_exploration
from rdkit import Chem
from chempath_core import (create_connection, search_compounds, insert_compound, 
                           get_therapeutic_areas, predict_therapeutic_areas, 
                           optimize_structure, perform_retrosynthesis, 
                           train_ml_model, get_compound_by_smiles, 
                           get_retrosynthesis_data, retrosynthesis_informed_optimization, 
                           store_retrosynthesis_informed_optimization, create_tables, create_indexes)

from chempath_ml_models import (train_retrosynthesis_model, train_reaction_prediction_model,
                                predict_retrosynthesis, predict_reaction, smiles_to_fingerprint)
import numpy as np
from ai_optimization import train_optimization_model, prioritize_optimization_strategies
import sqlite3

class ChemPathAPI:
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

    def table_exists(self, table_name):
        cursor = self.conn.cursor()
        cursor.execute(f"SELECT name FROM sqlite_master WHERE type='table' AND name='{table_name}'")
        return cursor.fetchone() is not None    

    def search_compounds(self, query):
        cursor = self.conn.cursor()
        cursor.execute("""
            SELECT name, smiles, molecular_weight, logp 
            FROM plant_compounds 
            WHERE name LIKE ? OR smiles LIKE ?
    """, ('%'+query+'%', '%'+query+'%'))
        results = cursor.fetchall()
        print(f"Search results: {results}")
        return results
    
    
    def get_all_compounds(self):
        cursor = self.conn.cursor()
        cursor.execute("SELECT name, molecular_weight, logp, h_bond_donors, h_bond_acceptors FROM plant_compounds")
        compounds = [{'name': row[0], 'molecular_weight': row[1], 'logp': row[2], 'h_bond_donors': row[3], 'h_bond_acceptors': row[4]} for row in cursor.fetchall()]
        print(f"Retrieved {len(compounds)} compounds from the database.")
        return compounds

    def add_compound(self, compound_data):
        try:
            cursor = self.conn.cursor()
            query = """
                INSERT INTO plant_compounds (
                    name, smiles, molecular_weight, logp, plant_source, biological_activity,
                    traditional_use, h_bond_donors, h_bond_acceptors,
                    polar_surface_area, rotatable_bonds
                ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
            """
            values = (
                compound_data['name'], compound_data['smiles'], compound_data['molecular_weight'],
                compound_data['logp'], compound_data['plant_source'], compound_data['biological_activity'],
                compound_data['traditional_use'], compound_data['h_bond_donors'], compound_data['h_bond_acceptors'],
                compound_data['polar_surface_area'], compound_data['rotatable_bonds']
            )
            cursor.execute(query, values)
            self.conn.commit()
            return True
        except sqlite3.Error as e:
            print(f"An error occurred: {e}")
            return False

    def get_optimization_data(self):
        return [
            {'molecular_weight': 300, 'logp': 2.5, 'h_bond_donors': 2, 'h_bond_acceptors': 5, 'optimization_score': 0.8},
            {'molecular_weight': 400, 'logp': 3.5, 'h_bond_donors': 1, 'h_bond_acceptors': 7, 'optimization_score': 0.6}
        ]

    def search(self, query=None, filters=None, page=1, per_page=10):
        return search_compounds(self.conn, query, filters, page, per_page)

    def get_therapeutic_areas(self):
        return get_therapeutic_areas(self.conn)

    def predict_therapeutic_areas(self, smiles):
        return predict_therapeutic_areas(smiles, self.model, self.scaler)

    def optimize_structure(self, smiles, optimization_type='functional_group', params={'target': '[OH]', 'replacement': '[NH2]'}):
        return optimize_structure(smiles, optimization_type, params)

    def explore_chemical_space(self, smiles, num_iterations=10):
        return chemical_space_exploration(smiles, num_iterations)

    def retrosynthesis_informed_optimization(self, smiles):
        compound = get_compound_by_smiles(self.conn, smiles)
        if compound:
            retrosynthesis_data = get_retrosynthesis_data(self.conn, compound[0])
            optimized_smiles = retrosynthesis_informed_optimization(smiles, retrosynthesis_data)
            store_retrosynthesis_informed_optimization(self.conn, compound[0], optimized_smiles)
            return optimized_smiles
        else:
            return None

    def train_ml_models(self):
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

    def optimize_compounds(self, compounds):
        optimization_data = self.get_optimization_data()
        optimization_model = train_optimization_model(optimization_data)
        prioritized_compounds = prioritize_optimization_strategies(compounds, optimization_model)
        return prioritized_compounds

    def perform_retrosynthesis(self, smiles):
        return perform_retrosynthesis(smiles)


    def close_connection(self):
        self.conn.close()
