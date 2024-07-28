import sqlite3
import numpy as np
from rdkit import Chem
from chempath_core import optimize_structure, chemical_space_exploration, get_retrosynthesis_data, retrosynthesis_informed_optimization, store_retrosynthesis_informed_optimization
from sklearn.model_selection import train_test_split
from sklearn.ensemble import RandomForestClassifier
from sklearn.preprocessing import StandardScaler
from chempath_utils import validate_smiles
import pandas as pd
import joblib
from rdkit.Chem import Descriptors
import requests

def calculate_molecular_descriptors(smiles):
    mol = Chem.MolFromSmiles(smiles)
    return {
        'molecular_weight': Descriptors.MolWt(mol),
        'logp': Descriptors.MolLogP(mol),
        'h_bond_donors': Descriptors.NumHDonors(mol),
        'h_bond_acceptors': Descriptors.NumHAcceptors(mol),
        'polar_surface_area': Descriptors.TPSA(mol),
        'rotatable_bonds': Descriptors.NumRotatableBonds(mol)
    }

def predict_therapeutic_areas_helper(smiles, all_therapeutic_areas, ph=None, temperature=None):
    model = joblib.load('therapeutic_area_predictor.joblib')
    scaler = joblib.load('feature_scaler.joblib')
    mlb = joblib.load('multilabel_binarizer.joblib')

    descriptors = calculate_molecular_descriptors(smiles)
    features = pd.DataFrame([descriptors])
    
    # Remove features not present during training
    features = features[['molecular_weight', 'logp', 'h_bond_donors', 'h_bond_acceptors', 'polar_surface_area', 'rotatable_bonds']]

    features_scaled = scaler.transform(features)

    predictions = model.predict_proba(features_scaled)
    threshold = 0.1
    
    predicted_areas = [area for area, prob in zip(mlb.classes_, predictions[0]) if prob.max() > threshold]

    return predicted_areas


def is_stable_at_ph(area, ph):
    # Implement logic to determine if the therapeutic area is stable at the given pH
    return True

def is_stable_at_temperature(area, temperature):
    # Implement logic to determine if the therapeutic area is stable at the given temperature
    return True

class ChemPathAPI:
    def __init__(self, db_path):
        self.conn = sqlite3.connect(db_path)
        self.cursor = self.conn.cursor()
        self.create_tables()
        self.create_indexes()

    def create_tables(self):
        cursor = self.conn.cursor()
        cursor.execute("""
            CREATE TABLE IF NOT EXISTS plant_compounds (
                id INTEGER PRIMARY KEY,
                name TEXT,
                smiles TEXT,
                molecular_weight REAL,
                logp REAL,
                plant_source TEXT,
                biological_activities TEXT,
                traditional_use TEXT,
                h_bond_donors INTEGER,
                h_bond_acceptors INTEGER,
                polar_surface_area REAL,
                rotatable_bonds INTEGER
            )
        """)
        cursor.execute("""
            CREATE TABLE IF NOT EXISTS predicted_therapeutic_areas (
                id INTEGER PRIMARY KEY,
                compound_id INTEGER,
                therapeutic_area TEXT
            )
        """)
        self.conn.commit()

    def create_indexes(self):
        cursor = self.conn.cursor()
        cursor.execute("CREATE INDEX IF NOT EXISTS idx_smiles ON plant_compounds (smiles)")
        self.conn.commit()

    def close_connection(self):
        if self.conn:
            self.conn.close()

    def insert_compound(self, compound_data):
        try:
            mol = Chem.MolFromSmiles(compound_data['smiles'])
            if mol is None:
                raise ValueError("Invalid SMILES structure")
            
            query = """
                INSERT INTO plant_compounds (
                    name, smiles, molecular_weight, logp, plant_source, biological_activities,
                    traditional_use, h_bond_donors, h_bond_acceptors,
                    polar_surface_area, rotatable_bonds
                ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
            """
            values = (
                compound_data['name'], compound_data['smiles'], compound_data['molecular_weight'],
                compound_data['logp'], compound_data['plant_source'], compound_data['biological_activities'],
                compound_data['traditional_use'], compound_data['h_bond_donors'], compound_data['h_bond_acceptors'],
                compound_data['polar_surface_area'], compound_data['rotatable_bonds']
            )
            self.cursor.execute(query, values)
            self.conn.commit()
            return self.cursor.lastrowid
        except Exception as e:
            print(f"Error inserting compound: {e}")
            raise
        
    def get_external_compound_info(self, compound_name):
        try:
            response = requests.get(f"https://api.example.com/compounds/{compound_name}")
            response.raise_for_status()
            return response.json()
        except requests.RequestException as e:
            raise ValueError(f"Error fetching compound info: {str(e)}")


    def get_compound(self, compound_id):
        self.cursor.execute("SELECT * FROM plant_compounds WHERE id = ?", (compound_id,))
        result = self.cursor.fetchone()
        if result:
            return {
                'name': result[1],
                'smiles': result[2],
                'molecular_weight': result[3],
                'logp': result[4],
                'plant_source': result[5],
                'biological_activities': result[6],
                'traditional_use': result[7]
        }
        return None

    def predict_therapeutic_areas(self, smiles, all_therapeutic_areas=None, ph=None, temperature=None):
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            raise ValueError("Invalid SMILES string")
    
        predicted_areas = all_therapeutic_areas or self.all_areas
    
        if ph is not None:
            if ph < 0 or ph > 14:
                predicted_areas.extend(["Antimicrobial", "Antifungal", "Enzyme Inhibitor"])
    
        if temperature is not None:
            if temperature < -273.15 or temperature > 100:
                predicted_areas.extend(["Protein Stabilizer", "Cryoprotectant", "Thermostable Enzyme"])
    
        return list(set(predicted_areas))




    def verify_compound_insertion(self, compound_id):
        cursor = self.conn.cursor()
        cursor.execute("SELECT * FROM plant_compounds WHERE id = ?", (compound_id,))
        result = cursor.fetchone()
        if result:
            print(f"Compound {compound_id} successfully inserted")
            return True
        else:
            print(f"Compound {compound_id} not found in database")
            return False

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
        from chempath_utils import search_compounds
        return search_compounds(self.conn, query, filters, page, per_page)

    def get_therapeutic_areas(self):
        from chempath_utils import get_therapeutic_areas
        return get_therapeutic_areas(self.conn)

    def optimize_structure(self, smiles, optimization_type='functional_group', params={'target': '[OH]', 'replacement': '[NH2]'}):
        if not Chem.MolFromSmiles(smiles):
            raise ValueError("Invalid SMILES string")

        return optimize_structure(smiles, optimization_type, params)

    def get_compound_by_smiles(self, smiles):
        cursor = self.conn.cursor()
        cursor.execute("SELECT * FROM plant_compounds WHERE smiles = ?", (smiles,))
        return cursor.fetchone()

    def validate_smiles(smiles):
        try:
            mol = Chem.MolFromSmiles(smiles)
            return mol is not None
        except Exception:
            return False
    
    def delete_database(self):
        self.conn.close()
        try:
            os.remove(self.db_path)
        except PermissionError:
            print(f"Unable to delete {self.db_path}. It may be in use by another process.")

    def explore_chemical_space(self, smiles, num_iterations=10):
        return chemical_space_exploration(smiles, num_iterations)

    def retrosynthesis_informed_optimization(self, smiles):
        compound = self.get_compound_by_smiles(smiles)
        if compound:
            retrosynthesis_data = get_retrosynthesis_data(self.conn, compound[0])
            optimized_smiles = retrosynthesis_informed_optimization(smiles, retrosynthesis_data)
            store_retrosynthesis_informed_optimization(self.conn, compound[0], optimized_smiles)
            return optimized_smiles
        else:
            return None

        def train_ml_models(self):
            from rdkit import Chem
            from rdkit.Chem import AllChem

            def smiles_to_fingerprint(smiles):
                mol = Chem.MolFromSmiles(smiles)
                return AllChem.GetMorganFingerprintAsBitVect(mol, 2, nBits=1024)

            cursor = self.conn.cursor()
            cursor.execute("SELECT smiles, retrosynthesis_feasibility, reaction_class FROM compounds")
            data = cursor.fetchall()

            X = np.array([smiles_to_fingerprint(row[0]) for row in data])
            y_retro = np.array([row[1] for row in data])
            y_reaction = np.array([row[2] for row in data])

            self.retro_model, self.retro_scaler = train_retrosynthesis_model(X, y_retro)
            self.reaction_model, self.reaction_scaler = train_reaction_classification_model(X, y_reaction)