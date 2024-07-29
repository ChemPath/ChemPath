import sqlite3
import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors
from chempath_core import optimize_structure, chemical_space_exploration, get_retrosynthesis_data, retrosynthesis_informed_optimization, store_retrosynthesis_informed_optimization
from sklearn.model_selection import train_test_split
from sklearn.ensemble import RandomForestClassifier
from sklearn.preprocessing import StandardScaler
from chempath_utils import validate_smiles, search_compounds, get_therapeutic_areas
import pandas as pd
import joblib
import requests
import logging

class ChemPathAPI:
    def __init__(self, db_path):
        self.conn = sqlite3.connect(db_path)
        self.cursor = self.conn.cursor()
        self.create_tables()
        self.create_indexes()
        logging.basicConfig(level=logging.INFO)

    def create_tables(self):
        self.cursor.execute("""
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
        self.cursor.execute("""
            CREATE TABLE IF NOT EXISTS predicted_therapeutic_areas (
                id INTEGER PRIMARY KEY,
                compound_id INTEGER,
                therapeutic_area TEXT
            )
        """)
        self.conn.commit()

    def create_indexes(self):
        self.cursor.execute("CREATE INDEX IF NOT EXISTS idx_smiles ON plant_compounds (smiles)")
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
            logging.error(f"Error inserting compound: {e}")
            self.conn.rollback()
            raise

    def get_external_compound_info(self, compound_name):
        try:
            response = requests.get(f"https://api.example.com/compounds/{compound_name}")
            response.raise_for_status()
            return response.json()
        except requests.RequestException as e:
            logging.error(f"Error fetching compound info: {str(e)}")
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
    
        predicted_areas = all_therapeutic_areas or self.get_therapeutic_areas()
    
        if ph is not None:
            if ph < 0 or ph > 14:
                predicted_areas.extend(["Antimicrobial", "Antifungal", "Enzyme Inhibitor"])
    
        if temperature is not None:
            if temperature < -273.15 or temperature > 100:
                predicted_areas.extend(["Protein Stabilizer", "Cryoprotectant", "Thermostable Enzyme"])
    
        return list(set(predicted_areas))

    def search_compounds(self, query=None, filters=None, page=1, per_page=10):
        return search_compounds(self.conn, query, filters, page, per_page)

    def perform_retrosynthesis(self, smiles):
        from retrosynthesis import perform_retrosynthesis as retro_synthesis
        compound = self.get_compound_by_smiles(smiles)
        if compound:
            compound_id = compound[0]  # Assuming the first element is the compound ID
            return retro_synthesis(self.conn, compound_id, smiles)
        else:
            return None  # Or handle the case where the compound is not found

    def get_compound_by_smiles(self, smiles):
        self.cursor.execute("SELECT * FROM plant_compounds WHERE smiles = ?", (smiles,))
        return self.cursor.fetchone()

    def retrosynthesis_informed_optimization(self, smiles):
        compound = self.get_compound_by_smiles(smiles)
        if compound:
            retrosynthesis_data = get_retrosynthesis_data(self.conn, compound[0])
            optimized_smiles = self._perform_retrosynthesis_optimization(smiles, retrosynthesis_data)
            store_retrosynthesis_informed_optimization(self.conn, compound[0], optimized_smiles)
            return optimized_smiles
        else:
            return None

    def _perform_retrosynthesis_optimization(self, smiles, retrosynthesis_data):
        mol = Chem.MolFromSmiles(smiles)
        if not mol:
            return None
        
        for step in retrosynthesis_data:
            reaction = AllChem.ReactionFromSmarts(step['reaction_smarts'])
            products = reaction.RunReactants((mol,))
            if products:
                best_product = max(products, key=lambda p: Chem.Crippen.MolLogP(p[0]))
                mol = best_product[0]
        
        mol = Chem.AddHs(mol)
        AllChem.EmbedMolecule(mol, randomSeed=42)
        AllChem.UFFOptimizeMolecule(mol, maxIters=200)
        mol = Chem.RemoveHs(mol)
        
        return Chem.MolToSmiles(mol)

    def verify_compound_insertion(self, compound_id):
        self.cursor.execute("SELECT * FROM plant_compounds WHERE id = ?", (compound_id,))
        result = self.cursor.fetchone()
        if result:
            logging.info(f"Compound {compound_id} successfully inserted")
            return True
        else:
            logging.warning(f"Compound {compound_id} not found in database")
            return False

    def table_exists(self, table_name):
        self.cursor.execute(f"SELECT name FROM sqlite_master WHERE type='table' AND name='{table_name}'")
        return self.cursor.fetchone() is not None

    def get_all_compounds(self):
        self.cursor.execute("SELECT name, molecular_weight, logp, h_bond_donors, h_bond_acceptors FROM plant_compounds")
        compounds = [{'name': row[0], 'molecular_weight': row[1], 'logp': row[2], 'h_bond_donors': row[3], 'h_bond_acceptors': row[4]} for row in self.cursor.fetchall()]
        logging.info(f"Retrieved {len(compounds)} compounds from the database.")
        return compounds

    def add_compound(self, compound_data):
        try:
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
            self.cursor.execute(query, values)
            self.conn.commit()
            return True
        except sqlite3.Error as e:
            logging.error(f"An error occurred: {e}")
            return False

    def get_optimization_data(self):
        return [
            {'molecular_weight': 300, 'logp': 2.5, 'h_bond_donors': 2, 'h_bond_acceptors': 5, 'optimization_score': 0.8},
            {'molecular_weight': 400, 'logp': 3.5, 'h_bond_donors': 1, 'h_bond_acceptors': 7, 'optimization_score': 0.6}
        ]

    def optimize_structure(self, smiles, optimization_type='functional_group', params={'target': '[OH]', 'replacement': '[NH2]'}):
        if not Chem.MolFromSmiles(smiles):
            raise ValueError("Invalid SMILES string")
        return optimize_structure(smiles, optimization_type, params)

    def validate_smiles(self, smiles):
        return validate_smiles(smiles)

    def delete_database(self):
        import os
        self.conn.close()
        try:
            os.remove(self.db_path)
        except PermissionError:
            logging.error(f"Unable to delete {self.db_path}. It may be in use by another process.")

    def explore_chemical_space(self, smiles, num_iterations=10):
        return chemical_space_exploration(smiles, num_iterations)

    def train_ml_models(self):
        def smiles_to_fingerprint(smiles):
            mol = Chem.MolFromSmiles(smiles)
            return AllChem.GetMorganFingerprintAsBitVect(mol, 2, nBits=1024)

        def train_model(X, y):
            scaler = StandardScaler()
            X_scaled = scaler.fit_transform(X)
            model = RandomForestClassifier(n_estimators=100, random_state=42)
            model.fit(X_scaled, y)
            return model, scaler

        self.cursor.execute("SELECT smiles, retrosynthesis_feasibility, reaction_class FROM compounds")
        data = self.cursor.fetchall()

        X = np.array([smiles_to_fingerprint(row[0]) for row in data])
        y_retro = np.array([row[1] for row in data])
        y_reaction = np.array([row[2] for row in data])

        self.retro_model, self.retro_scaler = train_model(X, y_retro)
        self.reaction_model, self.reaction_scaler = train_model(X, y_reaction)

    def calculate_molecular_descriptors(self, smiles):
        mol = Chem.MolFromSmiles(smiles)
        return {
            'molecular_weight': Descriptors.MolWt(mol),
            'logp': Descriptors.MolLogP(mol),
            'h_bond_donors': Descriptors.NumHDonors(mol),
            'h_bond_acceptors': Descriptors.NumHAcceptors(mol),
            'polar_surface_area': Descriptors.TPSA(mol),
            'rotatable_bonds': Descriptors.NumRotatableBonds(mol)
        }

    def predict_therapeutic_areas_helper(self, smiles, all_therapeutic_areas, ph=None, temperature=None):
        model = joblib.load('therapeutic_area_predictor.joblib')
        scaler = joblib.load('feature_scaler.joblib')
        mlb = joblib.load('multilabel_binarizer.joblib')

        descriptors = self.calculate_molecular_descriptors(smiles)
        features = pd.DataFrame([descriptors])
        
        features = features[['molecular_weight', 'logp', 'h_bond_donors', 'h_bond_acceptors', 'polar_surface_area', 'rotatable_bonds']]

        features_scaled = scaler.transform(features)

        predictions = model.predict_proba(features_scaled)
        threshold = 0.1
        
        predicted_areas = [area for area, prob in zip(mlb.classes_, predictions[0]) if prob.max() > threshold]

        return predicted_areas
    def is_stable_at_ph(self, area, ph):
        # Implement logic to determine if the therapeutic area is stable at the given pH
        return True

    def is_stable_at_temperature(self, area, temperature):
        # Implement logic to determine if the therapeutic area is stable at the given temperature
        return True

    def update_compound(self, compound_id, updated_data):
        try:
            query = """
                UPDATE plant_compounds
                SET name = ?, smiles = ?, molecular_weight = ?, logp = ?, plant_source = ?,
                    biological_activities = ?, traditional_use = ?, h_bond_donors = ?,
                    h_bond_acceptors = ?, polar_surface_area = ?, rotatable_bonds = ?
                WHERE id = ?
            """
            values = (
                updated_data['name'], updated_data['smiles'], updated_data['molecular_weight'],
                updated_data['logp'], updated_data['plant_source'], updated_data['biological_activities'],
                updated_data['traditional_use'], updated_data['h_bond_donors'], updated_data['h_bond_acceptors'],
                updated_data['polar_surface_area'], updated_data['rotatable_bonds'], compound_id
            )
            self.cursor.execute(query, values)
            self.conn.commit()
            return True
        except sqlite3.Error as e:
            logging.error(f"Error updating compound: {e}")
            self.conn.rollback()
            return False

    def get_compound_statistics(self):
        self.cursor.execute("""
            SELECT 
                COUNT(*) as total_compounds,
                AVG(molecular_weight) as avg_molecular_weight,
                AVG(logp) as avg_logp,
                AVG(h_bond_donors) as avg_h_bond_donors,
                AVG(h_bond_acceptors) as avg_h_bond_acceptors
            FROM plant_compounds
        """)
        return self.cursor.fetchone()

    def export_compounds_to_csv(self, filename):
        import csv
        self.cursor.execute("SELECT * FROM plant_compounds")
        with open(filename, 'w', newline='') as csvfile:
            csv_writer = csv.writer(csvfile)
            csv_writer.writerow([i[0] for i in self.cursor.description])  # write headers
            csv_writer.writerows(self.cursor)
        logging.info(f"Compounds exported to {filename}")

    def import_compounds_from_csv(self, filename):
        import csv
        with open(filename, 'r') as csvfile:
            csv_reader = csv.DictReader(csvfile)
            for row in csv_reader:
                self.add_compound(row)
        logging.info(f"Compounds imported from {filename}")

    def get_similar_compounds(self, smiles, threshold=0.7):
        from rdkit import Chem, DataStructs
        from rdkit.Chem import AllChem

        target_mol = Chem.MolFromSmiles(smiles)
        target_fp = AllChem.GetMorganFingerprintAsBitVect(target_mol, 2, nBits=1024)
        
        similar_compounds = []
        for compound in self.get_all_compounds():
            mol = Chem.MolFromSmiles(compound['smiles'])
            fp = AllChem.GetMorganFingerprintAsBitVect(mol, 2, nBits=1024)
            similarity = DataStructs.TanimotoSimilarity(target_fp, fp)
            if similarity >= threshold:
                similar_compounds.append((compound, similarity))
        
        return sorted(similar_compounds, key=lambda x: x[1], reverse=True)

# End of ChemPathAPI class
