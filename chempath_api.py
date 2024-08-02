import sqlite3
import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors
from src.core.chempath_core import optimize_structure, chemical_space_exploration, get_retrosynthesis_data, retrosynthesis_informed_optimization, store_retrosynthesis_informed_optimization
from sklearn.model_selection import train_test_split
from sklearn.ensemble import RandomForestClassifier
from sklearn.preprocessing import StandardScaler
from chempath_utils import validate_smiles, search_compounds, get_therapeutic_areas
import pandas as pd
import joblib
import requests
import logging
from sqlalchemy import create_engine, Column, Integer, String, Float, Text
from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy.orm import sessionmaker
from sqlalchemy.pool import QueuePool
import threading
import concurrent.futures
from threading import Lock
from sqlalchemy import or_

Base = declarative_base()

class PlantCompound(Base):
    __tablename__ = 'plant_compounds'

    id = Column(Integer, primary_key=True)
    name = Column(String)
    smiles = Column(String)
    molecular_weight = Column(Float)
    logp = Column(Float)
    plant_source = Column(String)
    biological_activities = Column(Text)
    traditional_use = Column(Text)
    h_bond_donors = Column(Integer)
    h_bond_acceptors = Column(Integer)
    polar_surface_area = Column(Float)
    rotatable_bonds = Column(Integer)

    
    def to_dict(self):
        return {
            'id': self.id,
            'name': self.name,
            'smiles': self.smiles,
            'molecular_weight': self.molecular_weight,
            'logp': self.logp,
            'plant_source': self.plant_source,
            'biological_activity': self.biological_activity,
            'traditional_use': self.traditional_use
        }

class PredictedTherapeuticArea(Base):
    __tablename__ = 'predicted_therapeutic_areas'

    id = Column(Integer, primary_key=True)
    compound_id = Column(Integer)
    therapeutic_area = Column(String)

class ChemPathAPI:
    def __init__(self, db_path):
        self.engine = create_engine(f'sqlite:///{db_path}',
                                    poolclass=QueuePool,
                                    pool_size=10,
                                    max_overflow=20,
                                    pool_timeout=30)
        Base.metadata.create_all(self.engine)
        self.Session = sessionmaker(bind=self.engine)
        self.local = threading.local()
        self.thread_pool = concurrent.futures.ThreadPoolExecutor(max_workers=10)
        self.lock = Lock()
        logging.basicConfig(level=logging.INFO)

    def process_user_request(self, request_type, *args, **kwargs):
        return self.thread_pool.submit(self._handle_request, request_type, *args, **kwargs)

    def _handle_request(self, request_type, *args, **kwargs):
        with self.lock:
            if request_type == 'search':
                return self.search_compounds(*args, **kwargs)
            elif request_type == 'predict':
                return self.predict_therapeutic_areas(*args, **kwargs)
            elif request_type == 'add':
                return self.add_compound(*args, **kwargs)

    def get_session(self):
        if not hasattr(self.local, 'session'):
            self.local.session = self.Session()
        return self.local.session

    def close_session(self):
        if hasattr(self.local, 'session'):
            self.local.session.close()
            del self.local.session

    def insert_compound(self, compound_data):
        try:
            session = self.get_session()
            mol = Chem.MolFromSmiles(compound_data['smiles'])
            if mol is None:
                raise ValueError("Invalid SMILES structure")
            
            new_compound = PlantCompound(
                name=compound_data['name'],
                smiles=compound_data['smiles'],
                molecular_weight=compound_data['molecular_weight'],
                logp=compound_data['logp'],
                plant_source=compound_data['plant_source'],
                biological_activities=compound_data['biological_activities'],
                traditional_use=compound_data['traditional_use'],
                h_bond_donors=compound_data['h_bond_donors'],
                h_bond_acceptors=compound_data['h_bond_acceptors'],
                polar_surface_area=compound_data['polar_surface_area'],
                rotatable_bonds=compound_data['rotatable_bonds']
            )
            session.add(new_compound)
            session.commit()
            return new_compound.id
        except Exception as e:
            logging.error(f"Error inserting compound: {e}")
            session.rollback()
            raise
        finally:
            self.close_session()

    def get_external_compound_info(self, compound_name):
        try:
            response = requests.get(f"https://api.example.com/compounds/{compound_name}")
            response.raise_for_status()
            return response.json()
        except requests.RequestException as e:
            logging.error(f"Error fetching compound info: {str(e)}")
            raise ValueError(f"Error fetching compound info: {str(e)}")
    def get_compound(self, compound_id):
        session = self.get_session()
        try:
            compound = session.query(PlantCompound).filter(PlantCompound.id == compound_id).first()
            if compound:
                return {
                    'name': compound.name,
                    'smiles': compound.smiles,
                    'molecular_weight': compound.molecular_weight,
                    'logp': compound.logp,
                    'plant_source': compound.plant_source,
                    'biological_activities': compound.biological_activities,
                    'traditional_use': compound.traditional_use
                }
            return None
        finally:
            self.close_session()    

    def predict_therapeutic_areas(self, smiles, all_therapeutic_areas=None, ph=None, temperature=None):
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            raise ValueError("Invalid SMILES string")

        predicted_areas = all_therapeutic_areas or self.get_therapeutic_areas()

        # Use machine learning model for prediction
        model = joblib.load('therapeutic_area_predictor.joblib')
        scaler = joblib.load('feature_scaler.joblib')
        mlb = joblib.load('multilabel_binarizer.joblib')

        descriptors = self.calculate_molecular_descriptors(smiles)
        features = pd.DataFrame([descriptors])
        features = features[['molecular_weight', 'logp', 'h_bond_donors', 'h_bond_acceptors', 'polar_surface_area', 'rotatable_bonds']]

        features_scaled = scaler.transform(features)
        predictions = model.predict_proba(features_scaled)
        threshold = 0.1
    
        ml_predicted_areas = [area for area, prob in zip(mlb.classes_, predictions[0]) if prob.max() > threshold]
        predicted_areas.extend(ml_predicted_areas)

        if ph is not None:
            if ph < 0 or ph > 14:
                predicted_areas.extend(["Antimicrobial", "Antifungal", "Enzyme Inhibitor"])

        if temperature is not None:
            if temperature < -273.15 or temperature > 100:
                predicted_areas.extend(["Protein Stabilizer", "Cryoprotectant", "Thermostable Enzyme"])

        return list(set(predicted_areas))

    def search_compounds(self, query=None, filters=None, page=1, per_page=10):
        session = self.get_session()
        try:
            compounds = session.query(PlantCompound)
            if query:
                compounds = compounds.filter(PlantCompound.name.ilike(f"%{query}%"))
            if filters:
                for key, value in filters.items():
                    compounds = compounds.filter(getattr(PlantCompound, key) == value)
            compounds = compounds.offset((page - 1) * per_page).limit(per_page)
            return [c.to_dict() for c in compounds.all()]
        except Exception as e:
            # Handle any exceptions that might occur during the database query
            print(f"An error occurred: {str(e)}")
            return []
        finally:
            self.close_session()

    def perform_retrosynthesis(self, smiles):
        from retrosynthesis import perform_retrosynthesis as retro_synthesis
        compound = self.get_compound_by_smiles(smiles)
        if compound:
            compound_id = compound['id']
            session = self.get_session()
            try:
                return retro_synthesis(session, compound_id, smiles)
            finally:
                self.close_session()
        else:
            return None
    def get_therapeutic_areas(self):
        session = self.get_session()
        try:
            areas = session.query(PredictedTherapeuticArea.therapeutic_area).distinct().all()
            return [area[0] for area in areas]
        finally:
            self.close_session()    

    def add_compound(self, compound_data):
        session = self.get_session()
        try:
            new_compound = PlantCompound(**compound_data)
            session.add(new_compound)
            session.commit()
            return True
        except Exception as e:
            session.rollback()
            logging.error(f"Error adding compound: {e}")
            return False
        finally:
            self.close_session()

    def get_compound_by_smiles(self, smiles):
        session = self.get_session()
        try:
            compound = session.query(PlantCompound).filter(PlantCompound.smiles == smiles).first()
            if compound:
                return {
                    'id': compound.id,
                    'name': compound.name,
                    'smiles': compound.smiles,
                    'molecular_weight': compound.molecular_weight,
                    'logp': compound.logp,
                    'plant_source': compound.plant_source,
                    'biological_activities': compound.biological_activities,
                    'traditional_use': compound.traditional_use
                }
            return None
        finally:
            self.close_session()

    def retrosynthesis_informed_optimization(self, smiles):
        compound = self.get_compound_by_smiles(smiles)
        if compound:
            session = self.get_session()
            try:
                retrosynthesis_data = get_retrosynthesis_data(session, compound['id'])
                optimized_smiles = self._perform_retrosynthesis_optimization(smiles, retrosynthesis_data)
                store_retrosynthesis_informed_optimization(session, compound['id'], optimized_smiles)
                return optimized_smiles
            finally:
                self.close_session()
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
        session = self.get_session()
        try:
            compound = session.query(PlantCompound).filter(PlantCompound.id == compound_id).first()
            if compound:
                logging.info(f"Compound {compound_id} successfully inserted")
                return True
            else:
                logging.warning(f"Compound {compound_id} not found in database")
                return False
        finally:
            self.close_session()

    def table_exists(self, table_name):
        return self.engine.has_table(table_name)

    def get_all_compounds(self):
        session = self.get_session()
        try:
            compounds = session.query(PlantCompound.name, PlantCompound.molecular_weight, PlantCompound.logp, 
                                  PlantCompound.h_bond_donors, PlantCompound.h_bond_acceptors).all()
            result = [{'name': c.name, 'molecular_weight': c.molecular_weight, 'logp': c.logp, 
                   'h_bond_donors': c.h_bond_donors, 'h_bond_acceptors': c.h_bond_acceptors} for c in compounds]
            logging.info(f"Retrieved {len(result)} compounds from the database.")
            return result
        finally:
            self.close_session()

    def calculate_molecular_descriptors(self, smiles):
        mol = Chem.MolFromSmiles(smiles)
        if not mol:
            raise ValueError("Invalid SMILES string")
    
        molecular_weight = Descriptors.MolWt(mol)
        logp = Descriptors.MolLogP(mol)
        h_bond_donors = Descriptors.NumHDonors(mol)
        h_bond_acceptors = Descriptors.NumHAcceptors(mol)
        polar_surface_area = Descriptors.TPSA(mol)
        rotatable_bonds = Descriptors.NumRotatableBonds(mol)
    
    def calculate_molecular_descriptors(self, smiles):
        mol = Chem.MolFromSmiles(smiles)
        if not mol:
            raise ValueError("Invalid SMILES string")

        molecular_weight = Descriptors.MolWt(mol)
        logp = Descriptors.MolLogP(mol)
        h_bond_donors = Descriptors.NumHDonors(mol)
        h_bond_acceptors = Descriptors.NumHAcceptors(mol)
        polar_surface_area = Descriptors.TPSA(mol)
        rotatable_bonds = Descriptors.NumRotatableBonds(mol)

        return {
            'molecular_weight': molecular_weight,
            'logp': logp,
            'h_bond_donors': h_bond_donors,
            'h_bond_acceptors': h_bond_acceptors,
            'polar_surface_area': polar_surface_area,
            'rotatable_bonds': rotatable_bonds
        }
    
    def optimize_structure(self, smiles, optimization_type='functional_group', params={'target': '[OH]', 'replacement': '[NH2]'}):
        if not Chem.MolFromSmiles(smiles):
            raise ValueError("Invalid SMILES string")
        return optimize_structure(smiles, optimization_type, params)

    def validate_smiles(self, smiles):
        return validate_smiles(smiles)

    def delete_database(self):
        import os
        self.engine.dispose()
        try:
            os.remove(self.engine.url.database)
        except FileNotFoundError:
            logging.warning("Database file not found.")    

if __name__ == "__main__":
    api = ChemPathAPI("chempath.db")# Example usage:
# smiles = "CC(=O)Nc1ccc(cc1)S(=O)(=O)N"
# print(api.predict_therapeutic_areas(smiles))    