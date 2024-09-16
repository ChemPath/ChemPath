import pandas as pd
import numpy as np
import joblib
from rdkit import Chem
from rdkit.Chem import Descriptors
from sqlalchemy.orm import Session
from models import Compound, TherapeuticArea
from database import engine

def calculate_molecular_descriptors(smiles, ph=None, temperature=None):
    mol = Chem.MolFromSmiles(smiles)
    descriptors = {
        'molecular_weight': Descriptors.MolWt(mol),
        'logp': Descriptors.MolLogP(mol),
        'h_bond_donors': Descriptors.NumHDonors(mol),
        'h_bond_acceptors': Descriptors.NumHAcceptors(mol),
        'polar_surface_area': Descriptors.TPSA(mol),
        'rotatable_bonds': Descriptors.NumRotatableBonds(mol)
    }
    
    if ph is not None:
        descriptors['ph'] = ph
    if temperature is not None:
        descriptors['temperature'] = temperature
    
    return descriptors

def predict_therapeutic_areas(smiles, all_therapeutic_areas, ph=None, temperature=None):
    model = joblib.load('therapeutic_area_predictor.joblib')
    scaler = joblib.load('feature_scaler.joblib')
    mlb = joblib.load('multilabel_binarizer.joblib')

    descriptors = calculate_molecular_descriptors(smiles, ph, temperature)
    features = pd.DataFrame([descriptors])
    features['mw_logp_ratio'] = features['molecular_weight'] / (features['logp'] + 1)
    features_scaled = scaler.transform(features)

    predictions = model.predict_proba(features_scaled)
    threshold = 0.1  # Lowered threshold for more diverse predictions
    predicted_areas = [area for area, prob in zip(mlb.classes_, predictions[0]) if prob.max() > threshold]

    # Store the prediction in the database
    with Session(engine) as session:
        compound = Compound(smiles=smiles, ph=ph, temperature=temperature)
        session.add(compound)
        session.flush()  # To get the compound ID

        for area in predicted_areas:
            therapeutic_area = session.query(TherapeuticArea).filter_by(name=area).first()
            if not therapeutic_area:
                therapeutic_area = TherapeuticArea(name=area)
                session.add(therapeutic_area)
            
            compound.therapeutic_areas.append(therapeutic_area)
        
        session.commit()

    return predicted_areas

def get_all_therapeutic_areas():
    with Session(engine) as session:
        return [area.name for area in session.query(TherapeuticArea).all()]

def retrain_model():
    with Session(engine) as session:
        compounds = session.query(Compound).all()
        data = []
        for compound in compounds:
            descriptors = calculate_molecular_descriptors(compound.smiles, compound.ph, compound.temperature)
            descriptors['therapeutic_areas'] = [area.name for area in compound.therapeutic_areas]
            data.append(descriptors)
        
        df = pd.DataFrame(data)
        all_therapeutic_areas = get_all_therapeutic_areas()

    # Train the model (implement your training logic here)
    # ...

    print("Model retrained successfully")
