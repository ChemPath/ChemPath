import pandas as pd
import numpy as np
import joblib
from rdkit import Chem
from rdkit.Chem import Descriptors

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

def predict_therapeutic_areas(smiles, all_therapeutic_areas):
    model = joblib.load('therapeutic_area_predictor.joblib')
    scaler = joblib.load('feature_scaler.joblib')
    mlb = joblib.load('multilabel_binarizer.joblib')

    descriptors = calculate_molecular_descriptors(smiles)
    features = pd.DataFrame([descriptors])
    features['mw_logp_ratio'] = features['molecular_weight'] / (features['logp'] + 1)
    features_scaled = scaler.transform(features)

    predictions = model.predict_proba(features_scaled)
    threshold = 0.1  # Lowered threshold for more diverse predictions
    predicted_areas = [area for area, prob in zip(mlb.classes_, predictions[0]) if prob.max() > threshold]

    return predicted_areas


