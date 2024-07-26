# File: ml_model.py
import warnings
warnings.filterwarnings("ignore", category=RuntimeWarning)
import pandas as pd
import numpy as np
from rdkit import Chem
from rdkit.Chem import Descriptors
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import StandardScaler
from sklearn.multioutput import MultiOutputClassifier
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import classification_report
import joblib

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

def prepare_data(conn):
    # Fetch data from the database
    query = "SELECT * FROM plant_compounds"
    df = pd.read_sql_query(query, conn)
    
    # Prepare features and target
    features = df[['molecular_weight', 'logp', 'h_bond_donors', 'h_bond_acceptors', 'polar_surface_area', 'rotatable_bonds']]
    target = df['biological_activity'].str.get_dummies(sep=', ')
    
    return features, target

def train_model(conn):
    features, target = prepare_data(conn)
    
    # Split the data
    X_train, X_test, y_train, y_test = train_test_split(features, target, test_size=0.2, random_state=42)
    
    # Scale the features
    scaler = StandardScaler()
    X_train_scaled = scaler.fit_transform(X_train)
    X_test_scaled = scaler.transform(X_test)
    
    # Train the model
    rf = RandomForestClassifier(n_estimators=100, random_state=42)
    model = MultiOutputClassifier(rf)
    model.fit(X_train_scaled, y_train)
    
    # Evaluate the model
    y_pred = model.predict(X_test_scaled)
    print(classification_report(y_test, y_pred, target_names=target.columns))
    
    # Save the model and scaler
    joblib.dump(model, 'therapeutic_area_predictor.joblib')
    joblib.dump(scaler, 'feature_scaler.joblib')
    
    return model, scaler

def predict_therapeutic_areas(smiles, model, scaler):
    descriptors = calculate_molecular_descriptors(smiles)
    features = pd.DataFrame([descriptors])
    features_scaled = scaler.transform(features)
    predictions = model.predict(features_scaled)
    therapeutic_areas = [model.classes_[i] for i, pred in enumerate(predictions[0]) if pred]
    return therapeutic_areas

