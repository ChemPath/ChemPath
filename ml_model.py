import warnings
warnings.filterwarnings("ignore", category=RuntimeWarning)
import pandas as pd
import numpy as np
from rdkit import Chem
from rdkit.Chem import Descriptors, AllChem
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import StandardScaler
from sklearn.multioutput import MultiOutputClassifier
from sklearn.ensemble import RandomForestClassifier, GradientBoostingClassifier
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
        'rotatable_bonds': Descriptors.NumRotatableBonds(mol),
        'aromatic_rings': Descriptors.NumAromaticRings(mol),
        'heavy_atom_count': Descriptors.HeavyAtomCount(mol),
        'complexity': Descriptors.BertzCT(mol)
    }

def prepare_data(conn):
    query = "SELECT * FROM plant_compounds"
    df = pd.read_sql_query(query, conn)
    print(f"Retrieved {len(df)} rows from plant_compounds")
    
    features = df.apply(lambda row: pd.Series(calculate_molecular_descriptors(row['smiles'])), axis=1)
    target = df['biological_activity'].str.get_dummies(sep=', ')
    
    return features, target

def train_model(conn):
    features, target = prepare_data(conn)
    if features is None or target is None:
        print("Failed to prepare data. Aborting model training.")
        return None, None
    
    X_train, X_test, y_train, y_test = train_test_split(features, target, test_size=0.2, random_state=42)
    
    scaler = StandardScaler()
    X_train_scaled = scaler.fit_transform(X_train)
    X_test_scaled = scaler.transform(X_test)
    
    # Random Forest Model
    rf_model = MultiOutputClassifier(RandomForestClassifier(n_estimators=100, random_state=42))
    rf_model.fit(X_train_scaled, y_train)
    
    rf_predictions = rf_model.predict(X_test_scaled)
    print("Random Forest Model Performance:")
    print(classification_report(y_test, rf_predictions, target_names=target.columns))
    
    # Gradient Boosting Model
    gb_model = MultiOutputClassifier(GradientBoostingClassifier(n_estimators=100, random_state=42))
    gb_model.fit(X_train_scaled, y_train)
    
    gb_predictions = gb_model.predict(X_test_scaled)
    print("Gradient Boosting Model Performance:")
    print(classification_report(y_test, gb_predictions, target_names=target.columns))
    
    # Save models
    joblib.dump(rf_model, 'random_forest_model.joblib')
    joblib.dump(gb_model, 'gradient_boosting_model.joblib')
    joblib.dump(scaler, 'feature_scaler.joblib')
    
    return rf_model, gb_model, scaler

def predict_therapeutic_areas(smiles, rf_model, gb_model, scaler):
    descriptors = calculate_molecular_descriptors(smiles)
    features = pd.DataFrame([descriptors])
    features_scaled = scaler.transform(features)
    
    rf_predictions = rf_model.predict(features_scaled)
    gb_predictions = gb_model.predict(features_scaled)
    
    rf_therapeutic_areas = [rf_model.classes_[i] for i, pred in enumerate(rf_predictions[0]) if pred]
    gb_therapeutic_areas = [gb_model.classes_[i] for i, pred in enumerate(gb_predictions[0]) if pred]
    
    return {
        'random_forest': rf_therapeutic_areas,
        'gradient_boosting': gb_therapeutic_areas
    }

if __name__ == "__main__":
    from src.database.chempath_database import create_connection
    conn = create_connection("chempath_database.db")
    rf_model, gb_model, scaler = train_model(conn)
    conn.close()
    
    test_smiles = "CC(=O)OC1=CC=CC=C1C(=O)O"  # Aspirin
    predictions = predict_therapeutic_areas(test_smiles, rf_model, gb_model, scaler)
    print("Predicted Therapeutic Areas:")
    print("Random Forest:", predictions['random_forest'])
    print("Gradient Boosting:", predictions['gradient_boosting'])
