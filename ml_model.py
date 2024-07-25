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
        'MolWt': Descriptors.MolWt(mol),
        'LogP': Descriptors.MolLogP(mol),
        'NumHDonors': Descriptors.NumHDonors(mol),
        'NumHAcceptors': Descriptors.NumHAcceptors(mol),
        'NumRotatableBonds': Descriptors.NumRotatableBonds(mol),
        'NumAromaticRings': Descriptors.NumAromaticRings(mol),
        'TPSA': Descriptors.TPSA(mol)
    }

def train_model(df, all_therapeutic_areas):
    # Calculate molecular descriptors for each compound
    df['molecular_descriptors'] = df['smiles'].apply(calculate_molecular_descriptors)

    # Extract features from the molecular_descriptors column
    features = pd.DataFrame(df['molecular_descriptors'].tolist(), index=df.index)

    # Prepare target variables (therapeutic areas)
    target = df['therapeutic_areas'].apply(lambda x: [1 if area in x else 0 for area in all_therapeutic_areas])
    target = pd.DataFrame(target.tolist(), columns=all_therapeutic_areas, index=df.index)

    # Split the data into training and testing sets
    X_train, X_test, y_train, y_test = train_test_split(features, target, test_size=0.2, random_state=42)

    # Scale the features
    scaler = StandardScaler()
    X_train_scaled = scaler.fit_transform(X_train)
    X_test_scaled = scaler.transform(X_test)

    # Initialize and train the model
    model = MultiOutputClassifier(RandomForestClassifier(n_estimators=100, random_state=42))
    model.fit(X_train_scaled, y_train)

    # Make predictions on the test set
    y_pred = model.predict(X_test_scaled)

    # Evaluate the model
    print(classification_report(y_test, y_pred, target_names=all_therapeutic_areas))

    # Save the model and scaler for later use
    joblib.dump(model, 'therapeutic_area_predictor.joblib')
    joblib.dump(scaler, 'feature_scaler.joblib')

    return model, scaler

# This function should be called with your data to train the model
# Example: train_model(your_dataframe, your_list_of_therapeutic_areas)
