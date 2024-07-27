# File: chempath_ml_models.py

import tensorflow as tf
from rdkit import Chem
from rdkit.Chem import AllChem
import numpy as np
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import StandardScaler

def smiles_to_fingerprint(smiles, size=2048):
    mol = Chem.MolFromSmiles(smiles)
    return list(AllChem.GetMorganFingerprintAsBitVect(mol, 2, nBits=size))

def create_retrosynthesis_model(input_shape):
    model = tf.keras.Sequential([
        tf.keras.layers.Dense(1024, activation='relu', input_shape=(input_shape,)),
        tf.keras.layers.Dense(512, activation='relu'),
        tf.keras.layers.Dense(256, activation='relu'),
        tf.keras.layers.Dense(128, activation='relu'),
        tf.keras.layers.Dense(64, activation='relu'),
        tf.keras.layers.Dense(1, activation='sigmoid')
    ])
    model.compile(optimizer='adam', loss='binary_crossentropy', metrics=['accuracy'])
    return model

def create_reaction_prediction_model(input_shape):
    model = tf.keras.Sequential([
        tf.keras.layers.Dense(1024, activation='relu', input_shape=(input_shape,)),
        tf.keras.layers.Dense(512, activation='relu'),
        tf.keras.layers.Dense(256, activation='relu'),
        tf.keras.layers.Dense(128, activation='relu'),
        tf.keras.layers.Dense(64, activation='relu'),
        tf.keras.layers.Dense(32, activation='softmax')  # Assuming 32 reaction classes
    ])
    model.compile(optimizer='adam', loss='categorical_crossentropy', metrics=['accuracy'])
    return model

def train_retrosynthesis_model(X, y):
    X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=42)
    
    scaler = StandardScaler()
    X_train_scaled = scaler.fit_transform(X_train)
    X_test_scaled = scaler.transform(X_test)
    
    model = create_retrosynthesis_model(X_train.shape[1])
    model.fit(X_train_scaled, y_train, epochs=50, batch_size=32, validation_split=0.2)
    
    loss, accuracy = model.evaluate(X_test_scaled, y_test)
    print(f"Retrosynthesis Model - Test Loss: {loss}, Test Accuracy: {accuracy}")
    
    return model, scaler

def train_reaction_prediction_model(X, y):
    X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=42)
    
    scaler = StandardScaler()
    X_train_scaled = scaler.fit_transform(X_train)
    X_test_scaled = scaler.transform(X_test)
    
    model = create_reaction_prediction_model(X_train.shape[1])
    model.fit(X_train_scaled, y_train, epochs=50, batch_size=32, validation_split=0.2)
    
    loss, accuracy = model.evaluate(X_test_scaled, y_test)
    print(f"Reaction Prediction Model - Test Loss: {loss}, Test Accuracy: {accuracy}")
    
    return model, scaler

def predict_retrosynthesis(model, scaler, smiles):
    fingerprint = smiles_to_fingerprint(smiles)
    scaled_fingerprint = scaler.transform([fingerprint])
    prediction = model.predict(scaled_fingerprint)[0][0]
    return prediction

def predict_reaction(model, scaler, smiles):
    fingerprint = smiles_to_fingerprint(smiles)
    scaled_fingerprint = scaler.transform([fingerprint])
    prediction = model.predict(scaled_fingerprint)[0]
    return prediction
def prepare_data(conn):
    try:
        query = "SELECT * FROM plant_compounds"
        df = pd.read_sql_query(query, conn)
        print(f"Retrieved {len(df)} rows from plant_compounds")
        
        features = df[['molecular_weight', 'logp', 'h_bond_donors', 'h_bond_acceptors', 'polar_surface_area', 'rotatable_bonds']]
        target = df['biological_activity']
        
        return features, target
    except pd.io.sql.DatabaseError as e:
        print(f"Error reading from database: {e}")
        return None, None

def train_model(conn):
    features, target = prepare_data(conn)
    if features is None or target is None:
        print("Failed to prepare data. Aborting model training.")
        return None, None
    
    X_train, X_test, y_train, y_test = train_test_split(features, target, test_size=0.2, random_state=42)
    
    scaler = StandardScaler()
    X_train_scaled = scaler.fit_transform(X_train)
    X_test_scaled = scaler.transform(X_test)
    
    model = RandomForestClassifier(n_estimators=100, random_state=42)
    model.fit(X_train_scaled, y_train)
    
    y_pred = model.predict(X_test_scaled)
    accuracy = accuracy_score(y_test, y_pred)
    print(f"Model accuracy: {accuracy:.2f}")
    
    return model, scaler
