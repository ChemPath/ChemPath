import pandas as pd
import numpy as np
from sklearn.model_selection import train_test_split, cross_val_score
from sklearn.ensemble import RandomForestClassifier
from sklearn.multioutput import MultiOutputClassifier
from sklearn.preprocessing import StandardScaler, MultiLabelBinarizer
from sklearn.metrics import f1_score
import joblib
from chempath_database import create_connection, get_all_compounds, get_therapeutic_areas

def train_model():
    conn = create_connection('chempath_database.db')
    compounds = get_all_compounds(conn)
    all_areas = get_therapeutic_areas(conn)

    df = pd.DataFrame(compounds, columns=['id', 'name', 'smiles', 'molecular_weight', 'logp', 'plant_source', 'biological_activity', 'traditional_use', 'h_bond_donors', 'h_bond_acceptors', 'polar_surface_area', 'rotatable_bonds'])

    # Enhanced feature engineering
    X = df[['molecular_weight', 'logp', 'h_bond_donors', 'h_bond_acceptors', 'polar_surface_area', 'rotatable_bonds']]
    X['mw_logp_ratio'] = X['molecular_weight'] / (X['logp'] + 1)  # Adding a new feature
    
    mlb = MultiLabelBinarizer()
    y = mlb.fit_transform(df['biological_activity'].apply(lambda x: x.split(',')))

    X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=42)

    scaler = StandardScaler()
    X_train_scaled = scaler.fit_transform(X_train)
    X_test_scaled = scaler.transform(X_test)

    # Fine-tuned Random Forest Classifier
    rf = RandomForestClassifier(n_estimators=200, max_depth=10, min_samples_split=5, random_state=42)
    model = MultiOutputClassifier(rf)

    # Cross-validation
    cv_scores = cross_val_score(model, X_train_scaled, y_train, cv=5, scoring='f1_macro')
    print(f"Cross-validation F1 scores: {cv_scores}")
    print(f"Mean CV F1 score: {np.mean(cv_scores)}")

    model.fit(X_train_scaled, y_train)

    # Evaluate on test set
    y_pred = model.predict(X_test_scaled)
    test_f1 = f1_score(y_test, y_pred, average='macro')
    print(f"Test set F1 score: {test_f1}")
    # Feature importance analysis
    feature_importance = model.estimators_[0].feature_importances_
    feature_names = X.columns
    for name, importance in zip(feature_names, feature_importance):
        print(f"Feature '{name}': {importance}")

    joblib.dump(model, 'therapeutic_area_predictor.joblib')
    joblib.dump(scaler, 'feature_scaler.joblib')
    joblib.dump(mlb, 'multilabel_binarizer.joblib')

    print("Model trained and saved successfully.")

if __name__ == '__main__':
    train_model()
