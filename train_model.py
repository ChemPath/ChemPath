import pandas as pd
import numpy as np
from sklearn.model_selection import train_test_split, cross_val_score
from sklearn.ensemble import RandomForestClassifier
from sklearn.multioutput import MultiOutputClassifier
from sklearn.preprocessing import StandardScaler, MultiLabelBinarizer
from sklearn.metrics import f1_score
import joblib
from sqlalchemy.orm import Session
from database import engine
from models import Compound, TherapeuticArea, MLModelPerformance

def train_model():
    with Session(engine) as session:
        compounds = session.query(Compound).all()
        all_areas = [area.name for area in session.query(TherapeuticArea).all()]

        data = []
        for compound in compounds:
            row = {
                'smiles': compound.smiles,
                'molecular_weight': compound.molecular_weight,
                'logp': compound.logp,
                'h_bond_donors': compound.h_bond_donors,
                'h_bond_acceptors': compound.h_bond_acceptors,
                'ph': compound.ph,
                'temperature': compound.temperature,
                'therapeutic_areas': [area.name for area in compound.therapeutic_areas]
            }
            data.append(row)

        df = pd.DataFrame(data)

        X = df[['molecular_weight', 'logp', 'h_bond_donors', 'h_bond_acceptors', 'ph', 'temperature']]
        X['mw_logp_ratio'] = X['molecular_weight'] / (X['logp'] + 1)
        
        mlb = MultiLabelBinarizer()
        y = mlb.fit_transform(df['therapeutic_areas'])

        X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=42)

        scaler = StandardScaler()
        X_train_scaled = scaler.fit_transform(X_train)
        X_test_scaled = scaler.transform(X_test)

        rf = RandomForestClassifier(n_estimators=200, max_depth=10, min_samples_split=5, random_state=42)
        model = MultiOutputClassifier(rf)

        cv_scores = cross_val_score(model, X_train_scaled, y_train, cv=5, scoring='f1_macro')
        print(f"Cross-validation F1 scores: {cv_scores}")
        print(f"Mean CV F1 score: {np.mean(cv_scores)}")

        model.fit(X_train_scaled, y_train)

        y_pred = model.predict(X_test_scaled)
        test_f1 = f1_score(y_test, y_pred, average='macro')
        print(f"Test set F1 score: {test_f1}")

        feature_importance = model.estimators_[0].feature_importances_
        feature_names = X.columns
        for name, importance in zip(feature_names, feature_importance):
            print(f"Feature '{name}': {importance}")

        joblib.dump(model, 'therapeutic_area_predictor.joblib')
        joblib.dump(scaler, 'feature_scaler.joblib')
        joblib.dump(mlb, 'multilabel_binarizer.joblib')

        performance = MLModelPerformance(
            model_name="RandomForestClassifier",
            accuracy=np.mean(cv_scores),
            f1_score=test_f1
        )
        session.add(performance)
        session.commit()

    print("Model trained and saved successfully.")

if __name__ == '__main__':
    train_model()
