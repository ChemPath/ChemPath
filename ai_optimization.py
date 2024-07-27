import numpy as np
from sklearn.ensemble import RandomForestRegressor

def train_optimization_model(compounds_data):
    X = np.array([[c['molecular_weight'], c['logp'], c['h_bond_donors'], c['h_bond_acceptors']] for c in compounds_data])
    y = np.array([c['optimization_score'] for c in compounds_data])
    model = RandomForestRegressor(n_estimators=100, random_state=42)
    model.fit(X, y)
    return model

def predict_optimization_priority(model, compound):
    features = np.array([[compound['molecular_weight'], compound['logp'], compound['h_bond_donors'], compound['h_bond_acceptors']]])
    prediction = model.predict(features)
    print(f"Calculated priority score: {prediction[0]:.2f}")
    return prediction[0]


def prioritize_optimization_strategies(compounds, model):
    print("Compounds before sorting:")
    for i, compound in enumerate(compounds[:3], 1):
        print(f"{i}. {compound['name']}")

    for compound in compounds:
        priority = predict_optimization_priority(model, compound)
        compound['optimization_priority'] = priority

    sorted_compounds = sorted(compounds, key=lambda x: x['optimization_priority'], reverse=True)

    print("\nCompounds after sorting:")
    for i, compound in enumerate(sorted_compounds[:3], 1):
        print(f"{i}. {compound['name']} - Priority: {compound['optimization_priority']:.2f}")

    return sorted_compounds

