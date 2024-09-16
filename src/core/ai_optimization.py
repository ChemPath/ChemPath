import numpy as np
from sklearn.ensemble import RandomForestRegressor
from sqlalchemy.orm import Session
from models import Compound, OptimizationResult
from database import engine

def train_optimization_model(session):
    compounds = session.query(Compound).all()
    X = np.array([[c.molecular_weight, c.logp, c.h_bond_donors, c.h_bond_acceptors] for c in compounds])
    y = np.array([c.optimization_score for c in compounds])
    model = RandomForestRegressor(n_estimators=100, random_state=42)
    model.fit(X, y)
    return model

def predict_optimization_priority(model, compound):
    features = np.array([[compound.molecular_weight, compound.logp, compound.h_bond_donors, compound.h_bond_acceptors]])
    prediction = model.predict(features)
    return prediction[0]

def prioritize_optimization_strategies(compounds, model):
    for compound in compounds:
        priority = predict_optimization_priority(model, compound)
        compound.optimization_priority = priority

    return sorted(compounds, key=lambda x: x.optimization_priority, reverse=True)

def optimize_compounds(compound_ids):
    with Session(engine) as session:
        model = train_optimization_model(session)
        compounds = session.query(Compound).filter(Compound.id.in_(compound_ids)).all()
        
        optimized_compounds = prioritize_optimization_strategies(compounds, model)
        
        results = []
        for compound in optimized_compounds:
            result = OptimizationResult(
                compound_id=compound.id,
                optimization_priority=compound.optimization_priority
            )
            session.add(result)
            results.append(result)
        
        session.commit()
        
        return [{"compound_id": r.compound_id, "optimization_priority": r.optimization_priority} for r in results]

def get_optimization_results(result_ids):
    with Session(engine) as session:
        results = session.query(OptimizationResult).filter(OptimizationResult.id.in_(result_ids)).all()
        return [{"id": r.id, "compound_id": r.compound_id, "optimization_priority": r.optimization_priority} for r in results]

# Example usage
if __name__ == "__main__":
    compound_ids = [1, 2, 3, 4, 5]  # Example compound IDs
    optimization_results = optimize_compounds(compound_ids)
    print("Optimization Results:")
    for result in optimization_results:
        print(f"Compound ID: {result['compound_id']}, Priority: {result['optimization_priority']:.2f}")
