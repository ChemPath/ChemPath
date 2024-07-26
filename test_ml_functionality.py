from chempath_api import ChemPathAPI

def test_ml_functionality():
    api = ChemPathAPI('chempath.db')

    # Create ml_model_performance table if it doesn't exist
    cursor = api.conn.cursor()
    cursor.execute("""
        CREATE TABLE IF NOT EXISTS ml_model_performance (
            id INTEGER PRIMARY KEY,
            model_name TEXT,
            accuracy REAL,
            precision REAL,
            recall REAL,
            f1_score REAL,
            support INTEGER
        );
    """)
    api.conn.commit()
    
    # Test model training
    model, scaler = api.train_ml_model(api.conn)
    
    # Test prediction
    sample_smiles = "O=C1C(O)=C(O)C(=O)C2=C1C=C(O)C(O)=C2O"  # Quercetin
    predicted_areas = api.predict_therapeutic_areas(sample_smiles)
    print(f"Predicted therapeutic areas for {sample_smiles}: {predicted_areas}")
    
    # Test searching compounds
    results = api.search_compounds(query="Quercetin")
    print(f"Search results: {results}")
    
    # Verify ML model performance storage
    cursor = api.conn.cursor()
    cursor.execute("SELECT * FROM ml_model_performance")
    performance_data = cursor.fetchall()
    print(f"ML model performance data: {performance_data}")
    
    # Test adding a new compound and predicting its therapeutic areas
    new_compound = ("New Compound", "CC1=CC=C(C=C1)C2=CC(=O)C3=C(O)C=C(O)C=C3O2", 286.24, 3.2, "Test Plant", "Unknown", "Test Use")
    api.add_compound(new_compound)
    new_smiles = new_compound[1]
    new_predictions = api.predict_therapeutic_areas(new_smiles)
    print(f"Predicted therapeutic areas for new compound: {new_predictions}")
    
    api.conn.close()

if __name__ == "__main__":
    test_ml_functionality()
