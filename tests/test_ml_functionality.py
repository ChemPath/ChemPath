import unittest
from chempath_api import ChemPathAPI
from database import engine
from sqlalchemy.orm import Session
from models import Compound, TherapeuticArea, MLModelPerformance

class TestMLFunctionality(unittest.TestCase):
    def setUp(self):
        self.api = ChemPathAPI()
        self.populate_sample_data()

    def populate_sample_data(self):
        with Session(engine) as session:
            compounds = [
                Compound(smiles="CC(=O)Oc1ccccc1C(=O)O", name="Aspirin"),
                Compound(smiles="CC12CCC3C(C1CCC2O)CCC4=CC(=O)CCC34C", name="Testosterone"),
                Compound(smiles="CN1CCC[C@H]1c2cccnc2", name="Nicotine"),
                Compound(smiles="CC(C)(C)NCC(O)c1ccc(O)c(O)c1", name="Salbutamol")
            ]
            session.add_all(compounds)
            session.commit()

    def test_train_ml_model(self):
        model, scaler = self.api.train_ml_model()
        self.assertIsNotNone(model)
        self.assertIsNotNone(scaler)

    def test_predict_therapeutic_areas(self):
        sample_smiles = "O=C1C(O)=C(O)C(=O)C2=C1C=C(O)C(O)=C2O"  # Quercetin
        predicted_areas = self.api.predict_therapeutic_areas(sample_smiles)
        self.assertIsInstance(predicted_areas, list)
        self.assertTrue(len(predicted_areas) > 0)

    def test_search_compounds(self):
        results = self.api.search_compounds(query="Aspirin")
        self.assertIsInstance(results, list)
        self.assertTrue(len(results) > 0)

    def test_ml_model_performance_storage(self):
        self.api.train_ml_model()
        with Session(engine) as session:
            performance_data = session.query(MLModelPerformance).all()
            self.assertTrue(len(performance_data) > 0)

    def test_add_new_compound_and_predict(self):
        new_compound = {
            "name": "New Compound",
            "smiles": "CC1=CC=C(C=C1)C2=CC(=O)C3=C(O)C=C(O)C=C3O2",
            "molecular_weight": 286.24,
            "logp": 3.2,
            "h_bond_donors": 2,
            "h_bond_acceptors": 4,
            "ph": 7.0,
            "temperature": 25.0
        }
        compound_id = self.api.insert_compound(new_compound)
        self.assertIsNotNone(compound_id)

        new_predictions = self.api.predict_therapeutic_areas(new_compound['smiles'])
        self.assertIsInstance(new_predictions, list)
        self.assertTrue(len(new_predictions) > 0)

    def tearDown(self):
        with Session(engine) as session:
            session.query(Compound).delete()
            session.query(TherapeuticArea).delete()
            session.query(MLModelPerformance).delete()
            session.commit()

if __name__ == "__main__":
    unittest.main()
