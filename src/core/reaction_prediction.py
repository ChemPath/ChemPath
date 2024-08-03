from sqlalchemy.orm import Session
from chempath_ml_models import create_reaction_prediction_model
from models import Reaction
from rdkit import Chem
from rdkit.Chem import AllChem

def predict_reaction(reactants: list, session: Session):
    input_shape = 2048  # Assuming Morgan fingerprint with 2048 bits

    model = create_reaction_prediction_model(input_shape=input_shape)

    # Implement reaction prediction logic
    combined_fp = get_combined_fingerprint(reactants)
    predicted_products = model.predict(combined_fp)

    # Convert predicted products to SMILES
    product_smiles = convert_prediction_to_smiles(predicted_products)

    predicted_reaction = Reaction(reactants=reactants, products=product_smiles)
    session.add(predicted_reaction)
    session.commit()

    return predicted_reaction

def get_combined_fingerprint(reactants):
    fps = []
    for smiles in reactants:
        mol = Chem.MolFromSmiles(smiles)
        fp = AllChem.GetMorganFingerprintAsBitVect(mol, 2, nBits=2048)
        fps.append(fp)
    return sum(fps)

def convert_prediction_to_smiles(predicted_products):
    # This function would depend on how your model outputs predictions
    # For this example, let's assume it returns SMILES strings directly
    return predicted_products

