from rdkit import Chem
from sqlalchemy.orm import Session
from models import Molecule

def process_molecule(smiles: str, session: Session):
    mol = Chem.MolFromSmiles(smiles)
    if mol:
        molecule = Molecule(smiles=smiles)
        session.add(molecule)
        session.commit()
        return molecule
    return None
