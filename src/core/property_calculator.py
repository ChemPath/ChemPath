from rdkit import Chem
from rdkit.Chem import Descriptors
from sqlalchemy.orm import Session
from models import Molecule

def calculate_properties(molecule_id: int, session: Session):
    molecule = session.query(Molecule).get(molecule_id)
    mol = Chem.MolFromSmiles(molecule.smiles)
    properties = {
        'molecular_weight': Descriptors.MolWt(mol),
        'logp': Descriptors.MolLogP(mol),
        'h_bond_donors': Descriptors.NumHDonors(mol),
        'h_bond_acceptors': Descriptors.NumHAcceptors(mol),
    }
    for key, value in properties.items():
        setattr(molecule, key, value)
    session.commit()
    return properties
