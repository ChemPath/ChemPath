from rdkit import Chem
from rdkit.Chem import AllChem, DataStructs
from sqlalchemy.orm import Session
from models import Molecule

def similarity_search(query_smiles: str, session: Session, threshold: float = 0.7):
    query_mol = Chem.MolFromSmiles(query_smiles)
    query_fp = AllChem.GetMorganFingerprintAsBitVect(query_mol, 2, nBits=1024)
    
    similar_molecules = []
    for molecule in session.query(Molecule).all():
        mol = Chem.MolFromSmiles(molecule.smiles)
        fp = AllChem.GetMorganFingerprintAsBitVect(mol, 2, nBits=1024)
        similarity = DataStructs.TanimotoSimilarity(query_fp, fp)
        if similarity >= threshold:
            similar_molecules.append((molecule, similarity))
    
    return sorted(similar_molecules, key=lambda x: x[1], reverse=True)
