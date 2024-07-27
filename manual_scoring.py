# manual_scoring.py

from rdkit import Chem
from rdkit.Chem import Descriptors

def calculate_manual_score(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None

    mw = Descriptors.ExactMolWt(mol)
    logp = Descriptors.MolLogP(mol)
    hbd = Descriptors.NumHDonors(mol)
    hba = Descriptors.NumHAcceptors(mol)

    # Simple scoring function based on Lipinski's Rule of Five
    score = 0
    if 200 <= mw <= 500:
        score += 1
    if logp <= 5:
        score += 1
    if hbd <= 5:
        score += 1
    if hba <= 10:
        score += 1

    return score / 4  # Normalize to 0-1 range
