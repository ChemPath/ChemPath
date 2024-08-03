from rdkit import Chem
from rdkit.Chem import Descriptors, Crippen, rdMolDescriptors

def calculate_manual_score(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None

    mw = Descriptors.ExactMolWt(mol)
    logp = Crippen.MolLogP(mol)
    hbd = rdMolDescriptors.CalcNumHBD(mol)
    hba = rdMolDescriptors.CalcNumHBA(mol)
    psa = Descriptors.TPSA(mol)
    rotatable_bonds = rdMolDescriptors.CalcNumRotatableBonds(mol)

    # Enhanced scoring function based on Lipinski's Rule of Five and additional criteria
    score = 0
    if 200 <= mw <= 500:
        score += 1
    if -0.4 <= logp <= 5.6:
        score += 1
    if hbd <= 5:
        score += 1
    if hba <= 10:
        score += 1
    if 20 <= psa <= 130:
        score += 1
    if rotatable_bonds <= 10:
        score += 1

    # Additional criteria
    num_rings = rdMolDescriptors.CalcNumRings(mol)
    if 1 <= num_rings <= 4:
        score += 1

    num_aromatic_rings = rdMolDescriptors.CalcNumAromaticRings(mol)
    if num_aromatic_rings <= 3:
        score += 1

    return score / 8  # Normalize to 0-1 range

def get_property_contributions(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None

    properties = {
        'Molecular Weight': Descriptors.ExactMolWt(mol),
        'LogP': Crippen.MolLogP(mol),
        'H-Bond Donors': rdMolDescriptors.CalcNumHBD(mol),
        'H-Bond Acceptors': rdMolDescriptors.CalcNumHBA(mol),
        'Polar Surface Area': Descriptors.TPSA(mol),
        'Rotatable Bonds': rdMolDescriptors.CalcNumRotatableBonds(mol),
        'Number of Rings': rdMolDescriptors.CalcNumRings(mol),
        'Aromatic Rings': rdMolDescriptors.CalcNumAromaticRings(mol)
    }

    return properties

def explain_score(smiles):
    score = calculate_manual_score(smiles)
    properties = get_property_contributions(smiles)

    explanation = f"The compound's manual score is {score:.2f}\n\n"
    explanation += "Property breakdown:\n"
    for prop, value in properties.items():
        explanation += f"{prop}: {value:.2f}\n"

    return explanation

if __name__ == "__main__":
    test_smiles = "CC(=O)OC1=CC=CC=C1C(=O)O"  # Aspirin
    score = calculate_manual_score(test_smiles)
    print(f"Manual score for Aspirin: {score:.2f}")
    print(explain_score(test_smiles))
