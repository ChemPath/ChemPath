from rdkit import Chem
from rdkit.Chem import Descriptors, Crippen

def calculate_molecular_descriptors(mol):
    return {
        'MolWt': Descriptors.ExactMolWt(mol),
        'LogP': Crippen.MolLogP(mol),
        'NumHDonors': Descriptors.NumHDonors(mol),
        'NumHAcceptors': Descriptors.NumHAcceptors(mol),
        'NumRotatableBonds': Descriptors.NumRotatableBonds(mol),
        'TPSA': Descriptors.TPSA(mol),
        'NumAromaticRings': Descriptors.NumAromaticRings(mol),
        'NumAliphaticRings': Descriptors.NumAliphaticRings(mol),
        'FractionCSP3': Descriptors.FractionCSP3(mol),
        'NumHeteroatoms': Descriptors.NumHeteroatoms(mol),
        'NumRadicalElectrons': Descriptors.NumRadicalElectrons(mol),
        'NumValenceElectrons': Descriptors.NumValenceElectrons(mol)
    }
