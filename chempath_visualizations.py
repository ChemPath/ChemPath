# chempath_visualizations.py

import matplotlib.pyplot as plt
from rdkit import Chem
from rdkit.Chem import Draw

def draw_molecule(smiles):
    mol = Chem.MolFromSmiles(smiles)
    return Draw.MolToImage(mol)

def draw_optimization_comparison(original_smiles, optimized_smiles):
    mol1 = Chem.MolFromSmiles(original_smiles)
    mol2 = Chem.MolFromSmiles(optimized_smiles)
    
    img = Draw.MolsToGridImage([mol1, mol2], molsPerRow=2, subImgSize=(300, 300), legends=["Original", "Optimized"])
    return img

def draw_retrosynthesis_pathway(target_smiles, steps):
    mols = [Chem.MolFromSmiles(target_smiles)]
    for step in steps:
        mols.append(Chem.MolFromSmiles(step['product']))
    
    img = Draw.MolsToGridImage(mols, molsPerRow=3, subImgSize=(200, 200), legends=["Target"] + [f"Step {i+1}" for i in range(len(steps))])
    return img
