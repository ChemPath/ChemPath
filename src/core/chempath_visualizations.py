# chempath_visualizations.py

import matplotlib.pyplot as plt
from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem import AllChem
import io
import base64

def draw_molecule(smiles):
    mol = Chem.MolFromSmiles(smiles)
    img = Draw.MolToImage(mol)
    return img_to_base64(img)

def draw_optimization_comparison(original_smiles, optimized_smiles):
    mol1 = Chem.MolFromSmiles(original_smiles)
    mol2 = Chem.MolFromSmiles(optimized_smiles)
    
    img = Draw.MolsToGridImage([mol1, mol2], molsPerRow=2, subImgSize=(300, 300), legends=["Original", "Optimized"])
    return img_to_base64(img)

def draw_retrosynthesis_pathway(target_smiles, steps):
    mols = [Chem.MolFromSmiles(target_smiles)]
    for step in steps:
        mols.append(Chem.MolFromSmiles(step['product']))
    
    img = Draw.MolsToGridImage(mols, molsPerRow=3, subImgSize=(200, 200), legends=["Target"] + [f"Step {i+1}" for i in range(len(steps))])
    return img_to_base64(img)

def draw_3d_conformation(smiles):
    mol = Chem.MolFromSmiles(smiles)
    mol = Chem.AddHs(mol)
    AllChem.EmbedMolecule(mol, randomSeed=42)
    AllChem.MMFFOptimizeMolecule(mol)
    img = Draw.MolToImage(mol, size=(300, 300))
    return img_to_base64(img)

def draw_similarity_matrix(similarity_matrix, compound_names):
    plt.figure(figsize=(10, 8))
    plt.imshow(similarity_matrix, cmap='viridis', aspect='auto')
    plt.colorbar()
    plt.xticks(range(len(compound_names)), compound_names, rotation=90)
    plt.yticks(range(len(compound_names)), compound_names)
    plt.title("Compound Similarity Matrix")
    plt.tight_layout()
    
    img_buffer = io.BytesIO()
    plt.savefig(img_buffer, format='png')
    img_buffer.seek(0)
    return base64.b64encode(img_buffer.getvalue()).decode()

def img_to_base64(img):
    img_buffer = io.BytesIO()
    img.save(img_buffer, format='PNG')
    img_buffer.seek(0)
    return base64.b64encode(img_buffer.getvalue()).decode()

def draw_property_distribution(properties, property_name):
    plt.figure(figsize=(8, 6))
    plt.hist(properties, bins=20, edgecolor='black')
    plt.title(f"Distribution of {property_name}")
    plt.xlabel(property_name)
    plt.ylabel("Frequency")
    plt.tight_layout()
    
    img_buffer = io.BytesIO()
    plt.savefig(img_buffer, format='png')
    img_buffer.seek(0)
    return base64.b64encode(img_buffer.getvalue()).decode()

