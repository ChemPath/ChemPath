import sys
import os
sys.path.append(os.path.dirname(os.path.abspath(__file__)))
import rdkit
from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors, Crippen
import pandas as pd
import numpy as np
from sklearn.ensemble import RandomForestRegressor
import json
from .molecular_descriptors import calculate_molecular_descriptors

# Step 1: Load bioisosteric groups from a JSON file

def load_bioisosteric_groups():
    current_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    file_path = os.path.join(current_dir, 'data', 'bioisosteric_groups.json')
    with open(file_path, 'r') as f:
        return json.load(f)

bioisosteric_groups = load_bioisosteric_groups()
# Step 2: Develop a database of bioisosteric functional groups
bioisosteric_groups = {
    '[CX3](=O)[OX1H1]': ['[SX4](=O)(=O)[OX1H1]', '[PX4](=O)([OX1H1])[OX1H1]', '[CX3](=O)[NX3H1][OX1H1]', '[#6]1[#7][#7][#7][#7]1'],
    '[OX2H1]': ['[NX3H2]', '[F]', '[SX2H1]', '[CH3]'],
    '[NX3H2]': ['[OX2H1]', '[F]', '[CH3]', '[CX4F3]'],
    '[CX3]=O': ['[CX3]=S', '[CX3]=[NX2R]', '[SX4](=O)(=O)', '[CX4F2]'],
    '[CH2]': ['[OX2]', '[SX2]', '[NX3H1]', '[CX4F2]']
}

# Step 3: Create a function to identify functional groups in a molecule
def identify_functional_groups(mol):
    functional_groups = []
    for smart in bioisosteric_groups.keys():
        pattern = Chem.MolFromSmarts(smart)
        if pattern is not None and mol.HasSubstructMatch(pattern):
            functional_groups.append(smart)
    return functional_groups

# Step 4: Implement a function to replace functional groups
def replace_functional_group(mol, original_group, new_group, max_replacements=3):
    original_pattern = Chem.MolFromSmarts(original_group)
    new_pattern = Chem.MolFromSmarts(new_group)
    
    if original_pattern is None or new_pattern is None:
        return None
    
    # Find all matches of the original group
    matches = mol.GetSubstructMatches(original_pattern)
    
    if not matches:
        return None
    
    # Replace up to max_replacements occurrences of the original group with the new group
    new_mol = AllChem.ReplaceSubstructs(mol, original_pattern, new_pattern, replaceAll=False, 
                                        replacementConnectionPoint=0, useChirality=False)
    
    return new_mol[0] if new_mol else None


def score_replacement(original_mol, new_mol):
    original_descriptors = calculate_molecular_descriptors(original_mol)
    new_descriptors = calculate_molecular_descriptors(new_mol)
    
    # Calculate a simple score based on the changes in key properties
    score = 0
    score += abs(new_descriptors['LogP'] - original_descriptors['LogP']) * 0.5
    score += abs(new_descriptors['MolWt'] - original_descriptors['MolWt']) * 0.01
    score += abs(new_descriptors['TPSA'] - original_descriptors['TPSA']) * 0.02
    
    return score

# Step 5: Create a function to generate analogs using functional group substitution
def generate_analogs(smiles):
    mol = Chem.MolFromSmiles(smiles)
    analogs = []
    
    functional_groups = identify_functional_groups(mol)
    
    for fg in functional_groups:
        for replacement in bioisosteric_groups[fg]:
            new_mol = replace_functional_group(mol, fg, replacement)
            if new_mol:
                new_smiles = Chem.MolToSmiles(new_mol)
                score = score_replacement(mol, new_mol)
                analogs.append((new_smiles, score))
    
    # Sort analogs by score in descending order
    analogs.sort(key=lambda x: x[1], reverse=True)
    
    return [analog[0] for analog in analogs]

# Step 6: Implement a simple machine learning model to predict activity
from rdkit.Chem import AllChem

def calculate_molecular_descriptors(mol):
    if isinstance(mol, str):
        mol = Chem.MolFromSmiles(mol)
    
    Chem.SanitizeMol(mol)  # Sanitize the molecule in place
    
    return {
        'MolWt': Descriptors.ExactMolWt(mol),
        'LogP': Crippen.MolLogP(mol),
        'NumHDonors': Descriptors.NumHDonors(mol),
        'NumHAcceptors': Descriptors.NumHAcceptors(mol),
        'NumRotatableBonds': Descriptors.NumRotatableBonds(mol),
        'TPSA': Descriptors.TPSA(mol),
        'NumAromaticRings': Descriptors.NumAromaticRings(mol)
    }


def train_activity_model(training_data):
    X = pd.DataFrame([calculate_molecular_descriptors(smiles) for smiles in training_data['SMILES']])
    y = training_data['Activity']
    
    model = RandomForestRegressor(n_estimators=100, random_state=42)
    model.fit(X, y)
    
    return model

def predict_activity(model, smiles):
    descriptors = calculate_molecular_descriptors(smiles)
    X = pd.DataFrame([descriptors])
    return model.predict(X)[0]

# Step 7: Put it all together in a main function
def optimize_compound(smiles, activity_model):
    analogs = generate_analogs(smiles)
    
    results = []
    for analog in analogs:
        predicted_activity = predict_activity(activity_model, analog)
        results.append({
            'SMILES': analog,
            'Predicted Activity': predicted_activity
        })
    
    results_df = pd.DataFrame(results)
    results_df = results_df.sort_values('Predicted Activity', ascending=False)
    
    return results_df

# Example usage
if __name__ == "__main__":
    # Assuming we have a trained activity model
    training_data = pd.DataFrame({
        'SMILES': ['CC(=O)OC1=CC=CC=C1C(=O)O', 'CC1=CC=C(C=C1)NC(=O)C2=CC=CC=C2O', 'CC1=C(C=CC=C1)NC(=O)C2=CC=CC=C2O'],
        'Activity': [0.5, 0.7, 0.6]
    })
    activity_model = train_activity_model(training_data)
    
    compound_smiles = "CC(=O)OC1=CC=CC=C1C(=O)O"  # Aspirin
    optimized_analogs = optimize_compound(compound_smiles, activity_model)
    print(optimized_analogs)
