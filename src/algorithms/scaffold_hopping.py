import rdkit
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Descriptors
from rdkit.Chem import Crippen
from rdkit.Chem.Scaffolds import MurckoScaffold
from rdkit.Chem.Pharm3D import Pharmacophore
from rdkit.Chem import ChemicalFeatures
from rdkit import RDConfig
import pandas as pd
import numpy as np
from sklearn.ensemble import RandomForestRegressor
import json
import os

def load_novel_scaffolds():
    with open('src/data/novel_scaffolds.json', 'r') as f:
        return json.load(f)

novel_scaffolds = {
    'benzene': 'c1ccccc1',
    'pyridine': 'c1ccncc1',
    'piperidine': 'C1CCNCC1',
    'cyclohexane': 'C1CCCCC1',
    'furan': 'c1ccoc1',
    'thiophene': 'c1ccsc1',
    'imidazole': 'c1cnc[nH]1',
    'oxazole': 'c1cocn1',
    'pyrazole': 'c1cn[nH]c1',
    'morpholine': 'C1COCCN1',
    'naphthalene': 'c1ccc2ccccc2c1',
    'quinoline': 'c1ccc2ncccc2c1',
    'indole': 'c1ccc2[nH]ccc2c1',
    'biphenyl': 'c1ccc(-c2ccccc2)cc1',
    'adamantane': 'C1C2CC3CC1CC(C2)C3'
}

def identify_advanced_pharmacophore(mol):
    mol_3d = AllChem.AddHs(mol)
    AllChem.EmbedMolecule(mol_3d)
    AllChem.UFFOptimizeMolecule(mol_3d)
    
    fdefName = os.path.join(RDConfig.RDDataDir, 'BaseFeatures.fdef')
    factory = ChemicalFeatures.BuildFeatureFactory(fdefName)
    feats = factory.GetFeaturesForMol(mol_3d)
    
    features = {
        'aromatic': mol_3d.GetSubstructMatches(Chem.MolFromSmarts('a')),
        'hydrogen_donor': [f.GetAtomIds() for f in feats if f.GetFamily() == 'Donor'],
        'hydrogen_acceptor': [f.GetAtomIds() for f in feats if f.GetFamily() == 'Acceptor'],
        'positive': [f.GetAtomIds() for f in feats if f.GetFamily() == 'PosIonizable'],
        'negative': [f.GetAtomIds() for f in feats if f.GetFamily() == 'NegIonizable'],
        'hydrophobic': [f.GetAtomIds() for f in feats if f.GetFamily() == 'Hydrophobe']
    }
    return features

def reconnect_side_chains(scaffold, side_chains):
    rxn = AllChem.ReactionFromSmarts('[*:1][H].[*:2][H]>>[*:1]-[*:2]')
    products = rxn.RunReactants((scaffold, side_chains))
    if products:
        return products[0][0]
    return None

def perform_scaffold_hop(mol, new_scaffold_smiles):
    original_scaffold = MurckoScaffold.GetScaffoldForMol(mol)
    new_scaffold = Chem.MolFromSmiles(new_scaffold_smiles)
    side_chains = AllChem.ReplaceCore(mol, original_scaffold)
    
    if side_chains is None:
        return None
    
    new_mol = reconnect_side_chains(new_scaffold, side_chains)
    
    if new_mol is None:
        return None
    
    try:
        Chem.SanitizeMol(new_mol)
        return new_mol
    except:
        return None


def generate_scaffold_analogs(smiles):
    mol = Chem.MolFromSmiles(smiles)
    analogs = []
    
    for scaffold_name, scaffold_smiles in novel_scaffolds.items():
        new_mol = perform_scaffold_hop(mol, scaffold_smiles)
        if new_mol:
            try:
                new_smiles = Chem.MolToSmiles(new_mol)
                analogs.append(new_smiles)
            except:
                continue
    
    return analogs

def calculate_molecular_descriptors(smiles):
    mol = Chem.MolFromSmiles(smiles)
    return {
        'MolWt': Descriptors.ExactMolWt(mol),
        'LogP': Crippen.MolLogP(mol),
        'NumHDonors': Descriptors.NumHDonors(mol),
        'NumHAcceptors': Descriptors.NumHAcceptors(mol),
        'NumRotatableBonds': Descriptors.NumRotatableBonds(mol)
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

def optimize_compound_scaffold(smiles, activity_model):
    analogs = generate_scaffold_analogs(smiles)
    
    results = []
    for analog in analogs:
        predicted_activity = predict_activity(activity_model, analog)
        results.append({
            'SMILES': analog,
            'Predicted Activity': predicted_activity
        })
    
    results_df = pd.DataFrame(results)
    
    if not results_df.empty and 'Predicted Activity' in results_df.columns:
        results_df = results_df.sort_values('Predicted Activity', ascending=False)
    
    return results_df

if __name__ == "__main__":
    training_data = pd.DataFrame({
        'SMILES': ['CC(=O)OC1=CC=CC=C1C(=O)O', 'CC1=CC=C(C=C1)NC(=O)C2=CC=CC=C2O', 'CC1=C(C=CC=C1)NC(=O)C2=CC=CC=C2O'],
        'Activity': [0.5, 0.7, 0.6]
    })
    activity_model = train_activity_model(training_data)
    
    compound_smiles = compound_smiles = "CN1C=NC2=C1C(=O)N(C)C(=O)N2C"  # Caffeine
    optimized_analogs = optimize_compound_scaffold(compound_smiles, activity_model)
    print(optimized_analogs)
