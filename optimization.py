from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Descriptors
from rdkit.Chem.Crippen import MolLogP
from database_operations import store_optimization_result, get_compound_by_smiles

def optimize_compound(conn, smiles):
    compound = get_compound_by_smiles(conn, smiles)
    if not compound:
        return None

    mol = Chem.MolFromSmiles(smiles)
    AllChem.EmbedMolecule(mol, randomSeed=42)
    AllChem.UFFOptimizeMolecule(mol, maxIters=200)

    optimization_data = {
        'optimized_smiles': Chem.MolToSmiles(mol),
        'energy': AllChem.UFFGetMoleculeForceField(mol).CalcEnergy(),
        'molecular_weight': Descriptors.ExactMolWt(mol),
        'logp': MolLogP(mol),
        'h_bond_donors': Descriptors.NumHDonors(mol),
        'h_bond_acceptors': Descriptors.NumHAcceptors(mol),
        'rotatable_bonds': Descriptors.NumRotatableBonds(mol),
        'polar_surface_area': Descriptors.TPSA(mol)
    }

    store_optimization_result(conn, compound[0], optimization_data)
    return optimization_data

def retrosynthesis_informed_optimization(smiles, retrosynthesis_data):
    mol = Chem.MolFromSmiles(smiles)
    
    # Apply retrosynthesis-informed modifications
    for step in retrosynthesis_data:
        # Example: Add functional groups based on retrosynthesis steps
        if 'add_group' in step:
            mol = Chem.RWMol(mol)
            mol.AddAtom(Chem.Atom(6))  # Add a carbon atom
            mol.AddBond(step['position'], mol.GetNumAtoms()-1, Chem.BondType.SINGLE)
    
    # Optimize the modified structure
    AllChem.EmbedMolecule(mol, randomSeed=42)
    AllChem.UFFOptimizeMolecule(mol, maxIters=500)
    
    return Chem.MolToSmiles(mol)

def evaluate_optimization(original_smiles, optimized_smiles):
    original_mol = Chem.MolFromSmiles(original_smiles)
    optimized_mol = Chem.MolFromSmiles(optimized_smiles)
    
    evaluation = {
        'original_mw': Descriptors.ExactMolWt(original_mol),
        'optimized_mw': Descriptors.ExactMolWt(optimized_mol),
        'original_logp': MolLogP(original_mol),
        'optimized_logp': MolLogP(optimized_mol),
        'original_tpsa': Descriptors.TPSA(original_mol),
        'optimized_tpsa': Descriptors.TPSA(optimized_mol)
    }
    
    return evaluation

if __name__ == "__main__":
    from src.database.chempath_database import create_connection
    
    conn = create_connection("chempath_database.db")
    test_smiles = "CC(=O)OC1=CC=CC=C1C(=O)O"  # Aspirin
    
    optimized_data = optimize_compound(conn, test_smiles)
    if optimized_data:
        print("Optimization results:")
        for key, value in optimized_data.items():
            print(f"{key}: {value}")
        
        evaluation = evaluate_optimization(test_smiles, optimized_data['optimized_smiles'])
        print("\nOptimization evaluation:")
        for key, value in evaluation.items():
            print(f"{key}: {value}")
    else:
        print("Optimization failed or compound not found in database.")
    
    conn.close()
