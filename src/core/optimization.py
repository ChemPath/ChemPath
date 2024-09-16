from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors, rdMolDescriptors
from rdkit.Chem.Crippen import MolLogP
from rdkit.Chem.rdMolDescriptors import GetMorganFingerprintAsBitVect
from database_operations import store_optimization_result, get_compound_by_smiles
from src.plant_retrosynthesis import PlantRetrosynthesis, suggest_precursors
import numpy as np

def optimize_compound(conn, smiles, iterations=100):
    compound = get_compound_by_smiles(conn, smiles)
    if not compound:
        return None

    mol = Chem.MolFromSmiles(smiles)
    best_mol = mol
    best_score = calculate_compound_score(mol)
    retro = PlantRetrosynthesis()

    for _ in range(iterations):
        new_mol = apply_random_modification(mol)
        new_score = calculate_compound_score(new_mol)
        
        # Consider retrosynthetic feasibility
        pathways = retro.analyze(Chem.MolToSmiles(new_mol))
        if len(pathways) > 0:
            new_score *= (1 + 0.1 * len(pathways))  # Boost score for synthetically feasible compounds
        
        if new_score > best_score:
            best_mol = new_mol
            best_score = new_score

    optimization_data = get_optimization_data(best_mol)
    store_optimization_result(conn, compound[0], optimization_data)
    return optimization_data


def calculate_compound_score(mol, target_property):
    if target_property == 'logP':
        return Descriptors.MolLogP(mol)
    elif target_property == 'molecular_weight':
        return Descriptors.ExactMolWt(mol)
    elif target_property == 'polararea':
        return Descriptors.TPSA(mol)
    elif target_property == 'complexity':
        return rdMolDescriptors.CalcCrippenDescriptors(mol)[1]
    elif target_property == 'hbonddonor':
        return Descriptors.NumHDonors(mol)
    elif target_property == 'hbondacc':
        return Descriptors.NumHAcceptors(mol)
    elif target_property == 'rotbonds':
        return Descriptors.NumRotatableBonds(mol)
    else:
        return 0
def apply_random_modification(mol):
    modifications = [add_atom, remove_atom, change_bond_order, add_ring]
    modification = np.random.choice(modifications)
    return modification(mol)

def add_atom(mol):
    editable_mol = Chem.RWMol(mol)
    atom_types = ['C', 'N', 'O']
    new_atom = Chem.Atom(np.random.choice(atom_types))
    atom_idx = editable_mol.AddAtom(new_atom)
    if atom_idx > 0:
        editable_mol.AddBond(atom_idx - 1, atom_idx, Chem.BondType.SINGLE)
    return editable_mol.GetMol()

def remove_atom(mol):
    editable_mol = Chem.RWMol(mol)
    if editable_mol.GetNumAtoms() > 1:
        atom_idx = np.random.randint(editable_mol.GetNumAtoms())
        editable_mol.RemoveAtom(atom_idx)
    return editable_mol.GetMol()

def change_bond_order(mol):
    editable_mol = Chem.RWMol(mol)
    if editable_mol.GetNumBonds() > 0:
        bond_idx = np.random.randint(editable_mol.GetNumBonds())
        bond = editable_mol.GetBondWithIdx(bond_idx)
        new_order = Chem.BondType.SINGLE if bond.GetBondType() != Chem.BondType.SINGLE else Chem.BondType.DOUBLE
        editable_mol.ReplaceBond(bond_idx, Chem.Bond(new_order))
    return editable_mol.GetMol()

def add_ring(mol):
    AllChem.EmbedMolecule(mol)
    ring_size = np.random.randint(3, 7)
    conf = mol.GetConformer()
    editable_mol = Chem.RWMol(mol)
    start_atom = np.random.randint(mol.GetNumAtoms())
    for i in range(ring_size):
        new_atom = Chem.Atom('C')
        new_idx = editable_mol.AddAtom(new_atom)
        if i == 0:
            editable_mol.AddBond(start_atom, new_idx, Chem.BondType.SINGLE)
        else:
            editable_mol.AddBond(new_idx - 1, new_idx, Chem.BondType.SINGLE)
        if i == ring_size - 1:
            editable_mol.AddBond(new_idx, start_atom, Chem.BondType.SINGLE)
    return editable_mol.GetMol()   

def apply_random_modification(mol):
    # Placeholder for applying random modifications
    # This should be replaced with actual modification logic
    return mol

def get_optimization_data(mol):
    return {
        'optimized_smiles': Chem.MolToSmiles(mol),
        'energy': AllChem.UFFGetMoleculeForceField(mol).CalcEnergy(),
        'molecular_weight': Descriptors.ExactMolWt(mol),
        'logp': MolLogP(mol),
        'h_bond_donors': Descriptors.NumHDonors(mol),
        'h_bond_acceptors': Descriptors.NumHAcceptors(mol),
        'rotatable_bonds': Descriptors.NumRotatableBonds(mol),
        'polar_surface_area': Descriptors.TPSA(mol),
        'natural_product_likeness': calculate_np_likeness(mol),
        'fraction_sp3': Descriptors.FractionCSP3(mol),
        'aromatic_rings': rdMolDescriptors.CalcNumAromaticRings(mol),
        'stereo_centers': rdMolDescriptors.CalcNumAtomStereoCenters(mol),
        'complexity': Descriptors.BertzCT(mol),
        'flexibility': calculate_flexibility(mol),
        'biodegradability': estimate_biodegradability(mol),
    }

def calculate_np_likeness(mol):
    morgan_fp = GetMorganFingerprintAsBitVect(mol, 2, nBits=1024)
    return sum(morgan_fp) / 1024

def calculate_flexibility(mol):
    rotatable_bonds = Descriptors.NumRotatableBonds(mol)
    total_bonds = mol.GetNumBonds()
    return rotatable_bonds / total_bonds if total_bonds > 0 else 0

def estimate_biodegradability(mol):
    logp = MolLogP(mol)
    mw = Descriptors.ExactMolWt(mol)
    return 1 / (1 + abs(logp)) * (300 / (300 + mw))

def normalize_logp(logp):
    return 1 / (1 + np.exp(-0.5 * (logp - 2)))

def normalize_tpsa(tpsa):
    return 1 / (1 + np.exp(-0.05 * (tpsa - 100)))

def normalize_complexity(complexity):
    return 1 / (1 + np.exp(-0.001 * (complexity - 500)))

def normalize_aromatic_rings(rings):
    return 1 / (1 + np.exp(-(rings - 2)))

def normalize_stereo_centers(centers):
    return 1 / (1 + np.exp(-(centers - 3)))

def retrosynthesis_informed_optimization(smiles, retrosynthesis_data):
    mol = Chem.MolFromSmiles(smiles)
    
    for step in retrosynthesis_data:
        if 'add_group' in step:
            mol = Chem.RWMol(mol)
            mol.AddAtom(Chem.Atom(6))
            mol.AddBond(step['position'], mol.GetNumAtoms()-1, Chem.BondType.SINGLE)
    
    AllChem.EmbedMolecule(mol, randomSeed=42)
    AllChem.UFFOptimizeMolecule(mol, maxIters=500)
    
    return Chem.MolToSmiles(mol)

def evaluate_optimization(original_smiles, optimized_smiles):
    original_mol = Chem.MolFromSmiles(original_smiles)
    optimized_mol = Chem.MolFromSmiles(optimized_smiles)
    
    evaluation = {
        'original_score': calculate_compound_score(original_mol),
        'optimized_score': calculate_compound_score(optimized_mol),
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
