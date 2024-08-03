from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors, Crippen
from rdkit.Chem import DataStructs
from rdkit.Chem.AllChem import GetMorganFingerprintAsBitVect
import random
from sqlalchemy.orm import Session
from models import Compound, StructuralOptimizationResult
from database import engine

def chemical_space_exploration(smiles, num_iterations=10):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return []

    explored_molecules = []
    for _ in range(num_iterations):
        new_mol = Chem.RWMol(mol)
        
        modification = random.choice(['add_atom', 'remove_atom', 'change_bond'])
        
        if modification == 'add_atom':
            atom = random.choice(['C', 'N', 'O'])
            new_mol.AddAtom(Chem.Atom(atom))
            idx = new_mol.GetNumAtoms() - 1
            random_atom = random.randint(0, idx - 1)
            new_mol.AddBond(random_atom, idx, Chem.BondType.SINGLE)
        
        elif modification == 'remove_atom':
            if new_mol.GetNumAtoms() > 1:
                idx = random.randint(0, new_mol.GetNumAtoms() - 1)
                new_mol.RemoveAtom(idx)
        
        elif modification == 'change_bond':
            if new_mol.GetNumBonds() > 0:
                bond_idx = random.randint(0, new_mol.GetNumBonds() - 1)
                bond = new_mol.GetBondWithIdx(bond_idx)
                new_order = random.choice([Chem.BondType.SINGLE, Chem.BondType.DOUBLE, Chem.BondType.TRIPLE])
                bond.SetBondType(new_order)
        
        new_smiles = Chem.MolToSmiles(new_mol)
        if new_smiles not in explored_molecules:
            explored_molecules.append(new_smiles)
    
    return explored_molecules

def functional_group_substitution(mol, target_group, replacement_group):
    target_pattern = Chem.MolFromSmarts(target_group)
    replacement_mol = Chem.MolFromSmiles(replacement_group)
    new_mol = AllChem.ReplaceSubstructs(mol, target_pattern, replacement_mol)
    return Chem.MolToSmiles(new_mol[0]) if new_mol else None

def ring_system_alteration(mol, alteration_type):
    rings = mol.GetRingInfo().AtomRings()
    if not rings:
        return None
    
    rw_mol = Chem.RWMol(mol)
    ring = random.choice(rings)
    
    if alteration_type == 'expand':
        if all(rw_mol.GetAtomWithIdx(idx).GetTotalValence() < 4 for idx in ring):
            atom_idx = ring[0]
            new_atom_idx = rw_mol.AddAtom(Chem.Atom(6))  # Add a carbon atom
            rw_mol.AddBond(atom_idx, new_atom_idx, Chem.BondType.SINGLE)
            rw_mol.AddBond(new_atom_idx, ring[1], Chem.BondType.SINGLE)
    elif alteration_type == 'contract':
        if len(ring) > 3:
            for atom_idx in ring:
                rw_mol.GetAtomWithIdx(atom_idx).SetIsAromatic(False)
            rw_mol.RemoveBond(ring[0], ring[-1])
            rw_mol.RemoveAtom(ring[0])
    elif alteration_type == 'fuse':
        if len(rings) > 1:
            other_ring = random.choice([r for r in rings if r != ring])
            if rw_mol.GetAtomWithIdx(ring[0]).GetTotalValence() < 4 and rw_mol.GetAtomWithIdx(other_ring[0]).GetTotalValence() < 4:
                rw_mol.AddBond(ring[0], other_ring[0], Chem.BondType.SINGLE)
    
    Chem.SanitizeMol(rw_mol)
    return Chem.MolToSmiles(rw_mol)

def scaffold_hopping(mol, scaffold_library):
    fp_gen = AllChem.GetMorganGenerator(radius=2, fpSize=2048)
    fp = fp_gen.GetFingerprint(mol)
    similar_mols = get_similar_molecules(fp, scaffold_library)
    hopped_mols = []
    for similar_mol in similar_mols:
        combined_mol = combine_fragments(mol, similar_mol)
        if combined_mol and combined_mol.GetNumAtoms() > 0:
            hopped_mols.append(Chem.MolToSmiles(combined_mol))
    return hopped_mols

def get_similar_molecules(fp, scaffold_library):
    similar_mols = []
    for smiles in scaffold_library:
        mol = Chem.MolFromSmiles(smiles)
        if mol:
            mol_fp = AllChem.GetMorganFingerprintAsBitVect(mol, 2, nBits=2048)
            similarity = DataStructs.TanimotoSimilarity(fp, mol_fp)
            if similarity > 0.7:
                similar_mols.append(mol)
    return similar_mols[:5]

def combine_fragments(mol1, mol2):
    combined = Chem.CombineMols(mol1, mol2)
    return combined

def generate_analogs(compound_id, num_analogs=5):
    with Session(engine) as session:
        compound = session.query(Compound).get(compound_id)
        if not compound:
            return {"error": "Compound not found"}

        mol = Chem.MolFromSmiles(compound.smiles)
        analogs = []
        for _ in range(num_analogs):
            analog = Chem.Mol(mol)
            
            modification = random.choice(['functional_group', 'ring_system', 'scaffold'])
            
            if modification == 'functional_group':
                target_group = '[OH]'  # Example target group
                replacement_group = '[NH2]'  # Example replacement group
                analog = functional_group_substitution(analog, target_group, replacement_group)
            elif modification == 'ring_system':
                alteration_type = random.choice(['expand', 'contract', 'fuse'])
                analog = ring_system_alteration(analog, alteration_type)
            elif modification == 'scaffold':
                scaffold_library = ['c1ccccc1', 'C1CCCCC1', 'c1ccncc1']
                analog = scaffold_hopping(analog, scaffold_library)
            
            if analog:
                if isinstance(analog, str):
                    analogs.append(analog)
                else:
                    analogs.append(Chem.MolToSmiles(analog))
        
        result = StructuralOptimizationResult(
            compound_id=compound_id,
            analogs=",".join(analogs)
        )
        session.add(result)
        session.commit()

        return {"result_id": result.id, "analogs": analogs}

def get_optimization_result(result_id):
    with Session(engine) as session:
        result = session.query(StructuralOptimizationResult).get(result_id)
        if result:
            return {
                "compound_id": result.compound_id,
                "analogs": result.analogs.split(",")
            }
        return {"error": "Result not found"}

# Example usage
if __name__ == "__main__":
    result = generate_analogs(1, num_analogs=5)
    print(result)
