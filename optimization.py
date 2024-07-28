# optimization.py


from rdkit import Chem
from rdkit.Chem import AllChem

from database_operations import store_optimization_result
def optimize_compound(conn, smiles):
    from database_operations import get_compound_by_smiles
    
    compound = get_compound_by_smiles(conn, smiles)
    if not compound:
        return None

    mol = Chem.MolFromSmiles(smiles)
    AllChem.EmbedMolecule(mol)
    AllChem.UFFOptimizeMolecule(mol)

    optimization_data = {
        'optimized_smiles': Chem.MolToSmiles(mol),
        'energy': AllChem.UFFGetMoleculeForceField(mol).CalcEnergy()
    }

    store_optimization_result(conn, compound[0], optimization_data)
    return optimization_data
