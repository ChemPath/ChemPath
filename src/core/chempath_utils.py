from rdkit import Chem
from rdkit.Chem import Descriptors, AllChem
from rdkit.Chem import Crippen
from rdkit.Chem import rdMolDescriptors
import numpy as np
from rdkit import DataStructs


def validate_smiles(smiles):
    mol = Chem.MolFromSmiles(smiles)
    return mol is not None

def calculate_molecular_properties(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    
    properties = {
        'molecular_weight': Descriptors.ExactMolWt(mol),
        'logp': Crippen.MolLogP(mol),
        'h_bond_donors': rdMolDescriptors.CalcNumHBD(mol),
        'h_bond_acceptors': rdMolDescriptors.CalcNumHBA(mol),
        'rotatable_bonds': rdMolDescriptors.CalcNumRotatableBonds(mol),
        'polar_surface_area': Descriptors.TPSA(mol),
        'num_atoms': mol.GetNumAtoms(),
        'num_rings': rdMolDescriptors.CalcNumRings(mol),
        'num_aromatic_rings': rdMolDescriptors.CalcNumAromaticRings(mol)
    }
    return properties

def generate_morgan_fingerprint(smiles, radius=2, nBits=2048):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    return AllChem.GetMorganFingerprintAsBitVect(mol, radius, nBits=nBits)

def calculate_similarity(fp1, fp2):
    return DataStructs.TanimotoSimilarity(fp1, fp2)

def search_compounds(conn, query=None, filters=None, page=1, per_page=10):
    cursor = conn.cursor()
    base_query = "SELECT * FROM plant_compounds WHERE 1=1"
    params = []

    if query:
        base_query += " AND (name LIKE ? OR smiles LIKE ?)"
        params.extend([f"%{query}%", f"%{query}%"])

    if filters:
        for key, value in filters.items():
            if key in ['molecular_weight', 'logp', 'h_bond_donors', 'h_bond_acceptors', 'polar_surface_area', 'rotatable_bonds']:
                base_query += f" AND {key} BETWEEN ? AND ?"
                params.extend(value)
            elif key in ['plant_source', 'biological_activity', 'traditional_use']:
                base_query += f" AND {key} LIKE ?"
                params.append(f"%{value}%")

    base_query += " ORDER BY name"
    base_query += f" LIMIT {per_page} OFFSET {(page - 1) * per_page}"

    cursor.execute(base_query, params)
    compounds = cursor.fetchall()

    count_query = f"SELECT COUNT(*) FROM ({base_query})"
    cursor.execute(count_query, params)
    total_count = cursor.fetchone()[0]

    return compounds, total_count

def get_therapeutic_areas(conn):
    cursor = conn.cursor()
    cursor.execute("SELECT DISTINCT therapeutic_area FROM predicted_therapeutic_areas")
    return [row[0] for row in cursor.fetchall()]

def calculate_druglikeness_score(properties):
    mw = properties['molecular_weight']
    logp = properties['logp']
    hbd = properties['h_bond_donors']
    hba = properties['h_bond_acceptors']
    psa = properties['polar_surface_area']
    rotatable_bonds = properties['rotatable_bonds']

    score = 0
    if 150 <= mw <= 500:
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

    return score / 6  # Normalize to 0-1 range

def generate_scaffold(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    scaffold = AllChem.MurckoScaffold.GetScaffoldForMol(mol)
    return Chem.MolToSmiles(scaffold)

def calculate_molecular_complexity(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    return Descriptors.BertzCT(mol)

def calculate_natural_product_likeness(properties):
    # This is a simplified version and can be expanded with more sophisticated rules
    mw = properties['molecular_weight']
    logp = properties['logp']
    num_rings = properties['num_rings']
    num_aromatic_rings = properties['num_aromatic_rings']

    score = 0
    if 200 <= mw <= 800:
        score += 1
    if 0 <= logp <= 5:
        score += 1
    if 1 <= num_rings <= 6:
        score += 1
    if num_aromatic_rings <= 3:
        score += 1

    return score / 4  # Normalize to 0-1 range

