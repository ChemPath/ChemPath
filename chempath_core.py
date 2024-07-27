import sqlite3
from rdkit import Chem
from rdkit.Chem import Descriptors
from rdkit.Chem import AllChem
import random
import logging
from structural_optimization import functional_group_substitution, ring_system_alteration, scaffold_hopping
import requests
from ml_model import predict_therapeutic_areas as ml_predict_therapeutic_areas
from ml_model import train_model


def setup_logging():
    logging.basicConfig(level=logging.DEBUG, format='%(asctime)s - %(levelname)s - %(message)s')
    return logging.getLogger(__name__)

def create_connection(db_file):
    conn = None
    try:
        conn = sqlite3.connect(db_file)
        print(f"Connected to SQLite version: {sqlite3.version}")
        return conn
    except sqlite3.Error as e:
        print(e)
    return conn

def calculate_descriptors(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    return {
        'h_bond_donors': Descriptors.NumHDonors(mol),
        'h_bond_acceptors': Descriptors.NumHAcceptors(mol),
        'polar_surface_area': Descriptors.TPSA(mol),
        'rotatable_bonds': Descriptors.NumRotatableBonds(mol)
    }

def fetch_random_compound_name():
    common_compounds = ["Aspirin", "Caffeine", "Ibuprofen", "Acetaminophen", "Penicillin", "Morphine", "Quinine", "Insulin", "Dopamine", "Serotonin"]
    return random.choice(common_compounds)

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

def insert_compound(conn, compound):
    sql = '''INSERT OR REPLACE INTO plant_compounds(name, smiles, molecular_weight, logp, plant_source, biological_activity, traditional_use, h_bond_donors, h_bond_acceptors, polar_surface_area, rotatable_bonds)
             VALUES(?,?,?,?,?,?,?,?,?,?,?)'''
    try:
        cursor = conn.cursor()
        cursor.execute(sql, compound)
        conn.commit()
        return cursor.lastrowid
    except sqlite3.Error as e:
        print(e)
        return None

def get_therapeutic_areas(conn):
    cursor = conn.cursor()
    cursor.execute("SELECT DISTINCT biological_activity FROM plant_compounds")
    areas = cursor.fetchall()
    return [area[0] for area in areas if area[0]]

def predict_therapeutic_areas(smiles, model, scaler):
    return ml_predict_therapeutic_areas(smiles, model, scaler)

def create_tables(conn):
    try:
        cursor = conn.cursor()
        cursor.execute('''
            CREATE TABLE IF NOT EXISTS plant_compounds (
                id INTEGER PRIMARY KEY,
                name TEXT NOT NULL UNIQUE,
                smiles TEXT,
                molecular_weight REAL,
                logp REAL,
                plant_source TEXT,
                biological_activity TEXT,
                traditional_use TEXT,
                h_bond_donors INTEGER,
                h_bond_acceptors INTEGER,
                polar_surface_area REAL,
                rotatable_bonds INTEGER
            )
        ''')
        print("Table 'plant_compounds' created successfully")
        
        cursor.execute('''
        CREATE TABLE IF NOT EXISTS predicted_therapeutic_areas (
            compound_id INTEGER PRIMARY KEY,
            smiles TEXT NOT NULL,
            predicted_areas TEXT NOT NULL,
            prediction_date TIMESTAMP DEFAULT CURRENT_TIMESTAMP,
            FOREIGN KEY (compound_id) REFERENCES plant_compounds(id)
        )
        ''')
        print("Table 'predicted_therapeutic_areas' created successfully")
    except sqlite3.Error as e:
        print(e)

def fetch_pubchem_data(compound_name):
    base_url = "https://pubchem.ncbi.nlm.nih.gov/rest/pug"
    search_url = f"{base_url}/compound/name/{compound_name}/JSON"
    
    response = requests.get(search_url)
    if response.status_code == 200:
        data = response.json()
        cid = data['PC_Compounds'][0]['id']['id']['cid']
        
        property_url = f"{base_url}/compound/cid/{cid}/property/MolecularFormula,MolecularWeight,XLogP,CanonicalSMILES/JSON"
        prop_response = requests.get(property_url)
        
        if prop_response.status_code == 200:
            prop_data = prop_response.json()['PropertyTable']['Properties'][0]
            smiles = prop_data['CanonicalSMILES']
            descriptors = calculate_descriptors(smiles)
            
            return {
                'name': compound_name,
                'smiles': smiles,
                'molecular_weight': prop_data['MolecularWeight'],
                'logp': prop_data['XLogP'],
                'h_bond_donors': descriptors['h_bond_donors'],
                'h_bond_acceptors': descriptors['h_bond_acceptors'],
                'polar_surface_area': descriptors['polar_surface_area'],
                'rotatable_bonds': descriptors['rotatable_bonds']
            }
    
    return None

def create_indexes(conn):
    cursor = conn.cursor()
    cursor.execute("CREATE INDEX IF NOT EXISTS idx_compound_name ON plant_compounds(name)")
    cursor.execute("CREATE INDEX IF NOT EXISTS idx_compound_smiles ON plant_compounds(smiles)")
    conn.commit()
    print("Indexes created successfully")

def optimize_structure(smiles, optimization_type, params):
    mol = Chem.MolFromSmiles(smiles)
    if optimization_type == 'functional_group':
        return functional_group_substitution(mol, params['target'], params['replacement'])
    elif optimization_type == 'ring_system':
        return ring_system_alteration(mol, params['alteration_type'])
    elif optimization_type == 'scaffold':
        return scaffold_hopping(mol, params['scaffold_library'])
    else:
        raise ValueError("Invalid optimization type")
    
def retrosynthesis_informed_optimization(smiles, retrosynthesis_data):
    if retrosynthesis_data is None:
        print("No retrosynthesis data available. Performing standard optimization.")
        # Implement a standard optimization method here
        return smiles  # Return original SMILES if no optimization is performed

    optimized_smiles = smiles
    for step in retrosynthesis_data:
        # Existing optimization logic
        pass

    return optimized_smiles



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

from pathlib import Path

def train_ml_model(conn):
    model, scaler = train_model(conn)
    return model, scaler
def main():
    database = Path("chempath_database.db")
    conn = create_connection(database)

    if conn is not None:
        create_tables(conn)
        create_indexes(conn)
        conn.close()
    else:
        print("Error! Cannot create the database connection.")

if __name__ == '__main__':
    main()
