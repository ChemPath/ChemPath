import sqlite3
from rdkit import Chem
from rdkit.Chem import Descriptors

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
import random

def fetch_random_compound_name():
    # This is a simplified version. In a real scenario, you might want to fetch this from a more comprehensive list or API
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

def predict_therapeutic_areas(smiles, all_therapeutic_areas):
    # This is a placeholder function. In a real-world scenario, you would implement
    # a machine learning model to predict therapeutic areas based on SMILES.
    # For now, we'll return a random selection of therapeutic areas.
    import random
    return random.sample(all_therapeutic_areas, min(3, len(all_therapeutic_areas)))

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
