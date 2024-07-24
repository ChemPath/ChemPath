import sqlite3
import csv
from pathlib import Path
import requests
from rdkit import Chem
from rdkit.Chem import Descriptors
import pandas as pd

def create_connection(db_file):
    conn = None
    try:
        conn = sqlite3.connect(db_file)
        print(f"Connected to SQLite version: {sqlite3.version}")
        return conn
    except sqlite3.Error as e:
        print(e)
    return conn

def get_all_compounds(conn):
    cursor = conn.cursor()
    cursor.execute("SELECT * FROM plant_compounds")
    compounds = cursor.fetchall()
    return compounds

def get_therapeutic_areas(conn):
    cursor = conn.cursor()
    cursor.execute("SELECT DISTINCT biological_activity FROM plant_compounds")
    areas = cursor.fetchall()
    return [area[0] for area in areas if area[0]]

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
    except sqlite3.Error as e:
        print(e)

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

def add_external_compound(conn):
    compound_name = input("Enter compound name to fetch from PubChem: ")
    compound_data = fetch_pubchem_data(compound_name)
    
    if compound_data:
        plant_source = input("Enter plant source: ")
        biological_activity = input("Enter biological activity: ")
        traditional_use = input("Enter traditional use: ")
        
        compound = (
            compound_data['name'],
            compound_data['smiles'],
            compound_data['molecular_weight'],
            compound_data['logp'],
            plant_source,
            biological_activity,
            traditional_use,
            compound_data['h_bond_donors'],
            compound_data['h_bond_acceptors'],
            compound_data['polar_surface_area'],
            compound_data['rotatable_bonds']
        )
        
        compound_id = insert_compound(conn, compound)
        print(f"Compound {compound_data['name']} added successfully with ID {compound_id}")
    else:
        print(f"Could not fetch data for compound {compound_name}")

def display_all_compounds(conn):
    compounds = get_all_compounds(conn)
    df = pd.DataFrame(compounds, columns=['id', 'name', 'smiles', 'molecular_weight', 'logp', 'plant_source', 'biological_activity', 'traditional_use', 'h_bond_donors', 'h_bond_acceptors', 'polar_surface_area', 'rotatable_bonds'])
    print(df)

def main():
    database = Path("chempath_database.db")
    conn = create_connection(database)

    if conn is not None:
        create_tables(conn)

        while True:
            print("\nChemPath Menu:")
            print("1. Display all compounds")
            print("2. Add compound from PubChem")
            print("3. Exit")
            
            choice = input("Enter your choice (1-3): ")
            
            if choice == '1':
                display_all_compounds(conn)
            elif choice == '2':
                add_external_compound(conn)
            elif choice == '3':
                break
            else:
                print("Invalid choice. Please try again.")

        conn.close()
    else:
        print("Error! Cannot create the database connection.")

if __name__ == '__main__':
    main()
