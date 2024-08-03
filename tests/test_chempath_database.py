import sqlite3
from pathlib import Path
from rdkit import Chem
from rdkit.Chem import Descriptors

def create_connection(db_file):
    conn = sqlite3.connect(db_file)
    return conn

def create_tables(conn):
    cursor = conn.cursor()
    
    cursor.execute('''
    CREATE TABLE IF NOT EXISTS compounds (
        id INTEGER PRIMARY KEY,
        smiles TEXT NOT NULL,
        name TEXT,
        molecular_weight REAL,
        logp REAL,
        h_bond_donors INTEGER,
        h_bond_acceptors INTEGER,
        ph REAL,
        temperature REAL
    )
    ''')
    
    cursor.execute('''
    CREATE TABLE IF NOT EXISTS therapeutic_areas (
        id INTEGER PRIMARY KEY,
        name TEXT NOT NULL UNIQUE
    )
    ''')
    
    cursor.execute('''
    CREATE TABLE IF NOT EXISTS compound_therapeutic_areas (
        compound_id INTEGER,
        therapeutic_area_id INTEGER,
        FOREIGN KEY (compound_id) REFERENCES compounds (id),
        FOREIGN KEY (therapeutic_area_id) REFERENCES therapeutic_areas (id),
        PRIMARY KEY (compound_id, therapeutic_area_id)
    )
    ''')
    
    conn.commit()

def insert_sample_data(conn):
    cursor = conn.cursor()
    
    sample_compounds = [
        ("CC(=O)Oc1ccccc1C(=O)O", "Aspirin"),
        ("CC12CCC3C(C1CCC2O)CCC4=CC(=O)CCC34C", "Testosterone"),
        ("CN1CCC[C@H]1c2cccnc2", "Nicotine"),
        ("CC(C)(C)NCC(O)c1ccc(O)c(O)c1", "Salbutamol"),
    ]
    
    for smiles, name in sample_compounds:
        mol = Chem.MolFromSmiles(smiles)
        cursor.execute('''
        INSERT INTO compounds (smiles, name, molecular_weight, logp, h_bond_donors, h_bond_acceptors)
        VALUES (?, ?, ?, ?, ?, ?)
        ''', (smiles, name, Descriptors.ExactMolWt(mol), Descriptors.MolLogP(mol),
              Descriptors.NumHDonors(mol), Descriptors.NumHAcceptors(mol)))
    
    sample_areas = [
        "Cardiovascular",
        "Respiratory",
        "Neurological",
        "Oncology",
        "Immunology"
    ]
    
    for area in sample_areas:
        cursor.execute('INSERT INTO therapeutic_areas (name) VALUES (?)', (area,))
    
    conn.commit()

def main():
    db_file = Path("test_chempath.db")
    conn = create_connection(db_file)
    
    if conn is not None:
        create_tables(conn)
        insert_sample_data(conn)
        conn.close()
        print("Database updated successfully.")
    else:
        print("Error! Cannot create the database connection.")

if __name__ == '__main__':
    main()
