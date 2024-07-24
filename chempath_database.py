# File: chempath_database.py

import sqlite3
import csv
from pathlib import Path

def create_connection(db_file):
    """Create a database connection to a SQLite database"""
    conn = None
    try:
        conn = sqlite3.connect(db_file)
        print(f"Connected to SQLite version: {sqlite3.version}")
        return conn
    except sqlite3.Error as e:
        print(e)
    return conn

def create_table(conn):
    """Create the plant_compounds table"""
    try:
        cursor = conn.cursor()
        cursor.execute('''
            CREATE TABLE IF NOT EXISTS plant_compounds (
                id INTEGER PRIMARY KEY,
                name TEXT NOT NULL,
                smiles TEXT,
                molecular_weight REAL,
                logp REAL,
                plant_source TEXT,
                biological_activity TEXT,
                traditional_use TEXT
            )
        ''')
        print("Table 'plant_compounds' created successfully")
    except sqlite3.Error as e:
        print(e)

def insert_compound(conn, compound):
    """Insert a new compound into the plant_compounds table"""
    sql = '''INSERT INTO plant_compounds(name, smiles, molecular_weight, logp, plant_source, biological_activity, traditional_use)
             VALUES(?,?,?,?,?,?,?)'''
    try:
        cursor = conn.cursor()
        cursor.execute(sql, compound)
        conn.commit()
        return cursor.lastrowid
    except sqlite3.Error as e:
        print(e)
        return None

def load_sample_data(conn):
    """Load sample data"""
    sample_data = [
        ("Quercetin", "O=C1C(O)=C(O)C(=O)C2=C1C=C(O)C(O)=C2O", 302.24, 1.54, "Various fruits and vegetables", "Antioxidant, anti-inflammatory", "Traditional Chinese Medicine for cardiovascular health"),
        ("Curcumin", "COC1=CC(=CC(=C1O)OC)C=CC(=O)CC(=O)C=CC2=CC(=C(C=C2)O)OC", 368.38, 3.29, "Turmeric (Curcuma longa)", "Anti-inflammatory, antioxidant", "Ayurvedic medicine for various ailments"),
    ]
    
    for compound in sample_data:
        insert_compound(conn, compound)
    print(f"Inserted {len(sample_data)} sample compounds")

def main():
    database = Path("chempath_database.db")
    conn = create_connection(database)

    if conn is not None:
        create_table(conn)
        load_sample_data(conn)
        conn.close()
    else:
        print("Error! Cannot create the database connection.")

if __name__ == '__main__':
    main()

# File: requirements.txt

sqlite3