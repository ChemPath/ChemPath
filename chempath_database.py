# File: chempath_database_operations.py

import sqlite3

def create_connection(db_file):
    """Create a database connection to the SQLite database specified by db_file"""
    conn = None
    try:
        conn = sqlite3.connect(db_file)
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

def get_compound_by_name(conn, name):
    """Retrieve a compound by its name"""
    cursor = conn.cursor()
    cursor.execute("SELECT * FROM plant_compounds WHERE name=?", (name,))
    return cursor.fetchone()

def update_compound(conn, compound):
    """Update an existing compound in the database"""
    sql = '''UPDATE plant_compounds
             SET smiles = ?, molecular_weight = ?, logp = ?, plant_source = ?, biological_activity = ?, traditional_use = ?
             WHERE name = ?'''
    try:
        cursor = conn.cursor()
        cursor.execute(sql, compound[1:] + (compound[0],))
        conn.commit()
    except sqlite3.Error as e:
        print(e)

def delete_compound(conn, name):
    """Delete a compound from the database by its name"""
    sql = 'DELETE FROM plant_compounds WHERE name=?'
    try:
        cursor = conn.cursor()
        cursor.execute(sql, (name,))
        conn.commit()
    except sqlite3.Error as e:
        print(e)

def get_all_compounds(conn):
    """Retrieve all compounds from the database"""
    cursor = conn.cursor()
    cursor.execute("SELECT * FROM plant_compounds")
    return cursor.fetchall()

# Example usage
if __name__ == '__main__':
    conn = create_connection("chempath_database.db")
    if conn is not None:
        create_table(conn)
        
        # Insert a compound
        new_compound = ("Quercetin", "O=C1C(O)=C(O)C(=O)C2=C1C=C(O)C(O)=C2O", 302.24, 1.54, "Various fruits and vegetables", "Antioxidant, anti-inflammatory", "Traditional Chinese Medicine for cardiovascular health")
        insert_compound(conn, new_compound)
        
        # Retrieve and print all compounds
        compounds = get_all_compounds(conn)
        for compound in compounds:
            print(compound)
        
        conn.close()
    else:
        print("Error! Cannot create the database connection.")