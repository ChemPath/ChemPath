# File: chempath_database.py

import sqlite3
import os

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
    """Create the plant_compounds table with updated structure"""
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
                therapeutic_areas TEXT,
                biological_activities TEXT,
                traditional_use TEXT
            )
        ''')
        print("Table 'plant_compounds' created successfully")
    except sqlite3.Error as e:
        print(e)

def update_database_schema(conn):
    """Update the existing database schema to include new columns"""
    try:
        cursor = conn.cursor()
        # Check if the columns exist, if not, add them
        cursor.execute("PRAGMA table_info(plant_compounds)")
        columns = [column[1] for column in cursor.fetchall()]
        
        if "therapeutic_areas" not in columns:
            cursor.execute("ALTER TABLE plant_compounds ADD COLUMN therapeutic_areas TEXT")
        if "biological_activities" not in columns:
            cursor.execute("ALTER TABLE plant_compounds ADD COLUMN biological_activities TEXT")
        
        print("Database schema updated successfully")
    except sqlite3.Error as e:
        print(f"Error updating database schema: {e}")

def insert_compound(conn, compound):
    """Insert a new compound into the plant_compounds table"""
    sql = '''INSERT INTO plant_compounds(name, smiles, molecular_weight, logp, plant_source, therapeutic_areas, biological_activities, traditional_use)
             VALUES(?,?,?,?,?,?,?,?)'''
    try:
        cursor = conn.cursor()
        cursor.execute(sql, compound)
        conn.commit()
        return cursor.lastrowid
    except sqlite3.Error as e:
        print(e)
        return None

# ... (keep any other existing functions)

def main():
    database = os.path.join(os.getcwd(), "chempath_database.db")
    conn = create_connection(database)

    if conn is not None:
        create_table(conn)
        update_database_schema(conn)
        conn.close()
    else:
        print("Error! Cannot create the database connection.")

if __name__ == '__main__':
    main()