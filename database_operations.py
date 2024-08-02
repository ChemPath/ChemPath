import sqlite
import json

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

def get_compound_by_smiles(conn, smiles):
    cursor = conn.cursor()
    cursor.execute("SELECT * FROM plant_compounds WHERE smiles = ?", (smiles,))
    return cursor.fetchone()

def store_optimization_result(conn, compound_id, optimization_data):
    cursor = conn.cursor()
    cursor.execute("""
        INSERT INTO optimization_results (compound_id, optimization_data)
        VALUES (?, ?)
    """, (compound_id, json.dumps(optimization_data)))
    conn.commit()

def store_retrosynthesis_result(conn, compound_id, retrosynthesis_data):
    cursor = conn.cursor()
    cursor.execute("""
        CREATE TABLE IF NOT EXISTS retrosynthesis_results (
            id INTEGER PRIMARY KEY,
            compound_id INTEGER,
            retrosynthesis_data TEXT,
            FOREIGN KEY (compound_id) REFERENCES plant_compounds (id)
        );
    """)
    cursor.execute("""
        INSERT INTO retrosynthesis_results (compound_id, retrosynthesis_data)
        VALUES (?, ?);
    """, (compound_id, json.dumps(retrosynthesis_data)))
    conn.commit()
def get_retrosynthesis_data(conn, compound_id):
    cursor = conn.cursor()
    cursor.execute("SELECT retrosynthesis_data FROM retrosynthesis_results WHERE compound_id = ?", (compound_id,))
    result = cursor.fetchone()
    return json.loads(result[0]) if result else None

def store_retrosynthesis_informed_optimization(conn, compound_id, optimized_smiles):
    cursor = conn.cursor()
    cursor.execute("""
        CREATE TABLE IF NOT EXISTS retrosynthesis_informed_optimizations (
            id INTEGER PRIMARY KEY,
            compound_id INTEGER,
            optimized_smiles TEXT,
            FOREIGN KEY (compound_id) REFERENCES plant_compounds (id)
        );
    """)
    cursor.execute("""
        INSERT INTO retrosynthesis_informed_optimizations (compound_id, optimized_smiles)
        VALUES (?, ?);
    """, (compound_id, optimized_smiles))
    conn.commit()

