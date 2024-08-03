import sqlite3
from src.database.chempath_database import create_connection

def recreate_table(conn):
    cursor = conn.cursor()

    # Drop the existing table
    cursor.execute("DROP TABLE IF EXISTS plant_compounds")

    # Create the new table with the updated structure
    cursor.execute("""
    CREATE TABLE plant_compounds (
        id INTEGER PRIMARY KEY AUTOINCREMENT,
        name TEXT NOT NULL,
        smiles TEXT NOT NULL,
        molecular_weight REAL,
        logp REAL,
        plant_source TEXT,
        therapeutic_areas TEXT,
        biological_activities TEXT,
        traditional_use TEXT,
        h_bond_donors INTEGER,
        h_bond_acceptors INTEGER,
        polar_surface_area REAL,
        rotatable_bonds INTEGER,
        pubchem_cid INTEGER,
        inchi TEXT,
        iupac_name TEXT,
        synonyms TEXT,
        complexity REAL,
        heavy_atom_count INTEGER,
        exact_mass REAL,
        monoisotopic_mass REAL
    )
    """)

    conn.commit()
    print("Table 'plant_compounds' has been recreated with the updated structure.")

def main():
    conn = create_connection("chempath_database.db")
    if conn is not None:
        recreate_table(conn)
        conn.close()
    else:
        print("Error! Cannot create the database connection.")

if __name__ == "__main__":
    main()
