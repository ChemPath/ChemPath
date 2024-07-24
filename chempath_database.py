import sqlite3
import csv
from pathlib import Path
import re
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

def create_tables(conn):
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
                traditional_use TEXT,
                h_bond_donors INTEGER,
                h_bond_acceptors INTEGER,
                polar_surface_area REAL,
                rotatable_bonds INTEGER
            )
        ''')
        cursor.execute('''
            CREATE TABLE IF NOT EXISTS compound_classes (
                id INTEGER PRIMARY KEY,
                name TEXT NOT NULL,
                description TEXT
            )
        ''')
        cursor.execute('''
            CREATE TABLE IF NOT EXISTS therapeutic_areas (
                id INTEGER PRIMARY KEY,
                name TEXT NOT NULL,
                description TEXT
            )
        ''')
        cursor.execute('''
            CREATE TABLE IF NOT EXISTS compound_class_relationships (
                compound_id INTEGER,
                compound_class_id INTEGER,
                PRIMARY KEY (compound_id, compound_class_id),
                FOREIGN KEY (compound_id) REFERENCES plant_compounds (id),
                FOREIGN KEY (compound_class_id) REFERENCES compound_classes (id)
            )
        ''')
        cursor.execute('''
            CREATE TABLE IF NOT EXISTS therapeutic_area_relationships (
                compound_id INTEGER,
                therapeutic_area_id INTEGER,
                PRIMARY KEY (compound_id, therapeutic_area_id),
                FOREIGN KEY (compound_id) REFERENCES plant_compounds (id),
                FOREIGN KEY (therapeutic_area_id) REFERENCES therapeutic_areas (id)
            )
        ''')
        print("All tables created successfully")
    except sqlite3.Error as e:
        print(e)

def alter_plant_compounds_table(conn):
    cursor = conn.cursor()
    try:
        cursor.execute("ALTER TABLE plant_compounds ADD COLUMN h_bond_donors INTEGER")
        cursor.execute("ALTER TABLE plant_compounds ADD COLUMN h_bond_acceptors INTEGER")
        cursor.execute("ALTER TABLE plant_compounds ADD COLUMN polar_surface_area REAL")
        cursor.execute("ALTER TABLE plant_compounds ADD COLUMN rotatable_bonds INTEGER")
        conn.commit()
        print("Table plant_compounds altered successfully")
    except sqlite3.Error as e:
        print(f"Error altering table: {e}")

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

def insert_compound(conn, compound):
    sql = '''INSERT INTO plant_compounds(name, smiles, molecular_weight, logp, plant_source, biological_activity, traditional_use, h_bond_donors, h_bond_acceptors, polar_surface_area, rotatable_bonds)
             VALUES(?,?,?,?,?,?,?,?,?,?,?)'''
    try:
        cursor = conn.cursor()
        descriptors = calculate_descriptors(compound[1])  # SMILES is at index 1
        if descriptors is None:
            print(f"Error: Invalid SMILES string for compound {compound[0]}")
            return None
        compound_with_descriptors = compound + tuple(descriptors.values())
        cursor.execute(sql, compound_with_descriptors)
        conn.commit()
        return cursor.lastrowid
    except sqlite3.Error as e:
        print(e)
        return None

def insert_compound_class(conn, class_name, description):
    sql = '''INSERT INTO compound_classes(name, description) VALUES(?,?)'''
    try:
        cursor = conn.cursor()
        cursor.execute(sql, (class_name, description))
        conn.commit()
        return cursor.lastrowid
    except sqlite3.Error as e:
        print(e)
        return None

def insert_therapeutic_area(conn, area_name, description):
    sql = '''INSERT INTO therapeutic_areas(name, description) VALUES(?,?)'''
    try:
        cursor = conn.cursor()
        cursor.execute(sql, (area_name, description))
        conn.commit()
        return cursor.lastrowid
    except sqlite3.Error as e:
        print(e)
        return None

def create_compound_class_relationship(conn, compound_id, class_id):
    sql = '''INSERT INTO compound_class_relationships(compound_id, compound_class_id) VALUES(?,?)'''
    try:
        cursor = conn.cursor()
        cursor.execute(sql, (compound_id, class_id))
        conn.commit()
    except sqlite3.Error as e:
        print(e)

def create_therapeutic_area_relationship(conn, compound_id, area_id):
    sql = '''INSERT INTO therapeutic_area_relationships(compound_id, therapeutic_area_id) VALUES(?,?)'''
    try:
        cursor = conn.cursor()
        cursor.execute(sql, (compound_id, area_id))
        conn.commit()
    except sqlite3.Error as e:
        print(e)

def load_sample_data(conn):
    classes = [
        ("Alkaloids", "Nitrogen-containing organic compounds"),
        ("Flavonoids", "Polyphenolic compounds"),
        ("Terpenes", "Hydrocarbons derived from isoprene units"),
        ("Glycosides", "Compounds containing a sugar molecule"),
        ("Phenolic acids", "Aromatic secondary plant metabolites")
    ]
    for class_name, description in classes:
        insert_compound_class(conn, class_name, description)
    
    areas = [
        ("Antioxidant", "Prevents or slows damage to cells"),
        ("Anti-inflammatory", "Reduces inflammation"),
        ("Anticancer", "Prevents or treats cancer"),
        ("Antimalarial", "Prevents or treats malaria"),
        ("Cardiovascular health", "Promotes heart health")
    ]
    for area_name, description in areas:
        insert_therapeutic_area(conn, area_name, description)
    
    sample_data = [
        ("Quercetin", "O=C1C(O)=C(O)C(=O)C2=C1C=C(O)C(O)=C2O", 302.24, 1.54, "Various fruits and vegetables", "Antioxidant, anti-inflammatory", "Traditional Chinese Medicine for cardiovascular health"),
        ("Curcumin", "COC1=CC(=CC(=C1O)OC)C=CC(=O)CC(=O)C=CC2=CC(=C(C=C2)O)OC", 368.38, 3.29, "Turmeric (Curcuma longa)", "Anti-inflammatory, antioxidant", "Ayurvedic medicine for various ailments"),
    ]
    
    for compound in sample_data:
        compound_id = insert_compound(conn, compound)
        create_compound_class_relationship(conn, compound_id, 2)  # Assuming Flavonoids is id 2
        create_therapeutic_area_relationship(conn, compound_id, 1)  # Assuming Antioxidant is id 1
        create_therapeutic_area_relationship(conn, compound_id, 2)  # Assuming Anti-inflammatory is id 2
    
    print(f"Inserted {len(sample_data)} sample compounds with relationships")

def get_all_compounds(conn):
    cursor = conn.cursor()
    cursor.execute("SELECT * FROM plant_compounds")
    return cursor.fetchall()

def get_therapeutic_areas(conn, compound_id):
    cursor = conn.cursor()
    cursor.execute("""
        SELECT ta.name, ta.description
        FROM therapeutic_areas ta
        JOIN therapeutic_area_relationships tar ON ta.id = tar.therapeutic_area_id
        WHERE tar.compound_id = ?
    """, (compound_id,))
    return cursor.fetchall()

def validate_compound_data(name, smiles, molecular_weight, logp, plant_source, biological_activity, traditional_use):
    if not name or not isinstance(name, str):
        raise ValueError("Name must be a non-empty string")
    if not smiles or not isinstance(smiles, str):
        raise ValueError("SMILES must be a non-empty string")
    if not isinstance(molecular_weight, (int, float)) or molecular_weight <= 0:
        raise ValueError("Molecular weight must be a positive number")
    if not isinstance(logp, (int, float)):
        raise ValueError("LogP must be a number")
    if not plant_source or not isinstance(plant_source, str):
        raise ValueError("Plant source must be a non-empty string")
    if not biological_activity or not isinstance(biological_activity, str):
        raise ValueError("Biological activity must be a non-empty string")
    if not traditional_use or not isinstance(traditional_use, str):
        raise ValueError("Traditional use must be a non-empty string")

def add_new_compound_interactive(conn):
    print("Enter new compound details:")
    name = input("Name: ")
    smiles = input("SMILES: ")
    molecular_weight = float(input("Molecular Weight: "))
    logp = float(input("LogP: "))
    plant_source = input("Plant Source: ")
    biological_activity = input("Biological Activity: ")
    traditional_use = input("Traditional Use: ")
    
    try:
        validate_compound_data(name, smiles, molecular_weight, logp, plant_source, biological_activity, traditional_use)
        compound = (name, smiles, molecular_weight, logp, plant_source, biological_activity, traditional_use)
        compound_id = insert_compound(conn, compound)
        print(f"Compound {name} added successfully with ID {compound_id}")
    except ValueError as e:
        print(f"Error: {e}")

def add_new_class_interactive(conn):
    print("Enter new compound class details:")
    name = input("Name: ")
    description = input("Description: ")
    
    if not name or not description:
        print("Error: Name and description cannot be empty")
    else:
        class_id = insert_compound_class(conn, name, description)
        print(f"Compound class {name} added successfully with ID {class_id}")

def add_new_therapeutic_area_interactive(conn):
    print("Enter new therapeutic area details:")
    name = input("Name: ")
    description = input("Description: ")
    
    if not name or not description:
        print("Error: Name and description cannot be empty")
    else:
        area_id = insert_therapeutic_area(conn, name, description)
        print(f"Therapeutic area {name} added successfully with ID {area_id}")

def display_all_compounds(conn):
    cursor = conn.cursor()
    cursor.execute("SELECT * FROM plant_compounds")
    rows = cursor.fetchall()
    for row in rows:
        print(f"\nCompound: {row[1]}")
        print(f"SMILES: {row[2]}")
        print(f"Molecular Weight: {row[3]}")
        print(f"LogP: {row[4]}")
        print(f"Plant Source: {row[5]}")
        print(f"Biological Activity: {row[6]}")
        print(f"Traditional Use: {row[7]}")
        print(f"H-Bond Donors: {row[8]}")
        print(f"H-Bond Acceptors: {row[9]}")
        print(f"Polar Surface Area: {row[10]}")
        print(f"Rotatable Bonds: {row[11]}")
        print("--------------------")

def search_compounds(conn, search_type, search_term):
    cursor = conn.cursor()
    if search_type == 'name':
        cursor.execute("SELECT * FROM plant_compounds WHERE name LIKE ?", ('%' + search_term + '%',))
    elif search_type == 'property':
        cursor.execute(f"SELECT * FROM plant_compounds WHERE {search_term[0]} {search_term[1]} ?", (search_term[2],))
    elif search_type == 'therapeutic_area':
        cursor.execute("""
            SELECT pc.* FROM plant_compounds pc
            JOIN therapeutic_area_relationships tar ON pc.id = tar.compound_id
            JOIN therapeutic_areas ta ON tar.therapeutic_area_id = ta.id
            WHERE ta.name LIKE ?
        """, ('%' + search_term + '%',))
    else:
        print("Invalid search type")
        return []
    
    return cursor.fetchall()

def display_search_results(results):
    if not results:
        print("No compounds found matching the search criteria.")
    else:
        for row in results:
            print(f"\nCompound: {row[1]}")
            print(f"SMILES: {row[2]}")
            print(f"Molecular Weight: {row[3]}")
            print(f"LogP: {row[4]}")
            print(f"Plant Source: {row[5]}")
            print(f"Biological Activity: {row[6]}")
            print(f"Traditional Use: {row[7]}")
            print(f"H-Bond Donors: {row[8]}")
            print(f"H-Bond Acceptors: {row[9]}")
            print(f"Polar Surface Area: {row[10]}")
            print(f"Rotatable Bonds: {row[11]}")
            print("--------------------")

def search_compounds_interactive(conn):
    print("\nSearch Compounds:")
    print("1. Search by name")
    print("2. Search by property")
    print("3. Search by therapeutic area")
    
    choice = input("Enter your choice (1-3): ")
    
    if choice == '1':
        search_term = input("Enter compound name to search: ")
        results = search_compounds(conn, 'name', search_term)
    elif choice == '2':
        property_name = input("Enter property name (molecular_weight, logp, h_bond_donors, h_bond_acceptors, polar_surface_area, rotatable_bonds): ")
        operator = input("Enter operator (>, <, =, >=, <=): ")
        value = input("Enter value: ")
        results = search_compounds(conn, 'property', (property_name, operator, value))
    elif choice == '3':
        search_term = input("Enter therapeutic area to search: ")
        results = search_compounds(conn, 'therapeutic_area', search_term)
    else:
        print("Invalid choice")
        return
    
    display_search_results(results)

def main():
    database = Path("chempath_database.db")
    conn = create_connection(database)

    if conn is not None:
        create_tables(conn)
        alter_plant_compounds_table(conn)
        load_sample_data(conn)

        while True:
            print("\nChemPath Menu:")
            print("1. Display all compounds")
            print("2. Add new compound")
            print("3. Add new compound class")
            print("4. Add new therapeutic area")
            print("5. Search compounds")
            print("6. Exit")
            
            choice = input("Enter your choice (1-6): ")
            
            if choice == '1':
                display_all_compounds(conn)
            elif choice == '2':
                add_new_compound_interactive(conn)
            elif choice == '3':
                add_new_class_interactive(conn)
            elif choice == '4':
                add_new_therapeutic_area_interactive(conn)
            elif choice == '5':
                search_compounds_interactive(conn)
            elif choice == '6':
                break
            else:
                print("Invalid choice. Please try again.")

        conn.close()
    else:
        print("Error! Cannot create the database connection.")

if __name__ == '__main__':
    main()

