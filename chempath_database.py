import requests
from rdkit import Chem
from chempath_core import (
    setup_logging, fetch_random_compound_name, fetch_pubchem_data,
    insert_compound, create_tables, create_indexes, search_compounds,
    get_therapeutic_areas, optimize_structure, chemical_space_exploration
)
from advanced_retrosynthesis import advanced_retrosynthetic_analysis
from synthetic_feasibility import analyze_synthetic_feasibility
from reagent_availability import analyze_reagent_availability
from retrosynthesis import perform_retrosynthesis
from predict_therapeutic_areas import predict_therapeutic_areas
from database_operations import create_connection, get_compound_by_smiles, store_optimization_result
from database_operations import create_connection, get_compound_by_smiles
from optimization import optimize_compound
from retrosynthesis import retrosynthetic_analysis
from chempath_core import create_tables, create_indexes

print("Script started")


def validate_smiles(smiles):
    mol = Chem.MolFromSmiles(smiles)
    return mol is not None

def validate_smarts(smarts):
    pattern = Chem.MolFromSmarts(smarts)
    return pattern is not None
from rdkit.Chem import Descriptors
from chempath_api import ChemPathAPI
from chempath_core import fetch_pubchem_data
__all__ = ['create_connection', 'search_compounds', 'insert_compound', 'get_therapeutic_areas', 'predict_therapeutic_areas']
import random
from predict_therapeutic_areas import predict_therapeutic_areas
import sqlite3
import csv
from pathlib import Path
from chempath_api import ChemPathAPI
import requests
from rdkit import Chem
from rdkit.Chem import Descriptors
import pandas as pd



def get_all_compounds(conn):
    cursor = conn.cursor()
    cursor.execute("SELECT * FROM plant_compounds")
    compounds = cursor.fetchall()
    return compounds

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

    # Get total count for pagination
    count_query = f"SELECT COUNT(*) FROM ({base_query})"
    cursor.execute(count_query, params)
    total_count = cursor.fetchone()[0]

    return compounds, total_count
def display_compounds(compounds):
    if not compounds:
        print("No compounds found.")
        return
    
    headers = ["ID", "Name", "SMILES", "Mol Weight", "LogP", "Plant Source", "Bio Activity", "Trad Use"]
    row_format = "{:>5} {:<20} {:<20} {:>10} {:>8} {:<15} {:<15} {:<15}"
    
    print(row_format.format(*headers))
    print("-" * 120)
    
    for compound in compounds:
        print(row_format.format(
            compound[0],
            compound[1][:20],
            compound[2][:20],
            f"{compound[3]:.2f}",
            f"{compound[4]:.2f}",
            compound[5][:15],
            compound[6][:15],
            compound[7][:15]
        ))

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
                biological_activities TEXT,
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


        cursor.execute('''
            CREATE TABLE IF NOT EXISTS retrosynthesis_results (
                id INTEGER PRIMARY KEY,
                compound_id INTEGER,
                retrosynthesis_data TEXT,
                FOREIGN KEY (compound_id) REFERENCES plant_compounds (id)
            )
        ''')
        print("Table 'retrosynthesis_results' created successfully")

    except sqlite3.Error as e:
        print(f"Error creating table: {e}")

    conn.rollback()


def create_indexes(conn):
    cursor = conn.cursor()
    cursor.execute("CREATE INDEX IF NOT EXISTS idx_compound_name ON plant_compounds(name)")
    cursor.execute("CREATE INDEX IF NOT EXISTS idx_compound_smiles ON plant_compounds(smiles)")
    conn.commit()
    print("Indexes created successfully")

def add_predicted_therapeutic_areas_table(conn):
    try:
        cursor = conn.cursor()
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

def display_predicted_areas(conn):
    cursor = conn.cursor()
    cursor.execute("SELECT * FROM predicted_therapeutic_areas")
    predictions = cursor.fetchall()
    df = pd.DataFrame(predictions, columns=['compound_id', 'smiles', 'predicted_areas', 'prediction_date'])
    print(df)
   

from predict_therapeutic_areas import predict_therapeutic_areas

import random
import json

def expand_dataset_from_pubchem():
    conn = create_connection('chempath_database.db')
    
    # Load previously added compounds
    try:
        with open('added_compounds.json', 'r') as f:
            added_compounds = set(json.load(f))
    except FileNotFoundError:
        added_compounds = set()

    # Larger list of compound names
    all_compounds = [
        "Aspirin", "Ibuprofen", "Paracetamol", "Metformin", "Atorvastatin",
        "Omeprazole", "Amoxicillin", "Lisinopril", "Levothyroxine", "Amlodipine",
        "Curcumin", "Quercetin", "Resveratrol", "Epigallocatechin gallate", "Berberine",
        "Capsaicin", "Gingerol", "Allicin", "Lycopene", "Sulforaphane",
        "Caffeine", "Theobromine", "Theophylline", "Nicotine", "Morphine",
        "Codeine", "Quinine", "Artemisinin", "Taxol", "Vincristine"
    ]

    # Randomly select 10 compounds that haven't been added before
    compounds_to_add = random.sample([c for c in all_compounds if c not in added_compounds], min(10, len(all_compounds) - len(added_compounds)))

    therapeutic_areas = ["anti-inflammatory", "analgesic", "antipyretic", "antidiabetic", "lipid-lowering",
                         "antacid", "antibiotic", "antihypertensive", "thyroid hormone", "calcium channel blocker",
                         "antioxidant", "anticancer", "neuroprotective", "cardioprotective", "antimicrobial"]

    for compound_name in compounds_to_add:
        compound_data = fetch_pubchem_data(compound_name)
        if compound_data:
            area = random.choice(therapeutic_areas)
            compound = (
                compound_data['name'],
                compound_data['smiles'],
                compound_data['molecular_weight'],
                compound_data['logp'],
                "Various",  # plant_source
                area,  # biological_activity
                "Various",  # traditional_use
                compound_data['h_bond_donors'],
                compound_data['h_bond_acceptors'],
                compound_data['polar_surface_area'],
                compound_data['rotatable_bonds']
            )
            insert_compound(conn, compound)
            added_compounds.add(compound_name)
            print(f"Added {compound_data['name']} with therapeutic area: {area}")

    conn.close()

    # Save updated list of added compounds
    with open('added_compounds.json', 'w') as f:
        json.dump(list(added_compounds), f)

    print("Dataset expansion completed.")




from pathlib import Path
from chempath_core import create_connection, create_tables, create_indexes
from chempath_api import ChemPathAPI

def display_compounds(compounds):
    if not compounds:
        print("No compounds found.")
        return
    
    headers = ["ID", "Name", "SMILES", "Mol Weight", "LogP", "Plant Source", "Bio Activity", "Trad Use"]
    row_format = "{:>5} {:<20} {:<20} {:>10} {:>8} {:<15} {:<15} {:<15}"
    
    print(row_format.format(*headers))
    print("-" * 120)
    
    for compound in compounds:
        print(row_format.format(
            compound[0],
            compound[1][:20],
            compound[2][:20],
            f"{compound[3]:.2f}",
            f"{compound[4]:.2f}",
            compound[5][:15],
            compound[6][:15],
            compound[7][:15]
        ))

import sqlite3
import csv
from pathlib import Path
from chempath_core import create_connection, create_tables, create_indexes, setup_logging
from chempath_api import ChemPathAPI
from rdkit import Chem

def validate_smiles(smiles):
    mol = Chem.MolFromSmiles(smiles)
    return mol is not None

def validate_smarts(smarts):
    pattern = Chem.MolFromSmarts(smarts)
    return pattern is not None


print("Defining main function")
def main():
    print("Entering main function")
    from chempath_api import ChemPathAPI
    database = Path("chempath_database.db")
    api = ChemPathAPI(database)
    create_tables(api.conn)
    if api.conn is not None:
        create_tables(api.conn)
        create_indexes(api.conn)

        while True:
            print("\nChemPath Menu:")
            print("1. Display all compounds")
            print("2. Add external compound")
            print("3. Predict therapeutic areas")
            print("4. Display predicted areas")
            print("5. Expand dataset from PubChem")
            print("6. Search compounds")
            print("7. Optimize structure")
            print("8. Generate analogs")
            print("9. Explore chemical space")
            print("10. Perform retrosynthetic analysis")
            print("11. Perform advanced retrosynthetic analysis")
            print("12. Analyze synthetic feasibility")
            print("13. Analyze reagent availability")
            print("14. Perform comprehensive analysis")
            print("15. Exit")

            choice = input("Enter your choice (1-15): ")
            
            if choice == '1':
                display_all_compounds(api.conn)
            elif choice == '2':
                add_external_compound(api.conn)
            elif choice == '3':
                smiles = input("Enter SMILES string for the compound: ")
                predicted_areas = api.predict_therapeutic_areas(smiles)
                print("Predicted therapeutic areas:", predicted_areas)
                cursor = api.conn.cursor()
                cursor.execute("INSERT INTO predicted_therapeutic_areas (smiles, predicted_areas) VALUES (?, ?)",
                               (smiles, ','.join(predicted_areas)))
                api.conn.commit()
                print("Prediction stored in database.")
            elif choice == '4':
                display_predicted_areas(api.conn)
            elif choice == '5':
                expand_dataset_from_pubchem(api)
            elif choice == '6':
                query = input("Enter search query (compound name or SMILES): ")
                filters = {}
                if input("Apply filters? (y/n): ").lower() == 'y':
                    filters['molecular_weight'] = [float(input("Min molecular weight: ")), float(input("Max molecular weight: "))]
                    filters['logp'] = [float(input("Min logP: ")), float(input("Max logP: "))]
                    filters['biological_activity'] = input("Biological activity: ")
                
                page = 1
                while True:
                    compounds, total_count = api.search(query, filters, page)
                    for compound in compounds:
                        print(compound)
                    print(f"Page {page} of {(total_count - 1) // 10 + 1}")
                    action = input("Next page (n), previous page (p), or quit (q)? ")
                    if action.lower() == 'n':
                        page += 1
                    elif action.lower() == 'p':
                        page = max(1, page - 1)
                    else:
                        break
            elif choice == '7':
                smiles = input("Enter SMILES string for the molecule: ")
                print("Choose optimization type:")
                print("1. Functional group substitution")
                print("2. Ring system alteration")
                print("3. Scaffold hopping")
                opt_choice = input("Enter your choice (1-3): ")
                
                if opt_choice == '1':
                    target = input("Enter target functional group (SMARTS): ")
                    replacement = input("Enter replacement group (SMARTS): ")
                    result = api.optimize_structure(smiles, 'functional_group', {'target': target, 'replacement': replacement})
                elif opt_choice == '2':
                    alteration_type = input("Enter alteration type (expand/contract/fuse): ")
                    result = api.optimize_structure(smiles, 'ring_system', {'alteration_type': alteration_type})
                elif opt_choice == '3':
                    scaffold_library = input("Enter scaffold SMILES (comma-separated): ").split(',')
                    result = api.optimize_structure(smiles, 'scaffold', {'scaffold_library': scaffold_library})
                else:
                    print("Invalid choice")
                    continue
                
                print("Optimization result:", result)
            elif choice == '8':
                smiles = input("Enter SMILES string for the compound: ")
                num_analogs = int(input("Enter number of analogs to generate: "))
                analogs = api.generate_analogs(smiles, num_analogs)
                print("Generated analogs:")
                for i, analog in enumerate(analogs, 1):
                    print(f"{i}. {analog}")
            elif choice == '9':
                smiles = input("Enter SMILES of the starting molecule: ")
                num_iterations = int(input("Enter number of iterations for exploration: "))
                explored_molecules = api.explore_chemical_space(smiles, num_iterations)
                print(f"Explored {len(explored_molecules)} new molecules:")
                for i, mol in enumerate(explored_molecules, 1):
                    print(f"{i}. {mol}")
            elif choice == '10':
                perform_retrosynthetic_analysis()
            elif choice == '11':
                perform_advanced_retrosynthetic_analysis()
            elif choice == '12':
                smiles = input("Enter the SMILES of the molecule to analyze: ")
                target_smiles = input("Enter the SMILES of the target molecule: ")
                analyze_synthetic_feasibility(smiles, target_smiles)
            elif choice == '13':
                reagent_name = input("Enter the name of the reagent to analyze: ")
                analyze_reagent_availability(reagent_name)
            elif choice == '14':
                smiles = input("Enter SMILES string for comprehensive analysis: ")
                results = comprehensive_analysis(api, smiles)
                print("\nComprehensive Analysis Results:")
                print(f"Optimization: {results['optimization']}")
                print(f"Retrosynthesis: {results['retrosynthesis']}")
                print(f"Therapeutic Areas: {results['therapeutic_areas']}")
            elif choice == '15':
                print("Exiting ChemPath. Goodbye!")
                break
            else:
                print("Invalid choice. Please try again.")

def perform_retrosynthetic_analysis():
    smiles = input("Enter the SMILES string of the target molecule: ")
    while True:
        try:
            depth = int(input("Enter the depth of analysis (1-3): "))
            if 1 <= depth <= 3:
                break
            else:
                print("Please enter a number between 1 and 3.")
        except ValueError:
            print("Please enter a valid integer.")
    retrosynthetic_analysis(smiles, depth)

def perform_advanced_retrosynthetic_analysis():
    smiles = input("Enter the SMILES string of the target molecule: ")
    while True:
        try:
            depth = int(input("Enter the depth of analysis (1-5): "))
            if 1 <= depth <= 5:
                break
            else:
                print("Please enter a number between 1 and 5.")
        except ValueError:
            print("Please enter a valid integer.")
    advanced_retrosynthetic_analysis(smiles, depth)

def get_optimization_result(conn, compound_id):
    cursor = conn.cursor()
    cursor.execute("SELECT optimization_data FROM optimization_results WHERE compound_id = ?", (compound_id,))
    result = cursor.fetchone()
    return json.loads(result[0]) if result else None

def store_retrosynthesis_result(conn, compound_id, retrosynthesis_data):
    cursor = conn.cursor()
    cursor.execute("""
        INSERT INTO retrosynthesis_results (compound_id, retrosynthesis_data)
        VALUES (?, ?)
    """, (compound_id, json.dumps(retrosynthesis_data)))
    conn.commit()

def get_retrosynthesis_result(conn, compound_id):
    cursor = conn.cursor()
    cursor.execute("SELECT retrosynthesis_data FROM retrosynthesis_results WHERE compound_id = ?", (compound_id,))
    result = cursor.fetchone()
    return json.loads(result[0]) if result else None

print("About to check if __name__ == '__main__'")

def comprehensive_analysis(api, smiles):
    print(f"\nPerforming comprehensive analysis for compound: {smiles}")
    
    # Step 1: Optimization
    print("\nStep 1: Structure Optimization")
    optimized_result = api.optimize_structure(smiles, 'functional_group', {'target': 'C(=O)O', 'replacement': 'C(=O)N'})
    print(f"Optimization result: {optimized_result}")
    
    # Step 2: Retrosynthesis
    print("\nStep 2: Retrosynthetic Analysis")
    retrosynthesis_result = perform_retrosynthesis(api.conn, 1, smiles)  # Assuming compound_id 1
    print(f"Retrosynthesis result: {retrosynthesis_result}")
    
    # Step 3: Therapeutic Area Prediction
    print("\nStep 3: Therapeutic Area Prediction")
    therapeutic_areas = predict_therapeutic_areas(smiles, [])
    print(f"Predicted therapeutic areas: {therapeutic_areas}")
    
    return {
        'optimization': optimized_result,
        'retrosynthesis': retrosynthesis_result,
        'therapeutic_areas': therapeutic_areas
    }

def load_sample_data(conn):
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
        create_tables(conn)
        load_sample_data(conn)
        conn.close()
    else:
        print("Error! Cannot create the database connection.")

if __name__ == '__main__':
    main()




