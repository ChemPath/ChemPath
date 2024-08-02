# File: chempath_final_structure_and_data_correction.py

import sqlite3
from src.database.chempath_database import create_connection

def correct_database_structure(conn):
    """Correct the database structure"""
    cursor = conn.cursor()
    
    try:
        # Rename 'biological_activity' to 'biological_activities'
        cursor.execute("ALTER TABLE plant_compounds RENAME COLUMN biological_activity TO biological_activities")
        print("Renamed 'biological_activity' to 'biological_activities'")
    except sqlite3.OperationalError:
        print("Column 'biological_activities' already exists")

    # Note: SQLite doesn't support dropping columns, so we can't remove the extra 'biological_activities' column

def update_database(conn):
    """Update the database with correct information"""
    cursor = conn.cursor()
    
    # Correct data for each compound
    correct_data = {
        "Quercetin": ("O=C1C(O)=C(O)C(=O)C2=C1C=C(O)C(O)=C2O", 302.24, 1.54, 
                      "Various fruits and vegetables",
                      "Antioxidant effects, Anti-inflammatory properties",
                      "Free radical scavenging, Enzyme inhibition", 
                      "Traditional Chinese Medicine for cardiovascular health"),
        "Curcumin": ("COC1=CC(=CC(=C1O)OC)C=CC(=O)CC(=O)C=CC2=CC(=C(C=C2)O)OC", 368.38, 3.29,
                     "Turmeric (Curcuma longa)",
                     "Anti-inflammatory effects, Cancer prevention",
                     "NF-κB inhibition, Antioxidant activity",
                     "Ayurvedic medicine for various ailments"),
        "Resveratrol": ("OC1=CC(=CC(=C1)O)C=CC2=CC(=C(C=C2)O)O", 228.24, 3.1,
                        "Grapes, red wine, peanuts",
                        "Cardiovascular health promotion, Anti-aging effects",
                        "Sirtuin activation, Antioxidant activity",
                        "Traditional use in some Asian medicines"),
        "Caffeine": ("CN1C=NC2=C1C(=O)N(C(=O)N2C)C", 194.19, -0.07,
                     "Coffee beans, tea leaves, cacao beans",
                     "Central nervous system stimulation, Cognitive enhancement",
                     "Adenosine receptor antagonism, Phosphodiesterase inhibition",
                     "Traditional use as a stimulant in various cultures"),
        "Theobromine": ("CN1C=NC2=C1C(=O)NC(=O)N2", 180.16, -0.78,
                        "Cacao beans, tea leaves",
                        "Mild stimulant effects, Cardiovascular effects",
                        "Phosphodiesterase inhibition, Adenosine receptor antagonism",
                        "Traditional use in chocolate and cocoa-based beverages"),
        "Capsaicin": ("COC1=C(C=CC(=C1)CNC(=O)CCCC\\C=C/C(C)C)O", 305.41, 3.04,
                      "Chili peppers",
                      "Pain relief, Anti-inflammatory effects",
                      "TRPV1 receptor activation, Substance P depletion",
                      "Traditional use in various cuisines and folk medicine"),
        "Lycopene": ("CC(C)=CCC\\C(C)=C\\C\\C=C(\\C)CCC=C(C)C=CC=C(C)C=CC=C(C)C=CC=C(C)C=C", 536.87, 9.16,
                     "Tomatoes, watermelons, pink grapefruits",
                     "Antioxidant effects, Potential cancer prevention",
                     "Free radical scavenging, Modulation of cell signaling pathways",
                     "Traditional consumption as part of Mediterranean diet"),
        "Gingerol": ("CCCCC[C@H](O)CC(=O)CCc1ccc(O)c(OC)c1", 294.39, 3.85,
                     "Ginger root",
                     "Anti-inflammatory effects, Digestive aid",
                     "COX-2 inhibition, Antioxidant activity",
                     "Traditional use in Asian medicine for various ailments")
    }
    
    try:
        conn.execute("BEGIN TRANSACTION")
        for name, data in correct_data.items():
            smiles, molecular_weight, logp, plant_source, therapeutic_areas, biological_activities, traditional_use = data
            
            update_sql = """
            UPDATE plant_compounds
            SET smiles = ?,
                molecular_weight = ?,
                logp = ?,
                plant_source = ?,
                therapeutic_areas = ?,
                biological_activities = ?,
                traditional_use = ?
            WHERE name = ?
            """
            cursor.execute(update_sql, (smiles, molecular_weight, logp, plant_source, therapeutic_areas, biological_activities, traditional_use, name))
            print(f"Updated {name}")
        
        conn.commit()
        print("Database update completed successfully")
    except sqlite3.Error as e:
        conn.rollback()
        print(f"An error occurred: {e}")
    finally:
        cursor.close()

def inspect_database(conn):
    """Inspect the current structure and content of the database"""
    cursor = conn.cursor()
    
    # Get table structure
    cursor.execute("PRAGMA table_info(plant_compounds)")
    columns = cursor.fetchall()
    print("Database Structure:")
    for column in columns:
        print(f"Column: {column[1]}, Type: {column[2]}")
    print("\n")

    # Get current data
    cursor.execute("SELECT * FROM plant_compounds")
    rows = cursor.fetchall()
    print("Current Data:")
    for row in rows:
        print(row)
    print("\n")

def main():
    conn = create_connection("chempath_database.db")
    if conn is not None:
        print("Before Update:")
        inspect_database(conn)
        
        correct_database_structure(conn)
        update_database(conn)
        
        print("\nAfter Update:")
        inspect_database(conn)
        
        conn.close()
    else:
        print("Error! Cannot create the database connection.")

if __name__ == "__main__":
    main()