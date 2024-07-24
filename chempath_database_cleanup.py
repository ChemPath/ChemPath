# File: chempath_database_cleanup.py

import sqlite3
from chempath_database import create_connection

def remove_duplicates(conn):
    """Remove duplicate entries, keeping the most complete one"""
    cursor = conn.cursor()
    cursor.execute("""
    DELETE FROM plant_compounds
    WHERE id NOT IN (
        SELECT MIN(id)
        FROM plant_compounds
        GROUP BY name
    )
    """)
    conn.commit()
    print(f"Removed {cursor.rowcount} duplicate entries")

def update_or_insert_compound(conn, compound):
    """Update existing compound or insert new one if it doesn't exist"""
    cursor = conn.cursor()
    
    # Check if compound already exists
    cursor.execute("SELECT * FROM plant_compounds WHERE name = ?", (compound[0],))
    existing = cursor.fetchone()
    
    if existing:
        # Update existing compound, but only if new data is not empty
        update_sql = """
        UPDATE plant_compounds
        SET smiles = CASE WHEN ? != '' THEN ? ELSE smiles END,
            molecular_weight = CASE WHEN ? != 0 THEN ? ELSE molecular_weight END,
            logp = CASE WHEN ? != 0 THEN ? ELSE logp END,
            plant_source = CASE WHEN ? != '' THEN ? ELSE plant_source END,
            therapeutic_areas = CASE WHEN ? != '' THEN ? ELSE therapeutic_areas END,
            biological_activities = CASE WHEN ? != '' THEN ? ELSE biological_activities END,
            traditional_use = CASE WHEN ? != '' THEN ? ELSE traditional_use END
        WHERE name = ?
        """
        cursor.execute(update_sql, compound[1:] + compound[1:] + (compound[0],))
    else:
        # Insert new compound
        insert_sql = """
        INSERT INTO plant_compounds (name, smiles, molecular_weight, logp, plant_source, therapeutic_areas, biological_activities, traditional_use)
        VALUES (?, ?, ?, ?, ?, ?, ?, ?)
        """
        cursor.execute(insert_sql, compound)
    
    conn.commit()

def cleanup_database():
    conn = create_connection("chempath_database.db")
    if conn is not None:
        remove_duplicates(conn)
        
        # Known compound data
        compounds = [
            ("Quercetin", "O=C1C(O)=C(O)C(=O)C2=C1C=C(O)C(O)=C2O", 302.24, 1.54, 
             "Various fruits and vegetables", 
             "Antioxidant effects, Anti-inflammatory properties",
             "Free radical scavenging, Enzyme inhibition", 
             "Traditional Chinese Medicine for cardiovascular health"),
            ("Curcumin", "COC1=CC(=CC(=C1O)OC)C=CC(=O)CC(=O)C=CC2=CC(=C(C=C2)O)OC", 368.38, 3.29,
             "Turmeric (Curcuma longa)", 
             "Anti-inflammatory effects, Cancer prevention",
             "NF-κB inhibition, Antioxidant activity",
             "Ayurvedic medicine for various ailments"),
            ("Resveratrol", "OC1=CC(=CC(=C1)O)C=CC2=CC(=C(C=C2)O)O", 228.24, 3.1,
             "Grapes, red wine, peanuts", 
             "Cardiovascular health promotion, Anti-aging effects",
             "Sirtuin activation, Antioxidant activity",
             "Traditional use in some Asian medicines"),
            ("Caffeine", "CN1C=NC2=C1C(=O)N(C(=O)N2C)C", 194.19, -0.07,
             "Coffee beans, tea leaves, cacao beans",
             "Central nervous system stimulation, Cognitive enhancement",
             "Adenosine receptor antagonism, Phosphodiesterase inhibition",
             "Traditional use as a stimulant in various cultures"),
            ("Theobromine", "CN1C=NC2=C1C(=O)NC(=O)N2", 180.16, -0.78,
             "Cacao beans, tea leaves",
             "Mild stimulant effects, Cardiovascular effects",
             "Phosphodiesterase inhibition, Adenosine receptor antagonism",
             "Traditional use in chocolate and cocoa-based beverages"),
            ("Capsaicin", "COC1=C(C=CC(=C1)CNC(=O)CCCC\\C=C/C(C)C)O", 305.41, 3.04,
             "Chili peppers",
             "Pain relief, Anti-inflammatory effects",
             "TRPV1 receptor activation, Substance P depletion",
             "Traditional use in various cuisines and folk medicine"),
            ("Lycopene", "CC(C)=CCC\\C(C)=C\\C\\C=C(\\C)CCC=C(C)C=CC=C(C)C=CC=C(C)C=CC=C(C)C=C", 536.87, 9.16,
             "Tomatoes, watermelons, pink grapefruits",
             "Antioxidant effects, Potential cancer prevention",
             "Free radical scavenging, Modulation of cell signaling pathways",
             "Traditional consumption as part of Mediterranean diet"),
            ("Gingerol", "CCCCC[C@H](O)CC(=O)CCc1ccc(O)c(OC)c1", 294.39, 3.85,
             "Ginger root",
             "Anti-inflammatory effects, Digestive aid",
             "COX-2 inhibition, Antioxidant activity",
             "Traditional use in Asian medicine for various ailments")
        ]
        
        for compound in compounds:
            update_or_insert_compound(conn, compound)
        
        print("Database cleaned and updated successfully")
        conn.close()
    else:
        print("Error! Cannot create the database connection.")

if __name__ == "__main__":
    cleanup_database()