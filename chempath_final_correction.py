# File: chempath_final_correction.py

import sqlite3
from chempath_database import create_connection

def final_data_correction(conn):
    """Make final corrections to the database"""
    cursor = conn.cursor()
    
    # Correct data for each compound
    correct_data = {
        "Quercetin": ("Antioxidant effects, Anti-inflammatory properties",
                      "Free radical scavenging, Enzyme inhibition", 
                      "Traditional Chinese Medicine for cardiovascular health"),
        "Curcumin": ("Anti-inflammatory effects, Cancer prevention",
                     "NF-κB inhibition, Antioxidant activity",
                     "Ayurvedic medicine for various ailments"),
        "Resveratrol": ("Cardiovascular health promotion, Anti-aging effects",
                        "Sirtuin activation, Antioxidant activity",
                        "Traditional use in some Asian medicines"),
        "Caffeine": ("Central nervous system stimulation, Cognitive enhancement",
                     "Adenosine receptor antagonism, Phosphodiesterase inhibition",
                     "Traditional use as a stimulant in various cultures"),
        "Theobromine": ("Mild stimulant effects, Cardiovascular effects",
                        "Phosphodiesterase inhibition, Adenosine receptor antagonism",
                        "Traditional use in chocolate and cocoa-based beverages"),
        "Capsaicin": ("Pain relief, Anti-inflammatory effects",
                      "TRPV1 receptor activation, Substance P depletion",
                      "Traditional use in various cuisines and folk medicine"),
        "Lycopene": ("Antioxidant effects, Potential cancer prevention",
                     "Free radical scavenging, Modulation of cell signaling pathways",
                     "Traditional consumption as part of Mediterranean diet"),
        "Gingerol": ("Anti-inflammatory effects, Digestive aid",
                     "COX-2 inhibition, Antioxidant activity",
                     "Traditional use in Asian medicine for various ailments")
    }
    
    for name, data in correct_data.items():
        therapeutic_areas, biological_activities, traditional_use = data
        
        update_sql = """
        UPDATE plant_compounds
        SET therapeutic_areas = ?,
            biological_activities = ?,
            traditional_use = ?
        WHERE name = ?
        """
        cursor.execute(update_sql, (therapeutic_areas, biological_activities, traditional_use, name))
    
    conn.commit()
    print("Final data correction completed")

def main():
    conn = create_connection("chempath_database.db")
    if conn is not None:
        final_data_correction(conn)
        conn.close()
    else:
        print("Error! Cannot create the database connection.")

if __name__ == "__main__":
    main()