# File: chempath_data_reorganization.py

import sqlite3
from src.database.chempath_database import create_connection

def reorganize_data(conn):
    """Reorganize misplaced data in the database"""
    cursor = conn.cursor()
    
    # Fetch all rows
    cursor.execute("SELECT * FROM plant_compounds")
    rows = cursor.fetchall()
    
    for row in rows:
        # Unpack the row based on the actual structure
        id, name, smiles, molecular_weight, logp, plant_source, biological_activity, traditional_use, therapeutic_areas, biological_activities = row
        
        # Reorganize the data
        reorganized_data = (
            name,
            logp,  # This is the correct SMILES
            float(smiles) if smiles and smiles.replace('.', '').isdigit() else 0,  # This is the correct Molecular Weight
            float(therapeutic_areas) if therapeutic_areas and therapeutic_areas.replace('.', '', 1).replace('-', '').isdigit() else 0,  # This is the correct LogP
            molecular_weight,  # This is the correct Plant Source
            biological_activities,  # This is the correct Therapeutic Areas
            biological_activity,  # This is the correct Biological Activities
            traditional_use
        )
        
        print(f"Reorganized data for {name}: {reorganized_data}")
        
        # Update the row
        update_sql = """
        UPDATE plant_compounds
        SET smiles = ?, 
            molecular_weight = ?, 
            logp = ?,
            plant_source = ?,
            therapeutic_areas = ?,
            biological_activities = ?,
            traditional_use = ?
        WHERE id = ?
        """
        cursor.execute(update_sql, reorganized_data[1:] + (id,))
    
    conn.commit()
    print("Data reorganization completed")

def main():
    conn = create_connection("chempath_database.db")
    if conn is not None:
        reorganize_data(conn)
        conn.close()
    else:
        print("Error! Cannot create the database connection.")

if __name__ == "__main__":
    main()