# File: chempath_data_correction.py

import sqlite3
from chempath_database import create_connection

def get_table_info(conn):
    """Get information about the structure of the plant_compounds table"""
    cursor = conn.cursor()
    cursor.execute("PRAGMA table_info(plant_compounds)")
    return cursor.fetchall()

def print_table_structure(table_info):
    """Print the structure of the plant_compounds table"""
    print("Current table structure:")
    for column in table_info:
        print(f"Column {column[0]}: {column[1]} ({column[2]})")

def correct_data(conn):
    """Correct misplaced data in the database"""
    cursor = conn.cursor()
    
    # Get table information
    table_info = get_table_info(conn)
    print_table_structure(table_info)
    
    # Fetch all rows
    cursor.execute("SELECT * FROM plant_compounds")
    rows = cursor.fetchall()
    
    for row in rows:
        print(f"Original row: {row}")
        
        # Ensure we have the expected number of columns
        if len(row) != len(table_info):
            print(f"Warning: Row has {len(row)} columns, expected {len(table_info)}")
            continue
        
        # Unpack the row based on the actual number of columns
        id = row[0]
        name = row[1] if len(row) > 1 else ""
        smiles = row[2] if len(row) > 2 else ""
        molecular_weight = row[3] if len(row) > 3 else 0
        logp = row[4] if len(row) > 4 else 0
        plant_source = row[5] if len(row) > 5 else ""
        therapeutic_areas = row[6] if len(row) > 6 else ""
        biological_activities = row[7] if len(row) > 7 else ""
        traditional_use = row[8] if len(row) > 8 else ""
        
        # Correct the data
        corrected_data = (
            name,
            plant_source,  # This was in SMILES
            float(smiles) if smiles and smiles.replace('.', '').isdigit() else 0,  # This was in Molecular Weight
            float(traditional_use) if traditional_use and traditional_use.replace('.', '').isdigit() else 0,  # This was in LogP
            molecular_weight,  # This was in Plant Source
            therapeutic_areas,
            biological_activities,
            logp  # This was in Traditional Use
        )
        
        print(f"Corrected data: {corrected_data}")
        
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
        cursor.execute(update_sql, corrected_data[1:] + (id,))
    
    conn.commit()
    print("Data correction attempt completed")

def main():
    conn = create_connection("chempath_database.db")
    if conn is not None:
        correct_data(conn)
        conn.close()
    else:
        print("Error! Cannot create the database connection.")

if __name__ == "__main__":
    main()