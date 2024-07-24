# File: check_data.py

from chempath_database import create_connection

def get_compound_classes(cursor, compound_id):
    cursor.execute("""
        SELECT cc.name 
        FROM compound_classes cc
        JOIN compound_class_relationships ccr ON cc.id = ccr.compound_class_id
        WHERE ccr.compound_id = ?
    """, (compound_id,))
    return [row[0] for row in cursor.fetchall()]

def get_therapeutic_areas(cursor, compound_id):
    cursor.execute("""
        SELECT ta.name 
        FROM therapeutic_areas ta
        JOIN therapeutic_area_relationships tar ON ta.id = tar.therapeutic_area_id
        WHERE tar.compound_id = ?
    """, (compound_id,))
    return [row[0] for row in cursor.fetchall()]

def check_data():
    conn = create_connection("chempath_database.db")
    if conn is not None:
        cursor = conn.cursor()
        cursor.execute("SELECT * FROM plant_compounds")
        rows = cursor.fetchall()
        for row in rows:
            print(f"Name: {row[1]}")
            print(f"SMILES: {row[2]}")
            print(f"Molecular Weight: {row[3]}")
            print(f"LogP: {row[4]}")
            print(f"Plant Source: {row[5]}")
            print(f"Biological Activity: {row[6]}")
            print(f"Traditional Use: {row[7]}")
            
            compound_classes = get_compound_classes(cursor, row[0])
            print(f"Compound Classes: {', '.join(compound_classes)}")
            
            therapeutic_areas = get_therapeutic_areas(cursor, row[0])
            print(f"Therapeutic Areas: {', '.join(therapeutic_areas)}")
            
            print("--------------------")
        conn.close()
    else:
        print("Error! Cannot create the database connection.")

if __name__ == "__main__":
    check_data()
