# File: check_data.py

from chempath_database import create_connection

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
            print(f"Therapeutic Areas: {row[6]}")
            print(f"Biological Activities: {row[7]}")
            print(f"Traditional Use: {row[8]}")
            print("--------------------")
        conn.close()
    else:
        print("Error! Cannot create the database connection.")

if __name__ == "__main__":
    check_data()