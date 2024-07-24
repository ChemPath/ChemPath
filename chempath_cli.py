# File: chempath_cli.py

from chempath_database import create_connection, insert_compound, get_compound_by_name, update_compound

def display_menu():
    print("\nChemPath CLI")
    print("1. Add new compound")
    print("2. Update existing compound")
    print("3. View compound")
    print("4. Exit")

def add_compound(conn):
    name = input("Compound name: ")
    smiles = input("SMILES: ")
    molecular_weight = float(input("Molecular weight: "))
    logp = float(input("LogP: "))
    plant_source = input("Plant source: ")
    therapeutic_areas = input("Therapeutic areas (comma-separated): ")
    biological_activities = input("Biological activities (comma-separated): ")
    traditional_use = input("Traditional use: ")

    compound = (name, smiles, molecular_weight, logp, plant_source, therapeutic_areas, biological_activities, traditional_use)
    insert_compound(conn, compound)
    print(f"{name} added successfully.")

def update_compound(conn):
    name = input("Enter compound name to update: ")
    compound = get_compound_by_name(conn, name)
    if compound:
        print("Enter new values (press enter to keep current value):")
        smiles = input(f"SMILES [{compound[2]}]: ") or compound[2]
        molecular_weight = input(f"Molecular weight [{compound[3]}]: ") or compound[3]
        logp = input(f"LogP [{compound[4]}]: ") or compound[4]
        plant_source = input(f"Plant source [{compound[5]}]: ") or compound[5]
        therapeutic_areas = input(f"Therapeutic areas [{compound[6]}]: ") or compound[6]
        biological_activities = input(f"Biological activities [{compound[7]}]: ") or compound[7]
        traditional_use = input(f"Traditional use [{compound[8]}]: ") or compound[8]

        updated_compound = (name, smiles, molecular_weight, logp, plant_source, therapeutic_areas, biological_activities, traditional_use)
        update_compound(conn, updated_compound)
        print(f"{name} updated successfully.")
    else:
        print(f"Compound {name} not found.")

def view_compound(conn):
    name = input("Enter compound name to view: ")
    compound = get_compound_by_name(conn, name)
    if compound:
        print(f"\nName: {compound[1]}")
        print(f"SMILES: {compound[2]}")
        print(f"Molecular weight: {compound[3]}")
        print(f"LogP: {compound[4]}")
        print(f"Plant source: {compound[5]}")
        print(f"Therapeutic areas: {compound[6]}")
        print(f"Biological activities: {compound[7]}")
        print(f"Traditional use: {compound[8]}")
    else:
        print(f"Compound {name} not found.")

def main():
    conn = create_connection("chempath_database.db")
    if conn is None:
        return

    while True:
        display_menu()
        choice = input("Enter your choice: ")
        if choice == '1':
            add_compound(conn)
        elif choice == '2':
            update_compound(conn)
        elif choice == '3':
            view_compound(conn)
        elif choice == '4':
            break
        else:
            print("Invalid choice. Please try again.")

    conn.close()

if __name__ == "__main__":
    main()