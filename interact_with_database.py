from src.database.chempath_database import create_connection, create_table, insert_compound, get_compound_by_name, update_compound, delete_compound

def main():
    # Create a connection to the database
    conn = create_connection("chempath_database.db")

    # Create the table
    create_table(conn)

    # Insert a new compound
    new_compound = ("Quercetin", "O=C1C(O)=C(O)C(=O)C2=C1C=C(O)C(O)=C2O", 302.24, 1.54, "Various fruits and vegetables", "Antioxidant, Anti-inflammatory", "Traditional Chinese Medicine", 5, 7, 131.36, 1)
    insert_compound(conn, new_compound)
    print("New compound inserted.")

    # Retrieve a compound
    compound = get_compound_by_name(conn, "Quercetin")
    if compound:
        print(f"Retrieved compound: {compound}")
    else:
        print("Compound not found.")

    # Update a compound
    updated_compound = ("Quercetin", "O=C1C(O)=C(O)C(=O)C2=C1C=C(O)C(O)=C2O", 302.24, 1.54, "Various fruits and vegetables", "Antioxidant, Anti-inflammatory, Anticancer", "Traditional Chinese Medicine", 5, 7, 131.36, 1)
    update_compound(conn, updated_compound)
    print("Compound updated.")

    # Delete a compound
    delete_compound(conn, "Quercetin")
    print("Compound deleted.")

    # Don't forget to close the connection when you're done
    conn.close()

if __name__ == "__main__":
    main()
