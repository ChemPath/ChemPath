from src.database.chempath_database import create_connection, create_table, insert_compound

# Create a connection to the database
conn = create_connection("chempath_database.db")

# Create the table
create_table(conn)

# Insert a new compound
new_compound = ("Compound Name", "SMILES", 100.0, 1.0, "Plant Source", "Activity", "Traditional Use")
insert_compound(conn, new_compound)

# Don't forget to close the connection when you're done
conn.close()
