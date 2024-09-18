import os
import pandas as pd

# Absolute path to the cleaned data directory and the specific file
data_directory = 'C:/Users/Dr. Contessa Petrini/ChemPath/src/data/cleaned_data/'
file_path = os.path.join(data_directory, 'Aidicha_ingredients.xlsx')

# Read the Excel file into a DataFrame
df = pd.read_excel(file_path)

# Extract herb name from the file name (e.g., 'Aidicha' from 'Aidicha_ingredients.xlsx')
herb_name = os.path.basename(file_path).split('_')[0]

# Placeholder herb_id (replace this with the actual herb_id from your database)
herb_id = 1  # Change this to the actual herb_id of 'Aidicha' in your Supabase

# Generate SQL statements
sql_statements = []

for index, row in df.iterrows():
    # Extract values from the DataFrame row
    molecule_name = row['molecule_name'].replace("'", "''")  # Escape single quotes in molecule names

    # Construct the SQL INSERT statement using the correct column names from your table
    sql_statement = f"""
    INSERT INTO ingredients (name, herb_id, quantity)
    VALUES ('{molecule_name}', {herb_id}, NULL);
    """
    
    sql_statements.append(sql_statement)

# Write SQL statements to a .sql file using UTF-8 encoding
output_file = os.path.join(data_directory, 'insert_ingredients.sql')
with open(output_file, 'w', encoding='utf-8') as file:
    file.write('\n'.join(sql_statements))

print(f'SQL statements have been written to {output_file}')
