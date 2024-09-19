import os
import pandas as pd
import requests
from supabase import create_client, Client

# Supabase configuration
SUPABASE_URL = 'https://lqwxjiijnhefterkeety.supabase.co'
SUPABASE_KEY = 'eyJhbGciOiJIUzI1NiIsInR5cCI6IkpXVCJ9.eyJpc3MiOiJzdXBhYmFzZSIsInJlZiI6Imxxd3hqaWlqbmhlZnRlcmtlZXR5Iiwicm9sZSI6InNlcnZpY2Vfcm9sZSIsImlhdCI6MTcyNjU5NTk3OSwiZXhwIjoyMDQyMTcxOTc5fQ.bmEG-1VU5MgU9TAWsXS0FEwVpH3qgc49dyUqXxvK3Kk'
supabase_client: Client = create_client(SUPABASE_URL, SUPABASE_KEY)

# Set up headers for REST API requests
headers = {
    'Content-Type': 'application/json',
    'apikey': SUPABASE_KEY,
    'Authorization': f'Bearer {SUPABASE_KEY}'
}

# Directory where the new data files are stored
data_directory = 'C:/Users/Dr. Contessa Petrini/ChemPath/src/data/cleaned_data/'

# Function to insert a new herb into the herbs table if it doesn't already exist
def insert_herb(herb_name):
    # Check if the herb already exists
    response = supabase_client.table('herbs').select('id').eq('name', herb_name).execute()
    if response.data:
        print(f"Herb '{herb_name}' already exists in the herbs table.")
        return response.data[0]['id']
    
    # Insert the new herb
    data = {
        'name': herb_name,
        'description': f'Description for {herb_name}'  # Placeholder description
    }
    response = requests.post(f'{SUPABASE_URL}/rest/v1/herbs', json=data, headers=headers)
    if response.status_code == 201:
        print(f"Successfully inserted herb '{herb_name}' into herbs table.")
        # Fetch the newly inserted herb's ID
        herb_response = supabase_client.table('herbs').select('id').eq('name', herb_name).execute()
        return herb_response.data[0]['id']
    else:
        print(f"Error inserting herb '{herb_name}': {response.text}")
        return None

# Function to insert data into a specific table, avoiding duplicates
def insert_data(file_path, table_name, herb_id, herb_name):
    df = pd.read_excel(file_path)
    
    if table_name == 'diseases':
        for _, row in df.iterrows():
            disease_name = row['disease_name'].replace("'", "''")
            description = row.get('description', '').replace("'", "''")
            
            # Check for duplicates
            response = supabase_client.table('diseases').select('id').eq('disease_name', disease_name).eq('herb_id', herb_id).execute()
            if response.data:
                print(f"Disease '{disease_name}' already exists for this herb.")
                continue
            
            # Insert data into diseases table
            data = {
                'herb_id': herb_id,
                'disease_name': disease_name,
                'description': description
            }
            response = requests.post(f'{SUPABASE_URL}/rest/v1/diseases', json=data, headers=headers)
            if response.status_code != 201:
                print(f"Error inserting into diseases: {response.text}")
            else:
                print(f"Successfully inserted '{disease_name}' into diseases.")

    elif table_name == 'ingredients':
        for _, row in df.iterrows():
            ingredient_name = row['molecule_name'].replace("'", "''")
            
            # Check for duplicates
            response = supabase_client.table('ingredients').select('id').eq('name', ingredient_name).eq('herb_id', herb_id).execute()
            if response.data:
                print(f"Ingredient '{ingredient_name}' already exists for this herb.")
                continue
            
            # Insert data into ingredients table
            data = {
                'name': ingredient_name,
                'herb_id': herb_id,
                'quantity': None  # Assuming quantity is not available
            }
            response = requests.post(f'{SUPABASE_URL}/rest/v1/ingredients', json=data, headers=headers)
            if response.status_code != 201:
                print(f"Error inserting into ingredients: {response.text}")
            else:
                print(f"Successfully inserted '{ingredient_name}' into ingredients.")

    elif table_name == 'targets':
        for _, row in df.iterrows():
            target_name = row['target_name'].replace("'", "''")
            target_description = row.get('target_description', '').replace("'", "''")
            
            # Check for duplicates
            response = supabase_client.table('targets').select('id').eq('target_name', target_name).eq('herb_name', herb_name).execute()
            if response.data:
                print(f"Target '{target_name}' already exists for this herb.")
                continue
            
            # Insert data into targets table
            data = {
                'herb_name': herb_name,
                'target_name': target_name,
                'target_description': target_description
            }
            response = requests.post(f'{SUPABASE_URL}/rest/v1/targets', json=data, headers=headers)
            if response.status_code != 201:
                print(f"Error inserting into targets: {response.text}")
            else:
                print(f"Successfully inserted '{target_name}' into targets.")

# Main function to insert all information for each herb
def main():
    # List of new herb files
    herb_files = [
        'Mahuanggen_disease.xlsx', 'Mahuanggen_ingredients.xlsx', 'Mahuanggen_targets.xlsx',
        'Baizhu_disease.xlsx', 'Baizhu_ingredients.xlsx', 'Baizhu_targets.xlsx',
        'Chenpi_disease.xlsx', 'Chenpi_ingredients.xlsx', 'Chenpi_targets.xlsx',
        'Mahuang_disease.xlsx', 'Mahuang_ingredients.xlsx', 'Mahuang_targets.xlsx'
    ]
    
    for file_name in herb_files:
        # Extract herb name from file name
        herb_name = file_name.split('_')[0]
        
        # Insert herb and get its ID
        herb_id = insert_herb(herb_name)
        if not herb_id:
            continue  # Skip this herb if insertion fails

        # Determine the table name based on file name
        if 'disease' in file_name.lower():
            table_name = 'diseases'
        elif 'ingredients' in file_name.lower():
            table_name = 'ingredients'
        elif 'targets' in file_name.lower():
            table_name = 'targets'
        else:
            print(f"Unknown table type for file: {file_name}")
            continue
        
        # Insert data into the respective table
        file_path = os.path.join(data_directory, file_name)
        if os.path.exists(file_path):
            print(f"Inserting data from: {file_name} into {table_name} table.")
            insert_data(file_path, table_name, herb_id, herb_name)
        else:
            print(f"File not found: {file_path}")

if __name__ == '__main__':
    main()

supabase_client: Client = create_client(SUPABASE_URL, SUPABASE_KEY)

# Set up headers for REST API requests
headers = {
    'Content-Type': 'application/json',
    'apikey': SUPABASE_KEY,
    'Authorization': f'Bearer {SUPABASE_KEY}'
}

# Directory where the new data files are stored
data_directory = 'C:/Users/Dr. Contessa Petrini/ChemPath/src/data/cleaned_data/'

# Function to insert a new herb into the herbs table if it doesn't already exist
def insert_herb(herb_name):
    # Check if the herb already exists
    response = supabase_client.table('herbs').select('id').eq('name', herb_name).execute()
    if response.data:
        print(f"Herb '{herb_name}' already exists in the herbs table.")
        return response.data[0]['id']
    
    # Insert the new herb
    data = {
        'name': herb_name,
        'description': f'Description for {herb_name}'  # Placeholder description
    }
    response = requests.post(f'{SUPABASE_URL}/rest/v1/herbs', json=data, headers=headers)
    if response.status_code == 201:
        print(f"Successfully inserted herb '{herb_name}' into herbs table.")
        # Fetch the newly inserted herb's ID
        herb_response = supabase_client.table('herbs').select('id').eq('name', herb_name).execute()
        return herb_response.data[0]['id']
    else:
        print(f"Error inserting herb '{herb_name}': {response.text}")
        return None

# Function to insert data into a specific table, avoiding duplicates
def insert_data(file_path, table_name, herb_id, herb_name):
    df = pd.read_excel(file_path)
    
    if table_name == 'diseases':
        for _, row in df.iterrows():
            disease_name = row['disease_name'].replace("'", "''")
            description = row.get('description', '').replace("'", "''")
            
            # Check for duplicates
            response = supabase_client.table('diseases').select('id').eq('disease_name', disease_name).eq('herb_id', herb_id).execute()
            if response.data:
                print(f"Disease '{disease_name}' already exists for this herb. Skipping insertion.")
                continue
            
            # Insert data into diseases table
            data = {
                'herb_id': herb_id,
                'disease_name': disease_name,
                'description': description
            }
            response = requests.post(f'{SUPABASE_URL}/rest/v1/diseases', json=data, headers=headers)
            if response.status_code != 201:
                print(f"Error inserting into diseases: {response.text}")
            else:
                print(f"Successfully inserted '{disease_name}' into diseases.")

    elif table_name == 'ingredients':
        for _, row in df.iterrows():
            ingredient_name = row['molecule_name'].replace("'", "''")
            
            # Check for duplicates
            response = supabase_client.table('ingredients').select('id').eq('name', ingredient_name).eq('herb_id', herb_id).execute()
            if response.data:
                print(f"Ingredient '{ingredient_name}' already exists for this herb.")
                continue
            
            # Insert data into ingredients table
            data = {
                'name': ingredient_name,
                'herb_id': herb_id,
                'quantity': None  # Assuming quantity is not available
            }
            response = requests.post(f'{SUPABASE_URL}/rest/v1/ingredients', json=data, headers=headers)
            if response.status_code != 201:
                print(f"Error inserting into ingredients: {response.text}")
            else:
                print(f"Successfully inserted '{ingredient_name}' into ingredients.")

    elif table_name == 'targets':
        for _, row in df.iterrows():
            target_name = row['target_name'].replace("'", "''")
            target_description = row.get('target_description', '').replace("'", "''")
            
            # Check for duplicates
            response = supabase_client.table('targets').select('id').eq('target_name', target_name).eq('herb_name', herb_name).execute()
            if response.data:
                print(f"Target '{target_name}' already exists for this herb.")
                continue
            
            # Insert data into targets table
            data = {
                'herb_name': herb_name,
                'target_name': target_name,
                'target_description': target_description
            }
            response = requests.post(f'{SUPABASE_URL}/rest/v1/targets', json=data, headers=headers)
            if response.status_code != 201:
                print(f"Error inserting into targets: {response.text}")
            else:
                print(f"Successfully inserted '{target_name}' into targets.")

# Main function to insert all information for each herb
def main():
    # Automatically find all herb files in the cleaned_data directory
    herb_files = [f for f in os.listdir(data_directory) if f.endswith('.xlsx') and '_' in f]

    for file_name in herb_files:
        # Extract herb name from file name
        herb_name = file_name.split('_')[0]
        
        # Insert herb and get its ID
        herb_id = insert_herb(herb_name)
        if not herb_id:
            continue  # Skip this herb if insertion fails

        # Determine the table name based on file name
        if 'disease' in file_name.lower():
            table_name = 'diseases'
        elif 'ingredients' in file_name.lower():
            table_name = 'ingredients'
        elif 'targets' in file_name.lower():
            table_name = 'targets'
        else:
            print(f"Unknown table type for file: {file_name}")
            continue
        
        # Insert data into the respective table
        file_path = os.path.join(data_directory, file_name)
        if os.path.exists(file_path):
            print(f"Inserting data from: {file_name} into {table_name} table.")
            insert_data(file_path, table_name, herb_id, herb_name)
        else:
            print(f"File not found: {file_path}")


if __name__ == '__main__':
    main()
