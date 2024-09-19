import pandas as pd
from supabase import create_client, Client

# Supabase configuration
SUPABASE_URL = 'https://lqwxjiijnhefterkeety.supabase.co'
SUPABASE_KEY = 'eyJhbGciOiJIUzI1NiIsInR5cCI6IkpXVCJ9.eyJpc3MiOiJzdXBhYmFzZSIsInJlZiI6Imxxd3hqaWlqbmhlZnRlcmtlZXR5Iiwicm9sZSI6InNlcnZpY2Vfcm9sZSIsImlhdCI6MTcyNjU5NTk3OSwiZXhwIjoyMDQyMTcxOTc5fQ.bmEG-1VU5MgU9TAWsXS0FEwVpH3qgc49dyUqXxvK3Kk'
supabase_client: Client = create_client(SUPABASE_URL, SUPABASE_KEY)

# Path to the exported TSV file
tsv_file_path = 'path/to/PYTOHUB_9.tsv'

# Function to insert PhytoHub data into Supabase
def insert_phytohub_data_into_supabase():
    # Read TSV file
    try:
        df = pd.read_csv(tsv_file_path, sep='\t')
    except Exception as e:
        print(f"Error reading TSV file: {e}")
        return

    # Iterate through each row to insert data into Supabase
    for _, row in df.iterrows():
        compound_name = row['name']
        inchi = row['inchi']
        smiles = row['smiles']
        biological_source = row['biologicalsource']
        reference = row['reference']
        
        # Check for duplicate entry
        response = supabase_client.table('traditional_compounds').select('id').eq('name', compound_name).execute()
        if response.data:
            print(f"Compound '{compound_name}' already exists in the database.")
            continue

        # Insert compound data into Supabase
        data = {
            'name': compound_name,
            'inchi': inchi,
            'smiles': smiles,
            'biological_source': biological_source,
            'reference': reference
        }
        response = supabase_client.table('traditional_compounds').insert(data).execute()
        if response.status_code != 201:
            print(f"Error inserting compound '{compound_name}': {response.text}")
        else:
            print(f"Successfully inserted compound '{compound_name}' into traditional_compounds.")

# Main function
def main():
    insert_phytohub_data_into_supabase()

if __name__ == "__main__":
    main()
