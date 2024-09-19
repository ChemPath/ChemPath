from supabase import create_client, Client
import pandas as pd

# Supabase configuration
url = 'https://lqwxjiijnhefterkeety.supabase.co'  # Replace with your Supabase URL
key = 'eyJhbGciOiJIUzI1NiIsInR5cCI6IkpXVCJ9.eyJpc3MiOiJzdXBhYmFzZSIsInJlZiI6Imxxd3hqaWlqbmhlZnRlcmtlZXR5Iiwicm9sZSI6InNlcnZpY2Vfcm9sZSIsImlhdCI6MTcyNjU5NTk3OSwiZXhwIjoyMDQyMTcxOTc5fQ.bmEG-1VU5MgU9TAWsXS0FEwVpH3qgc49dyUqXxvK3Kk'
  # Replace with your Supabase API key
supabase: Client = create_client(url, key)

# Load the CSV file
file_path = r"C:\Users\Dr. Contessa Petrini\OneDrive\Documents\Duke-Source-CSV\ETHNOBOT.csv"
df = pd.read_csv(file_path)

# Rename columns to match the database schema
df = df.rename(columns={
    'ACTIVITY': 'activity',
    'GENUS': 'genus',
    'SPECIES': 'species',
    'FAMILY': 'family',
    'CNAME': 'common_name',
    'COUNTRY': 'country',
    'REFERENCE': 'reference',
    'LONGREF': 'long_reference',
    'TAXON': 'taxon',
    'TAXAUTHOR': 'tax_author',
    'EFFECTIVE': 'effective',
    'CREATED': 'created',
    'MODIFIED': 'modified'
})

# Ensure genus and species are of the correct data type (string)
df['genus'] = df['genus'].astype(str)
df['species'] = df['species'].astype(str)

# Fill missing values
df['common_name'] = df['common_name'].fillna('Unknown')
df['family'] = df['family'].fillna('Unknown')
df['country'] = df['country'].fillna('Not Specified')
df['reference'] = df['reference'].fillna('Not Available')
df['long_reference'] = df['long_reference'].fillna('Not Available')
df['effective'] = df['effective'].fillna('Not Available')
df['tax_author'] = df['tax_author'].fillna('Not Available')

# Handle date columns
df['created'] = pd.to_datetime(df['created'], errors='coerce')
df['modified'] = pd.to_datetime(df['modified'], errors='coerce').fillna(pd.Timestamp('1900-01-01'))

# Prepare data for batch upsert
batch_size = 100
rows = df.to_dict(orient='records')

for i in range(0, len(rows), batch_size):
    batch = rows[i:i + batch_size]
    batch_data = [
        {
            'genus': row['genus'],
            'species': row['species'],
            'common_name': row['common_name'],
            'family': row['family'],
            'country': row['country'],
            'reference': row['reference'],
            'long_reference': row['long_reference'],
            'taxon': row['taxon'],
            'tax_author': row['tax_author'],
            'effective': row['effective'],
            'created': row['created'].isoformat() if pd.notna(row['created']) else None,
            'modified': row['modified'].isoformat()
        }
        for row in batch
    ]

    # Upsert batch data into the table using genus and species for conflict resolution
    supabase.table('ethnobotanical_uses').upsert(batch_data, on_conflict=['genus', 'species']).execute()
