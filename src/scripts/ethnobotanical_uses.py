import pandas as pd
from supabase import create_client, Client

# Supabase configuration
SUPABASE_URL = 'https://lqwxjiijnhefterkeety.supabase.co'
SUPABASE_KEY = 'eyJhbGciOiJIUzI1NiIsInR5cCI6IkpXVCJ9.eyJpc3MiOiJzdXBhYmFzZSIsInJlZiI6Imxxd3hqaWlqbmhlZnRlcmtlZXR5Iiwicm9sZSI6InNlcnZpY2Vfcm9sZSIsImlhdCI6MTcyNjU5NTk3OSwiZXhwIjoyMDQyMTcxOTc5fQ.bmEG-1VU5MgU9TAWsXS0FEwVpH3qgc49dyUqXxvK3Kk'
supabase_client: Client = create_client(SUPABASE_URL, SUPABASE_KEY)

# Load the CSV file
file_path = r"C:\Users\Dr. Contessa Petrini\OneDrive\Documents\Duke-Source-CSV\ETHNOBOT.csv"
df = pd.read_csv(file_path)

# Select relevant columns for the ethnobotanical_uses table
df_filtered = df[['ACTIVITY', 'GENUS', 'SPECIES', 'FAMILY', 'CNAME', 'COUNTRY', 
                  'REFERENCE', 'LONGREF', 'TAXON', 'TAXAUTHOR', 'EFFECTIVE', 'CREATED', 'MODIFIED']].copy()

# Rename columns to match the database schema
df_filtered = df_filtered.rename(columns={
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

# Function to handle NaN values and escape single quotes
def clean_value(value):
    if pd.isna(value):  # Check if the value is NaN
        return None  # Use None to insert NULL into the database
    return str(value).replace("'", "''")  # Convert to string and escape single quotes

# Insert data into the Supabase database without specifying 'ethno_id'
for _, row in df_filtered.iterrows():
    data = {
        'activity': clean_value(row['activity']),
        'genus': clean_value(row['genus']),
        'species': clean_value(row['species']),
        'family': clean_value(row['family']),
        'common_name': clean_value(row['common_name']),
        'country': clean_value(row['country']),
        'reference': clean_value(row['reference']),
        'long_reference': clean_value(row['long_reference']),
        'taxon': clean_value(row['taxon']),
        'tax_author': clean_value(row['tax_author']),
        'effective': clean_value(row['effective']),
        'created': clean_value(row['created']),
        'modified': clean_value(row['modified'])
    }
    
    # Insert data using Supabase client
    supabase_client.table('ethnobotanical_uses').insert(data).execute()