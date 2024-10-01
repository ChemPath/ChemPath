import pandas as pd
from supabase import create_client, Client
import numpy as np

# Initialize Supabase client
url = 'https://lqwxjiijnhefterkeety.supabase.co'
key = 'eyJhbGciOiJIUzI1NiIsInR5cCI6IkpXVCJ9.eyJpc3MiOiJzdXBhYmFzZSIsInJlZiI6Imxxd3hqaWlqbmhlZnRlcmtlZXR5Iiwicm9sZSI6InNlcnZpY2Vfcm9sZSIsImlhdCI6MTcyNjU5NTk3OSwiZXhwIjoyMDQyMTcxOTc5fQ.bmEG-1VU5MgU9TAWsXS0FEwVpH3qgc49dyUqXxvK3Kk'
supabase: Client = create_client(url, key)

# Load your data
file_path = 'C:/Users/Dr. Contessa Petrini/OneDrive/Desktop/Research Private/PubChem_compound_cache_kdY38CyBST1-E8EKQ3KIIhWR7PGJHdN5qVzINbJN2jSyVOY.csv'
df = pd.read_csv(file_path, low_memory=False)

# Strip leading/trailing spaces from column names to avoid KeyError
df.columns = df.columns.str.strip()

# Replace 'NaN' values with None for non-numeric fields and 0 for numeric fields
df = df.apply(lambda x: x.fillna('') if x.dtype == 'object' else x.fillna(0))

# Ensure all numeric columns are properly cast
numeric_columns = ['mw', 'exactmass', 'monoisotopicmass']
df[numeric_columns] = df[numeric_columns].apply(pd.to_numeric, errors='coerce')

# Specify the columns to use
columns_to_use = [
    'cid', 'cmpdname', 'cmpdsynonym', 'mw', 'mf', 'polararea', 'complexity', 'xlogp', 
    'heavycnt', 'hbonddonor', 'hbondacc', 'rotbonds', 'inchi', 'isosmiles', 
    'canonicalsmiles', 'inchikey', 'iupacname', 'exactmass', 'monoisotopicmass', 
    'charge', 'covalentunitcnt', 'isotopeatomcnt', 'totalatomstereocnt', 
    'definedatomstereocnt', 'undefinedatomstereocnt', 'totalbondstereocnt', 
    'definedbondstereocnt', 'undefinedbondstereocnt', 'pclidcnt', 'gpidcnt', 
    'gpfamilycnt', 'neighbortype', 'meshheadings', 'annothits', 'annothitcnt', 
    'aids', 'cidcdate', 'sidsrcname', 'depcatg', 'annotation'
]

# Filter the DataFrame
df_filtered = df[columns_to_use]

# Convert the DataFrame to a list of dictionaries for upsert
batch_size = 1000

# Insert data into Supabase in batches
for i in range(0, len(df_filtered), batch_size):
    batch_data = df_filtered.iloc[i:i + batch_size].to_dict(orient="records")
    try:
        # Insert data into the 'compounds' table
        response = supabase.table('compounds').upsert(batch_data).execute()
        print(f"Inserted batch {i // batch_size + 1}: {response}")
    except Exception as e:
        print(f"Error inserting batch {i // batch_size + 1}: {e}")