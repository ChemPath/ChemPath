import time
import psycopg2
import requests
import logging
from psycopg2 import sql

# Logging configuration
logging.basicConfig(filename='missing_pubchem_ids.log', level=logging.INFO)

# Retry logic
MAX_RETRIES = 3
RETRY_DELAY = 5  # seconds

# Database connection details
DB_NAME = 'postgres'  # Replace with your actual database name
DB_USER = 'postgres.lqwxjiijnhefterkeety'       # Replace with your actual username
DB_PASSWORD = 'Iamabundant1228!'   # Replace with your actual password
DB_HOST = 'aws-0-us-west-1.pooler.supabase.com'           # Replace with your actual host (e.g., 'localhost' or Supabase host)
DB_PORT = '6543'           # Replace with your actual port (e.g., '5432')

def main():
    try:
        # Establish connection
        connection = psycopg2.connect(
            database=DB_NAME,
            user=DB_USER,
            password=DB_PASSWORD,
            host=DB_HOST,
            port=DB_PORT
        )
        cursor = connection.cursor()
        # Your database operations here
        print("Database connection established successfully.")

    except Exception as error:
        print(f"Database connection error: {error}")

    finally:
        if connection:
            cursor.close()
            connection.close()
            print("PostgreSQL connection is closed")

if __name__ == "__main__":
    main()

def fetch_pubchem_id(ingredient_name):
    try:
        # Using PubChem PUG REST API to search by compound name
        base_url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{ingredient_name}/cids/TXT"
        response = requests.get(base_url)
        
        if response.status_code == 200:
            pubchem_id = response.text.strip()
            if pubchem_id:
                print(f"PubChem ID for {ingredient_name}: {pubchem_id}")
                return pubchem_id
            else:
                print(f"No PubChem ID found for {ingredient_name}")
                return None
        elif response.status_code == 404:
            print(f"PubChem ID not found for {ingredient_name} (404 Not Found). Skipping.")
            return None
        else:
            print(f"Error fetching data from PubChem for {ingredient_name}: {response.status_code}")
            return None
    
    except Exception as e:
        print(f"Error fetching PubChem ID for {ingredient_name}: {e}")
        return None

def main():
    try:
        # Establish connection
        connection = psycopg2.connect(
            database=DB_NAME,
            user=DB_USER,
            password=DB_PASSWORD,
            host=DB_HOST,
            port=DB_PORT
        )
        cursor = connection.cursor()
        print("Database connection established successfully.")
        
        # Fetch the ingredient names that are missing PubChem IDs from the database
        cursor.execute("SELECT name FROM ingredients WHERE pubchem_id IS NULL")
        ingredients = cursor.fetchall()

        for ingredient_tuple in ingredients:
            ingredient_name = ingredient_tuple[0]
            print(f"Fetching PubChem ID for ingredient: {ingredient_name}")
            
            # Fetch PubChem ID
            pubchem_id = fetch_pubchem_id(ingredient_name)
            
            # If a PubChem ID is found, update the database
            if pubchem_id:
                cursor.execute(
                    "UPDATE ingredients SET pubchem_id = %s WHERE name = %s", 
                    (pubchem_id, ingredient_name)
                )
                connection.commit()
                print(f"PubChem ID {pubchem_id} updated for {ingredient_name} in the database.")
            else:
                print(f"Skipping {ingredient_name} due to missing PubChem ID.")

    except Exception as error:
        print(f"Database connection error: {error}")

    finally:
        if connection:
            cursor.close()
            connection.close()
            print("PostgreSQL connection is closed")

if __name__ == "__main__":
    main()