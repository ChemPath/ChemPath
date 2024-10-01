import google.generativeai as genai
from loguru import logger
from supabase import create_client, Client
import time
from tenacity import retry, stop_after_attempt, wait_exponential

# Supabase client setup
url = 'https://lqwxjiijnhefterkeety.supabase.co'
key = 'eyJhbGciOiJIUzI1NiIsInR5cCI6IkpXVCJ9.eyJpc3MiOiJzdXBhYmFzZSIsInJlZiI6Imxxd3hqaWlqbmhlZnRlcmtlZXR5Iiwicm9sZSI6InNlcnZpY2Vfcm9sZSIsImlhdCI6MTcyNjU5NTk3OSwiZXhwIjoyMDQyMTcxOTc5fQ.bmEG-1VU5MgU9TAWsXS0FEwVpH3qgc49dyUqXxvK3Kk'
supabase: Client = create_client(url, key)

# Gemini AI setup
genai.configure(api_key='AIzaSyC6rBa5G4jex1Uik0Uzhq0OYureqJW6PwY')

def fetch_compounds_without_descriptions(offset, limit):
    try:
        query = supabase.table("compounds").select("id", "cmpdname", "iupacname", "molecular_formula").is_("description", "null").range(offset, offset + limit - 1)
        return query.execute().data
    except Exception as e:
        logger.error(f"Error fetching data from Supabase: {e}")
        return []

@retry(stop=stop_after_attempt(3), wait=wait_exponential(multiplier=1, min=4, max=10))
def generate_description(compound_name, molecular_formula):
    try:
        model = genai.GenerativeModel('gemini-pro')
        prompt = f"Provide a concise description for the compound: {compound_name}. Its molecular formula is {molecular_formula}. Focus on its structure, properties, and potential uses."
        response = model.generate_content(prompt)
        if response.parts:
            return response.parts[0].text.strip()
        else:
            return None
    except Exception as e:
        logger.error(f"Error generating description: {e}")
        raise

def update_database_description(compound_id, new_description):
    try:
        response = supabase.table("compounds").update({"description": new_description}).eq("id", compound_id).execute()
        return len(response.data) > 0
    except Exception as e:
        logger.error(f"Error updating database for compound ID {compound_id}: {e}")
        return False

def update_descriptions(limit=500):
    offset = 0
    total_updated = 0

    while total_updated < limit:
        compounds = fetch_compounds_without_descriptions(offset, limit - total_updated)
        if not compounds:
            break

        for compound in compounds:
            compound_id = compound['id']
            compound_name = compound.get('cmpdname') or compound.get('iupacname') or "Unnamed compound"
            molecular_formula = compound.get('molecular_formula') or "Unknown"
            
            new_description = generate_description(compound_name, molecular_formula)

            if new_description and update_database_description(compound_id, new_description):
                total_updated += 1
                logger.info(f"Updated description for compound ID {compound_id}")
            else:
                logger.warning(f"Failed to update description for compound ID {compound_id}")
            
            if total_updated >= limit:
                break
            
            time.sleep(1)  # To avoid rate limiting

        offset += len(compounds)
        logger.info(f"Processed batch. Total updated so far: {total_updated}")

    logger.info(f"Completed updating {total_updated} descriptions.")

if __name__ == "__main__":
    update_descriptions()