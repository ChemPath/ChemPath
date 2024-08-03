import requests
from bs4 import BeautifulSoup
from src.database.chempath_database import create_connection, insert_compound
import logging
import asyncio
import aiohttp

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

async def fetch_phytochemdb_data(session):
    url = 'https://phytochemdb.org/search/'
    async with session.get(url, params={'query': ''}) as response:
        soup = BeautifulSoup(await response.text(), 'html.parser')
        compounds = []
        for row in soup.find_all('tr'):
            cols = row.find_all('td')
            if len(cols) > 0:
                compound = {
                    'name': cols[0].text.strip(),
                    'smiles': cols[1].text.strip(),
                    'molecular_weight': float(cols[2].text.strip()) if cols[2].text.strip() else None,
                    'logp': float(cols[3].text.strip()) if cols[3].text.strip() else None,
                    'plant_source': cols[4].text.strip(),
                    'biological_activity': cols[5].text.strip(),
                    'traditional_use': cols[6].text.strip()
                }
                compounds.append(compound)
        return compounds

async def fetch_phytohub_data(session):
    url = 'https://phytohub.net/search/'
    async with session.get(url, params={'query': ''}) as response:
        soup = BeautifulSoup(await response.text(), 'html.parser')
        compounds = []
        for row in soup.find_all('tr'):
            cols = row.find_all('td')
            if len(cols) > 0:
                compound = {
                    'name': cols[0].text.strip(),
                    'smiles': cols[1].text.strip(),
                    'molecular_weight': float(cols[2].text.strip()) if cols[2].text.strip() else None,
                    'logp': float(cols[3].text.strip()) if cols[3].text.strip() else None,
                    'plant_source': cols[4].text.strip(),
                    'biological_activity': cols[5].text.strip(),
                    'traditional_use': cols[6].text.strip()
                }
                compounds.append(compound)
        return compounds

async def fetch_psc_db_data(session):
    url = 'https://psc-db.org/api/compounds'
    async with session.get(url) as response:
        compounds = await response.json()
        return [
            {
                'name': compound['name'],
                'smiles': compound['smiles'],
                'molecular_weight': compound.get('molecular_weight'),
                'logp': compound.get('logp'),
                'plant_source': compound.get('plant_source'),
                'biological_activity': compound.get('biological_activity'),
                'traditional_use': compound.get('traditional_use')
            }
            for compound in compounds
        ]

async def fetch_dr_dukes_data(session):
    url = 'https://phytochem.nal.usda.gov/phytochem/search'
    async with session.get(url, params={'query': ''}) as response:
        soup = BeautifulSoup(await response.text(), 'html.parser')
        compounds = []
        for row in soup.find_all('tr'):
            cols = row.find_all('td')
            if len(cols) > 0:
                compound = {
                    'name': cols[0].text.strip(),
                    'smiles': cols[1].text.strip(),
                    'molecular_weight': float(cols[2].text.strip()) if cols[2].text.strip() else None,
                    'logp': float(cols[3].text.strip()) if cols[3].text.strip() else None,
                    'plant_source': cols[4].text.strip(),
                    'biological_activity': cols[5].text.strip(),
                    'traditional_use': cols[6].text.strip()
                }
                compounds.append(compound)
        return compounds

async def process_and_insert_data(conn, data):
    for compound in data:
        try:
            await insert_compound(conn, compound)
        except Exception as e:
            logger.error(f"Error inserting compound {compound['name']}: {str(e)}")

async def populate_database():
    conn = await create_connection('chempath_database.db')
    
    databases = [
        ('PhytochemDB', fetch_phytochemdb_data),
        ('PhytoHub', fetch_phytohub_data),
        ('PSC-DB', fetch_psc_db_data),
        ("Dr. Duke's", fetch_dr_dukes_data)
    ]
    
    async with aiohttp.ClientSession() as session:
        for db_name, fetch_function in databases:
            try:
                logger.info(f"Fetching data from {db_name}...")
                data = await fetch_function(session)
                logger.info(f"Processing and inserting data from {db_name}...")
                await process_and_insert_data(conn, data)
                logger.info(f"Completed processing {db_name} data.")
            except Exception as e:
                logger.error(f"Error processing {db_name} data: {str(e)}")
    
    await conn.close()
    logger.info("Database population completed.")

if __name__ == "__main__":
    asyncio.run(populate_database())
