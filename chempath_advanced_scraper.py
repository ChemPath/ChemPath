import asyncio
import aiohttp
from bs4 import BeautifulSoup
from database_operations import create_engine_and_session, get_compound, update_compound
import logging
from database_operations import add_compound
from database_operations import Compound


logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

async def fetch_pubchem_data(session, cid):
    url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid}/JSON"
    async with session.get(url) as response:
        if response.status == 200:
            return await response.json()
        return None

async def fetch_wikipedia_data(session, name):
    url = f"https://en.wikipedia.org/api/rest_v1/page/summary/{name}"
    async with session.get(url) as response:
        if response.status == 200:
            return await response.json()
        return None

async def fetch_chemspider_data(session, name):
    url = f"http://www.chemspider.com/Search.aspx?q={name}"
    async with session.get(url) as response:
        if response.status == 200:
            html = await response.text()
            soup = BeautifulSoup(html, 'html.parser')
            # Extract relevant data from ChemSpider
            # This is a placeholder and should be replaced with actual parsing logic
            return {"chemspider_data": "placeholder"}
        return None

async def enrich_compound_data(db_session, compounds):
    async with aiohttp.ClientSession() as session:
        for compound in compounds:
            pubchem_data = await fetch_pubchem_data(session, compound.pubchem_cid)
            wiki_data = await fetch_wikipedia_data(session, compound.name)
            chemspider_data = await fetch_chemspider_data(session, compound.name)
            
            if pubchem_data and wiki_data and chemspider_data:
                try:
                    compound_dict = pubchem_data['PC_Compounds'][0]
                    wiki_summary = wiki_data.get('extract', '')
                    
                    update_data = {
                        'synonyms': ';'.join(compound_dict.get('synonyms', [])),
                        'biological_activity': compound_dict.get('pharmacology', {}).get('description'),
                        'traditional_use': wiki_summary,
                        'molecular_weight': compound_dict.get('molecular_weight'),
                        'iupac_name': compound_dict.get('iupac_name'),
                        'inchi': compound_dict.get('inchi'),
                        'logp': compound_dict.get('xlogp'),
                        'polararea': compound_dict.get('topological_polar_surface_area'),
                        'complexity': compound_dict.get('complexity'),
                        'heavycnt': compound_dict.get('heavy_atom_count'),
                        'hbonddonor': compound_dict.get('h_bond_donor_count'),
                        'hbondacc': compound_dict.get('h_bond_acceptor_count'),
                        'rotbonds': compound_dict.get('rotatable_bond_count'),
                        'exactmass': compound_dict.get('exact_mass'),
                        'monoisotopicmass': compound_dict.get('monoisotopic_mass'),
                        # Add any additional fields from ChemSpider data
                        'chemspider_data': chemspider_data.get('chemspider_data')
                    }
                    
                    update_compound(db_session, compound.id, update_data)
                    logger.info(f"Updated data for compound {compound.name} (ID: {compound.id})")
                except Exception as e:
                    logger.error(f"Error updating compound {compound.name} (ID: {compound.id}): {str(e)}")
            else:
                logger.warning(f"Could not fetch data for compound {compound.name} (ID: {compound.id})")

async def scrape_new_compounds(db_session, source_url):
    async with aiohttp.ClientSession() as session:
        async with session.get(source_url) as response:
            if response.status == 200:
                html = await response.text()
                soup = BeautifulSoup(html, 'html.parser')
                # Extract new compound data from the source
                # This is a placeholder and should be replaced with actual parsing logic
                new_compounds = [
                    {"name": "New Compound 1", "smiles": "CCO", "pubchem_cid": 702},
                    {"name": "New Compound 2", "smiles": "CC(=O)O", "pubchem_cid": 176}
                ]
                
                for compound_data in new_compounds:
                    add_compound(db_session, compound_data)
                    logger.info(f"Added new compound: {compound_data['name']}")
            else:
                logger.error(f"Failed to fetch data from {source_url}")

async def main():
    engine, db_session = create_engine_and_session()
    
    # Enrich existing compounds
    compounds = db_session.query(Compound).all()
    await enrich_compound_data(db_session, compounds)
    
    # Scrape for new compounds
    await scrape_new_compounds(db_session, "https://example.com/compound_list")
    
    db_session.close()

if __name__ == "__main__":
    asyncio.run(main())
