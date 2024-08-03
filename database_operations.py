from sqlalchemy import create_engine, Column, Integer, String, Float, Text
from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy.orm import sessionmaker
from sqlalchemy.types import TypeDecorator
import asyncio
import aiohttp
import logging

Base = declarative_base()
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

class CustomSmilesType(TypeDecorator):
    impl = Text

    def process_bind_param(self, value, dialect):
        # Add SMILES validation logic here
        return value

    def process_result_value(self, value, dialect):
        return value

class Compound(Base):
    __tablename__ = 'compounds'

    id = Column(Integer, primary_key=True)
    name = Column(String, nullable=False)
    smiles = Column(CustomSmilesType, nullable=False)
    molecular_weight = Column(Float)
    logp = Column(Float)
    plant_source = Column(String)
    biological_activity = Column(Text)
    traditional_use = Column(Text)
    pubchem_cid = Column(Integer)
    inchi = Column(Text)
    iupac_name = Column(String)
    synonyms = Column(Text)
    polararea = Column(Float)
    complexity = Column(Float)
    heavycnt = Column(Integer)
    hbonddonor = Column(Integer)
    hbondacc = Column(Integer)
    rotbonds = Column(Integer)
    exactmass = Column(Float)
    monoisotopicmass = Column(Float)

def create_engine_and_session(db_url='sqlite:///chempath_database.db'):
    engine = create_engine(db_url)
    Session = sessionmaker(bind=engine)
    return engine, Session()

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

async def enrich_compound_data(session, compounds, batch_size=10):
    async with aiohttp.ClientSession() as http_session:
        for i in range(0, len(compounds), batch_size):
            batch = compounds[i:i+batch_size]
            tasks = []
            for compound in batch:
                tasks.append(asyncio.create_task(fetch_pubchem_data(http_session, compound.pubchem_cid)))
                tasks.append(asyncio.create_task(fetch_wikipedia_data(http_session, compound.name)))
            
            results = await asyncio.gather(*tasks)
            
            for j in range(0, len(results), 2):
                compound = batch[j//2]
                pubchem_data, wiki_data = results[j], results[j+1]
                
                if pubchem_data and wiki_data:
                    try:
                        compound_dict = pubchem_data['PC_Compounds'][0]
                        wiki_summary = wiki_data.get('extract', '')
                        
                        compound.synonyms = ';'.join(compound_dict.get('synonyms', []))
                        compound.biological_activity = compound_dict.get('pharmacology', {}).get('description')
                        compound.traditional_use = wiki_summary
                        compound.molecular_weight = compound_dict.get('molecular_weight')
                        compound.iupac_name = compound_dict.get('iupac_name')
                        compound.inchi = compound_dict.get('inchi')
                        compound.logp = compound_dict.get('xlogp')
                        compound.polararea = compound_dict.get('topological_polar_surface_area')
                        compound.complexity = compound_dict.get('complexity')
                        compound.heavycnt = compound_dict.get('heavy_atom_count')
                        compound.hbonddonor = compound_dict.get('h_bond_donor_count')
                        compound.hbondacc = compound_dict.get('h_bond_acceptor_count')
                        compound.rotbonds = compound_dict.get('rotatable_bond_count')
                        compound.exactmass = compound_dict.get('exact_mass')
                        compound.monoisotopicmass = compound_dict.get('monoisotopic_mass')
                        
                        logger.info(f"Updated data for compound {compound.name} (ID: {compound.id})")
                    except Exception as e:
                        logger.error(f"Error updating compound {compound.name} (ID: {compound.id}): {str(e)}")
                else:
                    logger.warning(f"Could not fetch data for compound {compound.name} (ID: {compound.id})")
            
            session.commit()

def add_compound(session, compound_data):
    new_compound = Compound(**compound_data)
    session.add(new_compound)
    session.commit()
    return new_compound.id

def get_compound(session, compound_id):
    return session.query(Compound).filter(Compound.id == compound_id).first()

def search_compounds(session, query, search_field='name'):
    return session.query(Compound).filter(getattr(Compound, search_field).like(f"%{query}%")).all()

def update_compound(session, compound_id, update_data):
    compound = session.query(Compound).filter(Compound.id == compound_id).first()
    if compound:
        for key, value in update_data.items():
            setattr(compound, key, value)
        session.commit()
        return True
    return False

def delete_compound(session, compound_id):
    compound = session.query(Compound).filter(Compound.id == compound_id).first()
    if compound:
        session.delete(compound)
        session.commit()
        return True
    return False

def main():
    engine, session = create_engine_and_session()
    Base.metadata.create_all(engine)
    
    # Example usage
    compound_data = {
        'name': 'Quercetin',
        'smiles': 'O=C1C(O)=C(O)C(=O)C2=C1C=C(O)C(O)=C2O',
        'molecular_weight': 302.24,
        'logp': 1.54,
        'plant_source': 'Various fruits and vegetables',
        'biological_activity': 'Antioxidant, anti-inflammatory',
        'traditional_use': 'Traditional Chinese Medicine for cardiovascular health',
        'pubchem_cid': 5280343
    }

    compound_id = add_compound(session, compound_data)
    print(f"Added compound with ID: {compound_id}")

    compounds = search_compounds(session, 'Quercetin')
    asyncio.run(enrich_compound_data(session, compounds))

    session.close()

if __name__ == "__main__":
    main()
