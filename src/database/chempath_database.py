import requests
from bs4 import BeautifulSoup
import time
from flask import Flask
from flask_sqlalchemy import SQLAlchemy
import os
import logging
from requests.adapters import HTTPAdapter
from urllib3.util.retry import Retry
from flask import Flask
from flask_sqlalchemy import SQLAlchemy
from flask_migrate import Migrate
import re


app = Flask(__name__)
app.config['SQLALCHEMY_DATABASE_URI'] = 'postgresql://postgres:Iamabundant28@localhost/chempath_database'
app.config['SQLALCHEMY_TRACK_MODIFICATIONS'] = False
db = SQLAlchemy(app)
migrate = Migrate(app, db)

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

class PlantCompound(db.Model):
    __tablename__ = 'plant_compounds'
    id = db.Column(db.Integer, primary_key=True)
    phytohub_id = db.Column(db.String, unique=True, nullable=False)
    name = db.Column(db.String)
    systematic_name = db.Column(db.String)
    synonyms = db.Column(db.String)
    cas_number = db.Column(db.String)
    average_mass = db.Column(db.Float)
    monoisotopic_mass = db.Column(db.Float)
    chemical_formula = db.Column(db.String)
    iupac_name = db.Column(db.String)
    inchi_key = db.Column(db.String)
    inchi_identifier = db.Column(db.String)
    smiles = db.Column(db.String)
    solubility = db.Column(db.Float)
    log_s = db.Column(db.Float)
    log_p = db.Column(db.Float)

def fetch_compound_data(phytohub_id):
    url = f"https://phytohub.eu/entries/{phytohub_id}"
    session = requests_retry_session()
    
    try:
        response = session.get(url, timeout=10)
        response.raise_for_status()
    except requests.exceptions.RequestException as e:
        logger.error(f"Failed to fetch data for {phytohub_id}: {str(e)}")
        return None

    soup = BeautifulSoup(response.content, 'html.parser')
    compound_data = {'phytohub_id': phytohub_id}

    # Extract name from the page title
    title_element = soup.find('title')
    if title_element:
        title_text = title_element.text.strip()
        name_match = re.search(r'Showing entry for (.+)', title_text)
        if name_match:
            compound_data['name'] = name_match.group(1).strip()

    # Extract other fields
    fields = {
        'systematic_name': 'Systematic Name',
        'synonyms': 'Synonyms',
        'cas_number': 'CAS Number',
        'average_mass': 'Average Mass',
        'monoisotopic_mass': 'Monoisotopic Mass',
        'chemical_formula': 'Chemical Formula',
        'iupac_name': 'IUPAC Name',
        'inchi_key': 'InChI Key',
        'inchi_identifier': 'InChI Identifier',
        'smiles': 'SMILES',
        'solubility': 'Solubility (ALOGPS)',
        'log_s': 'LogS (ALOGPS)',
        'log_p': 'LogP (ALOGPS)'
    }

    for field, label in fields.items():
        element = soup.find('dt', string=label)
        if element and element.find_next('dd'):
            value = element.find_next('dd').text.strip()
            if value.lower() == 'not available':
                value = None
            elif field in ['average_mass', 'monoisotopic_mass', 'solubility', 'log_s', 'log_p']:
                try:
                    value = float(value)
                except ValueError:
                    value = None
            compound_data[field] = value
        else:
            compound_data[field] = None

    logger.info(f"Compound data for {phytohub_id}: {compound_data}")
    return compound_data

def reset_database():
    with app.app_context():
        db.drop_all()
        db.create_all()
        logger.info("Database reset complete")

def requests_retry_session(retries=3, backoff_factor=0.3, status_forcelist=(500, 502, 504)):
    session = requests.Session()
    retry = Retry(total=retries, read=retries, connect=retries,
                  backoff_factor=backoff_factor, status_forcelist=status_forcelist)
    adapter = HTTPAdapter(max_retries=retry)
    session.mount('http://', adapter)
    session.mount('https://', adapter)
    return session

def populate_database():
    with app.app_context():
        db.create_all()
        
        session = db.session
        
        for i in range(1, 1001):
            phytohub_id = f"PHUB{i:06d}"
            compound_data = fetch_compound_data(phytohub_id)
            if compound_data and compound_data['name']:
                try:
                    compound = PlantCompound(**compound_data)
                    session.add(compound)
                    session.flush()
                    logger.info(f"Added compound {phytohub_id} to session")
                except Exception as e:
                    logger.error(f"Error adding compound {phytohub_id} to session: {str(e)}")
                    session.rollback()
                
                if i % 100 == 0:
                    try:
                        session.commit()
                        logger.info(f"Committed batch of compounds to database")
                    except Exception as e:
                        logger.error(f"Error committing batch to database: {str(e)}")
                        session.rollback()
            
            time.sleep(1)
        
        try:
            session.commit()
            logger.info("Final commit of compounds to database")
        except Exception as e:
            logger.error(f"Error in final commit to database: {str(e)}")
            session.rollback()

    verify_data_insertion()

def verify_data_insertion():
    with app.app_context():
        count = PlantCompound.query.count()
        logger.info(f"Total compounds in database: {count}")

def main():
    reset_database()
    populate_database()
    verify_data_insertion()

if __name__ == '__main__':
    main()
