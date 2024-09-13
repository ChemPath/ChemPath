import requests
from bs4 import BeautifulSoup
import aiohttp
import asyncio
import time
from flask import Flask
from flask_sqlalchemy import SQLAlchemy
import logging
from requests.adapters import HTTPAdapter
from urllib3.util.retry import Retry
import re
import xml.etree.ElementTree as ET

app = Flask(__name__)
app.config['SQLALCHEMY_DATABASE_URI'] = 'postgresql://postgres:Iamabundant28@localhost/chempath_database'
app.config['SQLALCHEMY_TRACK_MODIFICATIONS'] = False
db = SQLAlchemy(app)

logging.basicConfig(filename='app.log', level=logging.DEBUG)
logger = logging.getLogger(__name__)

# Print and log when script starts
print("Starting script...")
logger.info("Starting script...")

### Define the PlantCompound model correctly
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
    pubmed_citations = db.Column(db.Integer)
    google_scholar_citations = db.Column(db.Integer)
    related_articles = db.Column(db.Text)

### Test database connection function using SQLAlchemy ORM
def test_database_connection():
    try:
        with app.app_context():
            # Test database connection with a simple query
            if db.session.query(PlantCompound).first() is not None:
                print("Database connection successful.")
                logger.info("Database connection successful.")
            else:
                print("Database connection successful but no data found.")
                logger.info("Database connection successful but no data found.")
    except Exception as e:
        print(f"Database connection failed: {e}")
        logger.error(f"Database connection failed: {e}")

### PhytoHub Data Loading Functions

def fetch_phytohub_data(phytohub_id):
    print(f"Fetching data for PhytoHub ID: {phytohub_id}")
    logger.info(f"Fetching data for PhytoHub ID: {phytohub_id}")
    url = f"https://phytohub.eu/entries/{phytohub_id}"
    try:
        response = requests.get(url)
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

    return compound_data

def load_phytohub_data():
    with app.app_context():
        session = db.session

        for i in range(1, 1001):
            phytohub_id = f"PHUB{i:06d}"
            compound_data = fetch_phytohub_data(phytohub_id)
            if compound_data and compound_data['name']:
                existing_compound = PlantCompound.query.filter_by(phytohub_id=phytohub_id).first()
                if existing_compound:
                    logger.info(f"Compound {phytohub_id} already exists in the database.")
                    continue
                try:
                    compound = PlantCompound(**compound_data)
                    session.add(compound)
                    session.commit()
                    logger.info(f"Committed compound {phytohub_id} to database")
                except Exception as e:
                    logger.error(f"Error adding compound {phytohub_id} to database: {str(e)}")
                    session.rollback()
            else:
                logger.warning(f"No data found for compound {phytohub_id}")
            time.sleep(0.5)

        logger.info("Completed loading PhytoHub data")

### PubMed Data Handling Functions

def fetch_pubmed_data(compound_name):
    print(f"Fetching PubMed data for compound: {compound_name}")
    logger.info(f"Fetching PubMed data for compound: {compound_name}")
    try:
        base_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"
        search_url = f"{base_url}esearch.fcgi"
        params = {
            'db': 'pubmed',
            'term': f"{compound_name} AND plant",
            'retmax': 100,
            'usehistory': 'y',
        }
        response = requests.get(search_url, params=params)
        response.raise_for_status()
        root = ET.fromstring(response.content)
        
        count = int(root.find(".//Count").text)
        query_key = root.find(".//QueryKey").text
        web_env = root.find(".//WebEnv").text
        
        fetch_url = f"{base_url}efetch.fcgi"
        fetch_params = {
            'db': 'pubmed',
            'query_key': query_key,
            'WebEnv': web_env,
            'retmode': 'xml',
            'retmax': 10,
        }
        fetch_response = requests.get(fetch_url, params=fetch_params)
        fetch_response.raise_for_status()
        fetch_root = ET.fromstring(fetch_response.content)
        
        related_articles = []
        for article in fetch_root.findall(".//PubmedArticle"):
            title_elem = article.find(".//ArticleTitle")
            pmid_elem = article.find(".//PMID")
            title = title_elem.text if title_elem is not None else "No Title"
            pmid = pmid_elem.text if pmid_elem is not None else "No PMID"
            related_articles.append(f"{title} (PMID: {pmid})")
        
        return count, related_articles
    except Exception as e:
        logger.error(f"Error fetching PubMed data for {compound_name}: {str(e)}")
        return 0, []

def update_with_pubmed_data():
    print("Starting PubMed data update...")
    logger.info("Starting PubMed data update...")
    with app.app_context():
        # Update only if both pubmed_citations and related_articles are missing
        compounds = PlantCompound.query.filter(
            PlantCompound.pubmed_citations.is_(None), 
            PlantCompound.related_articles.is_(None)
        ).all()
        for compound in compounds:
            try:
                pubmed_count, related_articles = fetch_pubmed_data(compound.name)
                compound.pubmed_citations = pubmed_count
                compound.related_articles = "\n".join(related_articles[:10])
                db.session.commit()
                logger.info(f"Updated PubMed data for {compound.name}")
            except Exception as e:
                logger.error(f"Error updating PubMed data for {compound.name}: {str(e)}")
                db.session.rollback()

### Main Workflow

def main():
    # Test the database connection
    test_database_connection()

    print("Loading PhytoHub data...")
    logger.info("Loading PhytoHub data...")
    load_phytohub_data()

    print("Forcing reload of PubMed data...")
    logger.info("Forcing reload of PubMed data...")
    update_with_pubmed_data()

if __name__ == "__main__":
    main()
