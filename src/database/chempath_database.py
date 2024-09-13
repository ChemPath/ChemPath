import requests
from bs4 import BeautifulSoup
import time
from flask import Flask
from flask_sqlalchemy import SQLAlchemy
import logging
from requests.adapters import HTTPAdapter
from urllib3.util.retry import Retry
import re
import xml.etree.ElementTree as ET
from scholarly import scholarly

app = Flask(__name__)
app.config['SQLALCHEMY_DATABASE_URI'] = 'postgresql://postgres:Iamabundant28@localhost/chempath_database'
app.config['SQLALCHEMY_TRACK_MODIFICATIONS'] = False
db = SQLAlchemy(app)

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
    pubmed_citations = db.Column(db.Integer)
    google_scholar_citations = db.Column(db.Integer)
    related_articles = db.Column(db.Text)

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
                    value = float(value.replace(',', ''))
                except ValueError:
                    value = None
            compound_data[field] = value
        else:
            compound_data[field] = None

    logger.info(f"Compound data for {phytohub_id}: {compound_data}")
    return compound_data

def load_phytohub_data():
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
                    session.commit()
                    logger.info(f"Committed compound {phytohub_id} to database")
                except Exception as e:
                    logger.error(f"Error adding compound {phytohub_id} to database: {str(e)}")
                    session.rollback()
            time.sleep(0.5)

        logger.info("Completed loading PhytoHub data")

def fetch_pubmed_data(compound_name):
    try:
        base_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"
        search_url = f"{base_url}esearch.fcgi"
        params = {
            'db': 'pubmed',
            'term': f"{compound_name} AND plant",
            'retmax': 100,
            'usehistory': 'y',
            # 'api_key': 'YOUR_API_KEY',  # Optional: Use your NCBI API key
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
            'retmax': 100,
            # 'api_key': 'YOUR_API_KEY',  # Optional: Use your NCBI API key
        }
        fetch_response = requests.get(fetch_url, params=fetch_params)
        fetch_response.raise_for_status()
        fetch_root = ET.fromstring(fetch_response.content)
        
        related_articles = []
        for article in fetch_root.findall(".//PubmedArticle"):
            title = article.find(".//ArticleTitle").text
            pmid = article.find(".//PMID").text
            related_articles.append(f"{title} (PMID: {pmid})")
        
        return count, related_articles
    except Exception as e:
        logger.error(f"Error fetching PubMed data for {compound_name}: {str(e)}")
        return 0, []

def fetch_google_scholar_data(compound_name):
    try:
        # Search for the compound on Google Scholar
        search_query = scholarly.search_pubs(compound_name + " plant")
        total_citations = 0
        for _ in range(10):  # Get citations from the first 10 publications
            try:
                publication = next(search_query)
                if 'num_citations' in publication.bib:
                    total_citations += publication.bib['num_citations']
            except StopIteration:
                break
        return total_citations
    except Exception as e:
        logger.error(f"Error fetching Google Scholar data for {compound_name}: {str(e)}")
        return 0

def update_compound_with_literature_data(compound):
    pubmed_count, related_articles = fetch_pubmed_data(compound.name)
    google_scholar_count = fetch_google_scholar_data(compound.name)
    
    compound.pubmed_citations = pubmed_count
    compound.google_scholar_citations = google_scholar_count
    compound.related_articles = "\n".join(related_articles[:10])

def update_with_google_scholar_data():
    with app.app_context():
        compounds = PlantCompound.query.filter(PlantCompound.google_scholar_citations.is_(None)).all()
        for compound in compounds:
            try:
                google_scholar_count = fetch_google_scholar_data(compound.name)
                compound.google_scholar_citations = google_scholar_count
                db.session.commit()
                logger.info(f"Updated Google Scholar data for {compound.name}")
            except Exception as e:
                logger.error(f"Error updating Google Scholar data for {compound.name}: {str(e)}")
                db.session.rollback()
            time.sleep(2)  # Delay to avoid overwhelming the server

def check_google_scholar_data_loaded():
    with app.app_context():
        return PlantCompound.query.filter(PlantCompound.google_scholar_citations.isnot(None)).count() > 0

def requests_retry_session(retries=3, backoff_factor=0.3, status_forcelist=(500, 502, 504)):
    session = requests.Session()
    retry = Retry(total=retries, read=retries, connect=retries,
                  backoff_factor=backoff_factor, status_forcelist=status_forcelist)
    adapter = HTTPAdapter(max_retries=retry)
    session.mount('http://', adapter)
    session.mount('https://', adapter)
    return session

def verify_data_insertion():
    with app.app_context():
        count = PlantCompound.query.count()
        logger.info(f"Total compounds in database: {count}")

def check_phytohub_data_loaded():
    with app.app_context():
        return PlantCompound.query.count() > 0

def main():
    if not check_phytohub_data_loaded():
        # Do not reset the database
        load_phytohub_data()
        verify_data_insertion()
    else:
        logger.info("PhytoHub data already loaded. Skipping data loading.")
    
    if not check_google_scholar_data_loaded():
        update_with_google_scholar_data()
    else:
        logger.info("Google Scholar data already loaded. Skipping data update.")

if __name__ == "__main__":
    main()
