import requests
from bs4 import BeautifulSoup
from tqdm import tqdm
import logging
from flask import Flask
from flask_sqlalchemy import SQLAlchemy
import os

app = Flask(__name__)
db_uri = os.environ.get('DATABASE_URI', "postgresql://postgres:Iamabundant28@localhost/chempath_database")
app.config['SQLALCHEMY_DATABASE_URI'] = db_uri
app.config['SQLALCHEMY_TRACK_MODIFICATIONS'] = False
db = SQLAlchemy(app)

class PlantCompound(db.Model):
    __tablename__ = 'plant_compounds'

    id = db.Column(db.Integer, primary_key=True)
    phytohub_id = db.Column(db.String, unique=True, nullable=False)
    name = db.Column(db.String)
    systematic_name = db.Column(db.String)
    synonyms = db.Column(db.Text)
    cas_number = db.Column(db.String)
    average_mass = db.Column(db.Float)
    monoisotopic_mass = db.Column(db.Float)
    chemical_formula = db.Column(db.String)
    iupac_name = db.Column(db.String)
    inchi_key = db.Column(db.String)
    inchi_identifier = db.Column(db.Text)
    smiles = db.Column(db.String)
    solubility = db.Column(db.String)
    logs = db.Column(db.String)
    logp = db.Column(db.String)
    hydrogen_acceptors = db.Column(db.String)
    hydrogen_donors = db.Column(db.String)
    rotatable_bond_count = db.Column(db.String)
    polar_surface_area = db.Column(db.String)
    refractivity = db.Column(db.String)
    polarizability = db.Column(db.String)
    formal_charge = db.Column(db.String)
    physiological_charge = db.Column(db.String)
    pka_strongest_basic = db.Column(db.String)
    pka_strongest_acidic = db.Column(db.String)
    number_of_rings = db.Column(db.String)
    rule_of_five = db.Column(db.String)
    bioavailability = db.Column(db.String)
    ghose_filter = db.Column(db.String)
    veber_rule = db.Column(db.String)
    mddr_like_rule = db.Column(db.String)
    chebi = db.Column(db.String)
    knapsack = db.Column(db.String)
    pubchem = db.Column(db.String)
    taxonomy_family = db.Column(db.String)
    taxonomy_class = db.Column(db.String)
    taxonomy_subclass = db.Column(db.String)
    food_sources = db.Column(db.Text)
    biomarker_roles = db.Column(db.Text)
    metabolism = db.Column(db.Text)
    inter_individual_variations = db.Column(db.Text)

def extract_single_compound(url):
    try:
        response = requests.get(url, timeout=30)
        response.raise_for_status()
        soup = BeautifulSoup(response.text, 'html.parser')
        
        compound = {}
        
        # Extract basic information
        dl = soup.find('dl', class_='dl-horizontal')
        if dl:
            for dt, dd in zip(dl.find_all('dt'), dl.find_all('dd')):
                key = dt.text.strip()
                value = dd.text.strip()
                if value and value != 'Not Available':
                    compound[key.lower().replace(' ', '_')] = value
        
        # Extract calculated properties
        properties_section = soup.find('section', id='cs')
        if properties_section:
            for dt, dd in zip(properties_section.find_all('dt'), properties_section.find_all('dd')):
                key = dt.text.strip()
                value = dd.text.strip()
                if value and value != 'Not Available':
                    compound[key.lower().replace(' ', '_')] = value
        
        # Extract classification information
        classification_section = soup.find('section', id='classification')
        if classification_section:
            for dt, dd in zip(classification_section.find_all('dt'), classification_section.find_all('dd')):
                key = dt.text.strip()
                value = dd.text.strip()
                if value and value != 'Not Available':
                    compound[key.lower().replace(' ', '_')] = value
        
        # Extract additional information
        additional_info = soup.find('section', id='additional-info')
        if additional_info:
            for dt, dd in zip(additional_info.find_all('dt'), additional_info.find_all('dd')):
                key = dt.text.strip()
                value = dd.text.strip()
                if value and value != 'Not Available':
                    compound[key.lower().replace(' ', '_')] = value
        
        # Ensure phytohub_id is present
        if 'phytohub_id' not in compound:
            logging.warning(f"No PhytoHub ID found for compound at {url}")
            return None
        
        logging.info(f"Successfully extracted data for compound {compound.get('phytohub_id', 'Unknown')}")
        return compound
    except requests.RequestException as e:
        logging.error(f"Error fetching data from {url}: {str(e)}")
    except Exception as e:
        logging.error(f"Unexpected error processing {url}: {str(e)}")
    return None

def extract_compound_data(url):
    logging.info(f"Extracting compound data from {url}")
    try:
        response = requests.get(url, timeout=30)
        response.raise_for_status()
        soup = BeautifulSoup(response.text, 'html.parser')
        table = soup.find('table', id='entries')
        if not table:
            logging.error("No table with id 'entries' found on the page.")
            return []
        compounds = []
        rows = table.find_all('tr')[1:1001]  # Get the first 1000 compounds
        
        for row in tqdm(rows, desc="Extracting compounds"):
            try:
                cols = row.find_all('td')
                if len(cols) < 4:
                    logging.warning(f"Row doesn't have enough columns: {row}")
                    continue
                phytohub_id = cols[3].text.strip()
                compound_url = f"https://phytohub.eu/entries/{phytohub_id}"
                compound_data = extract_single_compound(compound_url)
                if compound_data:
                    compounds.append(compound_data)
                    logging.debug(f"Extracted data for compound {phytohub_id}")
                else:
                    logging.warning(f"No data extracted for compound at {compound_url}")
            except Exception as e:
                logging.error(f"Error extracting data for compound {phytohub_id}: {str(e)}")
        
        logging.info(f"Extracted data for {len(compounds)} compounds")
        return compounds
    except requests.RequestException as e:
        logging.error(f"Error fetching data from {url}: {str(e)}")
        return []

def insert_compounds(compounds):
    successful_inserts = 0
    with db.session.begin_nested():
        for compound in compounds:
            try:
                if not compound.get('phytohub_id'):
                    logging.warning(f"Skipping compound with null phytohub_id: {compound}")
                    continue
                
                existing_compound = PlantCompound.query.filter_by(phytohub_id=compound['phytohub_id']).first()
                if existing_compound:
                    for key, value in compound.items():
                        if hasattr(existing_compound, key):
                            setattr(existing_compound, key, value)
                else:
                    new_compound = PlantCompound(**{k: v for k, v in compound.items() if hasattr(PlantCompound, k)})
                    db.session.add(new_compound)
                logging.info(f"Inserted/Updated compound: {compound['phytohub_id']}")
                successful_inserts += 1
            except Exception as e:
                logging.error(f"Error inserting compound {compound.get('phytohub_id', 'Unknown')}: {str(e)}")
        db.session.commit()
    return successful_inserts

def populate_database(app):
    logging.info("Starting database population from PhytoHub")
    url = "https://phytohub.eu/entries"
    compounds = extract_compound_data(url)
   
    logging.info(f"Extracted {len(compounds)} compounds from PhytoHub")
   
    total_inserted = 0
    chunks = [compounds[i:i+10] for i in range(0, len(compounds), 10)]
   
    with app.app_context():
        for chunk in tqdm(chunks, desc="Inserting compound chunks"):
            total_inserted += insert_compounds(chunk)
   
    logging.info(f"Database population complete. Inserted {total_inserted} compounds.")
    return total_inserted

if __name__ == '__main__':
    logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
    try:
        logging.info("Starting ChemPath database initialization")
        with app.app_context():
            db.create_all()
        logging.info("Tables created successfully")
       
        logging.info("Starting automated data collection")
        populate_database(app)
        logging.info("Automated data collection completed")
       
        app.run(debug=True)
    except Exception as e:
        logging.error(f"Error: {str(e)}")
