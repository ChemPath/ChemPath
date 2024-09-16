import sys
import os

# Add the project root directory to PYTHONPATH
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '../../')))

import requests
import logging
import xml.etree.ElementTree as ET
from flask import Flask
from flask_sqlalchemy import SQLAlchemy
from sqlalchemy import Column, Integer, String, Float, Text, ForeignKey
from sqlalchemy.exc import SQLAlchemyError
from sqlalchemy.orm import relationship
from src.scrapers.herb_scraper import scrape_herb_details  # Use an absolute import

# Initialize Flask and SQLAlchemy
app = Flask(__name__)
app.config['SQLALCHEMY_DATABASE_URI'] = 'postgresql://postgres:Iamabundant28@localhost/chempath_database'
app.config['SQLALCHEMY_TRACK_MODIFICATIONS'] = False
db = SQLAlchemy(app)

# Your code continues here...

# Define database models
class TraditionalCompound(db.Model):
    __tablename__ = 'traditional_compounds'

    id = Column(Integer, primary_key=True)
    name = Column(String, nullable=False)
    molecular_weight = Column(Float)
    alogp = Column(Float)
    h_don = Column(Integer)
    h_acc = Column(Integer)
    ob_percentage = Column(Float)  # Oral Bioavailability
    caco_2 = Column(Float)  # Caco-2 Permeability
    bbb = Column(Float)  # Blood-brain barrier penetration
    dl = Column(Float)  # Drug-likeness
    fasa = Column(Float)  # Fractional negative surface area
    tpsa = Column(Float)  # Topological polar surface area
    rbn = Column(Integer)  # Rotatable bond number
    hl = Column(Float)  # Half-life
    inchi_key = Column(String)
    pubchem_cid = Column(String)
    cas_number = Column(String)
    synonyms = Column(Text)
    related_articles = Column(Text)
    pubmed_citations = Column(Integer)  # or Column(String) if storing text data

    # Relationships
    herbs = relationship('RelatedHerb', backref='compound', lazy=True)
    targets = relationship('RelatedTarget', backref='compound', lazy=True)
    diseases = relationship('RelatedDisease', backref='compound', lazy=True)

class RelatedHerb(db.Model):
    __tablename__ = 'related_herbs'

    id = Column(Integer, primary_key=True)
    compound_id = Column(Integer, ForeignKey('traditional_compounds.id'), nullable=False)
    chinese_name = Column(String)
    latin_name = Column(String)

class RelatedTarget(db.Model):
    __tablename__ = 'related_targets'

    id = Column(Integer, primary_key=True)
    compound_id = Column(Integer, ForeignKey('traditional_compounds.id'), nullable=False)
    target_name = Column(String)
    source = Column(String)

class RelatedDisease(db.Model):
    __tablename__ = 'related_diseases'

    id = Column(Integer, primary_key=True)
    compound_id = Column(Integer, ForeignKey('traditional_compounds.id'), nullable=False)
    disease_name = Column(String)

# Save data to the database
def save_herb_data(herb_name, compounds, related_targets, related_diseases):
    with app.app_context():
        try:
            # Save compounds
            for compound_data in compounds:
                compound = TraditionalCompound(name=compound_data['name'], molecular_weight=compound_data['molecular_weight'])
                db.session.add(compound)

            # Save related targets
            for target_data in related_targets:
                target = RelatedTarget(compound_name=target_data['compound_name'], target_name=target_data['target_name'], source=target_data['source'])
                db.session.add(target)

            # Save related diseases
            for disease_data in related_diseases:
                disease = RelatedDisease(target_name=disease_data['target_name'], disease_name=disease_data['disease_name'])
                db.session.add(disease)

            db.session.commit()
            print(f"Data for {herb_name} saved successfully.")

        except SQLAlchemyError as e:
            print(f"Error saving data for {herb_name}: {str(e)}")
            db.session.rollback()

# Fetch PubMed data
def fetch_pubmed_data(compound_name):
    print(f"Fetching PubMed data for compound: {compound_name}")
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
        logging.error(f"Error fetching PubMed data for {compound_name}: {str(e)}")
        return 0, []

# Update compounds with PubMed data
def update_with_pubmed_data():
    print("Starting PubMed data update...")
    with app.app_context():
        compounds = TraditionalCompound.query.filter(
            TraditionalCompound.pubmed_citations.is_(None), 
            TraditionalCompound.related_articles.is_(None)
        ).all()
        for compound in compounds:
            try:
                pubmed_count, related_articles = fetch_pubmed_data(compound.name)
                compound.pubmed_citations = pubmed_count
                compound.related_articles = "\n".join(related_articles[:10])
                db.session.commit()
                logging.info(f"Updated PubMed data for {compound.name}")
            except Exception as e:
                logging.error(f"Error updating PubMed data for {compound.name}: {str(e)}")
                db.session.rollback()

# Main workflow
def main():
    # Sample usage: scraping data for a specific herb
    herb_name = 'Aidicha'
    compounds, related_targets, related_diseases = scrape_herb_details(herb_name)
    
    # Check if any data was retrieved before attempting to add to the database
    if compounds or related_targets or related_diseases:
        save_herb_data(herb_name, compounds, related_targets, related_diseases)
    else:
        print(f"No data found for herb: {herb_name}")

    # Update with PubMed data
    update_with_pubmed_data()

if __name__ == "__main__":
    main()
