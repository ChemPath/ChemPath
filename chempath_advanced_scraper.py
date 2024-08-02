# File: chempath_advanced_scraper.py

import requests
from bs4 import BeautifulSoup
import re
from src.database.chempath_database import create_connection
from src.database.chempath_database_cleanup import update_or_insert_compound

def safe_find(soup, class_name, string):
    """Safely find an element, returning None if not found"""
    element = soup.find('div', class_=class_name, string=string)
    return element.find_next('div') if element else None

def safe_extract_float(text):
    """Safely extract a float from text, returning None if not found"""
    if text:
        match = re.search(r'-?\d+\.?\d*', text)
        return float(match.group()) if match else None
    return None

def scrape_pubchem(compound_name):
    """Scrape compound information from PubChem with error handling"""
    url = f"https://pubchem.ncbi.nlm.nih.gov/compound/{compound_name}"
    try:
        response = requests.get(url)
        response.raise_for_status()
        soup = BeautifulSoup(response.text, 'html.parser')
        
        smiles_div = safe_find(soup, 'p-sm-top', 'Canonical SMILES')
        mw_div = safe_find(soup, 'p-sm-top', 'Molecular Weight')
        logp_div = safe_find(soup, 'p-sm-top', 'LogP')
        
        smiles = smiles_div.text.strip() if smiles_div else ""
        molecular_weight = safe_extract_float(mw_div.text.strip() if mw_div else "") or 0
        logp = safe_extract_float(logp_div.text.strip() if logp_div else "") or 0
        
        return smiles, molecular_weight, logp
    except requests.RequestException as e:
        print(f"Error fetching data for {compound_name}: {e}")
        return "", 0, 0

def scrape_compound(compound_name):
    """Combine data from multiple sources"""
    smiles, molecular_weight, logp = scrape_pubchem(compound_name)
    
    # Placeholder for additional scraping functions
    plant_source = ""
    therapeutic_areas = ""
    biological_activities = ""
    traditional_use = ""
    
    return (compound_name, smiles, molecular_weight, logp, plant_source, therapeutic_areas, biological_activities, traditional_use)

def main():
    compounds_to_scrape = ["Quercetin", "Curcumin", "Resveratrol", "Caffeine", "Theobromine", "Capsaicin", "Lycopene", "Gingerol"]
    conn = create_connection("chempath_database.db")
    
    if conn is not None:
        for compound in compounds_to_scrape:
            try:
                print(f"Scraping data for {compound}...")
                compound_data = scrape_compound(compound)
                update_or_insert_compound(conn, compound_data)
                print(f"Data for {compound} updated/inserted successfully.")
            except Exception as e:
                print(f"Error processing {compound}: {str(e)}")
        conn.close()
    else:
        print("Error! Cannot create the database connection.")

if __name__ == "__main__":
    main()