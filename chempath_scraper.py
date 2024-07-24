# File: chempath_scraper.py
from bs4 import BeautifulSoup

import requests

# If the above import fails, you may need to install the requests library:
# pip install requests
# If the above import fails, you may need to install the requests library:
# pip install requestsfrom bs4 import BeautifulSoup
import re
from chempath_database import create_connection, insert_compound

def scrape_pubchem(compound_name):
    """Scrape compound information from PubChem"""
    url = f"https://pubchem.ncbi.nlm.nih.gov/compound/{compound_name}"
    response = requests.get(url)
    soup = BeautifulSoup(response.text, 'html.parser')
    
    # Extract information with error handling
    smiles = molecular_weight = logp = ""
    
    smiles_div = soup.find('div', string='SMILES')
    if smiles_div and smiles_div.find_next('div'):
        smiles = smiles_div.find_next('div').text.strip()
    
    mw_div = soup.find('div', string='Molecular Weight')
    if mw_div and mw_div.find_next('div'):
        mw_text = mw_div.find_next('div').text.strip()
        molecular_weight = float(re.search(r'\d+\.\d+', mw_text).group()) if re.search(r'\d+\.\d+', mw_text) else 0.0
    
    logp_div = soup.find('div', string='LogP')
    if logp_div and logp_div.find_next('div'):
        logp_text = logp_div.find_next('div').text.strip()
        logp = float(re.search(r'-?\d+\.\d+', logp_text).group()) if re.search(r'-?\d+\.\d+', logp_text) else 0.0
    
    plant_source = ""
    biological_activity = ""
    traditional_use = ""
    
    return (compound_name, smiles, molecular_weight, logp, plant_source, biological_activity, traditional_use)

def main():
    # List of compounds to scrape
    compounds_to_scrape = [
        "Quercetin", 
        "Curcumin", 
        "Resveratrol",
        "Caffeine",
        "Theobromine",
        "Capsaicin",
        "Lycopene",
        "Gingerol"
    ]
    
    # Rest of the function remains the same

    
    # Connect to the database
    conn = create_connection("chempath_database.db")
    
    if conn is not None:
        for compound in compounds_to_scrape:
            try:
                print(f"Scraping data for {compound}...")
                compound_data = scrape_pubchem(compound)
                insert_compound(conn, compound_data)
                print(f"Data for {compound} inserted successfully.")
            except Exception as e:
                print(f"Error scraping {compound}: {str(e)}")
        
        conn.close()
    else:
        print("Error! Cannot create the database connection.")

if __name__ == "__main__":
    main()