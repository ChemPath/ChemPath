import sqlite3
import csv
from pathlib import Path
import pubchempy as pcp
import requests
from bs4 import BeautifulSoup
import asyncio
import aiohttp
import logging
from functools import lru_cache
from sqlalchemy import Column, Integer, String, Float, Text
from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy.orm import sessionmaker
from sqlalchemy import types
from marshmallow import Schema, fields, ValidationError
from flask import Flask, request, jsonify
from sqlalchemy import create_engine
from sqlalchemy.orm import sessionmaker
from models import Base
from rdkit.Chem import Descriptors
from rdkit import Chem
from sqlalchemy import update
from sqlalchemy.ext.asyncio import AsyncSession, create_async_engine
from sqlalchemy.future import select
from sqlalchemy.orm import sessionmaker
from sqlalchemy.orm import declarative_base
from marshmallow import Schema, fields
from flask_marshmallow import Marshmallow
import re

tasks = []
async def fetch_pubchem_data(http_session, cid):
    url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid}/JSON"
    async with http_session.get(url) as response:
        if response.status == 200:
            return await response.json()
        return None

async def enrich_compounds(compounds):
    async with aiohttp.ClientSession() as http_session:
        tasks = []
        for compound in compounds:
            tasks.append(asyncio.create_task(fetch_pubchem_data(http_session, compound.cid)))
        results = await asyncio.gather(*tasks)
        return results
# Create the engine
engine = create_engine('sqlite:///chempath.db')  # Adjust the database URL as needed

# Make sure to export the engine
__all__ = ['engine', 'fetch_pubchem_data', 'enrich_compounds']   

class DatabaseError(Exception):
    pass
def create_tables():
    Base.metadata.create_all(engine)
    print("Database tables created successfully.")

Session = sessionmaker(bind=engine)
Base = declarative_base()
app = Flask(__name__) 

import logging
from flask import Flask
from sqlalchemy import Column, Integer, String, Float, Text
from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy.orm import sessionmaker
from sqlalchemy import create_engine

# Create the engine
engine = create_engine('sqlite:///chempath.db')  # Adjust the database URL as needed

# Create the base
Base = declarative_base()

# Create the session
Session = sessionmaker(bind=engine)

# Create the app
app = Flask(__name__)

# Create the logger
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

class DatabaseError(Exception):
    """Base class for database-related exceptions."""
    pass

def create_tables():
    """Create the database tables."""
    Base.metadata.create_all(engine)
    logger.info("Database tables created successfully.")

class Compound(Base):
    """Represents a compound in the database."""
    __tablename__ = 'compounds'

    cid = Column(Integer, primary_key=True)
    cmpdname = Column(String)
    cmpdsynonym = Column(Text)
    molecular_weight = Column(Float)
    molecular_formula = Column(String)
    polararea = Column(Float)
    complexity = Column(Float)
    xlogp = Column(Float)
    heavycnt = Column(Integer)
    hbonddonor = Column(Integer)
    hbondacc = Column(Integer)
    rotbonds = Column(Integer)
    inchi = Column(Text)
    isosmiles = Column(Text)
    canonicalsmiles = Column(Text)
    inchikey = Column(String)
    iupacname = Column(Text)
    exactmass = Column(Float)
    monoisotopicmass = Column(Float)
    charge = Column(Integer)
    covalentunitcnt = Column(Integer)
    isotopeatomcnt = Column(Integer)
    totalatomstereocnt = Column(Integer)
    pclidcnt = Column(Integer)
    gpidcnt = Column(Integer)
    gpfamilycnt = Column(Integer)
    neighbortype = Column(String)
    annothits = Column(Integer)
    annotation = Column(Text)
    plant_source = Column(Text)
    biological_activity = Column(Text)
    traditional_use = Column(Text)

    def __repr__(self):
        return f"<Compound(cid={self.cid}, cmpdname='{self.cmpdname}')>"

__all__ = ['engine', 'Session', 'Base', 'Compound', 'create_tables']

import aiohttp
import asyncio
from functools import lru_cache
from sqlalchemy import update
from sqlalchemy.ext.asyncio import AsyncSession

@lru_cache(maxsize=100)
async def fetch_pubchem_data(session, cid):
    """
    Fetch PubChem data for a given CID.

    Args:
        session (aiohttp.ClientSession): The HTTP session to use.
        cid (int): The PubChem CID.

    Returns:
        dict or None: The PubChem data or None if the request failed.
    """
    url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid}/JSON"
    try:
        async with session.get(url) as response:
            if response.status == 200:
                return await response.json()
            else:
                logger.error(f"Failed to fetch PubChem data for CID {cid}: {response.status}")
                return None
    except Exception as e:
        logger.error(f"Failed to fetch PubChem data for CID {cid}: {str(e)}")
        return None

@lru_cache(maxsize=100)
async def fetch_wikipedia_data(session, compound_name):
    """
    Fetch Wikipedia data for a given compound name.

    Args:
        session (aiohttp.ClientSession): The HTTP session to use.
        compound_name (str): The compound name.

    Returns:
        dict or None: The Wikipedia data or None if the request failed.
    """
    url = f"https://en.wikipedia.org/w/api.php?action=query&prop=extracts&exintro&titles={compound_name}&format=json"
    try:
        async with session.get(url) as response:
            if response.status == 200:
                return await response.json()
            else:
                logger.error(f"Failed to fetch Wikipedia data for compound {compound_name}: {response.status}")
                return None
    except Exception as e:
        logger.error(f"Failed to fetch Wikipedia data for compound {compound_name}: {str(e)}")
        return None

async def update_compound(session: AsyncSession, compound_id: int, update_data: dict):
    """
    Update a compound in the database.

    Args:
        session (AsyncSession): The database session to use.
        compound_id (int): The ID of the compound to update.
        update_data (dict): The data to update.
    """
    stmt = update(Compound).where(Compound.cid == compound_id).values(**update_data)
    await session.execute(stmt)

def normalize_data(data):
    """
    Normalize the data by replacing None values with "NULL".

    Args:
        data (dict): The data to normalize.

    Returns:
        dict: The normalized data.
    """
    return {k: (v if v is not None else "NULL") for k, v in data.items()}

async def enrich_compound_data(session: AsyncSession, compounds, batch_size=10):
    """
    Enrich the compound data by fetching PubChem and Wikipedia data.

    Args:
        session (AsyncSession): The database session to use.
        compounds (list): The list of compounds to enrich.
        batch_size (int): The batch size to use.
    """
    async with aiohttp.ClientSession() as http_session:
        for i in range(0, len(compounds), batch_size):
            batch = compounds[i:i+batch_size]
            tasks = []
            for compound in batch:
                tasks.append(fetch_pubchem_data(http_session, compound.pubchem_cid))
                tasks.append(fetch_wikipedia_data(http_session, compound.cmpdname))
            
            results = await asyncio.gather(*tasks)
            
            for j in range(0, len(results), 2):
                compound = batch[j//2]
                pubchem_data, wiki_data = results[j], results[j+1]
                
                if pubchem_data and wiki_data:
                    try:
                        compound_dict = pubchem_data['PC_Compounds'][0]
                        wiki_summary = next(iter(wiki_data['query']['pages'].values())).get('extract', '')
                        
                        update_data = {
                            'cmpdsynonym': ';'.join(compound_dict.get('synonyms', [])),
                            'biological_activity': compound_dict.get('pharmacology', {}).get('description'),
                            'traditional_use': wiki_summary,
                            'molecular_weight': compound_dict.get('molecular_weight'),
                            'iupacname': compound_dict.get('iupac_name'),
                            'inchi': compound_dict.get('inchi'),
                            'xlogp': compound_dict.get('xlogp')
                        }
                        
                        await update_compound(session, compound.cid, update_data)
                        logger.info(f"Updated data for compound {compound.cmpdname} (ID: {compound.cid})")
                    except Exception as e:
                        logger.error(f"Error updating compound {compound.cmpdname} (ID: {compound.cid}): {str(e)}")
                else:
                    logger.warning(f"Could not fetch data for compound {compound.cmpdname} (ID: {compound.cid})")
            
            await session.commit()

import sqlite3

def display_compound_data(conn, compound_id):
    """
    Display the compound data for a given compound ID.

    Args:
        conn (sqlite3.Connection): The database connection to use.
        compound_id (int): The ID of the compound to display.

    Returns:
        None
    """
    cursor = conn.cursor()
    cursor.execute("SELECT * FROM plant_compounds WHERE id = ?", (compound_id,))
    columns = [description[0] for description in cursor.description]
    row = cursor.fetchone()
    
    if row:
        print("\n" + "="*50)
        print(f"Compound ID: {compound_id}")
        print("-"*50)
        for col, value in zip(columns, row):
            print(f"{col}: {value}")
        print("="*50 + "\n")
    else:
        print(f"No data found for compound ID: {compound_id}")

def create_connection(db_file):
    """
    Create a database connection to a SQLite database.

    Args:
        db_file (str): The path to the database file.

    Returns:
        sqlite3.Connection or None: The database connection or None if an error occurred.
    """
    conn = None
    try:
        conn = sqlite3.connect(db_file)
        print(f"Connected to SQLite version: {sqlite3.sqlite_version}")
        return conn
    except sqlite3.Error as e:
        print(f"Error creating connection: {e}")
    return conn

def create_table(conn):
    try:
        cursor = conn.cursor()
        cursor.execute('''
            CREATE TABLE IF NOT EXISTS plant_compounds (
                id INTEGER PRIMARY KEY,
                cmpdname TEXT NOT NULL,
                smiles TEXT,
                molecular_weight REAL,
                logp REAL,
                plant_source TEXT,
                biological_activities TEXT,
                traditional_use TEXT,
                pubchem_cid INTEGER,
                inchi TEXT,
                iupac_name TEXT,
                synonyms TEXT,
                polararea REAL,
                complexity REAL,
                heavycnt INTEGER,
                hbonddonor INTEGER,
                hbondacc INTEGER,
                rotbonds INTEGER,
                exactmass REAL,
                plant_part TEXT,
                extraction_method TEXT,
                known_interactions TEXT,
                ecological_role TEXT,
                seasonal_variation TEXT,
                therapeutic_areas TEXT
            )
        ''')
        print("Table 'plant_compounds' created successfully")
    except sqlite3.Error as e:
        print(f"Error creating table: {e}")

def drop_table(conn):
    try:
        cursor = conn.cursor()
        cursor.execute('DROP TABLE IF EXISTS plant_compounds')
        print("Existing 'plant_compounds' table dropped")
    except sqlite3.Error as e:
        print(f"Error dropping table: {e}")

def update_database(conn):
    cursor = conn.cursor()
    
    correct_data = {
        "Quercetin": ("O=C1C(O)=C(O)C(=O)C2=C1C=C(O)C(O)=C2O", 302.24, 1.54, 
                      "Various fruits and vegetables",
                      "Antioxidant effects, Anti-inflammatory properties",
                      "Free radical scavenging, Enzyme inhibition", 
                      "Traditional Chinese Medicine for cardiovascular health"),
        "Curcumin": ("COC1=CC(=CC(=C1O)OC)C=CC(=O)CC(=O)C=CC2=CC(=C(C=C2)O)OC", 368.38, 3.29,
                     "Turmeric (Curcuma longa)",
                     "Anti-inflammatory effects, Cancer prevention",
                     "NF-κB inhibition, Antioxidant activity",
                     "Ayurvedic medicine for various ailments"),
        "Resveratrol": ("OC1=CC(=CC(=C1)O)C=CC2=CC(=C(C=C2)O)O", 228.24, 3.1,
                        "Grapes, red wine, peanuts",
                        "Cardiovascular health promotion, Anti-aging effects",
                        "Sirtuin activation, Antioxidant activity",
                        "Traditional use in some Asian medicines"),
        "Caffeine": ("CN1C=NC2=C1C(=O)N(C(=O)N2C)C", 194.19, -0.07,
                     "Coffee beans, tea leaves, cacao beans",
                     "Central nervous system stimulation, Cognitive enhancement",
                     "Adenosine receptor antagonism, Phosphodiesterase inhibition",
                     "Traditional use as a stimulant in various cultures"),
        "Theobromine": ("CN1C=NC2=C1C(=O)NC(=O)N2", 180.16, -0.78,
                        "Cacao beans, tea leaves",
                        "Mild stimulant effects, Cardiovascular effects",
                        "Phosphodiesterase inhibition, Adenosine receptor antagonism",
                        "Traditional use in chocolate and cocoa-based beverages"),
        "Capsaicin": ("COC1=C(C=CC(=C1)CNC(=O)CCCC\\C=C/C(C)C)O", 305.41, 3.04,
                      "Chili peppers",
                      "Pain relief, Anti-inflammatory effects",
                      "TRPV1 receptor activation, Substance P depletion",
                      "Traditional use in various cuisines and folk medicine"),
        "Lycopene": ("CC(C)=CCC\\C(C)=C\\C\\C=C(\\C)CCC=C(C)C=CC=C(C)C=CC=C(C)C=C", 536.87, 9.16,
                     "Tomatoes, watermelons, pink grapefruits",
                     "Antioxidant effects, Potential cancer prevention",
                     "Free radical scavenging, Modulation of cell signaling pathways",
                     "Traditional consumption as part of Mediterranean diet"),
        "Gingerol": ("CCCCC[C@H](O)CC(=O)CCc1ccc(O)c(OC)c1", 294.39, 3.85,
                     "Ginger root",
                     "Anti-inflammatory effects, Digestive aid",
                     "COX-2 inhibition, Antioxidant activity",
                     "Traditional use in Asian medicine for various ailments"),
        "Epigallocatechin gallate": ("O=C(O[C@H]1Cc2c(O)cc(O)cc2O[C@@H]1c1cc(O)c(O)c(O)c1)c1cc(O)c(O)c(O)c1", 458.37, 1.5,
                                     "Green tea",
                                     "Antioxidant, Anti-cancer, Weight loss",
                                     "Free radical scavenging, Enzyme inhibition",
                                     "Traditional use in Chinese and Japanese medicine"),
        "Allicin": ("CC(=O)S[S]C(C)=O", 162.27, 1.35,
                    "Garlic",
                    "Antimicrobial, Cardiovascular health",
                    "Thiol modification, Enzyme inhibition",
                    "Traditional use in many cultures for health and flavor"),
        "Berberine": ("COc1ccc2cc3[n+](cc2c1OC)CCc1cc2c(cc1-3)OCO2", 336.36, -1.3,
                      "Berberis plants, goldenseal",
                      "Anti-diabetic, Anti-microbial",
                      "AMPK activation, Gut microbiome modulation",
                      "Traditional use in Chinese and Ayurvedic medicine"),
        "Piperine": ("O=C(/C=C/C=C/c1ccc2c(c1)OCO2)N1CCC[C@H]1C", 285.34, 2.78,
                     "Black pepper",
                     "Bioavailability enhancer, Anti-inflammatory",
                     "Inhibition of drug metabolism enzymes",
                     "Traditional use in Ayurvedic medicine"),
        "Eugenol": ("COc1cc(CC=C)ccc1O", 164.20, 2.27,
                    "Cloves, cinnamon, basil",
                    "Analgesic, Antimicrobial",
                    "Calcium channel modulation, Free radical scavenging",
                    "Traditional use in dental care and as a flavoring agent"),
        "Catechin": ("O=C1c2c(O)cc(O)cc2OC(c2ccc(O)c(O)c2)C1O", 290.27, 0.51,
                     "Tea, cocoa, berries",
                     "Antioxidant, Cardioprotective",
                     "Free radical scavenging, Enzyme inhibition",
                     "Traditional use in various herbal medicines"),
        "Kaempferol": ("O=c1c(O)c(-c2ccc(O)cc2)oc2cc(O)cc(O)c12", 286.24, 1.9,
                       "Broccoli, kale, tea",
                       "Anti-inflammatory, Antioxidant",
                       "NF-κB inhibition, Free radical scavenging",
                       "Traditional use in various herbal remedies"),
        "Limonene": ("CC(=CCC[C@H](C)C=C)C", 136.23, 4.57,
                     "Citrus fruits",
                     "Anti-cancer, Anxiolytic",
                     "Apoptosis induction, GABA receptor modulation",
                     "Traditional use in aromatherapy and as a flavoring agent"),
        "Menthol": ("CC(C)[C@H]1CC[C@@H](C)C[C@H]1O", 156.27, 3.4,
                    "Peppermint, other mint species",
                    "Analgesic, Decongestant",
                    "TRPM8 activation, κ-opioid receptor agonism",
                    "Traditional use in folk medicine for various ailments"),
        "Rosmarinic acid": ("O=C(O[C@H](Cc1ccc(O)c(O)c1)[C@H](O)C(=O)O)c1ccc(O)c(O)c1", 360.31, 1.82,
                            "Rosemary, sage, basil",
                            "Antioxidant, Anti-inflammatory",
                            "Free radical scavenging, Enzyme inhibition",
                            "Traditional use in Mediterranean herbal medicine"),
        "Sulforaphane": ("CS(=O)CCCN=C=S", 177.29, 0.72,
                         "Broccoli, Brussels sprouts, cabbage",
                         "Anti-cancer, Neuroprotective",
                         "Nrf2 activation, HDAC inhibition",
                         "Modern use in cancer prevention diets"),
        "Vincristine": ("CCC1(O)C(=O)OCC1C1(C(=O)OC)C(OC(C)=O)CC2C1CN1CCC(CC3=C(OC)C=CC(=C3)C1)C2C(=O)OC", 824.96, 2.82,
                        "Catharanthus roseus (Madagascar periwinkle)",
                        "Anti-cancer",
                        "Microtubule destabilization",
                        "Modern use in chemotherapy")
    }
    
    try:
        conn.execute("BEGIN TRANSACTION")
        for cmpdname, data in correct_data.items():
            smiles, molecular_weight, logp, plant_source, therapeutic_areas, biological_activities, traditional_use = data
            
            insert_sql = """
            INSERT INTO plant_compounds
            (cmpdname, smiles, molecular_weight, logp, plant_source, therapeutic_areas, biological_activities, traditional_use)
            VALUES (?, ?, ?, ?, ?, ?, ?, ?)
            """
            cursor.execute(insert_sql, (cmpdname, smiles, molecular_weight, logp, plant_source, therapeutic_areas, biological_activities, traditional_use))
            print(f"Inserted {cmpdname}")
        
        conn.commit()
        print("Database update completed successfully")
    except sqlite3.Error as e:
        conn.rollback()
        print(f"An error occurred: {e}")
    finally:
        cursor.close()
        
def main():
    db_file = "chempath.db"
    conn = create_connection(db_file)
    if conn:
        try:
            print(f"Connected to SQLite version: {sqlite3.sqlite_version}")
        except AttributeError:
            print("Could not determine SQLite version")
        drop_table(conn)
        create_table(conn)
        update_database(conn)
        conn.close()

if __name__ == "__main__":
    main()

class DatabaseError(Exception):
    """Base class for database-related exceptions."""
    pass

def insert_compound(conn, compound):
    """
    Insert a compound into the database.

    Args:
        conn (sqlite3.Connection): The database connection to use.
        compound (dict): The compound data to insert.

    Returns:
        int: The ID of the inserted compound.

    Raises:
        ValueError: If the compound data is invalid.
        DatabaseError: If an error occurs during insertion.
    """
    errors = validate_compound_data(compound)
    if errors:
        raise ValueError(f"Invalid compound data: {', '.join(errors)}")

    sql = '''INSERT INTO plant_compounds(cmpdname, smiles, molecular_weight, logp, plant_source, biological_activity, traditional_use, pubchem_cid, inchi, iupac_name, synonyms)
             VALUES(?,?,?,?,?,?,?,?,?,?,?)'''
    try:
        cursor = conn.cursor()
        cursor.execute(sql, tuple(compound.values()))
        conn.commit()
        return cursor.lastrowid
    except sqlite3.Error as e:
        conn.rollback()
        raise DatabaseError(f"Error inserting compound: {str(e)}")

def validate_compound_data(compound):
    errors = []
    
    # Validate cmpdname
    if not compound['cmpdname']:
        errors.append("Compound name is required")
    elif len(compound['cmpdname']) < 2:
        errors.append("Compound name must be at least 2 characters long")
    
    # Validate SMILES
    if not compound['smiles']:
        errors.append("SMILES is required")
    else:
        mol = Chem.MolFromSmiles(compound['smiles'])
        if not mol:
            errors.append("Invalid SMILES string")
        else:
            # Calculate and validate molecular properties
            compound['molecular_weight'] = Descriptors.ExactMolWt(mol)
            compound['logp'] = Descriptors.MolLogP(mol)
            
            if compound['molecular_weight'] < 0 or compound['molecular_weight'] > 2000:
                errors.append("Molecular weight must be between 0 and 2000 g/mol")
            
            if compound['logp'] < -10 or compound['logp'] > 10:
                errors.append("LogP must be between -10 and 10")
    
    # Validate plant_source
    if not compound['plant_source']:
        errors.append("Plant source is required")
    elif len(compound['plant_source']) < 3:
        errors.append("Plant source must be at least 3 characters long")
    elif not re.match(r'^[A-Za-z\s]+$', compound['plant_source']):
        errors.append("Plant source should only contain letters and spaces")
    
    # Validate biological_activity
    if not compound['biological_activity']:
        errors.append("Biological activity is required")
    elif len(compound['biological_activity']) < 5:
        errors.append("Biological activity description must be at least 5 characters long")
    
    # Validate traditional_use
    if compound['traditional_use'] and len(compound['traditional_use']) < 5:
        errors.append("Traditional use description, if provided, must be at least 5 characters long")
    
    # Validate pubchem_cid
    if 'pubchem_cid' in compound and compound['pubchem_cid']:
        try:
            pubchem_cid = int(compound['pubchem_cid'])
            if pubchem_cid < 1:
                errors.append("PubChem CID must be a positive integer")
        except ValueError:
            errors.append("PubChem CID must be a valid integer")
    
    # Check for data consistency
    if 'inchi' in compound and compound['inchi']:
        inchi_mol = Chem.MolFromInchi(compound['inchi'])
        if inchi_mol:
            inchi_smiles = Chem.MolToSmiles(inchi_mol)
            if inchi_smiles != compound['smiles']:
                errors.append("SMILES and InChI are inconsistent")
        else:
            errors.append("Invalid InChI string")
    # Ensure plant_source contains meaningful data
    if compound['plant_source']:
        if not re.search(r'\b[A-Z][a-z]+\b', compound['plant_source']):
            errors.append("Plant source should contain at least one properly capitalized plant name")
    
    # Ensure biological_activity contains meaningful data
    if compound['biological_activity']:
        activity_keywords = ['inhibitor', 'agonist', 'antagonist', 'modulator', 'inducer']
        if not any(keyword in compound['biological_activity'].lower() for keyword in activity_keywords):
            errors.append("Biological activity should contain at least one recognized activity term")
    
    # Check consistency between molecular_weight and SMILES
    if 'molecular_weight' in compound and compound['smiles']:
        mol = Chem.MolFromSmiles(compound['smiles'])
        if mol:
            calculated_mw = Descriptors.ExactMolWt(mol)
            if abs(compound['molecular_weight'] - calculated_mw) > 0.1:
                errors.append("Provided molecular weight is inconsistent with SMILES")
    
    # Check consistency between logp and SMILES
    if 'logp' in compound and compound['smiles']:
        mol = Chem.MolFromSmiles(compound['smiles'])
        if mol:
            calculated_logp = Descriptors.MolLogP(mol)
            if abs(compound['logp'] - calculated_logp) > 0.5:
                errors.append("Provided LogP is inconsistent with SMILES")
    
    # Check consistency between iupac_name and SMILES
    if 'iupac_name' in compound and compound['iupac_name'] and compound['smiles']:
        iupac_mol = Chem.MolFromIupacName(compound['iupac_name'])
        if iupac_mol:
            iupac_smiles = Chem.MolToSmiles(iupac_mol)
            if iupac_smiles != compound['smiles']:
                errors.append("IUPAC name is inconsistent with SMILES")
    
    return errors

def import_compound_data(conn, csv_file):
    """
    Import compound data from a CSV file.

    Args:
        conn (sqlite3.Connection): The database connection to use.
        csv_file (str): The path to the CSV file.

    Raises:
        DatabaseError: If an error occurs during import.
    """
    cursor = conn.cursor()
    with open(csv_file, 'r') as file:
        csv_reader = csv.DictReader(file)
        for row in csv_reader:
            errors = validate_compound_data(row)
            if errors:
                print(f"Skipping row due to errors: {', '.join(errors)}")
                continue
            try:
                cursor.execute('''
                    INSERT INTO plant_compounds 
                    (cmpdname, smiles, molecular_weight, logp, plant_source, biological_activity, traditional_use, pubchem_cid, inchi, iupac_name, synonyms)
                    VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
                ''', (
                    row['cmpdname'], row['smiles'], float(row['molecular_weight']), float(row['logp']),
                    row['plant_source'], row['biological_activity'], row['traditional_use'],
                    int(row['pubchem_cid']), row['inchi'], row['iupac_name'], row['synonyms']
                ))
            except sqlite3.Error as e:
                conn.rollback()
                raise DatabaseError(f"Error importing compound data: {str(e)}")
    conn.commit()
    print("Data imported successfully")

def create_indexes(conn):
    cursor = conn.cursor()
    try:
        # Existing indexes
        cursor.execute('CREATE INDEX IF NOT EXISTS idx_cmpdname ON plant_compounds(cmpdname)')
        cursor.execute('CREATE INDEX IF NOT EXISTS idx_smiles ON plant_compounds(smiles)')
        cursor.execute('CREATE INDEX IF NOT EXISTS idx_plant_source ON plant_compounds(plant_source)')
        
        # New composite indexes
        cursor.execute('CREATE INDEX IF NOT EXISTS idx_name_mw ON plant_compounds(cmpdname, molecular_weight)')
        cursor.execute('CREATE INDEX IF NOT EXISTS idx_source_activity ON plant_compounds(plant_source, biological_activity)')
        
        conn.commit()
        print("Indexes created successfully")
    except sqlite3.Error as e:
        conn.rollback()
        raise DatabaseError(f"Error creating indexes: {str(e)}")
def enable_fts(conn):
    cursor = conn.cursor()
    try:
        cursor.execute('''
        CREATE VIRTUAL TABLE IF NOT EXISTS plant_compounds_fts USING fts5(
            cmpdname, plant_source, biological_activity, traditional_use,
            content='plant_compounds', content_rowid='id'
        )
        ''')
        cursor.execute('''
        INSERT INTO plant_compounds_fts(plant_compounds_fts) VALUES('rebuild')
        ''')
        conn.commit()
        print("Full-text search enabled successfully")
    except sqlite3.Error as e:
        conn.rollback()
        raise DatabaseError(f"Error enabling full-text search: {str(e)}")

def fts_search(conn, query):
    cursor = conn.cursor()
    cursor.execute('''
    SELECT c.* FROM plant_compounds c
    JOIN plant_compounds_fts fts ON c.id = fts.rowid
    WHERE plant_compounds_fts MATCH ?
    ORDER BY rank
    ''', (query,))
    return cursor.fetchall()

def partition_database(conn):
    cursor = conn.cursor()
    try:
        # Create partitions based on molecular weight ranges
        cursor.execute('''
        CREATE TABLE IF NOT EXISTS compounds_light (
            CHECK( molecular_weight < 300 )
        ) INHERITS (plant_compounds)
        ''')
        cursor.execute('''
        CREATE TABLE IF NOT EXISTS compounds_medium (
            CHECK( molecular_weight >= 300 AND molecular_weight < 500 )
        ) INHERITS (plant_compounds)
        ''')
        cursor.execute('''
        CREATE TABLE IF NOT EXISTS compounds_heavy (
            CHECK( molecular_weight >= 500 )
        ) INHERITS (plant_compounds)
        ''')
        
        # Create a trigger to automatically insert into the correct partition
        cursor.execute('''
        CREATE OR REPLACE FUNCTION compound_insert_trigger()
        RETURNS TRIGGER AS $$
        BEGIN
            IF ( NEW.molecular_weight < 300 ) THEN
                INSERT INTO compounds_light VALUES (NEW.*);
            ELSIF ( NEW.molecular_weight >= 300 AND NEW.molecular_weight < 500 ) THEN
                INSERT INTO compounds_medium VALUES (NEW.*);
            ELSE
                INSERT INTO compounds_heavy VALUES (NEW.*);
            END IF;
            RETURN NULL;
        END;
        $$ LANGUAGE plpgsql;
        ''')
        
        cursor.execute('''
        CREATE TRIGGER insert_compound_trigger
        BEFORE INSERT ON plant_compounds
        FOR EACH ROW EXECUTE FUNCTION compound_insert_trigger();
        ''')
        
        conn.commit()
        print("Database partitioning implemented successfully")
    except sqlite3.Error as e:
        conn.rollback()
        raise DatabaseError(f"Error implementing database partitioning: {str(e)}")

def analyze_queries(conn):
    cursor = conn.cursor()
    cursor.execute("ANALYZE")
    cursor.execute("PRAGMA optimize")
    conn.commit()
    print("Database analyzed and optimized")

class ChemPathAPI:
    def __init__(self, db_path):
        """
        Initialize the ChemPath API.

        Args:
            db_path (str): The path to the database file.
        """
        self.conn = sqlite3.connect(db_path)

    def get_compound(self, compound_id):
        """
        Get a compound by ID.

        Args:
            compound_id (int): The ID of the compound to get.

        Returns:
            tuple or None: The compound data or None if not found.
        """
        cursor = self.conn.cursor()
        cursor.execute('SELECT * FROM plant_compounds WHERE id = ?', (compound_id,))
        return cursor.fetchone()

    def search_compounds(self, query, search_field='cmpdname'):
        """
        Search compounds by query.

        Args:
            query (str): The query to search for.
            search_field (str): The field to search in (default: 'cmpdname').

        Returns:
            list: The list of compound data matching the query.
        """
        cursor = self.conn.cursor()
        cursor.execute(f'''
            SELECT * FROM plant_compounds 
            WHERE {search_field} LIKE ?
        ''', (f'%{query}%',))
        return cursor.fetchall()

    def add_compound(self, compound_data):
        """
        Add a new compound to the database.

        Args:
            compound_data (dict): The data for the new compound.

        Returns:
            int: The ID of the newly added compound.
        """
        cursor = self.conn.cursor()
        cursor.execute('''
            INSERT INTO plant_compounds 
            (cmpdname, smiles, molecular_weight, logp, plant_source, biological_activity, traditional_use)
            VALUES (?, ?, ?, ?, ?, ?, ?)
        ''', (
            compound_data['cmpdname'], compound_data['smiles'], compound_data['molecular_weight'],
            compound_data['logp'], compound_data['plant_source'], compound_data['biological_activity'],
            compound_data['traditional_use']
        ))
        self.conn.commit()
        return cursor.lastrowid

    def close(self):
        """
        Close the database connection.
        """
        self.conn.close()

def main():
    db_file = "chempath.db"
    conn = create_connection(db_file)
    if conn:
        create_table(conn)
        update_database(conn)
        conn.close()

if __name__ == "__main__":
    main()
app = Flask(__name__)
ma = Marshmallow(app)

# Define the database connection URL
database_url = 'sqlite:///chempath_database.db'

# Create the database engine
engine = create_engine(database_url)

# Create a configured "Session" class
Session = sessionmaker(bind=engine)

# Create a base class for declarative class definitions
Base = declarative_base()

from sqlalchemy import Table, Column, Integer, ForeignKey
from sqlalchemy.orm import relationship

compound_pathway = Table('compound_pathway', Base.metadata,
    Column('compound_id', Integer, ForeignKey('plant_compounds.id')),
    Column('pathway_id', Integer, ForeignKey('biosynthetic_pathways.id'))
)

class Compound(Base):
    __tablename__ = 'compounds'
    id = Column(Integer, primary_key=True)
    cmpdname = Column(String)
    smiles = Column(String)
    molecular_weight = Column(Float)
    logp = Column(Float)
    plant_source = Column(String)
    biological_activity = Column(String)
    traditional_use = Column(String)
    pubchem_cid = Column(Integer)
    inchi = Column(String)
    iupac_name = Column(String)
    synonyms = Column(String)

    pathways = relationship("BiosyntheticPathway", secondary=compound_pathway, back_populates="compounds")

    def __repr__(self):
        return f"<Compound(cmpdname='{self.cmpdname}', smiles='{self.smiles}')>"

class CompoundSchema(Schema):
    id = fields.Int(dump_only=True)
    cmpdname = fields.Str(required=True)
    smiles = fields.Str(required=True)
    molecular_weight = fields.Float()
    logp = fields.Float()
    plant_source = fields.Str()
    biological_activity = fields.Str()
    traditional_use = fields.Str()
    pubchem_cid = fields.Int()
    inchi = fields.Str()
    iupac_name = fields.Str()
    synonyms = fields.Str()

compound_schema = CompoundSchema()
compounds_schema = CompoundSchema(many=True)




class BiosyntheticPathway(Base):
    __tablename__ = 'biosynthetic_pathways'

    id = Column(Integer, primary_key=True)
    name = Column(String, nullable=False)
    description = Column(String)
    organisms = Column(String)
    compounds = relationship("Compound", secondary=compound_pathway, back_populates="pathways")


def create_engine_and_session():
    engine = create_engine(database_url)
    Session = sessionmaker(bind=engine)
    return engine, Session()

def add_compound(session, compound_data):
    new_compound = Compound(**compound_data)
    session.add(new_compound)
    session.commit()
    return new_compound.id

def get_compound(session, compound_id):
    return session.query(Compound).filter(Compound.id == compound_id).first()

def search_compounds(session, query, search_field='cmpdname'):
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

async def enrich_compound_data(session, compounds):
    # TO DO: implement the enrichment logic here
    pass

@app.route('/compound', methods=['POST'])
def add_compound_route():
    data = request.json
    try:
        compound_data = compound_schema.load(data)
    except ValidationError as err:
        return jsonify(err.messages), 400
    
    compound_id = add_compound(Session(), compound_data)
    return jsonify({'message': 'Compound added successfully', 'id': compound_id}), 201

@app.route('/compound/<int:compound_id>', methods=['GET'])
def get_compound_route(compound_id):
    compound = get_compound(Session(), compound_id)
    if compound:
        return jsonify(compound_schema.dump(compound))
    return jsonify({'message': 'Compound not found'}), 404

@app.route('/compounds/search', methods=['GET'])
def search_compounds_route():
    query = request.args.get('query', '')
    search_field = request.args.get('field', 'cmpdname')
    compounds = search_compounds(Session(), query, search_field)
    return jsonify(compounds_schema.dump(compounds))

@app.route('/compound/<int:compound_id>', methods=['PUT'])
def update_compound_route(compound_id):
    data = request.json
    try:
        update_data = compound_schema.load(data, partial=True)
    except ValidationError as err:
        return jsonify(err.messages), 400
    
    success = update_compound(Session(), compound_id, update_data)
    if success:
        return jsonify({'message': 'Compound updated successfully'})
    return jsonify({'message': 'Compound not found'}), 404
@app.route('/compound/<int:compound_id>', methods=['DELETE'])
def delete_compound_route(compound_id):
    success = delete_compound(Session(), compound_id)
    if success:
        return jsonify({'message': 'Compound deleted successfully'})
    return jsonify({'message': 'Compound not found'}), 404

@app.route('/compounds/enrich', methods=['POST'])
def enrich_compounds_route():
    compound_ids = request.json.get('compound_ids', [])
    compounds = Session().query(Compound).filter(Compound.id.in_(compound_ids)).all()
    asyncio.run(enrich_compound_data(Session(), compounds))
    return jsonify({'message': 'Compound enrichment process started'})

if __name__ == '__main__':
    engine, session = create_engine_and_session()
    Base.metadata.create_all(engine)
    app.run(debug=True)