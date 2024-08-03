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
from sqlalchemy import create_engine, Column, Integer, String, Float
from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy.orm import sessionmaker
from sqlalchemy import types
from marshmallow import Schema, fields, ValidationError
from flask import Flask, request, jsonify
from sqlalchemy import create_engine
from sqlalchemy.orm import sessionmaker
from models import Base

# Create the engine
engine = create_engine('sqlite:///chempath.db')  # Adjust the database URL as needed

# Make sure to export the engine
__all__ = ['engine']

def create_tables():
    Base.metadata.create_all(engine)
    print("Database tables created successfully.")

Session = sessionmaker(bind=engine)
Base = declarative_base()
app = Flask(__name__)


logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

@lru_cache(maxsize=100)
async def fetch_pubchem_data(session, cid):
    url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid}/JSON"
    async with session.get(url) as response:
        if response.status == 200:
            return await response.json()
        return None

@lru_cache(maxsize=100)
async def fetch_wikipedia_data(session, name):
    url = f"https://en.wikipedia.org/api/rest_v1/page/summary/{name}"
    async with session.get(url) as response:
        if response.status == 200:
            return await response.json()
        return None

def normalize_data(data):
    return {k: (v if v is not None else "NULL") for k, v in data.items()}

async def enrich_compound_data(conn, compounds, batch_size=10):
    cursor = conn.cursor()
    
    async with aiohttp.ClientSession() as session:
        for i in range(0, len(compounds), batch_size):
            batch = compounds[i:i+batch_size]
            tasks = []
            for cid, name in batch:
                tasks.append(asyncio.create_task(fetch_pubchem_data(session, cid)))
                tasks.append(asyncio.create_task(fetch_wikipedia_data(session, name)))
            
            results = await asyncio.gather(*tasks)
            
            for j in range(0, len(results), 2):
                cid, name = batch[j//2]
                pubchem_data, wiki_data = results[j], results[j+1]
                
                if pubchem_data and wiki_data:
                    try:
                        compound_dict = pubchem_data['PC_Compounds'][0]
                        wiki_summary = wiki_data.get('extract', '')
                        
                        new_data = normalize_data({
                            'cmpdsynonym': ';'.join(compound_dict.get('synonyms', [])),
                            'biological_activity': compound_dict.get('pharmacology', {}).get('description'),
                            'traditional_use': wiki_summary,
                            'molecular_weight': compound_dict.get('molecular_weight'),
                            'molecular_formula': compound_dict.get('molecular_formula'),
                            'iupacname': compound_dict.get('iupac_name'),
                            'polararea': compound_dict.get('topological_polar_surface_area'),
                            'complexity': compound_dict.get('complexity'),
                            'xlogp': compound_dict.get('xlogp'),
                            'heavycnt': compound_dict.get('heavy_atom_count'),
                            'hbonddonor': compound_dict.get('h_bond_donor_count'),
                            'hbondacc': compound_dict.get('h_bond_acceptor_count'),
                            'rotbonds': compound_dict.get('rotatable_bond_count'),
                            'inchi': compound_dict.get('inchi'),
                            'isosmiles': compound_dict.get('isomeric_smiles'),
                            'canonicalsmiles': compound_dict.get('canonical_smiles'),
                            'inchikey': compound_dict.get('inchikey'),
                            'exactmass': compound_dict.get('exact_mass'),
                            'monoisotopicmass': compound_dict.get('monoisotopic_mass'),
                            'charge': compound_dict.get('charge'),
                            'covalentunitcnt': compound_dict.get('covalent_unit_count'),
                            'isotopeatomcnt': compound_dict.get('isotope_atom_count'),
                            'totalatomstereocnt': compound_dict.get('defined_atom_stereo_count'),
                        })
                        
                        cursor.execute('''
                        UPDATE compounds SET
                        cmpdsynonym = :cmpdsynonym,
                        biological_activity = :biological_activity,
                        traditional_use = :traditional_use,
                        molecular_weight = :molecular_weight,
                        molecular_formula = :molecular_formula,
                        iupacname = :iupacname,
                        polararea = :polararea,
                        complexity = :complexity,
                        xlogp = :xlogp,
                        heavycnt = :heavycnt,
                        hbonddonor = :hbonddonor,
                        hbondacc = :hbondacc,
                        rotbonds = :rotbonds,
                        inchi = :inchi,
                        isosmiles = :isosmiles,
                        canonicalsmiles = :canonicalsmiles,
                        inchikey = :inchikey,
                        exactmass = :exactmass,
                        monoisotopicmass = :monoisotopicmass,
                        charge = :charge,
                        covalentunitcnt = :covalentunitcnt,
                        isotopeatomcnt = :isotopeatomcnt,
                        totalatomstereocnt = :totalatomstereocnt
                        WHERE cid = :cid
                        ''', {**new_data, 'cid': cid})
                        
                        logger.info(f"Updated data for compound {name} (CID: {cid})")
                    except Exception as e:
                        logger.error(f"Error updating compound {name} (CID: {cid}): {str(e)}")
                else:
                    logger.warning(f"Could not fetch data for compound {name} (CID: {cid})")
            
            conn.commit()

def display_compound_data(conn, cid):
    cursor = conn.cursor()
    cursor.execute("SELECT * FROM compounds WHERE cid = ?", (cid,))
    columns = [description[0] for description in cursor.description]
    row = cursor.fetchone()
    
    if row:
        print("\n" + "="*50)
        for col, value in zip(columns, row):
            print(f"{col}: {value}")
    else:
        print(f"No data found for CID: {cid}")

def create_connection(db_file):
    """Create a database connection to a SQLite database"""
    conn = None
    try:
        conn = sqlite3.connect(db_file)
        print(f"Connected to SQLite version: {sqlite3.version}")
        return conn
    except sqlite3.Error as e:
        print(e)
    return conn

def create_table(conn):
    """Create the plant_compounds table"""
    try:
        cursor = conn.cursor()
        cursor.execute('''
            CREATE TABLE IF NOT EXISTS plant_compounds (
                id INTEGER PRIMARY KEY,
                name TEXT NOT NULL,
                smiles TEXT,
                molecular_weight REAL,
                logp REAL,
                plant_source TEXT,
                biological_activity TEXT,
                traditional_use TEXT,
                pubchem_cid INTEGER,
                inchi TEXT,
                iupac_name TEXT,
                synonyms TEXT
            )
        ''')
        print("Table 'plant_compounds' created successfully")
    except sqlite3.Error as e:
        print(e)

def insert_compound(conn, compound):
    """Insert a new compound into the plant_compounds table"""
    sql = '''INSERT INTO plant_compounds(name, smiles, molecular_weight, logp, plant_source, biological_activity, traditional_use, pubchem_cid, inchi, iupac_name, synonyms)
             VALUES(?,?,?,?,?,?,?,?,?,?,?)'''
    try:
        cursor = conn.cursor()
        cursor.execute(sql, compound)
        conn.commit()
        return cursor.lastrowid
    except sqlite3.Error as e:
        print(e)
        return None

def validate_compound_data(row):
    errors = []
    if not row['name']:
        errors.append("Name is required")
    if not row['smiles']:
        errors.append("SMILES is required")
    try:
        float(row['molecular_weight'])
    except ValueError:
        errors.append("Invalid molecular weight")
    try:
        float(row['logp'])
    except ValueError:
        errors.append("Invalid LogP value")
    return errors

def import_compound_data(conn, csv_file):
    cursor = conn.cursor()
    with open(csv_file, 'r') as file:
        csv_reader = csv.DictReader(file)
        for row in csv_reader:
            errors = validate_compound_data(row)
            if errors:
                print(f"Skipping row due to errors: {', '.join(errors)}")
                continue
            cursor.execute('''
                INSERT INTO plant_compounds 
                (name, smiles, molecular_weight, logp, plant_source, biological_activity, traditional_use, pubchem_cid, inchi, iupac_name, synonyms)
                VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
            ''', (
                row['name'], row['smiles'], float(row['molecular_weight']), float(row['logp']),
                row['plant_source'], row['biological_activity'], row['traditional_use'],
                int(row['pubchem_cid']), row['inchi'], row['iupac_name'], row['synonyms']
            ))
    conn.commit()
    print("Data imported successfully")

def create_indexes(conn):
    cursor = conn.cursor()
    cursor.execute('CREATE INDEX IF NOT EXISTS idx_name ON plant_compounds(name)')
    cursor.execute('CREATE INDEX IF NOT EXISTS idx_smiles ON plant_compounds(smiles)')
    cursor.execute('CREATE INDEX IF NOT EXISTS idx_plant_source ON plant_compounds(plant_source)')
    conn.commit()
    print("Indexes created successfully")

def search_compounds(conn, query, search_field='name'):
    cursor = conn.cursor()
    cursor.execute(f'''
        SELECT * FROM plant_compounds 
        WHERE {search_field} LIKE ?
    ''', (f'%{query}%',))
    return cursor.fetchall()

class ChemPathAPI:
    def __init__(self, db_path):
        self.conn = sqlite3.connect(db_path)
    
    def get_compound(self, compound_id):
        cursor = self.conn.cursor()
        cursor.execute('SELECT * FROM plant_compounds WHERE id = ?', (compound_id,))
        return cursor.fetchone()
    
    def search_compounds(self, query, search_field='name'):
        return search_compounds(self.conn, query, search_field)
    
    def add_compound(self, compound_data):
        cursor = self.conn.cursor()
        cursor.execute('''
            INSERT INTO plant_compounds 
            (name, smiles, molecular_weight, logp, plant_source, biological_activity, traditional_use, pubchem_cid, inchi, iupac_name, synonyms)
            VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
        ''', (
            compound_data['name'], compound_data['smiles'], compound_data['molecular_weight'],
            compound_data['logp'], compound_data['plant_source'], compound_data['biological_activity'],
            compound_data['traditional_use'], compound_data['pubchem_cid'], compound_data['inchi'],
            compound_data['iupac_name'], compound_data['synonyms']
        ))
        self.conn.commit()
        return cursor.lastrowid
    
    def close(self):
        self.conn.close()

def update_database(conn):
    """Update the database with correct information"""
    cursor = conn.cursor()
    
    # Correct data for each compound
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
        "Lycopene": ("CC(C)=CCC\\C(C)=C\\C\\C=C(\\C)CCC=C(C)C=CC=C(C)C=CC=C(C)C=CC=C(C)C=C", 536.87, 9.16,
                     "Tomatoes, watermelons, pink grapefruits",
                     "Antioxidant effects, Potential cancer prevention",
                     "Free radical scavenging, Modulation of cell signaling pathways",
                     "Traditional consumption as part of Mediterranean diet"),
        "Gingerol": ("CCCCC[C@H](O)CC(=O)CCc1ccc(O)c(OC)c1", 294.39, 3.85,
                     "Ginger root",
                     "Anti-inflammatory effects, Digestive aid",
                     "COX-2 inhibition, Antioxidant activity",
                     "Traditional use in Asian medicine for various ailments")
    }
    
    try:
        conn.execute("BEGIN TRANSACTION")
        for name, data in correct_data.items():
            smiles, molecular_weight, logp, plant_source, therapeutic_areas, biological_activities, traditional_use = data
            
            update_sql = """
            UPDATE plant_compounds
            SET smiles = ?,
                molecular_weight = ?,
                logp = ?,
                plant_source = ?,
                therapeutic_areas = ?,
                biological_activities = ?,
                traditional_use = ?
            WHERE name = ?
            """
            cursor.execute(update_sql, (smiles, molecular_weight, logp, plant_source, therapeutic_areas, biological_activities, traditional_use, name))
            print(f"Updated {name}")
        
        conn.commit()
        print("Database update completed successfully")
    except sqlite3.Error as e:
        conn.rollback()
        print(f"An error occurred: {e}")
    finally:
        cursor.close()

async def main():
    conn = sqlite3.connect("chempath_database.db")
    cursor = conn.cursor()
    cursor.execute("SELECT cid, cmpdname FROM compounds")
    compounds = cursor.fetchall()
    
    await enrich_compound_data(conn, compounds)
    
    display_compound_data(conn, 5)  # Display updated data for Artemisinin
    conn.close()

if __name__ == "__main__":
    asyncio.run(main())

class CustomSmilesType(types.TypeDecorator):
    impl = types.String

    def process_bind_param(self, value, dialect):
        if value is not None:
            # Validate SMILES string here
            if not is_valid_smiles(value):
                raise ValueError("Invalid SMILES string")
        return value

    def process_result_value(self, value, dialect):
        return value

def is_valid_smiles(smiles):
    # Implement SMILES validation logic here
    # This is a placeholder and should be replaced with actual validation
    return bool(smiles and isinstance(smiles, str))

class Compound(Base):
    __tablename__ = 'compounds'

    id = Column(Integer, primary_key=True)
    name = Column(String, nullable=False)
    smiles = Column(CustomSmilesType, nullable=False)
    molecular_weight = Column(Float)
    logp = Column(Float)
    plant_source = Column(String)
    biological_activity = Column(String)
    traditional_use = Column(String)
    pubchem_cid = Column(Integer)
    inchi = Column(String)
    iupac_name = Column(String)
    synonyms = Column(String)

    def __repr__(self):
        return f"<Compound(name='{self.name}', smiles='{self.smiles}')>"

def create_engine_and_session():
    engine = create_engine('sqlite:///chempath_database.db')
    Session = sessionmaker(bind=engine)
    return engine, Session()

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
                        
                        update_data = {
                            'synonyms': ';'.join(compound_dict.get('synonyms', [])),
                            'biological_activity': compound_dict.get('pharmacology', {}).get('description'),
                            'traditional_use': wiki_summary,
                            'molecular_weight': compound_dict.get('molecular_weight'),
                            'iupac_name': compound_dict.get('iupac_name'),
                            'inchi': compound_dict.get('inchi'),
                            'logp': compound_dict.get('xlogp')
                        }
                        
                        update_compound(session, compound.id, update_data)
                        
                        logger.info(f"Updated data for compound {compound.name} (ID: {compound.id})")
                    except Exception as e:
                        logger.error(f"Error updating compound {compound.name} (ID: {compound.id}): {str(e)}")
                else:
                    logger.warning(f"Could not fetch data for compound {compound.name} (ID: {compound.id})")
            
            session.commit()

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
    for compound in compounds:
        print(compound)

    asyncio.run(enrich_compound_data(session, compounds))

    session.close()

if __name__ == "__main__":
    main()

class CompoundSchema(Schema):
    id = fields.Int(dump_only=True)
    name = fields.Str(required=True)
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

@app.route('/compound', methods=['POST'])
def add_compound_route():
    data = request.json
    try:
        compound_data = compound_schema.load(data)
    except ValidationError as err:
        return jsonify(err.messages), 400
    
    compound_id = add_compound(session, compound_data)
    return jsonify({'message': 'Compound added successfully', 'id': compound_id}), 201

@app.route('/compound/<int:compound_id>', methods=['GET'])
def get_compound_route(compound_id):
    compound = get_compound(session, compound_id)
    if compound:
        return jsonify(compound_schema.dump(compound))
    return jsonify({'message': 'Compound not found'}), 404

@app.route('/compounds/search', methods=['GET'])
def search_compounds_route():
    query = request.args.get('query', '')
    search_field = request.args.get('field', 'name')
    compounds = search_compounds(session, query, search_field)
    return jsonify(compounds_schema.dump(compounds))

@app.route('/compound/<int:compound_id>', methods=['PUT'])
def update_compound_route(compound_id):
    data = request.json
    try:
        update_data = compound_schema.load(data, partial=True)
    except ValidationError as err:
        return jsonify(err.messages), 400
    
    success = update_compound(session, compound_id, update_data)
    if success:
        return jsonify({'message': 'Compound updated successfully'})
    return jsonify({'message': 'Compound not found'}), 404

@app.route('/compound/<int:compound_id>', methods=['DELETE'])
def delete_compound_route(compound_id):
    success = delete_compound(session, compound_id)
    if success:
        return jsonify({'message': 'Compound deleted successfully'})
    return jsonify({'message': 'Compound not found'}), 404

@app.route('/compounds/enrich', methods=['POST'])
def enrich_compounds_route():
    compound_ids = request.json.get('compound_ids', [])
    compounds = session.query(Compound).filter(Compound.id.in_(compound_ids)).all()
    asyncio.run(enrich_compound_data(session, compounds))
    return jsonify({'message': 'Compound enrichment process started'})

if __name__ == '__main__':
    engine, session = create_engine_and_session()
    Base.metadata.create_all(engine)
    app.run(debug=True)
