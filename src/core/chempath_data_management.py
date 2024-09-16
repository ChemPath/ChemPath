# File: chempath_data_management.py

import sqlite3
import schedule
import time
import asyncio
import aiohttp
from src.database.chempath_database import create_connection, get_all_compounds, update_compound, insert_compound
import logging

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

def get_table_info(conn):
    cursor = conn.cursor()
    cursor.execute("PRAGMA table_info(plant_compounds)")
    return cursor.fetchall()

def print_table_structure(table_info):
    logger.info("Current table structure:")
    for column in table_info:
        logger.info(f"Column {column[0]}: {column[1]} ({column[2]})")

def correct_data(conn):
    cursor = conn.cursor()
    table_info = get_table_info(conn)
    print_table_structure(table_info)
    
    cursor.execute("SELECT * FROM plant_compounds")
    rows = cursor.fetchall()
    
    for row in rows:
        logger.info(f"Original row: {row}")
        
        if len(row) != len(table_info):
            logger.warning(f"Warning: Row has {len(row)} columns, expected {len(table_info)}")
            continue
        
        id, name, smiles, molecular_weight, logp, plant_source, therapeutic_areas, biological_activities, traditional_use = row[:9]
        
        corrected_data = (
            name,
            plant_source,
            float(smiles) if smiles and smiles.replace('.', '').isdigit() else 0,
            float(traditional_use) if traditional_use and traditional_use.replace('.', '').isdigit() else 0,
            molecular_weight,
            therapeutic_areas,
            biological_activities,
            logp
        )
        
        logger.info(f"Corrected data: {corrected_data}")
        
        update_sql = """
        UPDATE plant_compounds
        SET smiles = ?, 
            molecular_weight = ?, 
            logp = ?,
            plant_source = ?,
            therapeutic_areas = ?,
            biological_activities = ?,
            traditional_use = ?
        WHERE id = ?
        """
        cursor.execute(update_sql, corrected_data[1:] + (id,))
    
    conn.commit()
    logger.info("Data correction completed")

def reorganize_data(conn):
    cursor = conn.cursor()
    
    cursor.execute("SELECT * FROM plant_compounds")
    rows = cursor.fetchall()
    
    for row in rows:
        id, name, smiles, molecular_weight, logp, plant_source, biological_activity, traditional_use, therapeutic_areas, biological_activities = row
        
        reorganized_data = (
            name,
            logp,
            float(smiles) if smiles and smiles.replace('.', '').isdigit() else 0,
            float(therapeutic_areas) if therapeutic_areas and therapeutic_areas.replace('.', '', 1).replace('-', '').isdigit() else 0,
            molecular_weight,
            biological_activities,
            biological_activity,
            traditional_use
        )
        
        logger.info(f"Reorganized data for {name}: {reorganized_data}")
        
        update_sql = """
        UPDATE plant_compounds
        SET smiles = ?, 
            molecular_weight = ?, 
            logp = ?,
            plant_source = ?,
            therapeutic_areas = ?,
            biological_activities = ?,
            traditional_use = ?
        WHERE id = ?
        """
        cursor.execute(update_sql, reorganized_data[1:] + (id,))
    
    conn.commit()
    logger.info("Data reorganization completed")

def definitive_data_correction(conn):
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
    }
    
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
        logger.info(f"Updated {name}")
    
    conn.commit()
    logger.info("Definitive data correction completed")

def final_structure_correction(conn):
    cursor = conn.cursor()
    
    try:
        cursor.execute("ALTER TABLE plant_compounds RENAME COLUMN biological_activity TO biological_activities")
        logger.info("Renamed 'biological_activity' to 'biological_activities'")
    except sqlite3.OperationalError:
        logger.info("Column 'biological_activities' already exists")

def inspect_database(conn):
    cursor = conn.cursor()
    
    cursor.execute("PRAGMA table_info(plant_compounds)")
    columns = cursor.fetchall()
    logger.info("Database Structure:")
    for column in columns:
        logger.info(f"Column: {column[1]}, Type: {column[2]}")

    cursor.execute("SELECT * FROM plant_compounds")
    rows = cursor.fetchall()
    logger.info("Current Data:")
    for row in rows:
        logger.info(row)

def swap_and_correct_data(conn):
    cursor = conn.cursor()
    
    correct_data = {
        "Quercetin": ("Antioxidant effects, Anti-inflammatory properties",
                      "Free radical scavenging, Enzyme inhibition", 
                      "Traditional Chinese Medicine for cardiovascular health"),
        "Curcumin": ("Anti-inflammatory effects, Cancer prevention",
                     "NF-κB inhibition, Antioxidant activity",
                     "Ayurvedic medicine for various ailments"),
    }
    
    for name, data in correct_data.items():
        therapeutic_areas, biological_activities, traditional_use = data
        
        update_sql = """
        UPDATE plant_compounds
        SET therapeutic_areas = ?,
            biological_activities = ?,
            traditional_use = ?
        WHERE name = ?
        """
        cursor.execute(update_sql, (therapeutic_areas, biological_activities, traditional_use, name))
    
    conn.commit()
    logger.info("Data swap and correction completed")

def ultimate_data_correction(conn):
    cursor = conn.cursor()
    
    correct_data = {
        "Quercetin": ("Antioxidant effects, Anti-inflammatory properties",
                      "Free radical scavenging, Enzyme inhibition", 
                      "Traditional Chinese Medicine for cardiovascular health"),
        "Curcumin": ("Anti-inflammatory effects, Cancer prevention",
                     "NF-κB inhibition, Antioxidant activity",
                     "Ayurvedic medicine for various ailments"),
    }
    
    for name, data in correct_data.items():
        therapeutic_areas, biological_activities, traditional_use = data
        
        update_sql = """
        UPDATE plant_compounds
        SET therapeutic_areas = ?,
            biological_activities = ?,
            traditional_use = ?
        WHERE name = ?
        """
        cursor.execute(update_sql, (therapeutic_areas, biological_activities, traditional_use, name))
    
    conn.commit()
    logger.info("Ultimate data correction completed")

async def fetch_compound_data(session, compound_name):
    url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{compound_name}/JSON"
    async with session.get(url) as response:
        if response.status == 200:
            return await response.json()
        return None

async def update_database():
    conn = create_connection("chempath_database.db")
    if conn is None:
        return

    existing_compounds = get_all_compounds(conn)
    async with aiohttp.ClientSession() as session:
        for compound in existing_compounds:
            try:
                updated_data = await fetch_compound_data(session, compound[1])  # compound[1] is the name
                if updated_data:
                    # Process and update the compound data
                    update_compound(conn, updated_data)
                    logger.info(f"Updated {compound[1]}")
                else:
                    logger.warning(f"No data found for {compound[1]}")
            except Exception as e:
                logger.error(f"Error updating {compound[1]}: {str(e)}")

    conn.close()

def schedule_updates():
    schedule.every().week.do(asyncio.run, update_database())
    
    while True:
        schedule.run_pending()
        time.sleep(1)

def main():
    conn = create_connection("chempath_database.db")
    if conn is not None:
        logger.info("Before corrections:")
        inspect_database(conn)
        
        correct_data(conn)
        reorganize_data(conn)
        definitive_data_correction(conn)
        final_structure_correction(conn)
        swap_and_correct_data(conn)
        ultimate_data_correction(conn)
        
        logger.info("After corrections:")
        inspect_database(conn)
        
        conn.close()
        
        # Start the scheduled updates
        schedule_updates()
    else:
        logger.error("Error! Cannot create the database connection.")

if __name__ == "__main__":
    main()
