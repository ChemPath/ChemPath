# File: chempath_updater.py

import schedule
import time
from chempath_advanced_scraper import scrape_compound
from chempath_database import create_connection, get_all_compounds, update_compound, insert_compound

def validate_compound(compound):
    """Perform basic validation checks on compound data"""
    if not compound[1] or not compound[2] or compound[3] <= 0 or not isinstance(compound[3], (int, float)):
        return False
    return True

def update_database():
    conn = create_connection("chempath_database.db")
    if conn is None:
        return

    # Update existing compounds
    existing_compounds = get_all_compounds(conn)
    for compound in existing_compounds:
        try:
            updated_data = scrape_compound(compound[1])  # compound[1] is the name
            if validate_compound(updated_data):
                update_compound(conn, updated_data)
                print(f"Updated {compound[1]}")
            else:
                print(f"Validation failed for {compound[1]}")
        except Exception as e:
            print(f"Error updating {compound[1]}: {str(e)}")

    # Add new compounds (you'd need to maintain a list of compounds to add)
    new_compounds = ["New Compound 1", "New Compound 2"]  # Example list
    for compound_name in new_compounds:
        try:
            new_data = scrape_compound(compound_name)
            if validate_compound(new_data):
                insert_compound(conn, new_data)
                print(f"Added {compound_name}")
            else:
                print(f"Validation failed for {compound_name}")
        except Exception as e:
            print(f"Error adding {compound_name}: {str(e)}")

    conn.close()

def main():
    schedule.every().week.do(update_database)
    
    while True:
        schedule.run_pending()
        time.sleep(1)

if __name__ == "__main__":
    main()