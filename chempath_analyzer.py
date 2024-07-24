# File: chempath_analyzer.py

from chempath_database import create_connection

def search_compounds(conn, criteria):
    """Search compounds based on given criteria"""
    cursor = conn.cursor()
    query = "SELECT * FROM plant_compounds WHERE "
    conditions = []
    values = []
    for key, value in criteria.items():
        conditions.append(f"{key} LIKE ?")
        values.append(f"%{value}%")
    query += " AND ".join(conditions)
    cursor.execute(query, values)
    return cursor.fetchall()

def analyze_therapeutic_areas(conn):
    """Analyze frequency of therapeutic areas"""
    cursor = conn.cursor()
    cursor.execute("SELECT therapeutic_areas FROM plant_compounds")
    all_areas = cursor.fetchall()
    area_freq = {}
    for areas in all_areas:
        for area in areas[0].split(', '):
            area_freq[area] = area_freq.get(area, 0) + 1
    return area_freq

def main():
    conn = create_connection("chempath_database.db")
    if conn is None:
        return

    # Example search
    search_results = search_compounds(conn, {"therapeutic_areas": "cancer", "biological_activities": "antioxidant"})
    print("Compounds for cancer with antioxidant activity:")
    for result in search_results:
        print(result[1])  # Print compound names

    # Example analysis
    area_frequencies = analyze_therapeutic_areas(conn)
    print("\nTherapeutic Area Frequencies:")
    for area, freq in sorted(area_frequencies.items(), key=lambda x: x[1], reverse=True):
        print(f"{area}: {freq}")

    conn.close()

if __name__ == "__main__":
    main()