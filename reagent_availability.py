# reagent_availability.py

import csv
from pathlib import Path

class ReagentDatabase:
    def __init__(self):
        self.reagents = {}
        self.load_database()

    def load_database(self):
        database_path = Path(__file__).parent / "reagent_database.csv"
        with open(database_path, 'r') as file:
            reader = csv.DictReader(file)
            for row in reader:
                self.reagents[row['name']] = {
                    'availability': float(row['availability']),
                    'cost': float(row['cost']),
                    'handling': int(row['handling']),
                    'environmental_impact': float(row['environmental_impact'])
                }

    def get_reagent_info(self, reagent_name):
        return self.reagents.get(reagent_name, None)

def calculate_availability_score(reagent_name, database):
    reagent_info = database.get_reagent_info(reagent_name)
    if reagent_info:
        return reagent_info['availability']
    return 0

def estimate_cost(reagent_name, database):
    reagent_info = database.get_reagent_info(reagent_name)
    if reagent_info:
        return reagent_info['cost']
    return float('inf')

def assess_handling_requirements(reagent_name, database):
    reagent_info = database.get_reagent_info(reagent_name)
    if reagent_info:
        handling_levels = ['Low', 'Medium', 'High']
        return handling_levels[reagent_info['handling']]
    return 'Unknown'

def evaluate_environmental_impact(reagent_name, database):
    reagent_info = database.get_reagent_info(reagent_name)
    if reagent_info:
        return reagent_info['environmental_impact']
    return float('inf')

def analyze_reagent_availability(reagent_name):
    database = ReagentDatabase()
    availability = calculate_availability_score(reagent_name, database)
    cost = estimate_cost(reagent_name, database)
    handling = assess_handling_requirements(reagent_name, database)
    environmental_impact = evaluate_environmental_impact(reagent_name, database)

    print(f"Reagent Availability Analysis for {reagent_name}:")
    print(f"Availability Score: {availability:.2f}/10")
    print(f"Estimated Cost: ${cost:.2f}")
    print(f"Handling Requirements: {handling}")
    print(f"Environmental Impact Score: {environmental_impact:.2f}/10")

# Example usage
if __name__ == "__main__":
    analyze_reagent_availability("Sodium borohydride")
