from sqlalchemy.orm import Session
from models import Reagent
from database import engine
from pathlib import Path
import csv

class ReagentDatabase:
    def __init__(self):
        self.load_database()

    def load_database(self):
        database_path = Path(__file__).parent / "reagent_database.csv"
        with Session(engine) as session:
            with open(database_path, 'r') as file:
                reader = csv.DictReader(file)
                for row in reader:
                    reagent = Reagent(
                        name=row['name'],
                        availability=float(row['availability']),
                        cost=float(row['cost']),
                        handling=int(row['handling']),
                        environmental_impact=float(row['environmental_impact'])
                    )
                    session.merge(reagent)
            session.commit()

    def get_reagent_info(self, reagent_name):
        with Session(engine) as session:
            reagent = session.query(Reagent).filter_by(name=reagent_name).first()
            return reagent

def calculate_availability_score(reagent_name, database):
    reagent = database.get_reagent_info(reagent_name)
    return reagent.availability if reagent else 0

def estimate_cost(reagent_name, database):
    reagent = database.get_reagent_info(reagent_name)
    return reagent.cost if reagent else float('inf')

def assess_handling_requirements(reagent_name, database):
    reagent = database.get_reagent_info(reagent_name)
    if reagent:
        handling_levels = ['Low', 'Medium', 'High']
        return handling_levels[reagent.handling]
    return 'Unknown'

def evaluate_environmental_impact(reagent_name, database):
    reagent = database.get_reagent_info(reagent_name)
    return reagent.environmental_impact if reagent else float('inf')

def analyze_reagent_availability(reagent_name):
    database = ReagentDatabase()
    availability = calculate_availability_score(reagent_name, database)
    cost = estimate_cost(reagent_name, database)
    handling = assess_handling_requirements(reagent_name, database)
    environmental_impact = evaluate_environmental_impact(reagent_name, database)

    return {
        "reagent_name": reagent_name,
        "availability_score": round(availability, 2),
        "estimated_cost": round(cost, 2),
        "handling_requirements": handling,
        "environmental_impact_score": round(environmental_impact, 2)
    }

if __name__ == "__main__":
    result = analyze_reagent_availability("Sodium borohydride")
    print(f"Reagent Availability Analysis for {result['reagent_name']}:")
    print(f"Availability Score: {result['availability_score']}/10")
    print(f"Estimated Cost: ${result['estimated_cost']}")
    print(f"Handling Requirements: {result['handling_requirements']}")
    print(f"Environmental Impact Score: {result['environmental_impact_score']}/10")
