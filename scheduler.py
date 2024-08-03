import threading
import time
from sqlalchemy.orm import Session
from database import engine
from models import Compound, TherapeuticArea
from predict_therapeutic_areas import train_model, get_all_therapeutic_areas
import pandas as pd

def retrain_model():
    with Session(engine) as session:
        compounds = session.query(Compound).all()
        all_therapeutic_areas = get_all_therapeutic_areas()

        data = []
        for compound in compounds:
            row = {
                'smiles': compound.smiles,
                'molecular_weight': compound.molecular_weight,
                'logp': compound.logp,
                'h_bond_donors': compound.h_bond_donors,
                'h_bond_acceptors': compound.h_bond_acceptors,
                'therapeutic_areas': [area.name for area in compound.therapeutic_areas]
            }
            data.append(row)

        df = pd.DataFrame(data)
        train_model(df, all_therapeutic_areas)

    print("Model retrained successfully")

def scheduler_thread():
    while True:
        retrain_model()
        time.sleep(7 * 24 * 60 * 60)  # Sleep for 1 week

def start_scheduler():
    thread = threading.Thread(target=scheduler_thread)
    thread.daemon = True
    thread.start()

if __name__ == '__main__':
    start_scheduler()
