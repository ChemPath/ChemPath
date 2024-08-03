import pandas as pd
from sqlalchemy.orm import Session
from models import Molecule

def preprocess_data(session: Session):
    molecules = session.query(Molecule).all()
    df = pd.DataFrame([vars(m) for m in molecules])
    # Implement data cleaning and preparation steps
    return df
