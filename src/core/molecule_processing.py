from rdkit import Chem
from sqlalchemy.orm import Session
from models import Molecule
from sqlalchemy.exc import IntegrityError

def process_molecule(smiles, session):
    try:
        molecule = Molecule(smiles=smiles)
        session.add(molecule)
        session.commit()
        return molecule
    except IntegrityError:
        session.rollback()
        return session.query(Molecule).filter_by(smiles=smiles).first()
