import sys
import os
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import pytest
from core.molecule_processing import process_molecule
from core.property_calculator import calculate_properties
from core.similarity_search import similarity_search
from core.visualization_utils import visualize_molecule
from models import Base
from sqlalchemy import create_engine
from sqlalchemy.orm import sessionmaker

@pytest.fixture(scope="module")
def session():
    engine = create_engine('sqlite:///:memory:')
    Base.metadata.create_all(engine)
    Session = sessionmaker(bind=engine)
    return Session()

@pytest.fixture(autouse=True)
def reset_session(session):
    session.rollback()
    yield
    session.close()

def test_process_molecule(session):
    smiles = "CC(=O)OC1=CC=CC=C1C(=O)O"
    molecule = process_molecule(smiles, session)
    assert molecule is not None
    assert molecule.smiles == smiles

def test_calculate_properties(session):
    smiles = "CC(=O)OC1=CC=CC=C1C(=O)O"
    molecule = process_molecule(smiles, session)
    properties = calculate_properties(molecule.id, session)
    assert 'molecular_weight' in properties
    assert 'logp' in properties

def test_similarity_search(session):
    query_smiles = "CC(=O)OC1=CC=CC=C1C(=O)O"
    process_molecule(query_smiles, session)
    results = similarity_search(query_smiles, session, threshold=0.7)
    assert len(results) > 0
    assert all(similarity >= 0.7 for _, similarity in results)

def test_visualize_molecule():
    smiles = "CC(=O)OC1=CC=CC=C1C(=O)O"
    image_data = visualize_molecule(smiles)
    assert image_data.startswith('iVBORw0KGgoAAAANSUhEUgAAASwAAAEsCAIAAAD2HxkiAAAsX0lEQVR4nO3deVxU9f')