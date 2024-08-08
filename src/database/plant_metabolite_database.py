from sqlalchemy import create_engine, Column, Integer, String, Text, Float
from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy.orm import sessionmaker
from sqlalchemy import ForeignKey
from sqlalchemy.orm import relationship


Base = declarative_base()

class PlantMetabolite(Base):
    __tablename__ = 'plant_metabolites'

    id = Column(Integer, primary_key=True)
    name = Column(String, nullable=False)
    smiles = Column(String, nullable=False)
    molecular_weight = Column(Float)
    biosynthetic_pathway_id = Column(Integer, ForeignKey('biosynthetic_pathways.id'))
    biosynthetic_pathway = relationship("BiosyntheticPathway", back_populates="metabolites")

class BiosyntheticPathway(Base):
    __tablename__ = 'biosynthetic_pathways'

    id = Column(Integer, primary_key=True)
    name = Column(String, nullable=False)
    metabolites = relationship("PlantMetabolite", back_populates="biosynthetic_pathway")

def create_engine_and_session():
    engine = create_engine('sqlite:///plant_metabolites.db')
    Base.metadata.create_all(engine)
    Session = sessionmaker(bind=engine)
    return engine, Session()

def initialize_database():
    engine, _ = create_engine_and_session()
    Base.metadata.create_all(engine)

def add_metabolite(session, name, smiles, molecular_weight, pathway_name):
    pathway = session.query(BiosyntheticPathway).filter_by(name=pathway_name).first()
    if not pathway:
        pathway = BiosyntheticPathway(name=pathway_name)
        session.add(pathway)
    
    metabolite = PlantMetabolite(name=name, smiles=smiles, molecular_weight=molecular_weight, biosynthetic_pathway=pathway)
    session.add(metabolite)
    session.commit()

def get_metabolite_by_smiles(session, smiles):
    return session.query(PlantMetabolite).filter_by(smiles=smiles).first()

def insert_metabolite(session, metabolite):
    session.add(metabolite)
    session.commit()
    return metabolite.id


def search_metabolites(session, query, search_field='name'):
    return session.query(PlantMetabolite).filter(getattr(PlantMetabolite, search_field).like(f"%{query}%")).all()

    

if __name__ == "__main__":
    initialize_database()
    engine, session = create_engine_and_session()

    # Example usage
    sample_metabolite = {
        'name': 'Quercetin',
        'smiles': 'O=C1C(O)=C(O)C(=O)C2=C1C=C(O)C(O)=C2O',
        'biosynthetic_pathway': 'Flavonoid biosynthesis',
        'plant_source': 'Various fruits and vegetables',
        'molecular_weight': 302.24,
        'known_activities': 'Antioxidant, anti-inflammatory'
    }

    metabolite_id = insert_metabolite(session, sample_metabolite)
    print(f"Inserted metabolite with ID: {metabolite_id}")

    retrieved_metabolite = get_metabolite_by_smiles(session, sample_metabolite['smiles'])
    if retrieved_metabolite:
        print(f"Retrieved metabolite: {retrieved_metabolite.name}")

    session.close()
