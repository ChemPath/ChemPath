from sqlalchemy import Column, Integer, String, Float
from sqlalchemy.ext.declarative import declarative_base

Base = declarative_base()

class Compound(Base):
    __tablename__ = 'compounds'

    id = Column(Integer, primary_key=True)
    name = Column(String, nullable=False)
    smiles = Column(String, nullable=False)
    molecular_weight = Column(Float)
    logp = Column(Float)
    # Add other fields as needed
