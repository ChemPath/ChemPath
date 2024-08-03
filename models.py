from sqlalchemy import Column, Integer, String, Float, Table, ForeignKey
from sqlalchemy.orm import relationship
from sqlalchemy.ext.declarative import declarative_base

Base = declarative_base()

compound_therapeutic_area = Table('compound_therapeutic_area', Base.metadata,
    Column('compound_id', Integer, ForeignKey('compounds.id')),
    Column('therapeutic_area_id', Integer, ForeignKey('therapeutic_areas.id'))
)

class Compound(Base):
    __tablename__ = 'compounds'

    id = Column(Integer, primary_key=True)
    smiles = Column(String, nullable=False)
    ph = Column(Float)
    temperature = Column(Float)
    therapeutic_areas = relationship('TherapeuticArea', secondary=compound_therapeutic_area, back_populates='compounds')

class TherapeuticArea(Base):
    __tablename__ = 'therapeutic_areas'

    id = Column(Integer, primary_key=True)
    name = Column(String, nullable=False, unique=True)
    compounds = relationship('Compound', secondary=compound_therapeutic_area, back_populates='therapeutic_areas')
