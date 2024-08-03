from sqlalchemy import Column, Integer, String, Float, Table, ForeignKey
from sqlalchemy.orm import relationship
from sqlalchemy.ext.declarative import declarative_base

Base = declarative_base()

# Association table for many-to-many relationship between Compound and TherapeuticArea
compound_therapeutic_area = Table('compound_therapeutic_area', Base.metadata,
    Column('compound_id', Integer, ForeignKey('compounds.id')),
    Column('therapeutic_area_id', Integer, ForeignKey('therapeutic_areas.id'))
)

class Molecule(Base):
    __tablename__ = 'molecules'
    id = Column(Integer, primary_key=True)
    smiles = Column(String, unique=True, nullable=False)

class Compound(Base):
    __tablename__ = 'compounds'

    id = Column(Integer, primary_key=True)
    name = Column(String, nullable=False)
    smiles = Column(String)
    molecular_weight = Column(Float)
    logp = Column(Float)
    plant_source = Column(String)
    biological_activity = Column(String)
    traditional_use = Column(String)

    def __repr__(self):
        return f"<Compound(id={self.id}, name='{self.name}', smiles='{self.smiles}')>"

    # Relationships
    advanced_retrosynthesis_results = relationship("AdvancedRetrosynthesisResult", back_populates="compound")
    therapeutic_areas = relationship("TherapeuticArea", secondary=compound_therapeutic_area, back_populates="compounds")


class AdvancedRetrosynthesisResult(Base):
    __tablename__ = 'advanced_retrosynthesis_results'
    id = Column(Integer, primary_key=True)
    compound_id = Column(Integer, ForeignKey('compounds.id'))
    depth = Column(Integer)
    strategy = Column(String)
    num_nodes = Column(Integer)
    num_edges = Column(Integer)
    tree_image = Column(String)


    # Relationships
    compound = relationship("Compound", back_populates="advanced_retrosynthesis_results")


class TherapeuticArea(Base):
    __tablename__ = 'therapeutic_areas'

    id = Column(Integer, primary_key=True)
    name = Column(String, nullable=False)

    # Relationships
    compounds = relationship("Compound", secondary=compound_therapeutic_area, back_populates="therapeutic_areas")
