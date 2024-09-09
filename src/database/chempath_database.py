import sys
import os
from flask import Flask
from flask_sqlalchemy import SQLAlchemy
from sqlalchemy import or_, func, text, Index, cast, String, Column, Integer, String, Float, Text, ARRAY, Enum, ForeignKey
from sqlalchemy.orm import relationship
from sqlalchemy.ext.declarative import declared_attr
from sqlalchemy.exc import IntegrityError
from apscheduler.schedulers.background import BackgroundScheduler
from apscheduler.triggers.cron import CronTrigger
import requests
from sqlalchemy.schema import DDL
from sqlalchemy.event import listen
import logging
from enum import Enum as PythonEnum

app = Flask(__name__)
db_uri = os.environ.get('DATABASE_URI', "postgresql://postgres:Iamabundant28@localhost/chempath_database")
app.config['SQLALCHEMY_DATABASE_URI'] = db_uri
app.config['SQLALCHEMY_TRACK_MODIFICATIONS'] = False
db = SQLAlchemy(app)

def update_database():
    print("Updating database...")
    update_from_phytohub()
    update_from_usda()
    update_from_pubchem()
    print("Database update complete.")

def update_from_phytohub():
    url = "https://phytohub.eu/api/compounds"
    response = requests.get(url)
    if response.status_code == 200:
        compounds = response.json()
        for compound in compounds:
            try:
                plant_compound = PlantCompound(
                    name=compound['name'],
                    phytohub_id=compound['id'],
                    smiles=compound.get('smiles'),
                    inchi_key=compound.get('inchi_key'),
                    chemical_formula=compound.get('chemical_formula'),
                    average_mass=compound.get('average_mass'),
                    source=DatabaseSource.PHYTOHUB
                )
                db.session.add(plant_compound)
                db.session.commit()
            except IntegrityError:
                db.session.rollback()

def update_from_usda():
    # Implementation for USDA database update
    pass

def update_from_pubchem():
    # Implementation for PubChem database update
    pass

scheduler = BackgroundScheduler()
scheduler.add_job(
    update_database,
    trigger=CronTrigger(day_of_week="mon", hour="0", minute="0"),
    id="database_update",
    name="Update database with new compounds",
    replace_existing=True,
)
scheduler.start()

import atexit
atexit.register(lambda: scheduler.shutdown())

class DatabaseSource(PythonEnum):
    PHYTOHUB = 'phytohub'
    CUSTOM = 'custom'
    OTHER = 'other'

class CompoundType(db.Model):
    __tablename__ = 'compound_types'
    id = db.Column(db.Integer, primary_key=True)
    name = db.Column(db.String(50), unique=True, nullable=False)

class BaseCompound:
    @declared_attr
    def id(cls):
        return db.Column(db.Integer, primary_key=True)

    @declared_attr
    def name(cls):
        return db.Column(db.String, nullable=False, index=True)

    @declared_attr
    def systematic_name(cls):
        return db.Column(db.String)

    @declared_attr
    def synonyms(cls):
        return db.Column(db.ARRAY(db.String))

    @declared_attr
    def cas_number(cls):
        return db.Column(db.String, index=True)

    @declared_attr
    def chemical_formula(cls):
        return db.Column(db.String)

    @declared_attr
    def smiles(cls):
        return db.Column(db.String, index=True)

    @declared_attr
    def inchi_key(cls):
        return db.Column(db.String, index=True)

    @declared_attr
    def average_mass(cls):
        return db.Column(db.Float)

    @declared_attr
    def monoisotopic_mass(cls):
        return db.Column(db.Float)

    @declared_attr
    def logp(cls):
        return db.Column(db.Float)

    @declared_attr
    def hydrogen_bond_donors(cls):
        return db.Column(db.Integer)

    @declared_attr
    def hydrogen_bond_acceptors(cls):
        return db.Column(db.Integer)

    @declared_attr
    def polar_surface_area(cls):
        return db.Column(db.Float)

    @declared_attr
    def rotatable_bond_count(cls):
        return db.Column(db.Integer)

    @declared_attr
    def compound_type_id(cls):
        return db.Column(db.Integer, db.ForeignKey('compound_types.id'), nullable=False)

    @declared_attr
    def compound_type(cls):
        return db.relationship('CompoundType', backref=db.backref('compounds', lazy=True))

class PlantCompound(db.Model):
    __tablename__ = 'plant_compounds'

    id = db.Column(db.Integer, primary_key=True)
    name = db.Column(db.String, nullable=False, index=True)
    source = db.Column(Enum(DatabaseSource), nullable=False)
    systematic_name = db.Column(db.String)
    synonyms = db.Column(db.ARRAY(db.String))
    cas_number = db.Column(db.String, index=True)
    chemical_formula = db.Column(db.String)
    smiles = db.Column(db.String, index=True)
    inchi_key = db.Column(db.String, index=True)
    average_mass = db.Column(db.Float)
    monoisotopic_mass = db.Column(db.Float)
    logp = db.Column(db.Float)
    hydrogen_bond_donors = db.Column(db.Integer)
    hydrogen_bond_acceptors = db.Column(db.Integer)
    polar_surface_area = db.Column(db.Float)
    rotatable_bond_count = db.Column(db.Integer)
    compound_type_id = db.Column(db.Integer, db.ForeignKey('compound_types.id'), nullable=False)
    compound_type = db.relationship('CompoundType', backref=db.backref('compounds', lazy=True))

    phytohub_id = db.Column(db.String, unique=True)
    taxonomy_family = db.Column(db.String)
    taxonomy_class = db.Column(db.String)
    taxonomy_subclass = db.Column(db.String)
    food_sources = db.Column(db.Text)
    biomarker_roles = db.Column(db.Text)
    metabolism = db.Column(db.Text)
    inter_individual_variations = db.Column(db.Text)

    text_search = db.Column(db.Text)

    __table_args__ = (
        db.Index('idx_plant_compound_name_inchi', 'name', 'inchi_key'),
        db.Index('idx_plant_compound_fts', 'text_search', postgresql_using='gin', postgresql_ops={'text_search': 'gin_trgm_ops'}),
    )

    @classmethod
    def create_text_search_trigger(cls):
        trigger_sql = DDL('''
            CREATE FUNCTION update_plant_compound_text_search() RETURNS TRIGGER AS $$
            BEGIN
                NEW.text_search := COALESCE(NEW.name, '') || ' ' || COALESCE(NEW.systematic_name, '') || ' ' || COALESCE(array_to_string(NEW.synonyms, ' '), '');
                RETURN NEW;
            END;
            $$ LANGUAGE plpgsql;

            CREATE TRIGGER update_plant_compound_text_search_trigger
            BEFORE INSERT OR UPDATE ON plant_compounds
            FOR EACH ROW EXECUTE FUNCTION update_plant_compound_text_search();
        ''')
        listen(cls.__table__, 'after_create', trigger_sql)

    @classmethod
    def find_by_taxonomy(cls, taxonomy):
        return cls.query.filter(
            or_(
                cls.taxonomy_family.ilike(f"%{taxonomy}%"),
                cls.taxonomy_class.ilike(f"%{taxonomy}%"),
                cls.taxonomy_subclass.ilike(f"%{taxonomy}%")
            )
        ).all()

    @classmethod
    def find_by_food_source(cls, food_source):
        return cls.query.filter(cls.food_sources.ilike(f"%{food_source}%")).all()

class SyntheticCompound(db.Model):
    __tablename__ = 'synthetic_compounds'

    id = db.Column(db.Integer, primary_key=True)
    name = db.Column(db.String, nullable=False, index=True)
    chemical_formula = db.Column(db.String)
    smiles = db.Column(db.String, index=True)
    inchi_key = db.Column(db.String, index=True)
    average_mass = db.Column(db.Float)
    monoisotopic_mass = db.Column(db.Float)
    logp = db.Column(db.Float)
    hydrogen_bond_donors = db.Column(db.Integer)
    hydrogen_bond_acceptors = db.Column(db.Integer)
    polar_surface_area = db.Column(db.Float)
    rotatable_bond_count = db.Column(db.Integer)
    compound_type_id = db.Column(db.Integer, db.ForeignKey('compound_types.id'), nullable=False)
    compound_type = db.relationship('CompoundType', backref=db.backref('compounds', lazy=True))

    synthesis_method = db.Column(db.String)
    purity = db.Column(db.Float)

class PrescriptionCompound(BaseCompound, db.Model):
    __tablename__ = 'prescription_compounds'

    id = db.Column(db.Integer, primary_key=True)
    brand_name = db.Column(db.String)
    active_ingredient = db.Column(db.String)
    dosage_form = db.Column(db.String)

def create_compound_types():
    compound_types = ['Plant', 'Synthetic', 'Prescription']
    for type_name in compound_types:
        existing_type = CompoundType.query.filter_by(name=type_name).first()
        if not existing_type:
            compound_type = CompoundType(name=type_name)
            db.session.add(compound_type)
    db.session.commit()
    logging.info("Compound types created or updated successfully")

def insert_sample_data():
    plant_type = CompoundType.query.filter_by(name='Plant').first()
    plant_compound = PlantCompound(
        name='Quercetin',
                systematic_name='2-(3,4-dihydroxyphenyl)-3,5,7-trihydroxy-4H-chromen-4-one',
        synonyms=['3,3\',4\',5,7-Pentahydroxyflavone', 'Sophoretin'],
        phytohub_id='PH000001',
        chemical_formula='C15H10O7',
        smiles='O=c1c(O)c(-c2ccc(O)c(O)c2)oc2cc(O)cc(O)c12',
        inchi_key='REFJWTPEDVJJIY-UHFFFAOYSA-N',
        average_mass=302.236,
        monoisotopic_mass=302.0427,
        logp=1.5,
        hydrogen_bond_donors=5,
        hydrogen_bond_acceptors=7,
        polar_surface_area=131.36,
        rotatable_bond_count=1,
        compound_type=plant_type,
        taxonomy_family='Flavonoids',
        food_sources='Apples, onions, tea',
        source=DatabaseSource.OTHER
    )
    db.session.add(plant_compound)

    synthetic_type = CompoundType.query.filter_by(name='Synthetic').first()
    synthetic_compound = SyntheticCompound(
        name='Aspirin',
        chemical_formula='C9H8O4',
        smiles='CC(=O)OC1=CC=CC=C1C(=O)O',
        inchi_key='BSYNRYMUTXBXSQ-UHFFFAOYSA-N',
        average_mass=180.157,
        monoisotopic_mass=180.0423,
        logp=1.19,
        hydrogen_bond_donors=1,
        hydrogen_bond_acceptors=4,
        polar_surface_area=63.6,
        rotatable_bond_count=3,
        compound_type=synthetic_type,
        synthesis_method='Reaction of salicylic acid with acetic anhydride',
        purity=99.5
    )
    db.session.add(synthetic_compound)

    prescription_type = CompoundType.query.filter_by(name='Prescription').first()
    prescription_compound = PrescriptionCompound(
        name='Ibuprofen',
        chemical_formula='C13H18O2',
        smiles='CC(C)Cc1ccc(cc1)[C@H](C)C(=O)O',
        inchi_key='HEFNNWSXXWATRW-GPIVLXJGSA-N',
        average_mass=206.285,
        monoisotopic_mass=206.1307,
        logp=3.97,
        hydrogen_bond_donors=1,
        hydrogen_bond_acceptors=2,
        polar_surface_area=37.3,
        rotatable_bond_count=4,
        compound_type=prescription_type,
        brand_name='Advil',
        active_ingredient='Ibuprofen',
        dosage_form='Tablet'
    )
    db.session.add(prescription_compound)

    db.session.commit()
    logging.info("Sample data inserted successfully")

def main():
    logging.basicConfig(level=logging.DEBUG)
    with app.app_context():
        db.drop_all()  # Drop all existing tables
        db.create_all()  # Create tables based on current model definitions
        db.session.execute(text("""
        CREATE OR REPLACE FUNCTION to_tsvector_immutable(text) RETURNS tsvector AS $$
        BEGIN
            RETURN to_tsvector('english', $1);
        END;
        $$ LANGUAGE plpgsql IMMUTABLE;
    """))
        db.session.commit()
        create_compound_types()
        insert_sample_data()

        # Query and print details
        plant_compounds = PlantCompound.query.all()
        synthetic_compounds = SyntheticCompound.query.all()
        prescription_compounds = PrescriptionCompound.query.all()

        print(f"Number of plant compounds: {len(plant_compounds)}")
        print(f"Number of synthetic compounds: {len(synthetic_compounds)}")
        print(f"Number of prescription compounds: {len(prescription_compounds)}")

        for compound in plant_compounds:
            print(f"Plant Compound: {compound.name}, Type: {compound.compound_type.name}")
        for compound in synthetic_compounds:
            print(f"Synthetic Compound: {compound.name}, Type: {compound.compound_type.name}")
        for compound in prescription_compounds:
            print(f"Prescription Compound: {compound.name}, Type: {compound.compound_type.name}")

if __name__ == '__main__':
    main()

