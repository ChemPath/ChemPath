from flask_sqlalchemy import SQLAlchemy

db = SQLAlchemy()

class TraditionalCompound(db.Model):
    __tablename__ = 'traditional_compounds'
    id = db.Column(db.Integer, primary_key=True)
    name = db.Column(db.String, nullable=False)
    molecular_weight = db.Column(db.Float)
    # Add other columns as needed

class RelatedHerb(db.Model):
    __tablename__ = 'related_herbs'
    id = db.Column(db.Integer, primary_key=True)
    compound_name = db.Column(db.String)
    chinese_name = db.Column(db.String)
    latin_name = db.Column(db.String)

class RelatedTarget(db.Model):
    __tablename__ = 'related_targets'
    id = db.Column(db.Integer, primary_key=True)
    compound_name = db.Column(db.String)
    target_name = db.Column(db.String)
    source = db.Column(db.String)

class RelatedDisease(db.Model):
    __tablename__ = 'related_diseases'
    id = db.Column(db.Integer, primary_key=True)
    target_name = db.Column(db.String)
    disease_name = db.Column(db.String)
