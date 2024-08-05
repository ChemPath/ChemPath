from flask import Flask, request, jsonify
from marshmallow import Schema, fields, ValidationError
from database_operations import (
    create_engine_and_session, add_compound, get_compound,
    search_compounds, update_compound, delete_compound, enrich_compound_data
)
import asyncio

app = Flask(__name__)
engine, session = create_engine_and_session()
import sqlite3
from rdkit import Chem
from rdkit.Chem import AllChem
from src.core.chempath_core import optimize_structure, predict_properties, retrosynthesize

class ChemPathAPI:
    def __init__(self, db_path):
        self.db_path = db_path

    def search_compounds(self, query):
        conn = sqlite3.connect(self.db_path)
        cursor = conn.cursor()
        query = f"%{query}%"
        cursor.execute("""
            SELECT name, smiles, molecular_weight, logp 
            FROM plant_compounds 
            WHERE name LIKE ? OR smiles LIKE ?
        """, (query, query))
        results = cursor.fetchall()
        total_count = len(results)
        conn.close()
        return results, total_count

    def optimize_structure(self, smiles):
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return "Invalid SMILES"
        optimized_mol = optimize_structure(mol)
        return Chem.MolToSmiles(optimized_mol)

    def perform_retrosynthesis(self, smiles):
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return {"error": "Invalid SMILES"}
        retro_steps = retrosynthesize(mol)
        return {"steps": retro_steps}

    def predict_properties(self, smiles):
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return {"error": "Invalid SMILES"}
        properties = predict_properties(mol)
        return properties

class CompoundSchema(Schema):
    id = fields.Int(dump_only=True)
    name = fields.Str(required=True)
    smiles = fields.Str(required=True)
    molecular_weight = fields.Float()
    logp = fields.Float()
    plant_source = fields.Str()
    biological_activity = fields.Str()
    traditional_use = fields.Str()
    pubchem_cid = fields.Int()
    inchi = fields.Str()
    iupac_name = fields.Str()
    synonyms = fields.Str()
    polararea = fields.Float()
    complexity = fields.Float()
    heavycnt = fields.Int()
    hbonddonor = fields.Int()
    hbondacc = fields.Int()
    rotbonds = fields.Int()
    exactmass = fields.Float()
    monoisotopicmass = fields.Float()

compound_schema = CompoundSchema()
compounds_schema = CompoundSchema(many=True)

@app.route('/compound', methods=['POST'])
def add_compound_route():
    data = request.json
    try:
        compound_data = compound_schema.load(data)
    except ValidationError as err:
        return jsonify(err.messages), 400
    
    compound_id = add_compound(session, compound_data)
    return jsonify({'message': 'Compound added successfully', 'id': compound_id}), 201

@app.route('/compound/<int:compound_id>', methods=['GET'])
def get_compound_route(compound_id):
    compound = get_compound(session, compound_id)
    if compound:
        return jsonify(compound_schema.dump(compound))
    return jsonify({'message': 'Compound not found'}), 404

@app.route('/compounds/search', methods=['GET'])
def search_compounds_route():
    query = request.args.get('query', '')
    search_field = request.args.get('field', 'name')
    compounds = search_compounds(session, query, search_field)
    return jsonify(compounds_schema.dump(compounds))

@app.route('/compound/<int:compound_id>', methods=['PUT'])
def update_compound_route(compound_id):
    data = request.json
    try:
        update_data = compound_schema.load(data, partial=True)
    except ValidationError as err:
        return jsonify(err.messages), 400
    
    success = update_compound(session, compound_id, update_data)
    if success:
        return jsonify({'message': 'Compound updated successfully'})
    return jsonify({'message': 'Compound not found'}), 404

@app.route('/compound/<int:compound_id>', methods=['DELETE'])
def delete_compound_route(compound_id):
    success = delete_compound(session, compound_id)
    if success:
        return jsonify({'message': 'Compound deleted successfully'})
    return jsonify({'message': 'Compound not found'}), 404

@app.route('/compounds/enrich', methods=['POST'])
def enrich_compounds_route():
    compound_ids = request.json.get('compound_ids', [])
    compounds = [get_compound(session, cid) for cid in compound_ids]
    asyncio.run(enrich_compound_data(session, compounds))
    return jsonify({'message': 'Compound enrichment process started'})

if __name__ == '__main__':
    app.run(debug=True)
