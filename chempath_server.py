from flask import Flask, request, jsonify
from chempath_api import ChemPathAPI
import asyncio

app = Flask(__name__)
api = ChemPathAPI("chempath_database.db")

@app.route('/api/process_request', methods=['POST'])
async def process_request():
    data = request.json
    request_type = data.pop('request_type', None)
    
    if not request_type:
        return jsonify({"error": "Missing request_type"}), 400
    
    try:
        result = await api.process_user_request(request_type, **data)
        return jsonify(result)
    except Exception as e:
        return jsonify({"error": str(e)}), 400

@app.route('/api/search_compounds', methods=['GET'])
async def search_compounds():
    query = request.args.get('query', '')
    filters = request.args.get('filters', {})
    page = int(request.args.get('page', 1))
    per_page = int(request.args.get('per_page', 10))
    
    try:
        results = await api.search_compounds(query, filters, page, per_page)
        return jsonify(results)
    except Exception as e:
        return jsonify({"error": str(e)}), 400

@app.route('/api/add_compound', methods=['POST'])
async def add_compound():
    compound_data = request.json
    
    try:
        result = await api.add_compound(compound_data)
        return jsonify(result)
    except Exception as e:
        return jsonify({"error": str(e)}), 400

@app.route('/api/predict_therapeutic_areas', methods=['POST'])
async def predict_therapeutic_areas():
    smiles = request.json.get('smiles')
    
    if not smiles:
        return jsonify({"error": "Missing SMILES string"}), 400
    
    try:
        predictions = await api.predict_therapeutic_areas(smiles)
        return jsonify(predictions)
    except Exception as e:
        return jsonify({"error": str(e)}), 400

@app.route('/api/optimize_structure', methods=['POST'])
async def optimize_structure():
    smiles = request.json.get('smiles')
    
    if not smiles:
        return jsonify({"error": "Missing SMILES string"}), 400
    
    try:
        optimized = await api.optimize_structure(smiles)
        return jsonify({"optimized_smiles": optimized})
    except Exception as e:
        return jsonify({"error": str(e)}), 400

if __name__ == '__main__':
    app.run(debug=True)

