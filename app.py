from flask import Flask, request, jsonify
from sqlalchemy.orm import Session
from database import engine
from models import Compound, TherapeuticArea
from predict_therapeutic_areas import predict_therapeutic_areas, get_all_therapeutic_areas
from retrosynthesis import perform_retrosynthesis, get_retrosynthesis_result
from advanced_retrosynthesis import advanced_retrosynthetic_analysis, get_advanced_retrosynthesis_result
from ai_optimization import optimize_compounds, get_optimization_results
from reagent_availability import analyze_reagent_availability
from scheduler import start_scheduler

app = Flask(__name__)

# Start the scheduler when the app starts
start_scheduler()

@app.route('/predict_therapeutic_areas', methods=['POST'])
def predict_areas():
    data = request.json
    smiles = data.get('smiles')
    ph = data.get('ph')
    temperature = data.get('temperature')

    if not smiles:
        return jsonify({'error': 'SMILES string is required'}), 400

    with Session(engine) as session:
        existing_compound = session.query(Compound).filter_by(smiles=smiles).first()
        if existing_compound:
            return jsonify({'message': 'Compound already exists', 'compound_id': existing_compound.id})

        all_therapeutic_areas = get_all_therapeutic_areas()
        predicted_areas = predict_therapeutic_areas(smiles, all_therapeutic_areas, ph=ph, temperature=temperature)

        compound = Compound(smiles=smiles, ph=ph, temperature=temperature)
        session.add(compound)
        session.flush()

        for area in predicted_areas:
            therapeutic_area = session.query(TherapeuticArea).filter_by(name=area).first()
            if not therapeutic_area:
                therapeutic_area = TherapeuticArea(name=area)
                session.add(therapeutic_area)
            compound.therapeutic_areas.append(therapeutic_area)

        session.commit()

        return jsonify({
            'compound_id': compound.id,
            'smiles': smiles,
            'predicted_therapeutic_areas': predicted_areas,
            'ph': ph,
            'temperature': temperature
        })

@app.route('/retrosynthesis', methods=['POST'])
def retrosynthesis():
    data = request.json
    compound_id = data.get('compound_id')
    smiles = data.get('smiles')
    depth = data.get('depth', 1)

    if not compound_id or not smiles:
        return jsonify({'error': 'Compound ID and SMILES are required'}), 400

    result = perform_retrosynthesis(compound_id, smiles, depth)
    return jsonify(result)

@app.route('/advanced_retrosynthesis', methods=['POST'])
def advanced_retrosynthesis():
    data = request.json
    compound_id = data.get('compound_id')
    smiles = data.get('smiles')
    depth = data.get('depth', 3)

    if not compound_id or not smiles:
        return jsonify({'error': 'Compound ID and SMILES are required'}), 400

    result = advanced_retrosynthetic_analysis(compound_id, smiles, depth)
    return jsonify(result)

@app.route('/optimize_compounds', methods=['POST'])
def optimize():
    data = request.json
    compound_ids = data.get('compound_ids')

    if not compound_ids:
        return jsonify({'error': 'List of compound IDs is required'}), 400

    results = optimize_compounds(compound_ids)
    return jsonify(results)

@app.route('/analyze_reagent', methods=['POST'])
def analyze_reagent():
    data = request.json
    reagent_name = data.get('reagent_name')

    if not reagent_name:
        return jsonify({'error': 'Reagent name is required'}), 400

    result = analyze_reagent_availability(reagent_name)
    return jsonify(result)

if __name__ == '__main__':
    app.run(debug=True)
