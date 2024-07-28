from flask import Flask, request, jsonify
from your_database_module import add_compound, get_compound, add_predicted_areas
from predict_therapeutic_areas import predict_therapeutic_areas
from scheduler import start_scheduler

app = Flask(__name__)

# Start the scheduler when the app starts
start_scheduler()

@app.route('/predict_therapeutic_areas', methods=['POST'])
def predict_areas():
    data = request.json
    smiles = data.get('smiles')

    if not smiles:
        return jsonify({'error': 'SMILES string is required'}), 400

    # Check if the compound already exists in the database
    existing_compound = get_compound(smiles)
    if existing_compound:
        return jsonify({'message': 'Compound already exists', 'compound_id': existing_compound['id']})

    

from flask import request, jsonify
from chempath_database import get_compound, add_compound, add_predicted_areas
from predict_therapeutic_areas import predict_therapeutic_areas
from chempath_core import fetch_all_therapeutic_areas

def predict_areas():
    data = request.json
    smiles = data.get('smiles')
    ph = data.get('ph')
    temperature = data.get('temperature')

    if not smiles:
        return jsonify({'error': 'SMILES string is required'}), 400

    if ph is not None:
        try:
            ph = float(ph)
            if ph < 0 or ph > 14:
                return jsonify({'error': 'pH must be between 0 and 14'}), 400
        except ValueError:
            return jsonify({'error': 'Invalid pH value'}), 400

    if temperature is not None:
        try:
            temperature = float(temperature)
            if temperature < -273.15:
                return jsonify({'error': 'Temperature must be above absolute zero (-273.15°C)'}), 400
        except ValueError:
            return jsonify({'error': 'Invalid temperature value'}), 400

    # Check if the compound already exists in the database
    existing_compound = get_compound(smiles)
    if existing_compound:
        return jsonify({'message': 'Compound already exists', 'compound_id': existing_compound['id']})

    # Predict therapeutic areas
    all_therapeutic_areas = fetch_all_therapeutic_areas()
    predicted_areas = predict_therapeutic_areas(smiles, all_therapeutic_areas, ph=ph, temperature=temperature)

    # Add the new compound to the database
    compound_id = add_compound(smiles, ph=ph, temperature=temperature)

    # Add predicted areas to the database
    add_predicted_areas(compound_id, smiles, predicted_areas)

    return jsonify({
        'compound_id': compound_id,
        'smiles': smiles,
        'predicted_therapeutic_areas': predicted_areas,
        'ph': ph,
        'temperature': temperature
    })


def fetch_all_therapeutic_areas():
    # Implement this function to get all possible therapeutic areas
    # This could involve querying the database or loading from a predefined list
    pass





    # Add the new compound to the database
    compound_id = add_compound(smiles)

    # Add predicted areas to the database
    add_predicted_areas(compound_id, smiles, predicted_areas)

    return jsonify({
        'compound_id': compound_id,
        'smiles': smiles,
        'predicted_therapeutic_areas': predicted_areas
    })

if __name__ == '__main__':
    app.run(debug=True)
