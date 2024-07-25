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

    # Predict therapeutic areas
    all_therapeutic_areas = fetch_all_therapeutic_areas()  # Implement this function to get all possible areas
    predicted_areas = predict_therapeutic_areas(smiles, all_therapeutic_areas)

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
