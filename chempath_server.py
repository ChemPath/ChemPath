from flask import Flask, request, jsonify
from chempath_api import ChemPathAPI

app = Flask(__name__)
api = ChemPathAPI("chempath_database.db")

@app.route('/api/process_request', methods=['POST'])
def process_request():
    data = request.json
    request_type = data.pop('request_type', None)
    
    if not request_type:
        return jsonify({"error": "Missing request_type"}), 400
    
    try:
        future = api.process_user_request(request_type, **data)
        result = future.result()  # Wait for the Future to complete
        return jsonify(result)
    except Exception as e:
        return jsonify({"error": str(e)}), 400



if __name__ == '__main__':
    app.run(debug=True)
