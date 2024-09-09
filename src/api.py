from flask import Flask
from flask_restful import Api, Resource, reqparse
from database.chempath_database import db, PlantCompound, SyntheticCompound, PrescriptionCompound


app = Flask(__name__)
api = Api(app)

class CompoundList(Resource):
    def get(self, compound_type):
        parser = reqparse.RequestParser()
        parser.add_argument('name', type=str)
        args = parser.parse_args()

        if compound_type == 'plant':
            query = PlantCompound.query
        elif compound_type == 'synthetic':
            query = SyntheticCompound.query
        elif compound_type == 'prescription':
            query = PrescriptionCompound.query
        else:
            return {'message': 'Invalid compound type'}, 400

        if args['name']:
            query = query.filter(PlantCompound.name.ilike(f"%{args['name']}%"))

        compounds = query.all()
        return {'compounds': [{'id': c.id, 'name': c.name} for c in compounds]}

api.add_resource(CompoundList, '/compounds/<string:compound_type>')

if __name__ == '__main__':
    app.run(debug=True)
