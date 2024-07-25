from chempath_core import create_connection, search_compounds, insert_compound, get_therapeutic_areas, predict_therapeutic_areas

class ChemPathAPI:
    def __init__(self, db_path):
        self.conn = create_connection(db_path)

    def search(self, query=None, filters=None, page=1, per_page=10):
        return search_compounds(self.conn, query, filters, page, per_page)

    def add_compound(self, compound_data):
        return insert_compound(self.conn, compound_data)

    def get_therapeutic_areas(self):
        return get_therapeutic_areas(self.conn)

    def predict_therapeutic_areas(self, smiles):
        all_therapeutic_areas = self.get_therapeutic_areas()
        return predict_therapeutic_areas(smiles, all_therapeutic_areas)

    def close_connection(self):
        self.conn.close()
