# test_chempath_database.py
import unittest
import sqlite3
from chempath_database import create_connection, insert_compound, insert_compound_class, insert_therapeutic_area, validate_compound_data

class TestChempathDatabase(unittest.TestCase):
    def setUp(self):
        self.conn = create_connection(":memory:")
        self.cursor = self.conn.cursor()
        self.cursor.execute('''
            CREATE TABLE IF NOT EXISTS plant_compounds (
                id INTEGER PRIMARY KEY,
                name TEXT NOT NULL,
                smiles TEXT,
                molecular_weight REAL,
                logp REAL,
                plant_source TEXT,
                biological_activity TEXT,
                traditional_use TEXT
            )
        ''')
        self.cursor.execute('''
            CREATE TABLE IF NOT EXISTS compound_classes (
                id INTEGER PRIMARY KEY,
                name TEXT NOT NULL,
                description TEXT
            )
        ''')
        self.cursor.execute('''
            CREATE TABLE IF NOT EXISTS therapeutic_areas (
                id INTEGER PRIMARY KEY,
                name TEXT NOT NULL,
                description TEXT
            )
        ''')

    def test_insert_compound(self):
        compound = ("Test Compound", "C1=CC=CC=C1", 100.0, 1.5, "Test Plant", "Test Activity", "Test Use")
        compound_id = insert_compound(self.conn, compound)
        self.assertIsNotNone(compound_id)
        
        self.cursor.execute("SELECT * FROM plant_compounds WHERE id = ?", (compound_id,))
        result = self.cursor.fetchone()
        self.assertEqual(result[1], "Test Compound")

    def test_insert_compound_class(self):
        class_id = insert_compound_class(self.conn, "Test Class", "Test Description")
        self.assertIsNotNone(class_id)
        
        self.cursor.execute("SELECT * FROM compound_classes WHERE id = ?", (class_id,))
        result = self.cursor.fetchone()
        self.assertEqual(result[1], "Test Class")

    def test_insert_therapeutic_area(self):
        area_id = insert_therapeutic_area(self.conn, "Test Area", "Test Description")
        self.assertIsNotNone(area_id)
        
        self.cursor.execute("SELECT * FROM therapeutic_areas WHERE id = ?", (area_id,))
        result = self.cursor.fetchone()
        self.assertEqual(result[1], "Test Area")

    def test_validate_compound_data(self):
        valid_data = ("Test", "C1=CC=CC=C1", 100.0, 1.5, "Plant", "Activity", "Use")
        self.assertIsNone(validate_compound_data(*valid_data))
        
        with self.assertRaises(ValueError):
            validate_compound_data("", "C1=CC=CC=C1", 100.0, 1.5, "Plant", "Activity", "Use")
        
        with self.assertRaises(ValueError):
            validate_compound_data("Test", "C1=CC=CC=C1", -100.0, 1.5, "Plant", "Activity", "Use")
    import unittest
import sqlite3
from chempath_database import create_connection, create_table, insert_compound

class TestChempathDatabase(unittest.TestCase):
    def setUp(self):
        self.conn = create_connection(':memory:')
        create_table(self.conn)

    def tearDown(self):
        self.conn.close()

    def test_insert_compound_with_ph_and_temperature(self):
        compound = ("Test Compound", "C1=CC=CC=C1", 100.0, 1.5, "Test Plant", "Test Activity", "Test Use", 7.0, 25.0)
        compound_id = insert_compound(self.conn, compound)
        self.assertIsNotNone(compound_id)

        cursor = self.conn.cursor()
        cursor.execute("SELECT * FROM plant_compounds WHERE id = ?", (compound_id,))
        result = cursor.fetchone()
        self.assertEqual(result[1], "Test Compound")
        self.assertEqual(result[8], 7.0)  # pH
        self.assertEqual(result[9], 25.0)  # temperature

if __name__ == '__main__':
    unittest.main()

    def tearDown(self):
        self.conn.close()

if __name__ == '__main__':
    unittest.main()
