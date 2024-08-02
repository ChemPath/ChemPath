import sqlite3
import csv
import requests
from bs4 import BeautifulSoup
import pubchempy as pcp

def create_connection(db_file):
    conn = None
    try:
        conn = sqlite3.connect(db_file)
        print(f"Connected to SQLite database: {db_file}")
        return conn
    except sqlite3.Error as e:
        print(e)
    return conn

def create_tables(conn):
    cursor = conn.cursor()
    
    cursor.execute('''
    CREATE TABLE IF NOT EXISTS compounds (
        cid INTEGER PRIMARY KEY,
        cmpdname TEXT,
        cmpdsynonym TEXT,
        molecular_weight REAL,
        molecular_formula TEXT,
        polararea REAL,
        complexity INTEGER,
        xlogp REAL,
        heavycnt INTEGER,
        hbonddonor INTEGER,
        hbondacc INTEGER,
        rotbonds INTEGER,
        inchi TEXT,
        isosmiles TEXT,
        canonicalsmiles TEXT,
        inchikey TEXT,
        iupacname TEXT,
        exactmass REAL,
        monoisotopicmass REAL,
        charge INTEGER,
        covalentunitcnt INTEGER,
        isotopeatomcnt INTEGER,
        totalatomstereocnt INTEGER,
        pclidcnt INTEGER,
        gpidcnt INTEGER,
        gpfamilycnt INTEGER,
        neighbortype TEXT,
        meshheadings TEXT,
        annothits TEXT,
        annotation TEXT,
        plant_source TEXT,
        biological_activity TEXT,
        traditional_use TEXT
    )
    ''')
    
    cursor.execute('''
    CREATE TABLE IF NOT EXISTS plant_compounds (
        id INTEGER PRIMARY KEY AUTOINCREMENT,
        cid INTEGER,
        plant_name TEXT,
        additional_info TEXT,
        FOREIGN KEY (cid) REFERENCES compounds (cid)
    )
    ''')
    
    conn.commit()


def alter_table(conn):
    cursor = conn.cursor()
    try:
        cursor.execute("ALTER TABLE compounds ADD COLUMN plant_source TEXT")
        cursor.execute("ALTER TABLE compounds ADD COLUMN biological_activity TEXT")
        cursor.execute("ALTER TABLE compounds ADD COLUMN traditional_use TEXT")
        conn.commit()
        print("Table 'compounds' altered successfully")
    except sqlite3.OperationalError:
        print("Columns already exist")
def import_pubchem_data(conn, csv_file_path):
    cursor = conn.cursor()
    with open(csv_file_path, 'r', encoding='utf-8-sig') as csv_file:
        csv_reader = csv.DictReader(csv_file)
        print("CSV columns:", csv_reader.fieldnames)
        for row in csv_reader:
            print(f"CID: {row[' cid']}")
            print(f"pclidcnt: {row['pclidcnt']}")
            print(f"gpidcnt: {row['gpidcnt']}")
            print(f"gpfamilycnt: {row['gpfamilycnt']}")
            print(f"neighbortype: {row['neighbortype']}")
            print(f"meshheadings: {row['meshheadings']}")
            print(f"annothits: {row['annothits']}")
            
            cursor.execute('''
            INSERT OR REPLACE INTO compounds (
                cid, cmpdname, cmpdsynonym, molecular_weight, molecular_formula, polararea,
                complexity, xlogp, heavycnt, hbonddonor, hbondacc, rotbonds,
                inchi, isosmiles, canonicalsmiles, inchikey, iupacname,
                exactmass, monoisotopicmass, charge, covalentunitcnt,
                isotopeatomcnt, totalatomstereocnt, pclidcnt, gpidcnt, gpfamilycnt,
                neighbortype, meshheadings, annothits, annotation
            ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
            ''', (
                row[' cid'], row['cmpdname'], row['cmpdsynonym'], row['mw'],
                row['mf'], row['polararea'], row['complexity'], row['xlogp'],
                row['heavycnt'], row['hbonddonor'], row['hbondacc'], row['rotbonds'],
                row['inchi'], row['isosmiles'], row['canonicalsmiles'], row['inchikey'],
                row['iupacname'], row['exactmass'], row['monoisotopicmass'], row['charge'],
                row['covalentunitcnt'], row['isotopeatomcnt'], row['totalatomstereocnt'],
                row['pclidcnt'], row['gpidcnt'], row['gpfamilycnt'], row['neighbortype'],
                row['meshheadings'], row['annothits'], row['annotation']
            ))
    
    conn.commit()
    print(f"Data from {csv_file_path} has been imported to the database.")

def add_plant_compounds(conn):
    cursor = conn.cursor()
    plant_compounds = [
        (1, "Salicin", "OCC1OC(OC2=CC=CC=C2CO)C(O)C(O)C1O", 286.28, "C13H18O7", 119.61, 442, 0.7, 13, 5, 7, 3, "InChI=1S/C13H18O7/c14-6-9-10(16)11(17)12(18)13(19-9)20-8-5-3-1-2-4-7(5)15/h1-5,9-18H,6,8H2/t9-,10-,11+,12-,13-/m1/s1", "OCC1OC(OC2=CC=CC=C2CO)C(O)C(O)C1O", "OCC1OC(OC2=CC=CC=C2CO)C(O)C(O)C1O", "CKTSBUTUHBMZGZ-AQBUPJNZSA-N", "2-(hydroxymethyl)phenyl β-D-glucopyranoside", 286.2747, 286.2747, 0, 1, 0, 5, "Willow bark", "Analgesic, Anti-inflammatory", "Traditional pain relief"),
        (2, "Morphine", "CN1CCC23C4C1CC5=C2C(=C(C=C5)O)OC3C(C=C4)O", 285.34, "C17H19NO3", 52.93, 299, 0.9, 19, 2, 4, 0, "InChI=1S/C17H19NO3/c1-18-7-6-17-10-3-5-13(20)16(17)21-15-12(19)4-2-9(14(15)17)8-11(10)18/h2-5,10-11,13,16,19-20H,6-8H2,1H3/t10-,11+,13-,16-,17-/m0/s1", "CN1CCC23C4C1CC5=C2C(=C(C=C5)O)OC3C(C=C4)O", "CN1CCC23C4C1CC5=C2C(=C(C=C5)O)OC3C(C=C4)O", "BQJCRHHNABKAKU-KBQPJGBKSA-N", "(4R,4aR,7S,7aR,12bS)-3-methyl-2,4,4a,7,7a,13-hexahydro-1H-4,12-methanobenzofuro[3,2-e]isoquinoline-7,9-diol", 285.1365, 285.1365, 1, 1, 0, 5, "Opium poppy", "Analgesic", "Pain management"),
        (3, "Vincristine", "CCC1(CC2CC(C3=C(CCN(C2)C1)C4=CC=CC=C4N3)(C5=C(C=C6C(=C5)C78CCN9C7C(C(C(C9)(CC(C8(C6)OC)O)COC(=O)C(=CO)NC)O)(C(=O)OC)O)OC)O)O", 824.96, "C46H56N4O10", 189.37, 2077, 2.8, 60, 4, 13, 8, "InChI=1S/C46H56N4O10/c1-6-46(55)16-48-23(7-8-49-36-19-11-12-21-38(36)49)34-20-10-9-18-33(34)41-29-26(52)28(54)44(59-5)45(29)31-14-15-50-35(31)22-39(57)42(60-3)30(13-17-56-2)43(58-4)40(22)32(45)24(18)25(41)27(53)37(46)47-50/h7-12,19-21,23,26,28-31,35,37,39,41-42,44,52-55H,6,13-18,22,47H2,1-5H3/t26-,28+,29-,30+,31-,35+,37-,39-,41+,42+,44+,45+,46-/m1/s1", "CCC1(CC2CC(C3=C(CCN(C2)C1)C4=CC=CC=C4N3)(C5=C(C=C6C(=C5)C78CCN9C7C(C(C(C9)(CC(C8(C6)OC)O)COC(=O)C(=CO)NC)O)(C(=O)OC)O)OC)O)O", "CCC1(CC2CC(C3=C(CCN(C2)C1)C4=CC=CC=C4N3)(C5=C(C=C6C(=C5)C78CCN9C7C(C(C(C9)(CC(C8(C6)OC)O)COC(=O)C(=CO)NC)O)(C(=O)OC)O)OC)O)O", "AAOVKJBEBIDNHE-KBQPJGBKSA-N", "methyl (1S,2R,3R,4S,5S,7S,9S)-9-{[(2S)-4-ethyl-2-hydroxy-4-oxopiperidine-1-yl]carbonyl}-5-hydroxy-2-(hydroxymethyl)-1,2,3,4,5,7-hexahydro-3-(methoxycarbonyl)-2,6-methano-8H-azecino[5,4-b]indole-8-carboxylate", 824.3965, 824.3965, 2, 1, 0, 18, "Madagascar periwinkle", "Antineoplastic", "Cancer treatment"),
        (4, "Paclitaxel", "CC1=C2C(C(=O)C3(C(CC4C(C3C(C(C2(C)C)(CC1OC(=O)C(C(C5=CC=CC=C5)NC(=O)C6=CC=CC=C6)O)O)OC(=O)C7=CC=CC=C7)(CO4)OC(=O)C)O)C)OC(=O)C", 853.91, "C47H51NO14", 221.29, 2077, 3.2, 61, 4, 14, 14, "InChI=1S/C47H51NO14/c1-25-24-47(61-40(55)31-16-10-7-11-17-31)34(50)38-42(57)36(28(3)46(25,60-39(54)30-15-9-6-8-14-30)23-48-37(53)29-13-5-4-12-27(29)2)41(56)35(51)33(49)26(38)19-32(47)44(58-21(25)43(47,22-45(59)62)20-52)18-59/h4-17,26,32-36,38,41-42,44,49-52H,18-20,22-24H2,1-3H3,(H,48,53)/t26-,32+,33-,34-,35+,36-,38+,41-,42+,44-,47-/m0/s1", "CC1=C2C(C(=O)C3(C(CC4C(C3C(C(C2(C)C)(CC1OC(=O)C(C(C5=CC=CC=C5)NC(=O)C6=CC=CC=C6)O)O)OC(=O)C7=CC=CC=C7)(CO4)OC(=O)C)O)C)OC(=O)C", "CC1=C2C(C(=O)C3(C(CC4C(C3C(C(C2(C)C)(CC1OC(=O)C(C(C5=CC=CC=C5)NC(=O)C6=CC=CC=C6)O)O)OC(=O)C7=CC=CC=C7)(CO4)OC(=O)C)O)C)OC(=O)C", "RCINICONZNJXQF-MZXODVADSA-N", "[(2R,3S)-3-benzamido-2-hydroxy-3-phenylpropanoyl]oxy-1,7β,10β-trihydroxy-9-oxo-5β,20-epoxytax-11-en-2α-yl benzoate", 853.3322, 853.3322, 0, 1, 0, 15, "Pacific yew tree", "Antineoplastic", "Cancer treatment"),
        (5, "Artemisinin", "CC1CCC2C(C(C3C(C12C)C(=O)OC3(C)OO)C)OC(=O)C4CCC(=CC4)C", 282.33, "C15H22O5", 53.99, 514, 2.8, 20, 0, 5, 0, "InChI=1S/C15H22O5/c1-8-5-6-12-10(3)15(19-7-4-9(2)16)13(17)11(8)14(12,18-20-15)9-7-4-9/h9,11-13,17H,4-7H2,1-3H3/t9-,11-,12+,13+,15-/m0/s1", "CC1CCC2C(C(C3C(C12C)C(=O)OC3(C)OO)C)OC(=O)C4CCC(=CC4)C", "CC1CCC2C(C(C3C(C12C)C(=O)OC3(C)OO)C)OC(=O)C4CCC(=CC4)C", "XJGCQEAXAAWNLV-UHFFFAOYSA-N", "(3R,5aS,6R,8aS,9R,12S,12aR)-octahydro-3,6,9-trimethyl-3,12-epoxy-12H-pyrano[4,3-j]-1,2-benzodioxepin-10(3H)-one", 282.1467, 282.1467, 0, 1, 0, 7, "Sweet wormwood", "Antimalarial", "Malaria treatment"),
        (6, "Capsaicin", "COC1=C(C=CC(=C1)CNC(=O)CCCC/C=C/C(C)C)O", 305.41, "C18H27NO3", 58.56, 450, 3.8, 21, 2, 3, 9, "InChI=1S/C18H27NO3/c1-13(2)5-4-6-7-8-16(20)19-12-14-9-10-15(21)17(11-14)22-18(3)/h4-5,9-11,13,21H,6-8,12H2,1-3H3,(H,19,20)/b5-4+", "COC1=C(C=CC(=C1)CNC(=O)CCCC/C=C/C(C)C)O", "COC1=C(C=CC(=C1)CNC(=O)CCCC/C=C/C(C)C)O", "YKPUWZUDDOIDPM-JTQLQIEISA-N", "(E)-N-[(4-hydroxy-3-methoxyphenyl)methyl]-8-methylnon-6-enamide", 305.1991, 305.1991, 0, 1, 0, 4, "Chili peppers", "Analgesic", "Pain relief"),
        (7, "Ephedrine", "CNC(C)C(C1=CC=CC=C1)O", 165.23, "C10H15NO", 12.03, 145, 1.1, 11, 2, 2, 2, "InChI=1S/C10H15NO/c1-8(11-2)10(12)9-6-4-3-5-7-9/h3-8,10-12H,1-2H3/t8-,10-/m0/s1", "CNC(C)C(C1=CC=CC=C1)O", "CNC(C)C(C1=CC=CC=C1)O", "KWGRBVOPPLSCSI-WAAGHKOSSA-N", "(1R,2S)-2-(methylamino)-1-phenylpropan-1-ol", 165.1154, 165.1154, 1, 1, 0, 2, "Ephedra sinica", "Stimulant", "Asthma treatment"),
        (8, "Galantamine", "COC1=C(C=C2C(=C1)C34CCN(CC3C=CC(C4O2)O)C)OC", 287.35, "C17H21NO3", 41.93, 340, 1.8, 20, 1, 4, 1, "InChI=1S/C17H21NO3/c1-18-8-7-17-6-5-14(19)16(20-17)9-11-3-4-12(21-2)13(15(11)17)10-18/h3-6,9,14,16,19H,7-8,10H2,1-2H3/t14-,16+,17+/m0/s1", "COC1=C(C=C2C(=C1)C34CCN(CC3C=CC(C4O2)O)C)OC", "COC1=C(C=C2C(=C1)C34CCN(CC3C=CC(C4O2)O)C)OC", "BQCBDWQJXSXMJZ-UHFFFAOYSA-N", "(4aS,6R,8aS)-3-methoxy-11-methyl-5,6,9,10,11,12-hexahydro-4aH-[1]benzofuro[3a,3,2-ef][2]benzazepin-6-ol", 287.1521, 287.1521, 1, 1, 0, 4, "Snowdrops and daffodils", "Acetylcholinesterase inhibitor", "Alzheimer's treatment"),
        (9, "Pilocarpine", "CCC1C(COC(=O)N2CCCC2C=C1)OC", 208.26, "C11H16N2O2", 39.73, 253, 1.1, 14, 0, 3, 1, "InChI=1S/C11H16N2O2/c1-3-9-8(6-14-11(13)12-4-2-5-10(12)7-9)15-9/h7,10H,2-6H2,1H3/t9-,10+/m1/s1", "CCC1C(COC(=O)N2CCCC2C=C1)OC", "CCC1C(COC(=O)N2CCCC2C=C1)OC", "SHYZJGJDUPGXIY-VHSXEESVSA-N", "(3S,4R)-3-ethyl-4-[(1-methyl-1H-imidazol-5-yl)methyl]dihydrofuran-2(3H)-one", 208.1212, 208.1212, 0, 1, 0, 2, "Jaborandi plant", "Parasympathomimetic", "Glaucoma treatment"),
        (10, "Digitoxin", "CC1C(C(CC(O1)OC2C(C(C(OC2OC3C4CCC5(C(C4CCC5(C3)C)O)C)C)O)O)O)OC6C(C(C(C(O6)CO)O)O)O", 764.94, "C41H64O13", 203.06, 1311, 2.4, 54, 8, 13, 7, "InChI=1S/C41H64O13/c1-19-9-8-16-39(4)22(48)12-13-40(39,5)23(49)14-17-41(19,40)54-38-37(53-36-35(52-33-30(45)27(42)25(20(2)50)51-33)32(47)29(44)24(15-43)34-31(46)28(36)26(41-34)21(3)18-19)11-10-18-6-7-50/h16,20-30,32,34-38,42-49H,6-15,18H2,1-5H3/t20-,21+,22+,23-,24+,25+,26-,27-,28+,29-,30+,32+,34+,35-,36+,37+,38+,39-,40+,41-/m1/s1", "CC1C(C(CC(O1)OC2C(C(C(OC2OC3C4CCC5(C(C4CCC5(C3)C)O)C)C)O)O)O)OC6C(C(C(C(O6)CO)O)O)O", "CC1C(C(CC(O1)OC2C(C(C(OC2OC3C4CCC5(C(C4CCC5(C3)C)O)C)C)O)O)O)OC6C(C(C(C(O6)CO)O)O)O", "AMHBKBMVCFXYQA-UHFFFAOYSA-N", "(3β,5β)-3-[(O-2,6-dideoxy-β-D-ribo-hexopyranosyl-(1→4)-O-2,6-dideoxy-β-D-ribo-hexopyranosyl-(1→4)-2,6-dideoxy-β-D-ribo-hexopyranosyl)oxy]-14-hydroxycard-20(22)-enolide", 764.4347, 764.4347, 0, 1, 0, 14, "Foxglove", "Cardiac glycoside", "Heart failure treatment")
    ]
    
    cursor.executemany('''
    INSERT OR REPLACE INTO compounds (
        cid, cmpdname, canonicalsmiles, molecular_weight, molecular_formula, polararea,
        complexity, xlogp, heavycnt, hbonddonor, hbondacc, rotbonds,
        inchi, isosmiles, canonicalsmiles, inchikey, iupacname,
        exactmass, monoisotopicmass, charge, covalentunitcnt,
        isotopeatomcnt, totalatomstereocnt, plant_source, biological_activity, traditional_use
    ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
    ''', plant_compounds)
    
    conn.commit()
    print(f"Added {len(plant_compounds)} plant-based compounds to the database.")

def enrich_compound_data(conn):
    cursor = conn.cursor()
    cursor.execute("SELECT cid, cmpdname FROM compounds")
    compounds = cursor.fetchall()

    for cid, name in compounds:
        print(f"Enriching data for {name} (CID: {cid})")
        
        # PubChem API
        try:
            compound = pcp.Compound.from_cid(cid)
            
            # Update fields
            cursor.execute("""
                UPDATE compounds SET
                cmpdsynonym = ?,
                biological_activity = ?,
                traditional_use = ?,
                molecular_weight = ?,
                molecular_formula = ?,
                iupacname = ?
                WHERE cid = ?
            """, (
                ';'.join(compound.synonyms) if compound.synonyms else None,
                compound.pharmacology_description,
                compound.traditional_uses,
                compound.molecular_weight,
                compound.molecular_formula,
                compound.iupac_name,
                cid
            ))
        except pcp.PubChemHTTPError:
            print(f"Could not fetch PubChem data for {name}")

        # Additional web scraping can be added here
        # For example, scraping Wikipedia for traditional uses

    conn.commit()
    print("Data enrichment complete")

# Usage
if __name__ == "__main__":
    conn = create_connection("chempath_database.db")
    enrich_compound_data(conn)
    conn.close()    

if __name__ == "__main__":
    db_file = r"C:\Users\Dr. Contessa Petrini\ChemPath\chempath_database.db"
    csv_file = r"C:\Users\Dr. Contessa Petrini\ChemPath\PubChem.csv"
    conn = create_connection(db_file)
    if conn is not None:
        create_tables(conn)
        alter_table(conn)  # Add this line
        import_pubchem_data(conn, csv_file)
        add_plant_compounds(conn)
        conn.close()
    else:
        print("Error! Cannot create the database connection.")




