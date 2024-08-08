from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors
from src.database.plant_metabolite_database import create_engine_and_session, get_metabolites_by_pathway, PlantMetabolite
from optimization import optimize_compound


class PlantRetrosynthesis:
    def __init__(self):
        self.engine, self.session = create_engine_and_session()
        self.biosynthetic_rules = self.load_biosynthetic_rules()
        self.known_metabolites = self.load_known_metabolites()

    def load_biosynthetic_rules(self):
        return {
            "terpene_synthesis": self.terpene_synthesis,
            "shikimate_pathway": self.shikimate_pathway,
            "polyketide_synthesis": self.polyketide_synthesis,
            "alkaloid_synthesis": self.alkaloid_synthesis,
            "flavonoid_synthesis": self.flavonoid_synthesis
        }

    def load_known_metabolites(self):
        return self.session.query(PlantMetabolite).all()

    def analyze(self, smiles):
        mol = Chem.MolFromSmiles(smiles)
        possible_pathways = []
        for rule_name, rule_func in self.biosynthetic_rules.items():
            if rule_func(mol):
                possible_pathways.append(rule_name)
        return possible_pathways

    def terpene_synthesis(self, mol):
        terpene_metabolites = get_metabolites_by_pathway(self.session, "Terpene biosynthesis")
        return any(Chem.MolFromSmiles(m.smiles).HasSubstructMatch(mol) for m in terpene_metabolites)

    def shikimate_pathway(self, mol):
        shikimate_metabolites = get_metabolites_by_pathway(self.session, "Shikimate pathway")
        return any(Chem.MolFromSmiles(m.smiles).HasSubstructMatch(mol) for m in shikimate_metabolites)

    def polyketide_synthesis(self, mol):
        polyketide_metabolites = get_metabolites_by_pathway(self.session, "Polyketide biosynthesis")
        return any(Chem.MolFromSmiles(m.smiles).HasSubstructMatch(mol) for m in polyketide_metabolites)

    def alkaloid_synthesis(self, mol):
        alkaloid_metabolites = get_metabolites_by_pathway(self.session, "Alkaloid biosynthesis")
        return any(Chem.MolFromSmiles(m.smiles).HasSubstructMatch(mol) for m in alkaloid_metabolites)

    def flavonoid_synthesis(self, mol):
        flavonoid_metabolites = get_metabolites_by_pathway(self.session, "Flavonoid biosynthesis")
        return any(Chem.MolFromSmiles(m.smiles).HasSubstructMatch(mol) for m in flavonoid_metabolites)

    def suggest_precursors(self, smiles):
        pathways = self.analyze(smiles)
        precursors = []
        for pathway in pathways:
            precursors.extend(get_metabolites_by_pathway(self.session, pathway))
        return precursors

    def optimize_and_analyze(self, smiles, target_property, iterations=100):
        optimized_mol = optimize_compound(self.session, smiles, target_property, iterations)
        optimized_smiles = Chem.MolToSmiles(optimized_mol)
        pathways = self.analyze(optimized_smiles)
        precursors = self.suggest_precursors(optimized_smiles)
        
        return {
            'optimized_smiles': optimized_smiles,
            'pathways': pathways,
            'precursors': precursors
        }

    def score_compound(self, smiles, target_property):
        mol = Chem.MolFromSmiles(smiles)
        property_score = self.calculate_property_score(mol, target_property)
        synthetic_feasibility = len(self.analyze(smiles))
        similarity_score = self.calculate_similarity_to_known_compounds(mol)
        
        total_score = property_score * 0.4 + synthetic_feasibility * 0.3 + similarity_score * 0.3
        return total_score

    def calculate_property_score(self, mol, target_property):
        # Implement property calculation based on the target_property
        if target_property == 'logP':
            return Descriptors.MolLogP(mol)
        # Add more property calculations as needed
        return 0

    def calculate_similarity_to_known_compounds(self, mol):
        max_similarity = 0
        for known_metabolite in self.known_metabolites:
            known_mol = Chem.MolFromSmiles(known_metabolite.smiles)
            similarity = AllChem.GetMorganFingerprintAsBitVect(mol, 2, nBits=1024).TanimotoSimilarity(
                AllChem.GetMorganFingerprintAsBitVect(known_mol, 2, nBits=1024)
            )
            max_similarity = max(max_similarity, similarity)
        return max_similarity

def discover_novel_compounds(smiles, target_property):
    retro = PlantRetrosynthesis()
    result = retro.optimize_and_analyze(smiles, target_property)
    score = retro.score_compound(result['optimized_smiles'], target_property)
    
    return {
        'optimized_smiles': result['optimized_smiles'],
        'score': score,
        'pathways': result['pathways'],
        'precursors': [p.name for p in result['precursors']]
    }

if __name__ == "__main__":
    test_compound = "CC(C)C1=CC=C(C=C1)C(C)C(=O)O"  # Ibuprofen as an example
    result = discover_novel_compounds(test_compound, 'logP')
    print(f"Optimized compound: {result['optimized_smiles']}")
    print(f"Score: {result['score']}")
    print(f"Possible biosynthetic pathways: {', '.join(result['pathways'])}")
    print(f"Potential precursors: {', '.join(result['precursors'])}")
