from sqlalchemy.orm import Session
from models import Compound, RetrosynthesisResult, ReactionStep
from database import engine
from rdkit import Chem
from rdkit.Chem import AllChem

# Reaction database
reaction_database = [
    {
        "name": "Aspirin hydrolysis",
        "reactant": "CC(=O)Oc1ccccc1C(=O)O",
        "products": ["CC(=O)O", "Oc1ccccc1C(=O)O"]
    },
    {
        "name": "Ester hydrolysis",
        "reactant": "C(=O)OC",
        "products": ["C(=O)O", "CO"]
    },
    {
        "name": "Carboxylic acid reduction",
        "reactant": "CC(=O)O",
        "products": ["CCO"]
    },
    {
        "name": "Acetyl group removal",
        "reactant": "CC(=O)Oc1ccccc1",
        "products": ["Oc1ccccc1", "CC(=O)O"]
    },
    {
        "name": "Salicylic acid reduction",
        "reactant": "Oc1ccccc1C(=O)O",
        "products": ["Oc1ccccc1CO"]
    },
    {
        "name": "Phenol acetylation",
        "reactant": "Oc1ccccc1C(=O)O",
        "products": ["CC(=O)Oc1ccccc1C(=O)O"]
    },
    {
        "name": "Alcohol oxidation",
        "reactant": "CCO",
        "products": ["CC(=O)O"]
    }
]

def substructure_match(mol, substructure):
    pattern = Chem.MolFromSmiles(substructure)
    return mol.HasSubstructMatch(pattern)

def apply_transformation(mol, reaction):
    if substructure_match(mol, reaction["reactant"]):
        reactant_mol = Chem.MolFromSmiles(reaction["reactant"])
        product_mols = [Chem.MolFromSmiles(p) for p in reaction["products"]]
        
        rms = AllChem.ReplaceSubstructs(mol, reactant_mol, product_mols[0])
        if len(rms) > 0:
            products = [Chem.MolToSmiles(rms[0])] + reaction["products"][1:]
            return [p for p in products if Chem.MolFromSmiles(p) is not None]
    return []

def generate_smiles_variants(smiles):
    mol = Chem.MolFromSmiles(smiles)
    return [
        Chem.MolToSmiles(mol),
        Chem.MolToSmiles(mol, isomericSmiles=True),
        Chem.MolToSmiles(mol, kekuleSmiles=True),
        Chem.MolToSmiles(mol, allHsExplicit=True)
    ]

def perform_retrosynthesis(compound_id, smiles, depth=1):
    with Session(engine) as session:
        compound = session.query(Compound).get(compound_id)
        if not compound:
            return {"error": "Compound not found"}

        retrosynthesis_result = RetrosynthesisResult(compound_id=compound_id)
        session.add(retrosynthesis_result)
        session.flush()

        steps = analyze_compound(session, retrosynthesis_result.id, smiles, depth)
        
        retrosynthesis_result.steps = steps
        session.commit()

        return {"result_id": retrosynthesis_result.id, "steps": [step.to_dict() for step in steps]}

def analyze_compound(session, result_id, smiles, depth):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return []

    steps = []
    for reaction in reaction_database:
        products = apply_transformation(mol, reaction)
        if products:
            step = ReactionStep(
                result_id=result_id,
                reaction_name=reaction['name'],
                reactant=smiles,
                products=",".join(products)
            )
            session.add(step)
            steps.append(step)

            if depth > 1:
                for product in products:
                    for variant in generate_smiles_variants(product):
                        sub_steps = analyze_compound(session, result_id, variant, depth - 1)
                        steps.extend(sub_steps)

    return steps

def get_retrosynthesis_result(result_id):
    with Session(engine) as session:
        result = session.query(RetrosynthesisResult).get(result_id)
        if result:
            return {
                "compound_id": result.compound_id,
                "steps": [step.to_dict() for step in result.steps]
            }
        return {"error": "Result not found"}

# Example usage
if __name__ == "__main__":
    result = perform_retrosynthesis(1, "CC(=O)Oc1ccccc1C(=O)O", depth=2)
    print(result)
