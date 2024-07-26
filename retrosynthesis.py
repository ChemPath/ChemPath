# retrosynthesis.py

from database_operations import store_retrosynthesis_result
from rdkit import Chem
from rdkit.Chem import AllChem

# Updated reaction database with explicit SMILES-based reactions
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
        
        # Replace the substructure in the original molecule
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

def retrosynthetic_analysis(smiles, depth=1, indent=""):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        print(f"{indent}Invalid SMILES string")
        return

    print(f"{indent}Retrosynthetic analysis for {smiles}:")
   
    any_transformation_applied = False
    for reaction in reaction_database:
        products = apply_transformation(mol, reaction)
        if products:
            any_transformation_applied = True
            print(f"{indent}  {reaction['name']}:")
            for product in products:
                print(f"{indent}    - {product}")
                if depth > 1:
                    for variant in generate_smiles_variants(product):
                        retrosynthetic_analysis(variant, depth - 1, indent + "    ")
   
    if not any_transformation_applied:
        print(f"{indent}  No applicable transformations found")

    print(f"{indent}Retrosynthetic analysis complete")

def perform_retrosynthesis(conn, compound_id, smiles, depth=1):
    retrosynthesis_data = {}
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        retrosynthesis_data['error'] = "Invalid SMILES string"
    else:
        retrosynthesis_data['steps'] = []
        for reaction in reaction_database:
            products = apply_transformation(mol, reaction)
            if products:
                step = {
                    'reaction_name': reaction['name'],
                    'products': products
                }
                retrosynthesis_data['steps'].append(step)
                if depth > 1:
                    for product in products:
                        for variant in generate_smiles_variants(product):
                            sub_analysis = perform_retrosynthesis(conn, compound_id, variant, depth - 1)
                            if 'steps' in sub_analysis:
                                step['sub_steps'] = sub_analysis['steps']

    store_retrosynthesis_result(conn, compound_id, retrosynthesis_data)
    return retrosynthesis_data
