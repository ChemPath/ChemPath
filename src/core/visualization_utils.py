from rdkit import Chem
from rdkit.Chem import Draw
import io
import base64

def visualize_molecule(smiles: str):
    mol = Chem.MolFromSmiles(smiles)
    img = Draw.MolToImage(mol)
    buffered = io.BytesIO()
    img.save(buffered, format="PNG")
    return base64.b64encode(buffered.getvalue()).decode()

def visualize_reaction(reactants: list, products: list):
    reaction = Chem.rdChemReactions.ChemicalReaction()
    for r in reactants:
        reaction.AddReactantTemplate(Chem.MolFromSmiles(r))
    for p in products:
        reaction.AddProductTemplate(Chem.MolFromSmiles(p))
    img = Draw.ReactionToImage(reaction)
    buffered = io.BytesIO()
    img.save(buffered, format="PNG")
    return base64.b64encode(buffered.getvalue()).decode()
