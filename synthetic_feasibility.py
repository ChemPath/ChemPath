# synthetic_feasibility.py

from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Descriptors

def count_synthesis_steps(molecule, target):
    # Simplified step counting based on functional group differences
    return abs(len(molecule.GetSubstructMatches(Chem.MolFromSmarts('[!C;!H]'))) - 
               len(target.GetSubstructMatches(Chem.MolFromSmarts('[!C;!H]'))))

def assess_intermediate_complexity(molecule):
    return Descriptors.ExactMolWt(molecule) / Descriptors.NumHDonors(molecule)

def estimate_side_reactions(molecule):
    # Estimate based on number of reactive functional groups
    reactive_groups = ['[OH]', '[NH]', '[CH]=O', 'C(=O)[OH]', 'C(=O)O[CH3]']
    return sum(len(molecule.GetSubstructMatches(Chem.MolFromSmarts(smarts))) for smarts in reactive_groups)

def evaluate_stereochemistry(molecule):
    return len(Chem.FindMolChiralCenters(molecule, includeUnassigned=True))

def calculate_synthetic_feasibility_score(smiles, target_smiles):
    molecule = Chem.MolFromSmiles(smiles)
    target = Chem.MolFromSmiles(target_smiles)
    
    steps_score = 1 / (1 + count_synthesis_steps(molecule, target))
    complexity_score = 1 / (1 + assess_intermediate_complexity(molecule))
    side_reactions_score = 1 / (1 + estimate_side_reactions(molecule))
    stereochemistry_score = 1 / (1 + evaluate_stereochemistry(molecule))
    
    total_score = (steps_score + complexity_score + side_reactions_score + stereochemistry_score) / 4
    return total_score * 100  # Scale to 0-100

def analyze_synthetic_feasibility(smiles, target_smiles):
    score = calculate_synthetic_feasibility_score(smiles, target_smiles)
    print(f"Synthetic Feasibility Score: {score:.2f}/100")
    if score > 75:
        print("High synthetic feasibility")
    elif score > 50:
        print("Moderate synthetic feasibility")
    else:
        print("Low synthetic feasibility")
