import re
from rdkit import Chem
from rdkit.Chem import Descriptors
from database_operations import create_engine_and_session, get_compound, update_compound
from database_operations import Compound

def validate_smiles(smiles):
    mol = Chem.MolFromSmiles(smiles)
    return mol is not None

def validate_inchi(inchi):
    mol = Chem.MolFromInchi(inchi)
    return mol is not None

def validate_molecular_weight(mw):
    return isinstance(mw, (int, float)) and mw > 0

def validate_logp(logp):
    return isinstance(logp, (int, float))

def validate_pubchem_cid(cid):
    return isinstance(cid, int) and cid > 0

def validate_name(name):
    return isinstance(name, str) and len(name) > 0

def calculate_molecular_weight(smiles):
    mol = Chem.MolFromSmiles(smiles)
    return Descriptors.ExactMolWt(mol) if mol else None

def calculate_logp(smiles):
    mol = Chem.MolFromSmiles(smiles)
    return Descriptors.MolLogP(mol) if mol else None

def check_compound_data(compound):
    errors = []
    
    if not validate_name(compound.name):
        errors.append("Invalid compound name")
    
    if not validate_smiles(compound.smiles):
        errors.append("Invalid SMILES string")
    else:
        calculated_mw = calculate_molecular_weight(compound.smiles)
        if calculated_mw and abs(calculated_mw - compound.molecular_weight) > 0.1:
            errors.append(f"Molecular weight mismatch: calculated {calculated_mw}, stored {compound.molecular_weight}")
        
        calculated_logp = calculate_logp(compound.smiles)
        if calculated_logp and abs(calculated_logp - compound.logp) > 0.5:
            errors.append(f"LogP mismatch: calculated {calculated_logp}, stored {compound.logp}")
    
    if not validate_molecular_weight(compound.molecular_weight):
        errors.append("Invalid molecular weight")
    
    if not validate_logp(compound.logp):
        errors.append("Invalid LogP value")
    
    if not validate_pubchem_cid(compound.pubchem_cid):
        errors.append("Invalid PubChem CID")
    
    if compound.inchi and not validate_inchi(compound.inchi):
        errors.append("Invalid InChI string")
    
    return errors

def check_database_integrity():
    engine, session = create_engine_and_session()
    compounds = session.query(Compound).all()
    
    for compound in compounds:
        errors = check_compound_data(compound)
        if errors:
            print(f"Errors found for compound {compound.name} (ID: {compound.id}):")
            for error in errors:
                print(f"  - {error}")
        else:
            print(f"No errors found for compound {compound.name} (ID: {compound.id})")
    
    session.close()

if __name__ == "__main__":
    check_database_integrity()

