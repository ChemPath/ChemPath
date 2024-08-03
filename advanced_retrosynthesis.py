from rdkit import Chem
from rdkit.Chem import AllChem
import networkx as nx
from sqlalchemy.orm import Session
from models import Compound, AdvancedRetrosynthesisResult
from database import engine
import matplotlib.pyplot as plt

def identify_disconnection_points(mol):
    disconnection_points = []
    for bond in mol.GetBonds():
        if bond.GetBondType() == Chem.BondType.SINGLE:
            disconnection_points.append(bond.GetIdx())
    return disconnection_points

def perform_disconnection(mol, bond_idx):
    rwmol = Chem.RWMol(mol)
    rwmol.RemoveBond(rwmol.GetBondWithIdx(bond_idx).GetBeginAtomIdx(),
                     rwmol.GetBondWithIdx(bond_idx).GetEndAtomIdx())
    fragments = Chem.GetMolFrags(rwmol, asMols=True)
    return [Chem.MolToSmiles(frag) for frag in fragments]

def build_retrosynthesis_tree(smiles, depth=3):
    mol = Chem.MolFromSmiles(smiles)
    G = nx.DiGraph()
    G.add_node(smiles)
    
    def dfs(current_smiles, current_depth):
        if current_depth == 0:
            return
        
        current_mol = Chem.MolFromSmiles(current_smiles)
        disconnection_points = identify_disconnection_points(current_mol)
        
        for bond_idx in disconnection_points:
            fragments = perform_disconnection(current_mol, bond_idx)
            for fragment in fragments:
                G.add_edge(current_smiles, fragment)
                dfs(fragment, current_depth - 1)
    
    dfs(smiles, depth)
    return G

def visualize_retrosynthesis_tree(G, result_id):
    plt.figure(figsize=(12, 8))
    pos = nx.spring_layout(G)
    nx.draw(G, pos, with_labels=True, node_color='lightblue', 
            node_size=2000, font_size=8, font_weight='bold')
    plt.title(f"Retrosynthesis Tree (Result ID: {result_id})")
    plt.axis('off')
    plt.tight_layout()
    
    # Save the plot to a file
    filename = f"retrosynthesis_tree_{result_id}.png"
    plt.savefig(filename)
    plt.close()
    return filename

def advanced_retrosynthetic_analysis(compound_id, smiles, depth=3):
    with Session(engine) as session:
        compound = session.query(Compound).get(compound_id)
        if not compound:
            return {"error": "Compound not found"}

        G = build_retrosynthesis_tree(smiles, depth)
        
        result = AdvancedRetrosynthesisResult(
            compound_id=compound_id,
            depth=depth,
            num_nodes=G.number_of_nodes(),
            num_edges=G.number_of_edges()
        )
        session.add(result)
        session.flush()

        tree_image = visualize_retrosynthesis_tree(G, result.id)
        result.tree_image = tree_image
        
        session.commit()

        return {
            "result_id": result.id,
            "num_nodes": result.num_nodes,
            "num_edges": result.num_edges,
            "tree_image": result.tree_image
        }

def get_advanced_retrosynthesis_result(result_id):
    with Session(engine) as session:
        result = session.query(AdvancedRetrosynthesisResult).get(result_id)
        if result:
            return {
                "compound_id": result.compound_id,
                "depth": result.depth,
                "num_nodes": result.num_nodes,
                "num_edges": result.num_edges,
                "tree_image": result.tree_image
            }
        return {"error": "Result not found"}

# Example usage
if __name__ == "__main__":
    result = advanced_retrosynthetic_analysis(1, "CC(=O)Oc1ccccc1C(=O)O", depth=3)
    print(result)
