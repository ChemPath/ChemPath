from rdkit import Chem
from rdkit.Chem import AllChem
import networkx as nx
from sqlalchemy.orm import Session
from models import Compound, AdvancedRetrosynthesisResult
from src.database.chempath_database import engine
import matplotlib.pyplot as plt
from src.database.chempath_database import engine, create_tables, Session
from sqlalchemy.orm import sessionmaker
from models import Base, Compound, AdvancedRetrosynthesisResult

SessionLocal = sessionmaker(bind=engine)

class ChemicalTransformationRule:
    def __init__(self, reaction_type, conditions, constraints):
        self.reaction_type = reaction_type
        self.conditions = conditions
        self.constraints = constraints

    def apply(self, molecule):
        # Logic to apply the transformation rule
        pass

def load_external_reactions(database_name):
    # Placeholder for loading reactions from external databases
    # This function would connect to Reaxys, SciFinder, or other databases
    # and return a list of reaction dictionaries
    pass

def identify_disconnection_points(mol):
    disconnection_points = []
    for bond in mol.GetBonds():
        if bond.GetBondType() == Chem.BondType.SINGLE:
            disconnection_points.append(("bond", bond.GetIdx()))
        # Add functional group disconnections
        # Add ring disconnections
    return disconnection_points

def perform_disconnection(mol, disconnection):
    disconnection_type, idx = disconnection
    if disconnection_type == "bond":
        rwmol = Chem.RWMol(mol)
        rwmol.RemoveBond(rwmol.GetBondWithIdx(idx).GetBeginAtomIdx(),
                         rwmol.GetBondWithIdx(idx).GetEndAtomIdx())
        fragments = Chem.GetMolFrags(rwmol, asMols=True)
        return [Chem.MolToSmiles(frag) for frag in fragments]
    # Handle other disconnection types
    return []

def build_retrosynthesis_tree(smiles, depth=3, strategy="dfs"):
    mol = Chem.MolFromSmiles(smiles)
    G = nx.DiGraph()
    G.add_node(smiles)
    
    if strategy == "dfs":
        dfs(smiles, depth, G)
    elif strategy == "bfs":
        bfs(smiles, depth, G)
    # Add other search strategies (Dijkstra's, A*)
    
    return G

def dfs(current_smiles, current_depth, G):
    if current_depth == 0:
        return
    
    current_mol = Chem.MolFromSmiles(current_smiles)
    disconnection_points = identify_disconnection_points(current_mol)
    
    for disconnection in disconnection_points:
        fragments = perform_disconnection(current_mol, disconnection)
        for fragment in fragments:
            G.add_edge(current_smiles, fragment)
            dfs(fragment, current_depth - 1, G)

def bfs(start_smiles, max_depth, G):
    queue = [(start_smiles, 0)]
    while queue:
        current_smiles, current_depth = queue.pop(0)
        if current_depth == max_depth:
            continue
        current_mol = Chem.MolFromSmiles(current_smiles)
        disconnection_points = identify_disconnection_points(current_mol)
        for disconnection in disconnection_points:
            fragments = perform_disconnection(current_mol, disconnection)
            for fragment in fragments:
                G.add_edge(current_smiles, fragment)
                queue.append((fragment, current_depth + 1))

def visualize_retrosynthesis_tree(G, result_id):
    num_nodes = G.number_of_nodes()
    fig_size = (max(8, num_nodes * 0.5), max(6, num_nodes * 0.3))
    fig, ax = plt.subplots(figsize=fig_size)
    
    pos = nx.spring_layout(G, k=0.9, iterations=50)
    
    node_size = max(300, 3000 / num_nodes)
    font_size = max(6, 12 - num_nodes * 0.1)
    
    nx.draw(G, pos, ax=ax, with_labels=True, node_color='lightblue', 
            node_size=node_size, font_size=font_size, font_weight='bold',
            arrows=True, arrowsize=10)
    
    ax.set_title(f"Retrosynthesis Tree (Result ID: {result_id})")
    ax.axis('off')
    
    filename = f"retrosynthesis_tree_{result_id}.png"
    plt.savefig(filename, dpi=300, bbox_inches='tight')
    plt.close(fig)
    return filename


def insert_sample_compounds():
    with SessionLocal() as session:
        sample_compounds = [
            Compound(smiles="CC(=O)Oc1ccccc1C(=O)O", name="Aspirin"),
            Compound(smiles="CC12CCC3C(C1CCC2O)CCC4=CC(=O)CCC34C", name="Testosterone"),
            Compound(smiles="CN1CCC[C@H]1c2cccnc2", name="Nicotine"),
            Compound(smiles="CC(C)(C)NCC(O)c1ccc(O)c(O)c1", name="Salbutamol")
        ]
        session.add_all(sample_compounds)
        session.commit()

def create_tables():
    Base.metadata.create_all(engine)
    print("Database tables created successfully.")

def advanced_retrosynthetic_analysis(compound_id, smiles, depth=3, strategy="dfs"):
    session = Session()
    try:
        compound = session.get(Compound, compound_id)
        if not compound:
            return {"error": "Compound not found"}

        G = build_retrosynthesis_tree(smiles, depth, strategy)
        
        result = AdvancedRetrosynthesisResult(
            compound_id=compound_id,
            depth=depth,
            strategy=strategy,
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
    finally:
        session.close()

def get_advanced_retrosynthesis_result(result_id):
    with Session(engine) as session:
        result = session.query(AdvancedRetrosynthesisResult).get(result_id)
        if result:
            return {
                "compound_id": result.compound_id,
                "depth": result.depth,
                "strategy": result.strategy,
                "num_nodes": result.num_nodes,
                "num_edges": result.num_edges,
                "tree_image": result.tree_image
            }
        return {"error": "Result not found"}

if __name__ == "__main__":
    create_tables()
    insert_sample_compounds()
    
    with Session() as session:
        compound = session.query(Compound).filter_by(name="Aspirin").first()
        if compound:
            result = advanced_retrosynthetic_analysis(compound.id, compound.smiles, depth=3, strategy="bfs")
            print(result)
        else:
            print("Sample compound not found.")


