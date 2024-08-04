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
import heapq

SessionLocal = sessionmaker(bind=engine)

class ChemicalTransformationRule:
    def __init__(self, reaction_type, conditions, constraints):
        self.reaction_type = reaction_type
        self.conditions = conditions
        self.constraints = constraints

    def apply(self, molecule):
        # Logic to apply the transformation rule
        pass

def build_retrosynthesis_tree(smiles, depth=3, strategy="dfs"):
    mol = Chem.MolFromSmiles(smiles)
    G = nx.DiGraph()
    G.add_node(smiles)
    
    if strategy == "dfs":
        dfs(smiles, depth, G)
    elif strategy == "bfs":
        bfs(smiles, depth, G)
    elif strategy == "dijkstra":
        dijkstra(smiles, depth, G)
    # Add other search strategies (A*)
    
    return G


def dijkstra(start_smiles, max_depth, G):
    distances = {start_smiles: 0}
    pq = [(0, start_smiles)]
    while pq:
        current_distance, current_smiles = heapq.heappop(pq)
        
        if current_distance > max_depth:
            continue
        
        current_mol = Chem.MolFromSmiles(current_smiles)
        disconnection_points = identify_disconnection_points(current_mol)
        
        for disconnection in disconnection_points:
            fragments = perform_disconnection(current_mol, disconnection)
            for fragment in fragments:
                new_distance = current_distance + 1
                if fragment not in distances or new_distance < distances[fragment]:
                    distances[fragment] = new_distance
                    G.add_edge(current_smiles, fragment)
                    heapq.heappush(pq, (new_distance, fragment))
    
    return G
def heuristic(smiles):
    mol = Chem.MolFromSmiles(smiles)
    return mol.GetNumAtoms() + mol.GetNumBonds()

def astar(start_smiles, max_depth, G):
    open_set = [(0, 0, start_smiles)]
    g_score = {start_smiles: 0}
    f_score = {start_smiles: heuristic(start_smiles)}
    
    while open_set:
        current_f, current_g, current_smiles = heapq.heappop(open_set)
        
        if current_g > max_depth:
            continue
        
        current_mol = Chem.MolFromSmiles(current_smiles)
        disconnection_points = identify_disconnection_points(current_mol)
        
        for disconnection in disconnection_points:
            fragments = perform_disconnection(current_mol, disconnection)
            for fragment in fragments:
                tentative_g = current_g + 1
                
                if fragment not in g_score or tentative_g < g_score[fragment]:
                    g_score[fragment] = tentative_g
                    f_score[fragment] = tentative_g + heuristic(fragment)
                    G.add_edge(current_smiles, fragment)
                    heapq.heappush(open_set, (f_score[fragment], tentative_g, fragment))
    
    return G

def load_external_reactions(database_name):
    # Placeholder for loading reactions from external databases
    # This function would connect to Reaxys, SciFinder, or other databases
    # and return a list of reaction dictionaries
    pass

def identify_disconnection_points(mol):
    disconnection_points = []
    
    # Bond disconnections
    for bond in mol.GetBonds():
        if bond.GetBondType() == Chem.BondType.SINGLE:
            disconnection_points.append(("bond", bond.GetIdx()))
    
    # Functional group disconnections
    fg_patterns = {
        "ester": Chem.MolFromSmarts("C(=O)O[C,H]"),
        "amide": Chem.MolFromSmarts("C(=O)N[C,H]"),
        "ether": Chem.MolFromSmarts("[C,H]O[C,H]"),
    }
    for fg_name, fg_pattern in fg_patterns.items():
        matches = mol.GetSubstructMatches(fg_pattern)
        for match in matches:
            disconnection_points.append(("functional_group", (fg_name, match)))
    
    # Ring disconnections
    ri = mol.GetRingInfo()
    for ring in ri.AtomRings():
        if len(ring) >= 5:  # Only consider rings with 5 or more atoms
            disconnection_points.append(("ring", ring))
    
    return disconnection_points

def perform_disconnection(current_mol, disconnection):
    rwmol = Chem.RWMol(current_mol)
    print("Molecule:", Chem.MolToSmiles(rwmol))
    
    aromatic_atoms = [atom.GetIdx() for atom in rwmol.GetAromaticAtoms()]
    print("Aromatic atoms:", aromatic_atoms)
    
    disconnection_type, data = disconnection
    if disconnection_type == "bond":
        rwmol.RemoveBond(rwmol.GetBondWithIdx(data).GetBeginAtomIdx(),
                         rwmol.GetBondWithIdx(data).GetEndAtomIdx())
    
    print("After disconnection:", Chem.MolToSmiles(rwmol))
    
    try:
        Chem.SanitizeMol(rwmol)
        print("After sanitization:", Chem.MolToSmiles(rwmol))
        fragments = Chem.GetMolFrags(rwmol, asMols=True)
        return [Chem.MolToSmiles(frag) for frag in fragments]
    except Chem.AtomKekulizeException as e:
        print(f"Error: {e}")
        # Try alternative disconnection strategies
        try:
            fragments = Chem.FragmentOnBRICSBonds(rwmol)
            return [Chem.MolToSmiles(frag) for frag in fragments]
        except:
            print("Alternative disconnection strategy failed")
            return []
    
    if disconnection_type == "functional_group":
        fg_name, match = data
        rwmol = Chem.RWMol(current_mol)
        if fg_name == "ester":
            rwmol.RemoveBond(match[1], match[3])
        elif fg_name == "amide":
            rwmol.RemoveBond(match[1], match[2])
        elif fg_name == "ether":
            rwmol.RemoveBond(match[0], match[1])
        fragments = Chem.GetMolFrags(rwmol, asMols=True)
        return [Chem.MolToSmiles(frag) for frag in fragments]
    
    elif disconnection_type == "ring":
        rwmol = Chem.RWMol(current_mol)
        rwmol.RemoveBond(data[0], data[-1])  # Break the ring at an arbitrary point
        fragments = Chem.GetMolFrags(rwmol, asMols=True)
        return [Chem.MolToSmiles(frag) for frag in fragments]
    
    return []    # Handle other disconnection types...
    return []    # Handle other disconnection types...
      

def bfs(start_smiles, max_depth, G):
    queue = [(start_smiles, 0)]
    while queue:
        current_smiles, current_depth = queue.pop(0)
        if current_depth == max_depth:
            continue
        try:
            current_mol = Chem.MolFromSmiles(current_smiles, sanitize=False)
            Chem.SanitizeMol(current_mol, sanitizeOps=Chem.SanitizeFlags.SANITIZE_ALL^Chem.SanitizeFlags.SANITIZE_KEKULIZE)
            disconnection_points = identify_disconnection_points(current_mol)
            for disconnection in disconnection_points:
                fragments = perform_disconnection(current_mol, disconnection)
                for fragment in fragments:
                    G.add_edge(current_smiles, fragment)
                    queue.append((fragment, current_depth + 1))
        except Chem.AtomKekulizeException as e:
            print(f"Error in BFS: {e}")
            continue




def build_retrosynthesis_tree(smiles, depth=3, strategy="dfs"):
    mol = Chem.MolFromSmiles(smiles)
    G = nx.DiGraph()
    G.add_node(smiles)
    
    if strategy == "dfs":
        dfs(smiles, depth, G)
    elif strategy == "bfs":
        bfs(smiles, depth, G)
    elif strategy == "dijkstra":
        dijkstra(smiles, depth, G)
    elif strategy == "astar":
        astar(smiles, depth, G)
    
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


