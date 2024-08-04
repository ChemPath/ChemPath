import heapq
import sqlite3
from collections import defaultdict

import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors
from sqlalchemy.orm import Session, sessionmaker

from models import Base, Compound, AdvancedRetrosynthesisResult
from src.database.chempath_database import engine, Session

Session = sessionmaker(bind=engine)

import logging

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

def get_db_connection():
    return sqlite3.connect('chempath.db')

def get_reaction_data(reaction_smiles):
    with get_db_connection() as conn:
        cursor = conn.cursor()
        cursor.execute("SELECT * FROM reactions WHERE reaction_smiles = ?", (reaction_smiles,))
        return cursor.fetchone()

def get_chemical_data(smiles):
    with get_db_connection() as conn:
        cursor = conn.cursor()
        cursor.execute("SELECT * FROM compounds WHERE smiles = ?", (smiles,))
        return cursor.fetchone()

def assess_reagent_availability(reagents):
    availability_scores = {}
    for reagent in reagents:
        commercial_score = check_commercial_availability(reagent)
        synthetic_score = assess_synthetic_accessibility(reagent)
        cost_score = evaluate_cost_and_scalability(reagent)
        overall_score = commercial_score * 0.4 + synthetic_score * 0.3 + cost_score * 0.3
        availability_scores[reagent] = overall_score
    return availability_scores

def check_commercial_availability(reagent):
    with get_db_connection() as conn:
        cursor = conn.cursor()
        cursor.execute("SELECT availability FROM plant_compounds WHERE smiles = ?", (reagent,))
        result = cursor.fetchone()
    return result[0] / 100 if result else 0.1

def assess_synthetic_accessibility(reagent):
    mol = Chem.MolFromSmiles(reagent)
    sa_score = Descriptors.MolLogP(mol)
    return 1 / (1 + np.exp(-0.5 * (sa_score - 2)))

def evaluate_cost_and_scalability(reagent):
    with get_db_connection() as conn:
        cursor = conn.cursor()
        cursor.execute("SELECT molecular_weight FROM plant_compounds WHERE smiles = ?", (reagent,))
        result = cursor.fetchone()
    if result:
        mol_weight = result[0]
        cost_factor = 1 / (1 + np.exp(-0.01 * (mol_weight - 200)))
        scalability_factor = 1 - cost_factor
        return (cost_factor + scalability_factor) / 2
    return 0.5

class ChemicalTransformationRule:
    def __init__(self, reaction_type, conditions, constraints):
        self.reaction_type = reaction_type
        self.conditions = conditions
        self.constraints = constraints

    def apply(self, molecule):
        # Logic to apply the transformation rule
        pass

class Node:
    def __init__(self, smiles, depth):
        self.smiles = smiles
        self.children = []
        self.depth = depth
        self.feasibility_score = None
        self.availability_score = None

def build_retrosynthesis_tree(start_smiles, max_depth):
    root = Node(start_smiles, 0)
    queue = [root]

    while queue:
        current_node = queue.pop(0)
        if current_node.depth < max_depth:
            mol = Chem.MolFromSmiles(current_node.smiles)
            disconnection_points = identify_disconnection_points(mol)
            for disconnection in disconnection_points:
                fragments = perform_disconnection(mol, disconnection)
                for fragment in fragments:
                    child_node = Node(fragment, current_node.depth + 1)
                    current_node.children.append(child_node)
                    queue.append(child_node)
            
            pathway = get_pathway_to_root(current_node)
            current_node.feasibility_score = assess_synthetic_feasibility(pathway)
            current_node.availability_score = assess_reagent_availability([current_node.smiles])

    return root

def get_pathway_to_root(node):
    pathway = []
    current = node
    while current:
        pathway.append(current.smiles)
        current = current.parent
    return pathway[::-1]

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
    pass

def identify_disconnection_points(mol):
    disconnection_points = []
    
    for bond in mol.GetBonds():
        if bond.GetBondType() == Chem.BondType.SINGLE:
            disconnection_points.append(("bond", bond.GetIdx()))
    
    fg_patterns = {
        "ester": Chem.MolFromSmarts("C(=O)O[C,H]"),
        "amide": Chem.MolFromSmarts("C(=O)N[C,H]"),
        "ether": Chem.MolFromSmarts("[C,H]O[C,H]"),
    }
    for fg_name, fg_pattern in fg_patterns.items():
        matches = mol.GetSubstructMatches(fg_pattern)
        for match in matches:
            disconnection_points.append(("functional_group", (fg_name, match)))
    
    ri = mol.GetRingInfo()
    for ring in ri.AtomRings():
        if len(ring) >= 5:
            disconnection_points.append(("ring", ring))
    
    return disconnection_points

def perform_disconnection(current_mol, disconnection):
    rwmol = Chem.RWMol(current_mol)
    logger.info("Molecule: %s", Chem.MolToSmiles(rwmol))
    
    aromatic_atoms = [atom.GetIdx() for atom in rwmol.GetAromaticAtoms()]
    logger.info("Aromatic atoms: %s", aromatic_atoms)
    
    disconnection_type, data = disconnection
    if disconnection_type == "bond":
        rwmol.RemoveBond(rwmol.GetBondWithIdx(data).GetBeginAtomIdx(),
                         rwmol.GetBondWithIdx(data).GetEndAtomIdx())
    
    logger.info("After disconnection: %s", Chem.MolToSmiles(rwmol))
    
    try:
        Chem.SanitizeMol(rwmol)
        logger.info("After sanitization: %s", Chem.MolToSmiles(rwmol))
        fragments = Chem.GetMolFrags(rwmol, asMols=True)
        return [Chem.MolToSmiles(frag) for frag in fragments]
    except Chem.AtomKekulizeException as e:
        logger.info(f"Error: {e}")
        try:
            fragments = Chem.FragmentOnBRICSBonds(rwmol)
            return [Chem.MolToSmiles(frag) for frag in fragments]
        except:
            logger.info("Alternative disconnection strategy failed")
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
        rwmol.RemoveBond(data[0], data[-1])
        fragments = Chem.GetMolFrags(rwmol, asMols=True)
        return [Chem.MolToSmiles(frag) for frag in fragments]
    
    return []

def assess_synthetic_feasibility(pathway):
    complexity_score = evaluate_reaction_complexity(pathway)
    availability_score = evaluate_reactant_availability(pathway)
    conditions_score = evaluate_reaction_conditions(pathway)
    yield_score = evaluate_yield_and_selectivity(pathway)
    
    weights = [0.3, 0.3, 0.2, 0.2]
    feasibility_score = np.dot([complexity_score, availability_score, conditions_score, yield_score], weights)
    
    return feasibility_score

def evaluate_reaction_complexity(pathway):
    complexity_scores = []
    for reaction in pathway:
        reaction_data = get_reaction_data(reaction)
        num_reactants = len(reaction_data['reactants'])
        num_products = len(reaction_data['products'])
        num_steps = reaction_data['num_steps']
        
        complexity = (num_reactants + num_products) * num_steps / 10
        complexity_scores.append(min(complexity, 1))
    
    return np.mean(complexity_scores)

def evaluate_reactant_availability(pathway):
    availability_scores = []
    for reaction in pathway:
        reaction_data = get_reaction_data(reaction)
        for reactant in reaction_data['reactants']:
            chemical_data = get_chemical_data(reactant)
            availability = chemical_data['availability_score']
            availability_scores.append(availability)
    
    return np.mean(availability_scores)

def evaluate_reaction_conditions(pathway):
    condition_scores = []
    for reaction in pathway:
        reaction_data = get_reaction_data(reaction)
        temperature = reaction_data['temperature']
        pressure = reaction_data['pressure']
        
        temp_score = 1 - (abs(temperature - 25) / 200)
        pressure_score = 1 - (abs(pressure - 1) / 100)
        
        condition_scores.append((temp_score + pressure_score) / 2)
    
    return np.mean(condition_scores)

def evaluate_yield_and_selectivity(pathway):
    yield_scores = []
    for reaction in pathway:
        reaction_data = get_reaction_data(reaction)
        yield_value = reaction_data['yield']
        selectivity = reaction_data['selectivity']
        
        yield_scores.append((yield_value * selectivity) / 100)
    
    return np.mean(yield_scores)

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
            logger.info(f"Error in BFS: {e}")
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

def visualize_retrosynthesis_tree(G, result_id):
    pos = nx.spring_layout(G)
    nx.draw(G, pos, with_labels=True, node_color='lightblue', node_size=500, font_size=8)
    labels = {}
    for node in G.nodes():
        feasibility = G.nodes[node].get('feasibility_score', 'N/A')
        availability = G.nodes[node].get('availability_score', 'N/A')
        labels[node] = f"{node}\nFeasibility: {feasibility}\nAvailability: {availability}"
    
    nx.draw_networkx_labels(G, pos, labels, font_size=6)
    plt.axis('off')
    
    filename = f'retrosynthesis_tree_{result_id}.png'
    plt.savefig(filename)
    plt.close()
    
    return filename




def add_nodes_and_edges(G, node):
    G.add_node(node.smiles, feasibility_score=node.feasibility_score, availability_score=node.availability_score)
    for child in node.children:
        G.add_edge(node.smiles, child.smiles)
        add_nodes_and_edges(G, child)

def insert_sample_compounds():
    with Session() as session:
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
    logger.info("Database tables created successfully.")

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

