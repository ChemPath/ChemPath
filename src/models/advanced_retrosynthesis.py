import heapq
import sqlite3
from collections import defaultdict
import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors
from sqlalchemy.orm import Session, sessionmaker
from .base import Base
from .compound import Compound
from src.database.chempath_database import engine
from sqlalchemy import Column, Integer, String, Float, JSON
from networkx.readwrite import json_graph
Session = sessionmaker(bind=engine)

import logging

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

def get_db_connection():
    """Returns a connection to the SQLite database."""
    return sqlite3.connect('chempath.db')

def get_reaction_data(reaction_smiles):
    """Retrieves reaction data from the database."""
    with get_db_connection() as conn:
        cursor = conn.cursor()
        cursor.execute("SELECT * FROM reactions WHERE reaction_smiles = ?", (reaction_smiles,))
        return cursor.fetchone()

def get_chemical_data(smiles):
    """Retrieves chemical data from the database."""
    with get_db_connection() as conn:
        cursor = conn.cursor()
        cursor.execute("SELECT * FROM compounds WHERE smiles = ?", (smiles,))
        return cursor.fetchone()

def assess_reagent_availability(reagents):
    """Assesses the availability of reagents."""
    availability_scores = {}
    for reagent in reagents:
        commercial_score = check_commercial_availability(reagent)
        synthetic_score = assess_synthetic_accessibility(reagent)
        cost_score = evaluate_cost_and_scalability(reagent)
        overall_score = commercial_score * 0.4 + synthetic_score * 0.3 + cost_score * 0.3
        availability_scores[reagent] = overall_score
    return availability_scores

def check_commercial_availability(reagent):
    """Checks commercial availability of a reagent."""
    with get_db_connection() as conn:
        cursor = conn.cursor()
        cursor.execute("SELECT availability FROM plant_compounds WHERE smiles = ?", (reagent,))
        result = cursor.fetchone()
    return result[0] / 100 if result else 0.1

def assess_synthetic_accessibility(reagent):
    """Assesses synthetic accessibility of a reagent."""
    mol = Chem.MolFromSmiles(reagent)
    sa_score = Descriptors.MolLogP(mol)
    return 1 / (1 + np.exp(-0.5 * (sa_score - 2)))

def evaluate_cost_and_scalability(reagent):
    """Evaluates cost and scalability of a reagent."""
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
    """Represents a chemical transformation rule."""
    def __init__(self, reaction_type, conditions, constraints):
        self.reaction_type = reaction_type
        self.conditions = conditions
        self.constraints = constraints

    def apply(self, molecule):
        # Logic to apply the transformation rule
        pass

class Node:
    """Represents a node in the retrosynthesis tree."""
    def __init__(self, smiles, depth):
        self.smiles = smiles
        self.children = []
        self.depth = depth
        self.feasibility_score = None
        self.availability_score = None

def build_retrosynthesis_tree(start_smiles, max_depth, strategy="dfs"):
    """Builds a retrosynthesis tree using the specified strategy."""
    mol = Chem.MolFromSmiles(start_smiles)
    G = nx.DiGraph()
    G.add_node(start_smiles)
    
    if strategy == "dfs":
        dfs(start_smiles, max_depth, G)
    elif strategy == "bfs":
        bfs(start_smiles, max_depth, G)
    elif strategy == "dijkstra":
        dijkstra(start_smiles, max_depth, G)
    elif strategy == "astar":
        astar(start_smiles, max_depth, G)
    
    return G

def dfs(current_smiles, current_depth, G):
    """Performs a depth-first search to build the retrosynthesis tree."""
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
    """Performs a breadth-first search to build the retrosynthesis tree."""
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

def dijkstra(start_smiles, max_depth, G):
    """Performs Dijkstra's algorithm to build the retrosynthesis tree."""
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

def astar(start_smiles, max_depth, G):
    r"""Performs A* search to build the retrosynthesis tree."""
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

def heuristic(smiles):
    """Calculates the heuristic score for a molecule."""
    mol = Chem.MolFromSmiles(smiles)
    return mol.GetNumAtoms() + mol.GetNumBonds()

def identify_disconnection_points(mol):
    """Identifies disconnection points in a molecule."""
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
    """Performs disconnection on a molecule."""
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
    """Assesses synthetic feasibility of a pathway."""
    complexity_score = evaluate_reaction_complexity(pathway)
    availability_score = evaluate_reactant_availability(pathway)
    conditions_score = evaluate_reaction_conditions(pathway)
    yield_score = evaluate_yield_and_selectivity(pathway)
    
    weights = [0.3, 0.3, 0.2, 0.2]
    feasibility_score = np.dot([complexity_score, availability_score, conditions_score, yield_score], weights)
    
    return feasibility_score

def evaluate_reaction_complexity(pathway):
    """Evaluates reaction complexity of a pathway."""
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
    """Evaluates reactant availability of a pathway."""
    availability_scores = []
    for reaction in pathway:
        reaction_data = get_reaction_data(reaction)
        for reactant in reaction_data['reactants']:
            chemical_data = get_chemical_data(reactant)
            availability = chemical_data['availability_score']
            availability_scores.append(availability)
    
    return np.mean(availability_scores)

def evaluate_reaction_conditions(pathway):
    """Evaluates reaction conditions of a pathway."""
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
    """Evaluates yield and selectivity of a pathway."""
    yield_scores = []
    for reaction in pathway:
        reaction_data = get_reaction_data(reaction)
        yield_value = reaction_data['yield']
        selectivity = reaction_data['selectivity']
        
        yield_scores.append((yield_value * selectivity) / 100)
    
    return np.mean(yield_scores)

def visualize_retrosynthesis_tree(G, result_id):
    """Visualizes the retrosynthesis tree."""
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
    """Adds nodes and edges to the graph."""
    G.add_node(node.smiles, feasibility_score=node.feasibility_score, availability_score=node.availability_score)
    for child in node.children:
        G.add_edge(node.smiles, child.smiles)
        add_nodes_and_edges(G, child)

def insert_sample_compounds():
    """Inserts sample compounds into the database."""
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
    """Creates tables in the database."""
    Base.metadata.create_all(engine)
    logger.info("Database tables created successfully.")

def advanced_retrosynthetic_analysis(compound_id, smiles, depth=3, strategy="dfs"):
    """Performs advanced retrosynthetic analysis."""
    session = Session()
    try:
        compound = session.get(Compound, compound_id)
        if not compound:
            return {"error": "Compound not found"}

        G = build_retrosynthesis_tree(smiles, depth, strategy)

        result = {
            'graph': G,
            'num_nodes': G.number_of_nodes(),
            'num_edges': G.number_of_edges(),
            'num_pathways': len(list(nx.all_simple_paths(G, smiles, [node for node in G.nodes() if G.out_degree(node) == 0]))),
            'tree_image': visualize_retrosynthesis_tree(G, compound_id)
        }

        graph_data = json_graph.node_link_data(result['graph'])
        result['graph'] = graph_data
    
        return result
    finally:
        session.close()

from networkx.readwrite import json_graph
def get_advanced_retrosynthesis_result(result_id):
    result = AdvancedRetrosynthesisResult.query.get(result_id)
    if result:
        graph = json_graph.node_link_graph(result.graph)
        return {
            'result_id': result.id,
            'compound_id': result.compound_id,
            'depth': result.depth,
            'strategy': result.strategy,
            'num_nodes': result.num_nodes,
            'num_edges': result.num_edges,
            'num_pathways': result.num_pathways,
            'tree_image': result.tree_image,
            'graph': graph
        }
    return None

        


class AdvancedRetrosynthesisResult(Base):
    __tablename__ = 'advanced_retrosynthesis_results'

    id = Column(Integer, primary_key=True)
    compound_id = Column(Integer, nullable=False)
    depth = Column(Integer, nullable=False)
    strategy = Column(String, nullable=False)
    num_nodes = Column(Integer)
    num_edges = Column(Integer)
    num_pathways = Column(Integer)
    tree_image = Column(String)
    graph = Column(JSON)  # New field to store graph data

    def __repr__(self):
        return f"<AdvancedRetrosynthesisResult(id={self.id}, compound_id={self.compound_id})>"

def main():
    create_tables()
    insert_sample_compounds()
    
    compound_id = 1
    smiles = "CC(=O)Oc1ccccc1C(=O)O"
    depth = 3
    strategy = "dfs"
    
    result = advanced_retrosynthetic_analysis(compound_id, smiles, depth, strategy)
    print(result)

if __name__ == "__main__":
    main()   
    