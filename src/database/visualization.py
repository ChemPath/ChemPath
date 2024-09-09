import sys
import os
sys.path.append(os.path.dirname(os.path.abspath(__file__)))

import networkx as nx
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA
from sklearn.manifold import TSNE
import plotly.graph_objects as go
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit import DataStructs
from database.chempath_database import PlantCompound, SyntheticCompound, PrescriptionCompound, db
from flask import current_app


def create_similarity_network(threshold=0.7):
    G = nx.Graph()
    compounds = PlantCompound.query.all() + SyntheticCompound.query.all() + PrescriptionCompound.query.all()
    
    for i, comp1 in enumerate(compounds):
        mol1 = Chem.MolFromSmiles(comp1.smiles)
        fp1 = AllChem.GetMorganFingerprintAsBitVect(mol1, 2, nBits=1024)
        G.add_node(comp1.name)
        
        for comp2 in compounds[i+1:]:
            mol2 = Chem.MolFromSmiles(comp2.smiles)
            fp2 = AllChem.GetMorganFingerprintAsBitVect(mol2, 2, nBits=1024)
            similarity = DataStructs.TanimotoSimilarity(fp1, fp2)
            
            if similarity >= threshold:
                G.add_edge(comp1.name, comp2.name, weight=similarity)
    
    return G

def visualize_network(G):
    pos = nx.spring_layout(G)
    plt.figure(figsize=(12, 8))
    nx.draw(G, pos, with_labels=True, node_color='lightblue', node_size=500, font_size=8)
    plt.title("Compound Similarity Network")
    plt.axis('off')
    plt.tight_layout()
    plt.savefig('compound_network.png', dpi=300, bbox_inches='tight')
    plt.close()

def perform_dimensionality_reduction():
    compounds = PlantCompound.query.all() + SyntheticCompound.query.all() + PrescriptionCompound.query.all()
    features = [[c.average_mass, c.logp, c.hydrogen_bond_donors, c.hydrogen_bond_acceptors] for c in compounds]
    
    pca = PCA(n_components=2)
    pca_result = pca.fit_transform(features)
    
    tsne = TSNE(n_components=2, random_state=42)
    tsne_result = tsne.fit_transform(features)
    
    return compounds, pca_result, tsne_result

def visualize_clusters(compounds, pca_result, tsne_result):
    fig = go.Figure()
    
    # PCA plot
    fig.add_trace(go.Scatter(
        x=pca_result[:, 0], y=pca_result[:, 1],
        mode='markers',
        marker=dict(size=8, color=[c.__class__.__name__ for c in compounds], colorscale='Viridis'),
        text=[c.name for c in compounds],
        name='PCA'
    ))
    
    # t-SNE plot
    fig.add_trace(go.Scatter(
        x=tsne_result[:, 0], y=tsne_result[:, 1],
        mode='markers',
        marker=dict(size=8, color=[c.__class__.__name__ for c in compounds], colorscale='Viridis'),
        text=[c.name for c in compounds],
        name='t-SNE',
        visible=False
    ))
    
    fig.update_layout(
        title='Compound Clustering',
        xaxis_title='Dimension 1',
        yaxis_title='Dimension 2',
        updatemenus=[dict(
            type="buttons",
            direction="right",
            x=0.7,
            y=1.2,
            showactive=True,
            buttons=list([
                dict(label="PCA",
                     method="update",
                     args=[{"visible": [True, False]}]),
                dict(label="t-SNE",
                     method="update",
                     args=[{"visible": [False, True]}]),
            ]),
        )]
    )
    
    fig.write_html('compound_clusters.html')

def generate_visualizations():
    G = create_similarity_network()
    visualize_network(G)
    compounds, pca_result, tsne_result = perform_dimensionality_reduction()
    visualize_clusters(compounds, pca_result, tsne_result)

if __name__ == "__main__":
    generate_visualizations()
