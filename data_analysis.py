from src.database.chempath_database import create_connection, get_all_compounds
from pathlib import Path
import pandas as pd
import numpy as np
from sklearn.cluster import KMeans
from sklearn.preprocessing import StandardScaler
import matplotlib.pyplot as plt
import networkx as nx
from src.database.chempath_database import create_connection, get_all_compounds, get_therapeutic_areas

def load_compound_data():
    database = Path("chempath_database.db")
    conn = create_connection(database)
    compounds = get_all_compounds(conn)
    conn.close()
    
    df = pd.DataFrame(compounds, columns=['id', 'name', 'smiles', 'molecular_weight', 'logp', 'plant_source', 'biological_activity', 'traditional_use', 'h_bond_donors', 'h_bond_acceptors', 'polar_surface_area', 'rotatable_bonds'])
    return df


def cluster_compounds(df, n_clusters=3):
    features = df[['molecular_weight', 'logp']].dropna()
    scaler = StandardScaler()
    features_normalized = scaler.fit_transform(features)
    kmeans = KMeans(n_clusters=n_clusters, random_state=42)
    df.loc[features.index, 'cluster'] = kmeans.fit_predict(features_normalized)
    return df

def visualize_clusters(df):
    plt.figure(figsize=(10, 6))
    scatter = plt.scatter(df['molecular_weight'], df['logp'], c=df['cluster'], cmap='viridis')
    plt.xlabel('Molecular Weight')
    plt.ylabel('LogP')
    plt.title('Compound Clusters based on Molecular Properties')
    plt.colorbar(scatter, label='Cluster')
    plt.savefig('compound_clusters.png')
    plt.close()

def create_compound_therapeutic_area_network(df):
    G = nx.Graph()
    for _, compound in df.iterrows():
        G.add_node(compound['name'], type='compound')
        areas = compound['biological_activity'].split(', ')
        for area in areas:
            G.add_node(area, type='therapeutic_area')
            G.add_edge(compound['name'], area)
    return G

def visualize_network(G):
    plt.figure(figsize=(12, 8))
    pos = nx.spring_layout(G)
    compound_nodes = [node for node, data in G.nodes(data=True) if data['type'] == 'compound']
    area_nodes = [node for node, data in G.nodes(data=True) if data['type'] == 'therapeutic_area']
    
    nx.draw_networkx_nodes(G, pos, nodelist=compound_nodes, node_color='lightblue', node_size=500, alpha=0.8)
    nx.draw_networkx_nodes(G, pos, nodelist=area_nodes, node_color='lightgreen', node_size=700, alpha=0.8)
    nx.draw_networkx_edges(G, pos, width=1, alpha=0.5)
    nx.draw_networkx_labels(G, pos, font_size=8)
    
    plt.title('Compound-Therapeutic Area Network')
    plt.axis('off')
    plt.tight_layout()
    plt.savefig('compound_therapeutic_area_network.png')
    plt.close()
import matplotlib.pyplot as plt

def plot_molecular_weight_distribution(df):
    plt.figure(figsize=(10, 6))
    plt.hist(df['molecular_weight'], bins=20, edgecolor='black')
    plt.title('Distribution of Molecular Weights')
    plt.xlabel('Molecular Weight')
    plt.ylabel('Frequency')
    plt.savefig('molecular_weight_distribution.png')
    plt.close()

def analyze_compounds():
    df = load_compound_data()
    
    print("Basic statistics of molecular properties:")
    print(df[['molecular_weight', 'logp']].describe())
    
    G = create_compound_therapeutic_area_network(df)
    print(f"Network created with {G.number_of_nodes()} nodes and {G.number_of_edges()} edges")
    
    plot_molecular_weight_distribution(df)


if __name__ == "__main__":
    analyze_compounds()
