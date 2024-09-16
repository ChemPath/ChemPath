from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors
from database_operations import create_engine_and_session, get_compound, update_compound
import numpy as np
from sklearn.decomposition import PCA
from sklearn.cluster import KMeans
import matplotlib.pyplot as plt
import seaborn as sns
from rdkit import DataStructs
from database_operations import Compound

def calculate_molecular_descriptors(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    
    descriptors = {
        'MolWt': Descriptors.ExactMolWt(mol),
        'LogP': Descriptors.MolLogP(mol),
        'NumHDonors': Descriptors.NumHDonors(mol),
        'NumHAcceptors': Descriptors.NumHAcceptors(mol),
        'NumRotatableBonds': Descriptors.NumRotatableBonds(mol),
        'TPSA': Descriptors.TPSA(mol),
        'NumAromaticRings': Descriptors.NumAromaticRings(mol)
    }
    return descriptors

def perform_pca(descriptors_list):
    X = np.array(descriptors_list)
    pca = PCA(n_components=2)
    pca_result = pca.fit_transform(X)
    return pca_result

def perform_clustering(pca_result, n_clusters=5):
    kmeans = KMeans(n_clusters=n_clusters, random_state=42)
    cluster_labels = kmeans.fit_predict(pca_result)
    return cluster_labels

def visualize_clusters(pca_result, cluster_labels, compound_names):
    plt.figure(figsize=(12, 8))
    scatter = plt.scatter(pca_result[:, 0], pca_result[:, 1], c=cluster_labels, cmap='viridis')
    plt.colorbar(scatter)
    
    for i, name in enumerate(compound_names):
        plt.annotate(name, (pca_result[i, 0], pca_result[i, 1]))
    
    plt.title('Compound Clustering based on Molecular Descriptors')
    plt.xlabel('First Principal Component')
    plt.ylabel('Second Principal Component')
    plt.tight_layout()
    plt.savefig('compound_clusters.png')
    plt.close()

def analyze_structural_similarity(compounds):
    fingerprints = [AllChem.GetMorganFingerprintAsBitVect(Chem.MolFromSmiles(c.smiles), 2, nBits=1024) for c in compounds]
    similarity_matrix = []
    for fp1 in fingerprints:
        row = []
        for fp2 in fingerprints:
            similarity = DataStructs.TanimotoSimilarity(fp1, fp2)
            row.append(similarity)
        similarity_matrix.append(row)
    
    plt.figure(figsize=(12, 10))
    sns.heatmap(similarity_matrix, annot=True, cmap='YlGnBu', xticklabels=[c.name for c in compounds], yticklabels=[c.name for c in compounds])
    plt.title('Structural Similarity Heatmap')
    plt.tight_layout()
    plt.savefig('structural_similarity_heatmap.png')
    plt.close()

def main():
    engine, session = create_engine_and_session()
    compounds = session.query(Compound).all()
    
    descriptors_list = []
    compound_names = []
    valid_compounds = []
    
    for compound in compounds:
        descriptors = calculate_molecular_descriptors(compound.smiles)
        if descriptors:
            descriptors_list.append(list(descriptors.values()))
            compound_names.append(compound.name)
            valid_compounds.append(compound)
    
    pca_result = perform_pca(descriptors_list)
    cluster_labels = perform_clustering(pca_result)
    
    visualize_clusters(pca_result, cluster_labels, compound_names)
    analyze_structural_similarity(valid_compounds)
    
    print("Analysis complete. Check 'compound_clusters.png' and 'structural_similarity_heatmap.png' for visualizations.")
    
    session.close()

if __name__ == "__main__":
    main()
