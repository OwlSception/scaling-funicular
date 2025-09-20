# Clustering Methods in Python

import numpy as np
import pandas as pd
from sklearn.cluster import KMeans, AgglomerativeClustering
from sklearn.metrics import silhouette_score
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.cluster.hierarchy import dendrogram, linkage
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
from sklearn.manifold import TSNE
import umap
import hdbscan


# Function to perform K-Means Clustering
def kmeans_clustering(data, n_clusters):
    scaler = StandardScaler()
    data_scaled = scaler.fit_transform(data)
    kmeans = KMeans(n_clusters=n_clusters, random_state=42)
    labels = kmeans.fit_predict(data_scaled)
    silhouette_avg = silhouette_score(data_scaled, labels)
    return labels, silhouette_avg


# Function to perform Hierarchical Clustering
def hierarchical_clustering(data, n_clusters):
    scaler = StandardScaler()
    data_scaled = scaler.fit_transform(data)
    hierarchical = AgglomerativeClustering(n_clusters=n_clusters)
    labels = hierarchical.fit_predict(data_scaled)
    silhouette_avg = silhouette_score(data_scaled, labels)
    return labels, silhouette_avg


# Function to plot dendrogram
def plot_dendrogram(data):
    scaler = StandardScaler()
    data_scaled = scaler.fit_transform(data)
    linked = linkage(data_scaled, "ward")
    plt.figure(figsize=(10, 7))
    dendrogram(linked)
    plt.show()


# Function to perform PCA
def pca_reduction(data, n_components=2):
    scaler = StandardScaler()
    data_scaled = scaler.fit_transform(data)
    pca = PCA(n_components=n_components)
    pca_result = pca.fit_transform(data_scaled)
    return pca_result


# Function to perform t-SNE
def tsne_reduction(data, n_components=2, perplexity=30):
    scaler = StandardScaler()
    data_scaled = scaler.fit_transform(data)
    tsne = TSNE(n_components=n_components, perplexity=perplexity, random_state=42)
    tsne_result = tsne.fit_transform(data_scaled)
    return tsne_result


# Function to perform UMAP
def umap_reduction(data, n_components=2, n_neighbors=15, min_dist=0.1):
    scaler = StandardScaler()
    data_scaled = scaler.fit_transform(data)
    umap_model = umap.UMAP(
        n_components=n_components,
        n_neighbors=n_neighbors,
        min_dist=min_dist,
        random_state=42,
    )
    umap_result = umap_model.fit_transform(data_scaled)
    return umap_result


# Function to perform HDBSCAN Clustering
def hdbscan_clustering(data, min_cluster_size=5):
    scaler = StandardScaler()
    data_scaled = scaler.fit_transform(data)
    clusterer = hdbscan.HDBSCAN(min_cluster_size=min_cluster_size)
    labels = clusterer.fit_predict(data_scaled)
    return labels, clusterer.probabilities_


# Function to visualize clusters
def plot_clusters(reduced_data, labels, title="Cluster Plot"):
    plt.figure(figsize=(10, 7))
    sns.scatterplot(
        x=reduced_data[:, 0],
        y=reduced_data[:, 1],
        hue=labels,
        palette="viridis",
        legend="full",
    )
    plt.title(title)
    plt.show()


# Example usage
# if __name__ == "__main__":
#     # Load your data
#     data = pd.read_csv('your_data.csv')  # Replace with your data file

#     # K-Means Clustering
#     kmeans_labels, kmeans_silhouette = kmeans_clustering(data, n_clusters=3)
#     print(f'K-Means Silhouette Score: {kmeans_silhouette}')

#     # Hierarchical Clustering
#     hierarchical_labels, hierarchical_silhouette = hierarchical_clustering(data, n_clusters=3)
#     print(f'Hierarchical Silhouette Score: {hierarchical_silhouette}')
#     plot_dendrogram(data)

#     # Dimensionality Reduction
#     pca_result = pca_reduction(data)
#     tsne_result = tsne_reduction(data)
#     umap_result = umap_reduction(data)

#     # Plotting results
#     plot_clusters(pca_result, kmeans_labels, title='K-Means Clusters (PCA)')
#     plot_clusters(tsne_result, hierarchical_labels, title='Hierarchical Clusters (t-SNE)')
#     plot_clusters(umap_result, kmeans_labels, title='K-Means Clusters (UMAP)')

#     # HDBSCAN Clustering
#     hdbscan_labels, hdbscan_probs = hdbscan_clustering(data)
#     plot_clusters(umap_result, hdbscan_labels, title='HDBSCAN Clusters (UMAP)')
