import numpy as np
import pandas as pd
import networkx as nx

def load_data(csv_file, graphml_file):
    """Load gene expression data from CSV and graph structure from GraphML."""
    df = pd.read_csv(csv_file)  # Load CSV
    G = nx.read_graphml(graphml_file)  # Load GraphML

    # Map node IDs (n1, n2, ...) to gene names
    node_name_mapping = {node: G.nodes[node].get("name", node) for node in G.nodes()}
    
    # Relabel graph with gene names
    G = nx.relabel_nodes(G, node_name_mapping, copy=True)




        # Extract the largest connected component
    if nx.is_directed(G):
        largest_cc = max(nx.strongly_connected_components(G), key=len)
    else:
        largest_cc = max(nx.connected_components(G), key=len)
    G_LCC = G.subgraph(largest_cc).copy()

    return df, G, G_LCC

def graph_fourier_transform(G, signal):
    """
    Compute the Graph Fourier Transform (GFT) for a given signal on a graph.
    
    Parameters:
    - G (networkx.Graph): Input graph
    - signal (np.array): Signal values for each node
    
    Returns:
    - F_hat (np.array): Graph Fourier Transform of the input signal
    - U (np.array): Eigenvectors of the Laplacian (Fourier basis)
    - lambda_ (np.array): Eigenvalues of the Laplacian
    """
    # Compute the combinatorial Laplacian
    L = nx.laplacian_matrix(G).toarray()
    
    # Compute eigenvalues and eigenvectors
    lambda_, U = np.linalg.eigh(L)

    # Compute Graph Fourier Transform
    F_hat = U.T @ signal
    
    return F_hat, U, lambda_


def save_results_to_csv(results, filename="gft_results.csv"):
    df = pd.DataFrame.from_dict(results, orient="index")
    df.index.name = "sampleid"
    df.to_csv(filename)


def compute_gft_for_samples(csv_file, graphml_file):
    """Compute GFT for each sample in the CSV file."""
    df, G, GLCC = load_data(csv_file, graphml_file)
    
    # Ensure the graph nodes match the genes in the CSV
    genes_in_graph = list(G.nodes())
    genes_in_csv = df.columns[1:]  # Skip 'sampleid' column
    
    # Get the intersection of genes present in both graph and CSV
    common_genes = list(set(genes_in_graph) & set(genes_in_csv))
    common_genes.sort()  # Ensure order consistency

    if not common_genes:
        raise ValueError("No matching genes found between the graph and CSV file.")

    # Extract expression data only for common genes
    df_filtered = df[['sampleid'] + common_genes]

    # Compute GFT for each sample
    results = {}
    for _, row in df_filtered.iterrows():
        sample_id = row['sampleid']
        signal = row[1:].values.astype(float)  # Convert to numeric array
        F_hat, U, lambda_ = graph_fourier_transform(G.subgraph(common_genes), signal)
        results[sample_id] = F_hat
    
    return results

# Example usage:
expression_data_path = "../../data/expression_data.csv"
cosmic_graph_path = "../g_COSMIC.graphml"
gft_results = compute_gft_for_samples(expression_data_path, cosmic_graph_path)
save_results_to_csv(gft_results, "gft_results.csv")

# Print the transformed signal for the first sample
for sample, F_hat in gft_results.items():
    print(f"Sample {sample}: {F_hat}")
    break  # Print only the first sample

