import os
import networkx as nx

def load_graph(graphml_path):
    """Load a graph from a GraphML file."""
    return nx.read_graphml(graphml_path)

def get_largest_connected_component(G):
    """Extract the largest connected component from a graph."""
    if not nx.is_connected(G):
        largest_cc = max(nx.connected_components(G), key=len)
        return G.subgraph(largest_cc).copy()
    return G.copy()

def create_personalization_vector(G, gene_list):
    """Create a personalization vector for personalized PageRank."""
    return {node: (1 if G.nodes[node].get('name', node) in gene_list else 0) for node in G.nodes}

def compute_personalized_pagerank(G, P0, alpha=0.85):
    """Compute personalized PageRank."""
    return nx.pagerank(G, alpha=alpha, personalization=P0)

def process_gene_list(G, gene_list_path):
    """Process a single gene list file and return the top PageRank results."""
    with open(gene_list_path, "r") as f:
        gene_list = {line.strip() for line in f}

    P0 = create_personalization_vector(G, gene_list)
    G_simple = nx.Graph(G)  # Convert to simple graph if needed
    pagerank_scores = compute_personalized_pagerank(G_simple, P0)

    # Sort genes by PageRank score
    return sorted(
        ((G.nodes[node].get('name', node), score) for node, score in pagerank_scores.items()),
        key=lambda x: x[1], reverse=True
    )[:50]  # Get top 50 genes

def save_top_genes(pagerank_scores, file_name):
    # Sort genes by their PageRank scores (highest first)
    sorted_genes = sorted(pagerank_scores.items(), key=lambda x: x[1], reverse=True)

    # Get the top 10, 20, and 50 genes (no scores, just names)
    top_10_genes = [gene for gene, score in sorted_genes[:10]]
    top_20_genes = [gene for gene, score in sorted_genes[:20]]
    top_50_genes = [gene for gene, score in sorted_genes[:50]]

    # Define file paths
    dir_path = os.path.dirname(file_name)  # Get the directory from the file path
    if not os.path.exists(dir_path):
        os.makedirs(dir_path)  # Create directory if it doesn't exist

    # Create and save the pagerank_{file}_10.txt file
    with open(f"{file_name}_10.txt", "w") as f:
        for gene in top_10_genes:
            f.write(f"{gene}\n")

    # Create and save the pagerank_{file}_20.txt file
    with open(f"{file_name}_20.txt", "w") as f:
        for gene in top_20_genes:
            f.write(f"{gene}\n")

    # Create and save the pagerank_{file}_50.txt file
    with open(f"{file_name}_50.txt", "w") as f:
        for gene in top_50_genes:
            f.write(f"{gene}\n")

    print("Files saved: ")
    print(f"{file_name}_10.txt")
    print(f"{file_name}_20.txt")
    print(f"{file_name}_50.txt")

def process_all_features(base_path, G, output_folder):
    """Iterate through all feature files and compute PageRank results."""
    for root, _, files in os.walk(base_path):
        if "features" in root:
            parent_folder = os.path.dirname(root)
            save_path = os.path.join(parent_folder, "pagerank_genes")
            os.makedirs(save_path, exist_ok=True)

            for file in files:
                if file.endswith(".txt"):
                    input_path = os.path.join(root, file)
                    base_name = os.path.splitext(file)[0]  # Strip the file extension from the file name
                    output_file = os.path.join(save_path, base_name)  # Final output file for pagerank results

                    print(f"Processing {input_path}...")

                    results = process_gene_list(G, input_path)
                    top_10_genes = [gene for gene, score in results[:10]]
                    top_20_genes = [gene for gene, score in results[:20]]
                    top_50_genes = [gene for gene, score in results[:50]]

                    with open(f"{output_file}.txt", "w") as out_f:
                        for gene, score in results:
                            out_f.write(f"{gene}\t{score}\n")

                    print(f"Saved results to {output_file}")

                    # Create and save the pagerank_{file}_10.txt file
                    with open(f"{output_file}_10.txt", "w") as f:
                        for gene in top_10_genes:
                            f.write(f"{gene}\n")

                    # Create and save the pagerank_{file}_20.txt file
                    with open(f"{output_file}_20.txt", "w") as f:
                        for gene in top_20_genes:
                            f.write(f"{gene}\n")

                    # Create and save the pagerank_{file}_50.txt file
                    with open(f"{output_file}_50.txt", "w") as f:
                        for gene in top_50_genes:
                            f.write(f"{gene}\n")

                    print("Files saved: ")
                    print(f"{output_file}_10.txt")
                    print(f"{output_file}_20.txt")
                    print(f"{output_file}_50.txt")

# Paths
graphml_path = "data/g_KEGG_UNION_LINCS.graphml"  # Change as needed
base_results_path = "results"  # Base path for results

# Load the graph and extract the largest connected component
G = load_graph(graphml_path)
G_largest = get_largest_connected_component(G)

# Iterate through datasets (gCSI, GDSC2), conditions (negative, positive), and methods (en, ridge, lasso)
for dataset in ["gCSI", "GDSC2"]:
    for condition in ["negative", "positive"]:
        for method in ["en", "ridge", "lasso"]:
            features_path = os.path.join(base_results_path, dataset, condition, method, "features")
            if os.path.exists(features_path):
                process_all_features(features_path, G_largest, "pagerank_genes")
