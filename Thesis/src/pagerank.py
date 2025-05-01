import os
import networkx as nx
import pandas as pd

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

def process_gene_list(G, gene_list):
    """Process a single gene list file and return the top PageRank results."""

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
    #sorted_genes = sorted(pagerank_scores.items(), key=lambda x: x[1], reverse=True)

    # Get the top 10, 20, and 50 genes (no scores, just names)
    top_10_genes = [gene for gene, score in pagerank_scores[:10]]
    top_20_genes = [gene for gene, score in pagerank_scores[:20]]
    top_50_genes = [gene for gene, score in pagerank_scores[:50]]

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

def process_all_features(base_path, dataset, condition, method, G, output_folder):
    """Iterate through all feature files and compute PageRank results."""
    for root, _, files in os.walk(base_path):
        if "features" in root:
            parent_folder = os.path.dirname(root)
            save_path = os.path.join(parent_folder, "pagerank_genes")
            os.makedirs(save_path, exist_ok=True)

            for file in files:
                if file.endswith(".txt"):
                    input_path = os.path.join(root, file)
                    base_name = os.path.basename(input_path)  # Get the parent folder name --> drugname
                    drug_name = os.path.splitext(base_name)[0]  # Strip the file extension from the file name ---> basename is the drug
                    output_file = os.path.join(output_folder, method, dataset, condition, drug_name)  # Final output file for pagerank results

                    print(f"Processing {input_path}...")

                    with open(input_path, "r") as f:
                        gene_list = {line.strip() for line in f}

                    results = process_gene_list(G, gene_list)
                    save_top_genes(results, output_file)

def run_pagerank_from_drugbank(path_drugbank_genes, G):

    # Load DrugBank data
    drugbank_data = pd.read_csv(path_drugbank_genes, sep="\t", header=0)
    drug_gene_map = {
        row['DRUG_NAME']: [gene.strip() for gene in row['GENES_AFFECTED'].split(",")]
        for _, row in drugbank_data.iterrows()
    }


    # Run PageRank with the filtered gene list
    for drug in drug_gene_map:
        gene_list = drug_gene_map[drug]
        print(f"Processing drug: {drug} with genes: {gene_list}")

        # Compute PageRank scores
        results = process_gene_list(G, gene_list)
        output_path = os.path.join(base_results_path, "drugbank", f"{drug}")
        save_top_genes(results, output_path)

        # Save results
        #save_top_genes(pagerank_scores, os.path.join("results", "pagerank_genes", f"pagerank_{drug}"))

    


# Paths
graphml_path = "data/g_KEGG_UNION_LINCS.graphml"  # Change as needed
base_results_path = os.path.join("results", "pagerank_output")  # Base path for results
results_path = "results"  # Base path for results

# Load the graph and extract the largest connected component
G = load_graph(graphml_path)
G_largest = get_largest_connected_component(G)
G_simple = nx.Graph(G_largest)

# Compute basic graph theory measures
deg = nx.degree_centrality(G_simple)
print("Degree centrality computed")
btw = nx.betweenness_centrality(G_simple)
print("Betweenness centrality computed")
eig = nx.eigenvector_centrality(G_simple, max_iter=1000)
print("Eigenvector centrality computed")
pr = nx.pagerank(G_simple, personalization=None)  # or just omit it

centrality_dir = "results/centrality_measures"
os.makedirs(centrality_dir, exist_ok=True)


for (centrality, name) in zip([deg, btw, eig, pr], ["degree", "betweenness", "eigenvector", "pagerank"]):
    sorted_nodes = sorted(centrality.items(), key=lambda x: x[1], reverse=True)
    top_genes_named = []
    for node_id, score in sorted_nodes:
        gene_name = G_simple.nodes[node_id].get("name", node_id)
        top_genes_named.append((gene_name, score))
    output_base = os.path.join(centrality_dir, f"{name}_top")
    save_top_genes(top_genes_named, output_base)



################################################ RUN PAGERANK in features from lasso/en/rige #########################################

# Iterate through datasets (gCSI, GDSC2), conditions (negative, positive), and methods (en, ridge, lasso)
#######for dataset in ["gCSI", "GDSC2"]:
#######    for condition in ["negative", "positive"]:
#######        for method in ["en", "ridge", "lasso"]:
#######            features_path = os.path.join(results_path, dataset, condition, method, "features")
#######            if os.path.exists(features_path):
#######                process_all_features(features_path, dataset, condition, method, G_largest, base_results_path)


###################################################### RUN PAGERANK in drugbank genes ################################################

path_drugbank_info = "results/drugbank/AffectedGenesByDrug.txt"
drugs_path = "data/common_drugs.txt"

# Load the list of drugs from the file
with open(drugs_path, "r") as f:
    drug_list = [line.strip() for line in f]

#pagerank_scores = run_pagerank_from_drugbank(path_drugbank_info, G)
