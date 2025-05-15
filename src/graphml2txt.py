import networkx as nx
import numpy as np
import os

def convert_graphml_to_txt(graphml_file, txt_file):
    # Load the GraphML file
    G = nx.read_graphml(graphml_file)
    
    # Open the output text file
    with open(txt_file, 'w') as f:
        for edge in G.edges():
            node1 = G.nodes[edge[0]].get("name", edge[0])  # Get name attribute, fallback to ID
            node2 = G.nodes[edge[1]].get("name", edge[1])  # Get name attribute, fallback to ID
            f.write(f"{node1}\t{node2}\n")

def convert_graphml_to_txt_LCC(graphml_file, txt_file):
    # Load the GraphML file
    G = nx.read_graphml(graphml_file)
    
    # Extract the largest connected component
    if nx.is_directed(G):
        largest_cc = max(nx.strongly_connected_components(G), key=len)
    else:
        largest_cc = max(nx.connected_components(G), key=len)
    G_lcc = G.subgraph(largest_cc).copy()
    
    # Open the output text file
    with open(txt_file, 'w') as f:
        for edge in G_lcc.edges():
            node1 = G_lcc.nodes[edge[0]].get("name", edge[0])  # Get name attribute, fallback to ID
            node2 = G_lcc.nodes[edge[1]].get("name", edge[1])  # Get name attribute, fallback to ID
            f.write(f"{node1}\t{node2}\n")


graph_file = "data/g_KEGG_UNION_LINCS.graphml"


# Convert graphml to txt

output_file = f"data/GGI/g_KEGG_UNION_LINCS.txt"
convert_graphml_to_txt(graph_file, output_file)
output_file_lcc = f"data/GGI/g_KEGG_UNION_LINCS_LCC.txt"
convert_graphml_to_txt_LCC(graph_file, output_file_lcc)






############################################ the following part of the script obtains the neighbours for each extracted feature ###############


def read_gene_list(gene_list_file):
    with open(gene_list_file, 'r') as f:
        genes = [line.strip() for line in f.readlines()]
    return genes


def all_neighbors(graph, node):
    """
    Get all neighbors (direct and indirect) of a node.
    This will return a set of all nodes that are connected to the input node.
    """
    neighbors = set(graph.neighbors(node))  # Direct neighbors
    for neighbor in list(neighbors):
        neighbors.update(graph.neighbors(neighbor))  # Adding second-order neighbors
    neighbors.remove(node)  # Remove the node itself if present
    return neighbors

def get_all_neighbors_for_genes(graphml_file, gene_list_file, output_file):
    # Load the GraphML file
    G = nx.read_graphml(graphml_file)
    
    # Read the list of genes
    genes = read_gene_list(gene_list_file)
    
    # Create a dictionary to map gene names to node IDs
    gene_to_node = {data["name"]: node for node, data in G.nodes(data=True) if "name" in data}
    
    # Open the output text file
    with open(output_file, 'w') as f:
        for gene in genes:
            # Check if the gene exists in the graph
            if gene in gene_to_node:
                node = gene_to_node[gene]  # Get the node ID corresponding to the gene
                neighbors = all_neighbors(G, node)  # Get all neighbors (direct and indirect)
                for neighbor in neighbors:
                    # Write the gene and its neighbor to the file
                    node1 = gene  # Gene name
                    node2 = G.nodes[neighbor].get("name", neighbor)  # Get the name attribute of the neighbor
                    f.write(f"{node1}\t{node2}\n")
            else:
                print(f"Gene {gene} not found in the graph.")


# Example usage
graphml_file = "../g_COSMIC.graphml"  # Example GraphML file
gene_list_file = "../../Results/Single_Drug/COSMIC/CLASSIFICATION/Bortezomib.txt"  # File containing the list of genes
output_file = "../../Results/Single_Drug/COSMIC/CLASSIFICATION/Bortezomib_neighbours"  # Output file to store the results

def find_gene_list_files(root_folder):
    gene_list_files = []
    for dirpath, _, filenames in os.walk(root_folder):
        for filename in filenames:
            if filename.endswith(".txt") and not filename.endswith("_metrics.txt"):
                gene_list_files.append(os.path.join(dirpath, filename))
    return gene_list_files

# Find all gene list files in the specified directory
root_folder = "../../Results2/Single_Drug/"
gene_list_files = find_gene_list_files(root_folder)

# Run the function for each gene list file
for gene_list_file in gene_list_files:
    # Determine the appropriate graphml file based on the second to last folder name
    second_last_folder = os.path.basename(os.path.dirname(os.path.dirname(gene_list_file)))
    graphml_file = f"../../Results2/Single_Drug/{second_last_folder}/g_{second_last_folder}.graphml"
    
    # Set the output file name
    output_file = gene_list_file.replace(".txt", "_neighbours.txt")
    get_all_neighbors_for_genes(graphml_file, gene_list_file, output_file)
# Run the function
get_all_neighbors_for_genes(graphml_file, gene_list_file, output_file)