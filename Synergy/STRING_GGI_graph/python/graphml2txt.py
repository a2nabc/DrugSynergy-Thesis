import networkx as nx

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

graph_inputs = {
    "protein_coding": "../g_protein_coding.graphml",
    "COSMIC": "../g_COSMIC.graphml",
    "KEGG": "../g_kegg.graphml",
    "LINCS": "../g_LINCS.graphml",
    "MSig_Hallmark": "../g_mdsig_hallmarks.graphml"
}


# Process each subset
def process_all_graphs(graph_inputs):
    for name, graph_file in graph_inputs.items():
        output_file = f"{name}_output.txt"
        convert_graphml_to_txt(graph_file, output_file)
        output_file_lcc = f"{name}_output_LCC.txt"
        convert_graphml_to_txt_LCC(graph_file, output_file_lcc)

# Run the processing function
process_all_graphs(graph_inputs)