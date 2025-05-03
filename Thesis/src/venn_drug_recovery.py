import matplotlib.pyplot as plt
from matplotlib_venn import venn2

# File paths
drugbank_file = "data/drugbank_gene_set.txt"
betweenness_file = "results/centrality_measures/betweenness_top_all.txt"
degree_file = "results/centrality_measures/degree_top_all.txt"
eigenvector_file = "results/centrality_measures/eigenvector_top_all.txt"
pagerank_file = "results/centrality_measures/pagerank_top_all.txt"

# Function to read and filter gene sets from files
def read_filtered_gene_set(file_path):
    with open(file_path, 'r') as file:
        return set(line.split('\t')[0] for line in file if float(line.split('\t')[1]) != 0)

# Replace the read_gene_set function with the filtered version
read_gene_set = read_filtered_gene_set

# Function to read gene sets from files
def read_gene_set(file_path):
    with open(file_path, 'r') as file:
        return set(line.strip() for line in file)

# Read gene sets
drugbank_genes = read_gene_set(drugbank_file)
betweenness_genes = read_gene_set(betweenness_file)
degree_genes = read_gene_set(degree_file)
eigenvector_genes = read_gene_set(eigenvector_file)
pagerank_genes = read_gene_set(pagerank_file)

# Function to create and save Venn diagrams
def create_venn(set1, set2, label1, label2, output_file):
    venn = venn2([set1, set2], (label1, label2))
    plt.title(f"Venn Diagram: {label1} vs {label2}")
    plt.savefig(output_file)
    plt.close()

# Generate Venn diagrams
create_venn(drugbank_genes, betweenness_genes, "DrugBank", "Betweenness", "figs/DRUGBANK_CENTRALITY/venn_drugbank_betweenness.png")
create_venn(drugbank_genes, degree_genes, "DrugBank", "Degree", "figs/DRUGBANK_CENTRALITY/venn_drugbank_degree.png")
create_venn(drugbank_genes, eigenvector_genes, "DrugBank", "Eigenvector", "figs/DRUGBANK_CENTRALITY/venn_drugbank_eigenvector.png")
create_venn(drugbank_genes, pagerank_genes, "DrugBank", "PageRank", "figs/DRUGBANK_CENTRALITY/venn_drugbank_pagerank.png")

print("Venn diagrams have been saved as PNG files.")