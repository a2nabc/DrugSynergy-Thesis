import matplotlib.pyplot as plt
from matplotlib_venn import venn2
from matplotlib_venn import venn3

# File paths
drugbank_file = "data/drugbank_gene_set.txt"
dtc_file = "data/dtc_gene_set.txt"
betweenness_file = "results/centrality_measures/betweenness_top_50.txt"
degree_file = "results/centrality_measures/degree_top_50.txt"
eigenvector_file = "results/centrality_measures/eigenvector_top_50.txt"
pagerank_file = "results/centrality_measures/pagerank_top_50.txt"

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
dtc_genes = read_gene_set(dtc_file)
betweenness_genes = read_gene_set(betweenness_file)
degree_genes = read_gene_set(degree_file)
eigenvector_genes = read_gene_set(eigenvector_file)
pagerank_genes = read_gene_set(pagerank_file)

# Function to create and save Venn diagrams for 2 sets
def create_venn(set1, set2, label1, label2, output_file):
    venn = venn2([set1, set2], (label1, label2))
    plt.title(f"Venn Diagram: {label1} vs {label2}")
    plt.savefig(output_file)
    plt.close()

# Function to create and save Venn diagrams for 3 sets
def create_venn_three(set1, set2, set3, label1, label2, label3, output_file):
    venn = venn3([set1, set2, set3], (label1, label2, label3))
    plt.title(f"Venn Diagram: {label1} vs {label2} vs {label3}")
    plt.savefig(output_file)
    plt.close()

# Generate Venn diagrams for 2 sets
create_venn(drugbank_genes, betweenness_genes, "DrugBank", "Betweenness", "figs/DRUGBANK_CENTRALITY/venn_drugbank_betweenness.png")
create_venn(drugbank_genes, degree_genes, "DrugBank", "Degree", "figs/DRUGBANK_CENTRALITY/venn_drugbank_degree.png")
create_venn(drugbank_genes, eigenvector_genes, "DrugBank", "Eigenvector", "figs/DRUGBANK_CENTRALITY/venn_drugbank_eigenvector.png")
create_venn(drugbank_genes, pagerank_genes, "DrugBank", "PageRank", "figs/DRUGBANK_CENTRALITY/venn_drugbank_pagerank.png")

create_venn(dtc_genes, betweenness_genes, "DTC", "Betweenness", "figs/DRUGBANK_CENTRALITY/venn_dtc_betweenness.png")
create_venn(dtc_genes, degree_genes, "DTC", "Degree", "figs/DRUGBANK_CENTRALITY/venn_dtc_degree.png")
create_venn(dtc_genes, eigenvector_genes, "DTC", "Eigenvector", "figs/DRUGBANK_CENTRALITY/venn_dtc_eigenvector.png")
create_venn(dtc_genes, pagerank_genes, "DTC", "PageRank", "figs/DRUGBANK_CENTRALITY/venn_dtc_pagerank.png")

# Generate Venn diagrams for 3 sets
create_venn_three(drugbank_genes, dtc_genes, betweenness_genes, "DrugBank", "DTC", "Betweenness", "figs/DRUGBANK_CENTRALITY/venn_drugbank_dtc_betweenness.png")
create_venn_three(drugbank_genes, dtc_genes, degree_genes, "DrugBank", "DTC", "Degree", "figs/DRUGBANK_CENTRALITY/venn_drugbank_dtc_degree.png")
create_venn_three(drugbank_genes, dtc_genes, eigenvector_genes, "DrugBank", "DTC", "Eigenvector", "figs/DRUGBANK_CENTRALITY/venn_drugbank_dtc_eigenvector.png")
create_venn_three(drugbank_genes, dtc_genes, pagerank_genes, "DrugBank", "DTC", "PageRank", "figs/DRUGBANK_CENTRALITY/venn_drugbank_dtc_pagerank.png")


print("Venn diagrams have been saved as PNG files.")