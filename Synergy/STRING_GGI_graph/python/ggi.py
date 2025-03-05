import networkx as nx
import urllib.request
import os
import pandas as pd
import argparse
import tqdm
from typing import Dict

# Define gene subset files
gene_subset_files = {
    "protein_coding": "gene_list.txt",
    "COSMIC": "../../data/single_drug/gene_subsets/COSMIC.txt",
    "KEGG": "../../data/single_drug/gene_subsets/kegg.txt",
    "LINCS": "../../data/single_drug/gene_subsets/LINCS.txt",
    "MSig_Hallmark": "../../data/single_drug/gene_subsets/mdsig_hallmarks.txt"
}

# Function to map Ensembl protein IDs to HUGO gene symbols
def map_ensemble_ids_to_hugo_symbols(aliases: pd.DataFrame) -> Dict[str, str]:
    id_to_hugo_symbol = {}
    print("\nMapping Ensemble IDs to HUGO Symbols...")
    for ens_id in tqdm.tqdm(pd.unique(aliases['#string_protein_id'])):
        temp = aliases[aliases['#string_protein_id'] == ens_id]
        unique_genes = pd.unique(temp['alias'])
        if len(unique_genes) == 1:
            id_to_hugo_symbol[ens_id] = unique_genes[0]
    return id_to_hugo_symbol

# Main function
def main(min_score: int):
    # Create output directories
    os.makedirs("./data/", exist_ok=True)
    os.makedirs("./edgelists/", exist_ok=True)

    # STRING database links
    link_file = "https://stringdb-downloads.org/download/protein.links.v12.0/9606.protein.links.v12.0.txt.gz"
    alias_file = "https://stringdb-downloads.org/download/protein.aliases.v12.0/9606.protein.aliases.v12.0.txt.gz"

    link_file_local = "./data/9606.protein.links.v12.0.txt.gz"
    link_file_out = "./data/STRING_links.txt"
    alias_file_local = "./data/9606.protein.aliases.v12.0.txt.gz"  # Fixed file name
    alias_file_out = "./data/STRING_aliases.txt"

    # Download STRING data if needed
    if not os.path.isfile(link_file_out):
        urllib.request.urlretrieve(link_file, link_file_local)
        os.system(f"gunzip -c {link_file_local} > {link_file_out}")

    if not os.path.isfile(alias_file_out):
        urllib.request.urlretrieve(alias_file, alias_file_local)
        os.system(f"gunzip -c {alias_file_local} > {alias_file_out}")

    # Load STRING interaction data
    string_links = pd.read_csv(link_file_out, sep=" ")
    string_links['combined_score'] = string_links['combined_score'].astype(int)

    # Load protein-to-gene alias mapping
    annotation = pd.read_csv(alias_file_out, sep="\t")
    hugo_aliases = annotation[annotation['source'].str.contains('HUGO', na=False)].reset_index(drop=True)

    # Create mapping from Ensembl protein IDs to HUGO gene symbols
    id_to_hugo_map = map_ensemble_ids_to_hugo_symbols(hugo_aliases)

    # Process each gene subset
    for subset_name, subset_path in gene_subset_files.items():
        print(f"\nProcessing subset: {subset_name}...")

        # Load gene subset
        if not os.path.isfile(subset_path):
            print(f"Skipping {subset_name} (File not found: {subset_path})")
            continue
        
        with open(subset_path, "r") as f:
            gene_subset = set(f.read().splitlines())  # Store genes in a set

        # Create a new graph for this subset
        G = nx.Graph()

        print(f"Building network for {subset_name}, please wait...")
        for idx, row in tqdm.tqdm(string_links.iterrows(), total=string_links.shape[0]):
            prot1, prot2, score = row['protein1'], row['protein2'], int(row['combined_score'])

            if score >= min_score:
                if prot1 in id_to_hugo_map and prot2 in id_to_hugo_map:
                    gene1, gene2 = id_to_hugo_map[prot1], id_to_hugo_map[prot2]

                    # Keep only edges where both genes are in the subset
                    if gene1 in gene_subset and gene2 in gene_subset:
                        G.add_edge(gene1, gene2)

        # Save full network for the subset
        edge_list = [f"{e[0]}\t{e[1]}\n" for e in G.edges()]
        full_output_path = f"./edgelists/full_{min_score}_{subset_name}.txt"
        with open(full_output_path, "w") as ostream:
            ostream.writelines(edge_list)

        print(f"Saved FULL network: {full_output_path} ({len(G.edges())} edges)")

        # Extract and save largest connected component (LCC)
        if G.number_of_nodes() > 0:
            LCC_genes = max(nx.connected_components(G), key=len)
            LCC = G.subgraph(LCC_genes)

            edge_list = [f"{e[0]}\t{e[1]}\n" for e in LCC.edges()]
            lcc_output_path = f"./edgelists/LCC_{min_score}_{subset_name}.txt"
            with open(lcc_output_path, "w") as ostream:
                ostream.writelines(edge_list)

            print(f"Saved LCC network: {lcc_output_path} ({len(LCC.edges())} edges)")

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        prog="Building String Networks for Gene Subsets",
        description="")

    parser.add_argument("-ms", help="The Minimum Score", required=True, type=int)

    args = parser.parse_args()
    main(args.ms)
