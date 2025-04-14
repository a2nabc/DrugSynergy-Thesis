import pandas as pd
from collections import defaultdict
import numpy as np
import scipy.stats as stat
from statsmodels.stats.multitest import multipletests
import matplotlib.pyplot as plt
import os


def perform_enrichment_analysis(lasso_genes, background_genes, pathway_dict):
    """
    Performs pathway enrichment analysis using the hypergeometric test.

    Parameters:
    - lasso_genes (set): Set of Lasso-selected gene symbols.
    - background_genes (set): Set of all background gene symbols.
    - pathway_dict (dict): Dictionary mapping pathway names/IDs to sets of gene symbols.

    Returns:
    - DataFrame: Enrichment results with pathway, p-value, FDR q-value, etc.
    """
    results = []
    M = len(background_genes)
    N = len(lasso_genes)

    for pathway_id, pathway_genes in pathway_dict.items():
        pathway_genes = set(pathway_genes)
        overlap_genes = lasso_genes & pathway_genes
        n = len(pathway_genes & background_genes)
        k = len(overlap_genes)

        if n == 0 or k == 0:
            continue  # Skip non-informative pathways

        p_value = stat.hypergeom.sf(k - 1, M, n, N)

        results.append({
            'Pathway ID': pathway_id,
            'Pathway Name': pathway_id,
            'Overlap': k,
            'Pathway Size': n,
            'P-Value': p_value,
            'Overlap Genes': ', '.join(overlap_genes)
        })

    results_df = pd.DataFrame(results)

    if results_df.empty:
        return results_df  # no enrichments found

    _, results_df['FDR Q-Value'], _, _ = multipletests(results_df['P-Value'])
    results_df = results_df.sort_values('FDR Q-Value')


    return results_df


def reactome_genes(filepath='./data/raw/ReactomePathways.gmt'):
    output = defaultdict(list)
    with open(filepath, 'r') as f:
        for line in f:
            parts = line.strip().split('\t')
            pathway = parts[0]
            for gene in parts[2:]:
                output[pathway].append(gene)
    return output


def get_lasso_genes(is_GDSC, is_positive, biomarker):
    target_set = "GDSC2" if is_GDSC else "gCSI"
    experiment = "positive" if is_positive else "negative"
    lasso_genes_path = f'./results/{target_set}/{experiment}/lasso/features/{biomarker}.txt'
    if not os.path.exists(lasso_genes_path):
        return set()
    with open(lasso_genes_path, 'r') as f:
        return set(line.strip() for line in f)

def get_pagerank_lasso_genes(is_GDSC, is_positive, biomarker):
    target_set = "GDSC2" if is_GDSC else "gCSI"
    experiment = "positive" if is_positive else "negative"
    pagerank_genes_path = f'./results/pagerank_output/lasso/{target_set}/{experiment}/{biomarker}_50.txt'
    if not os.path.exists(pagerank_genes_path):
        return set()
    with open(pagerank_genes_path, 'r') as f:
        return set(line.strip() for line in f)
    
def get_drugbank_genes(biomarker):
    path_drugbank_genes = f'./results/pagerank_output/drugbank/{biomarker}_50.txt'
    if not os.path.exists(path_drugbank_genes):
        return set()
    with open(path_drugbank_genes, 'r') as f:
        return set(line.strip() for line in f)

# Load data
reactome_pathways = reactome_genes()
with open('./data/background_genes.txt', 'r') as f:
    background_genes = set(line.strip() for line in f)
with open('./data/common_drugs.txt', 'r') as f:
    drugs = set(line.strip() for line in f)





def save_results_enrichment(enriched_df, results_path, biomarker, bp_dic):
    if not enriched_df.empty:
        bp_dic[biomarker] = enriched_df[enriched_df['FDR Q-Value'] < 0.2]['Pathway Name'].tolist()
        print(f"\nFiltered pathways with FDR q < 0.2 for {biomarker} saved at {results_path}:\n")

        if not os.path.exists(results_path):
            os.makedirs(results_path)  # Create directory if it doesn't exist

        file_name = os.path.join(results_path, biomarker)  # Final output file for pagerank results
        
        with open(f"{file_name}.txt", "w") as f:
            for pathway in bp_dic[biomarker]:
                f.write(f"{pathway}\n")



# Run enrichment for lasso features, and pagerank 
base_results_path = os.path.join("results", "pathway_enrichment")  # Base path for results
bp_dic_lasso = {}
bp_dic_pagerank_lasso = {}
bp_dic_pagerank_drugbank = {}
for biomarker in drugs:
    lasso_genes = get_lasso_genes(True, True, biomarker)
    pagerank_from_lasso = get_pagerank_lasso_genes(False, True, biomarker)
    pagerank_from_drugbank = get_drugbank_genes(biomarker)

    # Perform enrichment analysis for lasso_genes    
    if not lasso_genes:
        continue

    enriched_lasso_df = perform_enrichment_analysis(lasso_genes, background_genes, reactome_pathways)
    dir = os.path.join(base_results_path, "lasso_features")
    save_results_enrichment(enriched_lasso_df, dir, biomarker, bp_dic_lasso)


    # Perform enrichment analysis for pagerank_from_lasso
    if not pagerank_from_lasso:
        continue

    enriched_pagerank_lasso_df = perform_enrichment_analysis(pagerank_from_lasso, background_genes, reactome_pathways)
    dir = os.path.join(base_results_path, "pagerank_from_lasso")
    save_results_enrichment(enriched_pagerank_lasso_df, dir, biomarker, bp_dic_pagerank_lasso)


    # Perform enrichment analysis for pagerank_from_drugbank
    if not pagerank_from_drugbank:
        continue

    enriched_pagerank_drugbank_df = perform_enrichment_analysis(pagerank_from_drugbank, background_genes, reactome_pathways)
    dir = os.path.join(base_results_path, "pagerank_from_drugbank")
    save_results_enrichment(enriched_pagerank_drugbank_df, dir, biomarker, bp_dic_pagerank_drugbank)




# Optional: Plot for one biomarker
#if lasso_genes:
#    results_df = perform_enrichment_analysis(lasso_genes, background_genes, reactome_pathways)
#
#    if not results_df.empty:
#        plt.hist(results_df['P-Value'], bins=50, color='skyblue', edgecolor='black')
#        plt.title('Distribution of Enrichment P-Values')
#        plt.xlabel('P-Value')
#        plt.ylabel('Number of Pathways')
#        plt.grid(True)
#        plt.show()
#
#        top_pathways = results_df.nsmallest(10, 'FDR Q-Value')
#        plt.figure(figsize=(10, 6))
#        plt.barh(top_pathways['Pathway Name'], -np.log10(top_pathways['FDR Q-Value']), color='salmon')
#        plt.xlabel('-log10(FDR Q-Value)')
#        plt.title('Top 10 Enriched Pathways')
#        plt.gca().invert_yaxis()
#        plt.tight_layout()
#        plt.show()
#    else:
#        print("No enriched pathways for the selected Lasso genes.")



