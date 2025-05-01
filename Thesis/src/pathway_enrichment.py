import pandas as pd
from collections import defaultdict
import numpy as np
import scipy.stats as stat
from statsmodels.stats.multitest import multipletests
import matplotlib.pyplot as plt
import os
from matplotlib_venn import venn3
import matplotlib.pyplot as plt
from itertools import combinations


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


def get_feature_list(is_GDSC, is_positive, biomarker, feature_type): #feature_type = "lasso" or "en" or "ridge"
    target_set = "GDSC2" if is_GDSC else "gCSI"
    experiment = "positive" if is_positive else "negative"
    genes_path = f'./results/{target_set}/{experiment}/{feature_type}/features/{biomarker}.txt'
    if not os.path.exists(genes_path):
        return set()
    with open(genes_path, 'r') as f:
        return set(line.strip() for line in f)

def get_feature_list_correlations(is_GDSC, is_positive, biomarker, feature_type): #feature_type = "lasso" or "en" or "ridge"
    target_set = "GDSC2" if is_GDSC else "gCSI"
    experiment = "positive" if is_positive else "negative"
    genes_path = f'./results/{target_set}/{experiment}/{feature_type}/correlations/{biomarker}_7.txt'
    if not os.path.exists(genes_path):
        return set()
    with open(genes_path, 'r') as f:
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
    else:
        print(f"\nNot enrichment found for {biomarker}")

def venn_plot(dic1, dic2, dic3, labels, output_path):
    set1 = set([p for lst in dic1.values() for p in lst])
    set2 = set([p for lst in dic2.values() for p in lst])
    set3 = set([p for lst in dic3.values() for p in lst])

    plt.figure(figsize=(8, 6))
    venn3([set1, set2, set3], set_labels=labels)
    plt.title("Overlap of Enriched Pathways")
    plt.tight_layout()
    plt.savefig(output_path)
    plt.show()

# Load data
reactome_pathways = reactome_genes()
with open('./data/background_genes.txt', 'r') as f:
    background_genes = set(line.strip() for line in f)
with open('./data/common_drugs.txt', 'r') as f:
    drugs = set(line.strip() for line in f)



base_results_path = os.path.join("results", "pathway_enrichment")  # Base path for results

########################################### PATHWAY ENRICHMENT TO COMPARE DEGREE, BETWEENNESS, EIGENVECTOR AND PAGERANK results ##########################################
centrality_files = {
    "degree": "results/centrality_measures/degree_top_50.txt",
    "betweenness": "results/centrality_measures/betweenness_top_50.txt",
    "eigenvector": "results/centrality_measures/eigenvector_top_50.txt",
    "pagerank": "results/centrality_measures/pagerank_top_50.txt"
}

bp_dic_centrality = {}

for measure, filepath in centrality_files.items():
    if not os.path.exists(filepath):
        print(f"File not found: {filepath}")
        continue

    with open(filepath, 'r') as f:
        centrality_genes = set(line.strip() for line in f)

    if not centrality_genes:
        print(f"No genes found in {filepath}")
        continue

    enriched_centrality_df = perform_enrichment_analysis(centrality_genes, background_genes, reactome_pathways)
    dir = os.path.join(base_results_path, "centrality_measures")

    if not os.path.exists(dir):
        os.makedirs(dir)

    save_results_enrichment(enriched_centrality_df, dir, measure, bp_dic_centrality)

# Optional: Generate a Venn plot for centrality measures
fig_folder = "figs/CENTRALITY_MEASURES"
if not os.path.exists(fig_folder):
    os.makedirs(fig_folder)


# Generate all combinations of 3 measures out of the 4
measures = ["degree", "betweenness", "eigenvector", "pagerank"]
measure_combinations = list(combinations(measures, 3))

# Create Venn plots for each combination
for combo in measure_combinations:
    fig_folder = "figs/CENTRALITY_MEASURES_COMBINATIONS"
    if not os.path.exists(fig_folder):
        os.makedirs(fig_folder)

    venn_plot(
        {combo[0]: bp_dic_centrality.get(combo[0], [])},
        {combo[1]: bp_dic_centrality.get(combo[1], [])},
        {combo[2]: bp_dic_centrality.get(combo[2], [])},
        combo,
        os.path.join(fig_folder, f"venn_{'_'.join(combo)}.png")
    )


########################################### PATHWAY ENRICHMENT TO COMPARE LASSO, PAGERANK_LASSO AND PAGERANK_DRUGBANK ##########################################
bp_dic_lasso = {}
bp_dic_pagerank_lasso = {}
bp_dic_pagerank_drugbank = {}
for biomarker in drugs:

    is_GDSC = False
    is_positive = False
    target_set = "GDSC2" if is_GDSC else "gCSI"
    experiment = "positive" if is_positive else "negative"

    lasso_genes = get_feature_list(is_GDSC, is_positive, biomarker, "lasso")
    pagerank_from_lasso = get_pagerank_lasso_genes(is_GDSC, is_positive, biomarker)
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


################################################ Venn plot ################################################
fig_folder = "figs/LASSO_PAGERANK1_PAGERANK2"
if not os.path.exists(fig_folder):
            os.makedirs(fig_folder)  # Create directory if it doesn't exist
venn_plot(bp_dic_lasso, bp_dic_pagerank_lasso, bp_dic_pagerank_drugbank, ('Lasso', 'PageRank from Lasso', 'PageRank from DrugBank'), os.path.join(fig_folder, "pathway_enrichment_venn.png"))



########################################### PATHWAY ENRICHMENT TO COMPARE LASSO, EN AND RIDGE RESULTS ##########################################

bp_dic_en = {}
bp_dic_ridge = {}

for biomarker in drugs:
    # Perform enrichment analysis for Elastic Net (EN) features
   # en_genes = get_feature_list(True, True, biomarker, "en").union(get_feature_list(True, False, biomarker, "en"))  # Combine positive and negative EN genes
    is_GDSC = False
    is_positive = False
    target_set = "GDSC2" if is_GDSC else "gCSI"
    experiment = "positive" if is_positive else "negative"

    en_genes = get_feature_list(is_GDSC, is_positive, biomarker, "en") 

    if en_genes:
        enriched_en_df = perform_enrichment_analysis(en_genes, background_genes, reactome_pathways)
        dir = os.path.join(base_results_path, "en_features")
        save_results_enrichment(enriched_en_df, dir, biomarker, bp_dic_en)


    # Perform enrichment analysis for Ridge features
   # ridge_genes = get_feature_list(False, True, biomarker, "ridge").union(get_feature_list(False, False, biomarker, "ridge"))  # Combine positive and negative Ridge genes
    ridge_genes = get_feature_list(is_GDSC, is_positive, biomarker, "ridge")

    if ridge_genes:
        enriched_ridge_df = perform_enrichment_analysis(ridge_genes, background_genes, reactome_pathways)
        dir = os.path.join(base_results_path, "ridge_features")
        save_results_enrichment(enriched_ridge_df, dir, biomarker, bp_dic_ridge)



################################## PATHWAY ENRICHMENT IN HIGHEST CORRELATED GENES ##############################################

bp_dic_lasso_corr = {}
bp_dic_en_corr = {}
bp_dic_ridge_corr = {}

for biomarker in drugs:
    is_GDSC = False
    is_positive = False
    target_set = "GDSC2" if is_GDSC else "gCSI"
    experiment = "positive" if is_positive else "negative"

    lasso_genes_corr = get_feature_list_correlations(is_GDSC, is_positive, biomarker, "lasso")
    en_genes_corr = get_feature_list_correlations(is_GDSC, is_positive, biomarker, "en")
    ridge_genes_corr = get_feature_list_correlations(is_GDSC, is_positive, biomarker, "ridge")

    
    # Perform enrichment analysis for lasso_genes    
    if not lasso_genes_corr:
        continue
        
    enriched_lasso_df_corr = perform_enrichment_analysis(lasso_genes_corr, background_genes, reactome_pathways)
    dir = os.path.join(base_results_path, f"aac_correlations_lasso/{target_set}/{experiment}")
    save_results_enrichment(enriched_lasso_df_corr, dir, biomarker, bp_dic_lasso_corr)
    
    # Perform enrichment analysis for en_genes    
    if not en_genes_corr:
        continue
        
    enriched_en_df_corr = perform_enrichment_analysis(en_genes_corr, background_genes, reactome_pathways)
    dir = os.path.join(base_results_path, f"aac_correlations_en/{target_set}/{experiment}")
    save_results_enrichment(enriched_en_df_corr, dir, biomarker, bp_dic_en_corr)

    # Perform enrichment analysis for ridge_genes    
    if not ridge_genes_corr:
        continue
        
    enriched_ridge_df_corr = perform_enrichment_analysis(ridge_genes_corr, background_genes, reactome_pathways)
    dir = os.path.join(base_results_path, f"aac_correlations_ridge/{target_set}/{experiment}")
    save_results_enrichment(enriched_ridge_df_corr, dir, biomarker, bp_dic_ridge_corr)


#################################### Pathway enrichment for drugbank dt ######################################
bp_dic_drugbank = {}

drugbank_file = "results/drugbank/AffectedGenesByDrug.txt"
if os.path.exists(drugbank_file):
    with open(drugbank_file, 'r') as f:
        next(f)  # Skip header
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) < 3:
                continue
            drug_name = parts[1]
            affected_genes = set(parts[2].split(', '))

            enriched_drugbank_df = perform_enrichment_analysis(affected_genes, background_genes, reactome_pathways)
            dir = os.path.join(base_results_path, "drugbank")
            save_results_enrichment(enriched_drugbank_df, dir, drug_name, bp_dic_drugbank)
else:
    print(f"File not found: {drugbank_file}")


############################################## Pathway enrichment for DTC ##############################################
bp_dic_dtc = {}

dtc_file = "results/dtc/AffectedGenesByDrug.csv"
if os.path.exists(dtc_file):
    dtc_data = pd.read_csv(dtc_file)
    for _, row in dtc_data.iterrows():
        drug_name = row['compound_name']
        if pd.notna(row['gene_names']):
            affected_genes = set(row['gene_names'].split(', '))
        else:
            affected_genes = set()

        enriched_dtc_df = perform_enrichment_analysis(affected_genes, background_genes, reactome_pathways)
        dir = os.path.join(base_results_path, "dtc")
        save_results_enrichment(enriched_dtc_df, dir, drug_name, bp_dic_dtc)
else:
    print(f"File not found: {dtc_file}")

################################################ Venn plot ################################################
fig_folder = "figs/LASSO_EN_RIDGE"
if not os.path.exists(fig_folder):
            os.makedirs(fig_folder)
venn_plot(bp_dic_lasso, bp_dic_en, bp_dic_ridge, ('Lasso', 'Elastic Net', 'Ridge'), os.path.join(fig_folder, "pathway_enrichment_venn.png"))


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










############################################## HEATMAPS ##############################################  

import pandas as pd

def make_binary_matrix(bp_dic):
    all_pathways = set(p for lst in bp_dic.values() for p in lst)
    matrix = pd.DataFrame(0, index=bp_dic.keys(), columns=sorted(all_pathways))
    
    for drug, pathways in bp_dic.items():
        matrix.loc[drug, pathways] = 1
    
    return matrix

lasso_bin = make_binary_matrix(bp_dic_lasso)
pr_lasso_bin = make_binary_matrix(bp_dic_pagerank_lasso)
drugbank_bin = make_binary_matrix(bp_dic_pagerank_drugbank)



import seaborn as sns


def plot_heatmap(matrix, title, fname):
    plt.figure(figsize=(12, 6))
    sns.heatmap(matrix, cmap="Blues", cbar=True)
    plt.title(title)
    plt.ylabel("Drug")
    plt.xlabel("Pathway")
    plt.xticks(rotation=90)
    plt.tight_layout()
    plt.savefig(fname)
    plt.show()

#plot_heatmap(lasso_bin, "Enriched Pathways per Drug (Lasso)", "figs/lasso_pathways_heatmap.png")
#plot_heatmap(pr_lasso_bin, "Enriched Pathways per Drug (PR from Lasso)", "figs/pr_lasso_pathways_heatmap.png")
#plot_heatmap(drugbank_bin, "Enriched Pathways per Drug (DrugBank)", "figs/drugbank_pathways_heatmap.png")