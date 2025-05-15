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

def get_pagerank_dtc_genes(biomarker):
    path_dtc_genes = f'./results/pagerank_output/dtc/drugwise/{biomarker}_50.txt'
    if not os.path.exists(path_dtc_genes):
        return set()
    with open(path_dtc_genes, 'r') as f:
        return set(line.strip() for line in f)
    
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
    path_drugbank_genes = f'./results/pagerank_output/drugbank/drugwise/{biomarker}_50.txt'
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
    #plt.show()

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

# Get all files ending with .txt and _cancer.txt
txt_files = sorted([f for f in os.listdir(os.path.join(base_results_path, "centrality_measures")) if f.endswith(".txt") and not f.endswith("_cancer.txt")])
cancer_txt_files = sorted([f for f in os.listdir(os.path.join(base_results_path, "centrality_measures")) if f.endswith("_cancer.txt")])

# Generate combinations of 3 for each group
txt_combinations = list(combinations(txt_files, 3))
cancer_txt_combinations = list(combinations(cancer_txt_files, 3))

# Create Venn plots for each combination of .txt files
for combo in txt_combinations:
    fig_folder = "figs/CENTRALITY_MEASURES_COMBINATIONS"
    if not os.path.exists(fig_folder):
        os.makedirs(fig_folder)

    # Read pathways from saved files for each measure
    sets = []
    for measure in combo:
        file_path = os.path.join(base_results_path, "centrality_measures", measure)
        if os.path.exists(file_path):
            with open(file_path, 'r') as f:
                pathways = [line.strip() for line in f]
                sets.append(set(pathways))
        else:
            sets.append(set())  # Empty set if file doesn't exist

    venn_plot(
        {combo[0]: sets[0]},
        {combo[1]: sets[1]},
        {combo[2]: sets[2]},
        combo,
        os.path.join(fig_folder, f"venn_{'_'.join([os.path.splitext(c)[0] for c in combo])}.png")
    )

# Create Venn plots for each combination of _cancer.txt files
for combo in cancer_txt_combinations:
    fig_folder = "figs/CENTRALITY_MEASURES_COMBINATIONS_CANCER"
    if not os.path.exists(fig_folder):
        os.makedirs(fig_folder)

    # Read pathways from saved files for each measure
    sets = []
    for measure in combo:
        file_path = os.path.join(base_results_path, "centrality_measures", measure)
        if os.path.exists(file_path):
            with open(file_path, 'r') as f:
                pathways = [line.strip() for line in f]
                sets.append(set(pathways))
        else:
            sets.append(set())  # Empty set if file doesn't exist

    venn_plot(
        {combo[0]: sets[0]},
        {combo[1]: sets[1]},
        {combo[2]: sets[2]},
        combo,
        os.path.join(fig_folder, f"venn_{'_'.join([os.path.splitext(c)[0] for c in combo])}.png")
    )


########################################### PATHWAY ENRICHMENT TO COMPARE LASSO, PAGERANK_LASSO AND PAGERANK_DRUGBANK ##########################################


########################################### PATHWAY ENRICHMENT TO COMPARE LASSO, PAGERANK_LASSO, PAGERANK_DRUGBANK AND PAGERANK_DTC ##########################################
bp_dic_lasso = {}
bp_dic_pagerank_lasso = {}
bp_dic_pagerank_drugbank = {}
bp_dic_pagerank_dtc = {}

for biomarker in drugs:

    is_GDSC = False
    is_positive = True
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

    # Save cancer-specific enrichment results
    enriched_lasso_df_cancer = enriched_lasso_df[enriched_lasso_df['Pathway Name'].str.contains('cancer', case=False)]
    if not enriched_lasso_df_cancer.empty:
        enriched_lasso_df_cancer['Pathway Name'].to_csv(os.path.join(dir, f"{biomarker}_cancer.txt"), sep='\t', index=False, header=False)

    # Perform enrichment analysis for pagerank_from_lasso
    if not pagerank_from_lasso:
        continue

    enriched_pagerank_lasso_df = perform_enrichment_analysis(pagerank_from_lasso, background_genes, reactome_pathways)
    dir = os.path.join(base_results_path, "pagerank_from_lasso")
    save_results_enrichment(enriched_pagerank_lasso_df, dir, biomarker, bp_dic_pagerank_lasso)

    # Save cancer-specific enrichment results
    enriched_pagerank_lasso_df_cancer = enriched_pagerank_lasso_df[enriched_pagerank_lasso_df['Pathway Name'].str.contains('cancer', case=False)]
    if not enriched_pagerank_lasso_df_cancer.empty:
        enriched_pagerank_lasso_df_cancer['Pathway Name'].to_csv(os.path.join(dir, f"{biomarker}_cancer.txt"), sep='\t', index=False, header=False)

    # Perform enrichment analysis for pagerank_from_drugbank
    if not pagerank_from_drugbank:
        continue

    enriched_pagerank_drugbank_df = perform_enrichment_analysis(pagerank_from_drugbank, background_genes, reactome_pathways)
    dir = os.path.join(base_results_path, "pagerank_from_drugbank")
    save_results_enrichment(enriched_pagerank_drugbank_df, dir, biomarker, bp_dic_pagerank_drugbank)

    # Save cancer-specific enrichment results
    enriched_pagerank_drugbank_df_cancer = enriched_pagerank_drugbank_df[enriched_pagerank_drugbank_df['Pathway Name'].str.contains('cancer', case=False)]
    if not enriched_pagerank_drugbank_df_cancer.empty:
        enriched_pagerank_drugbank_df_cancer['Pathway Name'].to_csv(os.path.join(dir, f"{biomarker}_cancer.txt"), sep='\t', index=False, header=False)

    # Perform enrichment analysis for pagerank_from_dtc
    pagerank_from_dtc = get_pagerank_dtc_genes(biomarker)
    if not pagerank_from_dtc:
        continue

    enriched_pagerank_dtc_df = perform_enrichment_analysis(pagerank_from_dtc, background_genes, reactome_pathways)
    dir = os.path.join(base_results_path, "pagerank_from_dtc")
    save_results_enrichment(enriched_pagerank_dtc_df, dir, biomarker, bp_dic_pagerank_dtc)

    # Save cancer-specific enrichment results
    enriched_pagerank_dtc_df_cancer = enriched_pagerank_dtc_df[enriched_pagerank_dtc_df['Pathway Name'].str.contains('cancer', case=False)]
    if not enriched_pagerank_dtc_df_cancer.empty:
        enriched_pagerank_dtc_df_cancer['Pathway Name'].to_csv(os.path.join(dir, f"{biomarker}_cancer.txt"), sep='\t', index=False, header=False)


################################################ Venn plot ################################################
fig_folder = "figs/LASSO_PAGERANK1_PAGERANK2_CANCER"
if not os.path.exists(fig_folder):
    os.makedirs(fig_folder)  # Create directory if it doesn't exist

# Filter cancer-specific pathways
bp_dic_lasso_cancer = {k: [p for p in v if isinstance(p, str) and 'cancer' in p.lower()] for k, v in bp_dic_lasso.items()}
bp_dic_pagerank_lasso_cancer = {k: [p for p in v if isinstance(p, str) and 'cancer' in p.lower()] for k, v in bp_dic_pagerank_lasso.items()}
bp_dic_pagerank_drugbank_cancer = {k: [p for p in v if isinstance(p, str) and 'cancer' in p.lower()] for k, v in bp_dic_pagerank_drugbank.items()}
bp_dic_pagerank_dtc_cancer = {k: [p for p in v if isinstance(p, str) and 'cancer' in p.lower()] for k, v in bp_dic_pagerank_dtc.items()}


# Generate Venn plot for cancer-specific pathways for PPR-Lasso, PPR-DrugBank, and PPR-DTC
fig_folder = "figs/PPR_CANCER"
if not os.path.exists(fig_folder):
    os.makedirs(fig_folder)

venn_plot(bp_dic_pagerank_lasso_cancer, bp_dic_pagerank_drugbank_cancer, bp_dic_pagerank_dtc_cancer, 
          ('PPR-Lasso', 'PPR-DrugBank', 'PPR-DTC'), 
          os.path.join(fig_folder, "pathway_enrichment_venn_PPR_Cancer.png"))

# Generate Venn plot for cancer-specific pathways
venn_plot(bp_dic_lasso_cancer, bp_dic_pagerank_lasso_cancer, bp_dic_pagerank_drugbank_cancer, 
          ('Lasso', 'PageRank from Lasso', 'PageRank from DrugBank'), 
          os.path.join(fig_folder, "pathway_enrichment_venn_Lasso_Cancer.png"))

# Load DrugBank pathways into a dictionary
drugbank_pathways_dic = {}
drugbank_files = [f for f in os.listdir("results/pathway_enrichment/drugbank") if f.endswith("_cancer.txt")]

for file in drugbank_files:
    file_path = os.path.join("results/pathway_enrichment/drugbank", file)
    with open(file_path, 'r') as f:
        drug_name = os.path.splitext(file)[0].replace("_cancer", "")
        pathways = [line.strip() for line in f]
        drugbank_pathways_dic[drug_name] = pathways

# Define a helper function to read pathways from a file
def read_pathways_from_file(file_path):
    if os.path.exists(file_path):
        with open(file_path, 'r') as f:
            return set(line.strip() for line in f)
    return set()
# Generate Venn plot using drugbank_pathways_dic
# Read pathways from saved files for PageRank from Lasso and PageRank from DrugBank
pagerank_lasso_pathways = set()
pagerank_lasso_dir = os.path.join(base_results_path, "pagerank_from_lasso")
if os.path.exists(pagerank_lasso_dir):
    for file in os.listdir(pagerank_lasso_dir):
        if file.endswith(".txt") and "_cancer" in file:
            file_path = os.path.join(pagerank_lasso_dir, file)
            pagerank_lasso_pathways.update(read_pathways_from_file(file_path))

pagerank_drugbank_pathways = set()
pagerank_drugbank_dir = os.path.join(base_results_path, "pagerank_from_drugbank")
if os.path.exists(pagerank_drugbank_dir):
    for file in os.listdir(pagerank_drugbank_dir):
        if file.endswith(".txt") and "_cancer" in file:
            file_path = os.path.join(pagerank_drugbank_dir, file)
            pagerank_drugbank_pathways.update(read_pathways_from_file(file_path))

pagerank_dtc_pathways = set()
pagerank_dtc_dir = os.path.join(base_results_path, "pagerank_from_dtc")
if os.path.exists(pagerank_dtc_dir):
    for file in os.listdir(pagerank_dtc_dir):
        if file.endswith(".txt") and "_cancer" in file:
            file_path = os.path.join(pagerank_dtc_dir, file)
            pagerank_dtc_pathways.update(read_pathways_from_file(file_path))

# Generate Venn plot
# Convert drugbank_pathways_dic to a set of pathways
drugbank_pathways_set = set(p for pathways in drugbank_pathways_dic.values() for p in pathways)

print(pagerank_drugbank_pathways)
print(pagerank_lasso_pathways)
print(drugbank_pathways_set)
# Generate Venn plot

# Read pathways from the file data/drugbank_pathways_cancer.txt
drugbank_pathways_set = read_pathways_from_file("data/drugbank_pathways_cancer.txt")
venn_plot({'DrugBank': drugbank_pathways_set}, {'PageRank from Lasso': pagerank_lasso_pathways}, {'PageRank from DrugBank': pagerank_drugbank_pathways}, 
          ('DrugBank', 'PageRank from Lasso', 'PageRank from DrugBank'), 
          os.path.join(fig_folder, "pathway_enrichment_venn_DB_cancer.png"))

        # Read pathways from the file data/dtc_pathways_cancer.txt
dtc_pathways_set = read_pathways_from_file("data/dtc_pathways_cancer.txt")
venn_plot({'DTC': dtc_pathways_set}, {'PageRank from Lasso': pagerank_lasso_pathways}, {'PageRank from DrugBank': pagerank_drugbank_pathways}, 
            ('DTC', 'PageRank from Lasso', 'PageRank from DrugBank'), 
            os.path.join(fig_folder, "pathway_enrichment_venn_DTC_cancer.png"))

# Generate Venn plot for cancer-specific pathways for PPR-DrugBank, PPR-Lasso, and PPR-DTC
fig_folder = "figs/PPR_CANCER"
if not os.path.exists(fig_folder):
    os.makedirs(fig_folder)

# Filter cancer-specific pathways for PPR-DrugBank, PPR-Lasso, and PPR-DTC
bp_dic_pagerank_drugbank_cancer = {k: [p for p in v if 'cancer' in p.lower()] for k, v in bp_dic_pagerank_drugbank.items()}
bp_dic_pagerank_lasso_cancer = {k: [p for p in v if 'cancer' in p.lower()] for k, v in bp_dic_pagerank_lasso.items()}
bp_dic_pagerank_dtc_cancer = {k: [p for p in v if 'cancer' in p.lower()] for k, v in bp_dic_pagerank_dtc.items()}

# Generate Venn plot
venn_plot(bp_dic_pagerank_drugbank_cancer, bp_dic_pagerank_lasso_cancer, bp_dic_pagerank_dtc_cancer, 
            ('PPR-DrugBank', 'PPR-Lasso', 'PPR-DTC'), 
            os.path.join(fig_folder, "pathway_enrichment_venn_PPR_Cancer.png"))
########################################### PATHWAY ENRICHMENT TO COMPARE LASSO, EN AND RIDGE RESULTS ##########################################

bp_dic_en = {}
bp_dic_ridge = {}

for biomarker in drugs:
    # Perform enrichment analysis for Elastic Net (EN) features
   # en_genes = get_feature_list(True, True, biomarker, "en").union(get_feature_list(True, False, biomarker, "en"))  # Combine positive and negative EN genes
    is_GDSC = False
    is_positive = True
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


fig_folder = "figs/PAGERANK_COMPARISON_CANCER"
if not os.path.exists(fig_folder):
    os.makedirs(fig_folder)

# Filter cancer-specific pathways
bp_dic_pagerank_lasso_cancer = {k: [p for p in v if 'cancer' in p.lower()] for k, v in bp_dic_pagerank_lasso.items()}
bp_dic_pagerank_drugbank_cancer = {k: [p for p in v if 'cancer' in p.lower()] for k, v in bp_dic_pagerank_drugbank.items()}
bp_dic_pagerank_dtc_cancer = {k: [p for p in v if 'cancer' in p.lower()] for k, v in bp_dic_pagerank_dtc.items()}

venn_plot(bp_dic_pagerank_lasso_cancer, bp_dic_pagerank_drugbank_cancer, bp_dic_pagerank_dtc_cancer, 
            ('PageRank from Lasso', 'PageRank from DrugBank', 'PageRank from DTC'), 
            os.path.join(fig_folder, "pathway_enrichment_venn_pagerank_comparison_cancer.png"))

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
fig_folder = "figs/CENTRALITY_MEASURES_COMBINATIONS"



# File paths for DrugBank and DTC
drugbank_file_path = "results/pathway_enrichment/drugbank"
dtc_file_path = "results/pathway_enrichment/dtc"

# Centrality measure file paths
betweenness_file_path = "results/pathway_enrichment/centrality_measures/betweenness_cancer.txt"
degree_file_path = "results/pathway_enrichment/centrality_measures/degree_cancer.txt"
eigenvector_file_path = "results/pathway_enrichment/centrality_measures/eigenvector_cancer.txt"
pagerank_file_path = "results/pathway_enrichment/centrality_measures/pagerank_cancer.txt"

# Read pathways from files
# Read and merge pathways from all files ending in _cancer.txt in the DrugBank and DTC directories
drugbank_pathways_cancer = set()
if os.path.exists(drugbank_file_path):
    for file in os.listdir(drugbank_file_path):
        if file.endswith("_cancer.txt"):
            file_path = os.path.join(drugbank_file_path, file)
            drugbank_pathways_cancer.update(read_pathways_from_file(file_path))
else:
    print(f"Directory not found: {drugbank_file_path}")

dtc_pathways_cancer = set()
for file in os.listdir(dtc_file_path):
    if file.endswith("_cancer.txt"):
        file_path = os.path.join(dtc_file_path, file)
        dtc_pathways_cancer.update(read_pathways_from_file(file_path))

betweenness_pathways = read_pathways_from_file(betweenness_file_path)
degree_pathways = read_pathways_from_file(degree_file_path)
eigenvector_pathways = read_pathways_from_file(eigenvector_file_path)
pagerank_pathways = read_pathways_from_file(pagerank_file_path)

# Generate Venn plots
venn_plot({'DrugBank': drugbank_pathways_cancer}, {'DTC': dtc_pathways_cancer}, {'Betweenness': betweenness_pathways}, 
          ('DrugBank', 'DTC', 'Betweenness'), 
          os.path.join(fig_folder, "pathway_enrichment_venn_drugbank_dtc_betweenness_cancer.png"))

venn_plot({'DrugBank': drugbank_pathways_cancer}, {'DTC': dtc_pathways_cancer}, {'Degree': degree_pathways}, 
          ('DrugBank', 'DTC', 'Degree'), 
          os.path.join(fig_folder, "pathway_enrichment_venn_drugbank_dtc_degree_cancer.png"))

venn_plot({'DrugBank': drugbank_pathways_cancer}, {'DTC': dtc_pathways_cancer}, {'Eigenvector': eigenvector_pathways}, 
          ('DrugBank', 'DTC', 'Eigenvector'), 
          os.path.join(fig_folder, "pathway_enrichment_venn_drugbank_dtc_eigenvector_cancer.png"))

venn_plot({'DrugBank': drugbank_pathways_cancer}, {'DTC': dtc_pathways_cancer}, {'Pagerank': pagerank_pathways}, 
          ('DrugBank', 'DTC', 'Pagerank'), 
          os.path.join(fig_folder, "pathway_enrichment_venn_drugbank_dtc_pagerank_cancer.png"))




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
import glob


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

# Perform cancer pathway enrichment for DTC and DrugBank gene sets
dtc_gene_set_path = "data/dtc_gene_set.txt"
drugbank_gene_set_path = "data/drugbank_gene_set.txt"
dtc_output_path = "data/dtc_pathways_cancer.txt"
drugbank_output_path = "data/drugbank_pathways_cancer.txt"

# Helper function to read gene sets from a file
def read_gene_set(file_path):
    if os.path.exists(file_path):
        with open(file_path, 'r') as f:
            return set(line.strip() for line in f)
    return set()

# Read gene sets
dtc_genes = read_gene_set(dtc_gene_set_path)
drugbank_genes = read_gene_set(drugbank_gene_set_path)

# Perform enrichment analysis for DTC gene set
if dtc_genes:
    enriched_dtc_df = perform_enrichment_analysis(dtc_genes, background_genes, reactome_pathways)
    enriched_dtc_df_cancer = enriched_dtc_df[enriched_dtc_df['Pathway Name'].str.contains('cancer', case=False)]
    if not enriched_dtc_df_cancer.empty:
        enriched_dtc_df_cancer['Pathway Name'].to_csv(dtc_output_path, sep='\t', index=False, header=False)
    else:
        print(f"No cancer pathways found for DTC gene set.")

# Perform enrichment analysis for DrugBank gene set
if drugbank_genes:
    enriched_drugbank_df = perform_enrichment_analysis(drugbank_genes, background_genes, reactome_pathways)
    enriched_drugbank_df_cancer = enriched_drugbank_df[enriched_drugbank_df['Pathway Name'].str.contains('cancer', case=False)]
    if not enriched_drugbank_df_cancer.empty:
        enriched_drugbank_df_cancer['Pathway Name'].to_csv(drugbank_output_path, sep='\t', index=False, header=False)
    else:
        print(f"No cancer pathways found for DrugBank gene set.")

# Perform cancer-specific pathway enrichment for PPR-Lasso
ppr_lasso_cancer_pathways = {}
for biomarker in drugs:
    pagerank_from_lasso = get_drugbank_genes(biomarker)
    if not pagerank_from_lasso:
        continue

    enriched_ppr_lasso_df = perform_enrichment_analysis(pagerank_from_lasso, background_genes, reactome_pathways)

    # Save all pathway names
    output_dir_all = os.path.join(base_results_path, "ppr_db_all")
    if not os.path.exists(output_dir_all):
        os.makedirs(output_dir_all)

    if not enriched_ppr_lasso_df.empty:
        enriched_ppr_lasso_df['Pathway Name'].to_csv(os.path.join(output_dir_all, f"{biomarker}_all.txt"), sep='\t', index=False, header=False)

    enriched_ppr_lasso_df['Pathway Name'].to_csv(os.path.join(output_dir_all, f"{biomarker}_all.txt"), sep='\t', index=False)

    # Read cancer-related pathways from the file
    cancer_related_pathways = set()
    cancer_related_pathways_file = "data/cancer_related_pathways.txt"
    if os.path.exists(cancer_related_pathways_file):
        with open(cancer_related_pathways_file, 'r') as f:

            cancer_related_pathways = set(line.strip() for line in f)

    # Filter pathways based on the cancer-related pathways set
    enriched_ppr_lasso_df_cancer = enriched_ppr_lasso_df[enriched_ppr_lasso_df['Pathway Name'].isin(cancer_related_pathways)]

    # Save cancer-specific pathways
    if not enriched_ppr_lasso_df_cancer.empty:
        output_dir_cancer = os.path.join(base_results_path, "ppr_db_cancer")
        if not os.path.exists(output_dir_cancer):
            os.makedirs(output_dir_cancer)

        enriched_ppr_lasso_df_cancer['Pathway Name'].to_csv(os.path.join(output_dir_cancer, f"{biomarker}_cancer.txt"), sep='\t', index=False, header=False)