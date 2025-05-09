from itertools import combinations
import os
from matplotlib_venn import venn3

def jaccard_distance(set1, set2):
    """
    Compute the Jaccard distance between two sets.

    Jaccard distance is defined as 1 - Jaccard index, where
    Jaccard index = |Intersection(set1, set2)| / |Union(set1, set2)|

    :param set1: First set
    :param set2: Second set
    :return: Jaccard distance (float)
    """
    intersection = len(set1.intersection(set2))
    union = len(set1.union(set2))
    if union == 0:
        return 0.0  # Handle edge case where both sets are empty
    return 1 - (intersection / union)

def pairwise_jaccard_distances(sets):
    """
    Compute the Jaccard distances for all pairwise combinations of sets.

    :param sets: List of sets
    :return: List of tuples (set1, set2, jaccard_distance)
    """
    distances = []
    for set1, set2 in combinations(sets, 2):
        distance = jaccard_distance(set1, set2)
        distances.append((set1, set2, distance))
    return distances


def load_set_from_folder(folder_path, suffix="_cancer.txt"):
    """
    Load a set by taking the union of all elements in files with a specific suffix in a folder.

    :param folder_path: Path to the folder containing the files
    :param suffix: Suffix of the files to include
    :return: A set containing the union of all elements in the files
    """
    result_set = set()
    for filename in os.listdir(folder_path):
        if filename.endswith(suffix):
            file_path = os.path.join(folder_path, filename)
            with open(file_path, 'r') as file:
                for line in file:
                    result_set.add(line.strip())
    return result_set

def read_lines_from_file(file_path):
    """
    Read all lines from a file and return them as a set of strings.

    :param file_path: Path to the file
    :return: Set of lines (stripped of whitespace)
    """
    with open(file_path, 'r') as file:
        return {line.strip() for line in file}


paths_folder2 = {
    "lasso_features": "results/pathway_enrichment/lasso_features",
    "drugbank": "results/pathway_enrichment/drugbank",
    "pagerank_from_lasso": "results/pathway_enrichment/pagerank_from_lasso",
    "pagerank_from_drugbank": "results/pathway_enrichment/pagerank_from_drugbank",
}
paths_file2 = {
    "pagerank": "results/pathway_enrichment/centrality_measures/pagerank_cancer.txt"
}
# Define which experiment to run
experiment_to_run = 0

if experiment_to_run == 1:
    # Experiment 1
    paths_folder = {
        "lasso_features": "results/GDSC2/positive/lasso/features"
    }
    paths_file = {
        "DT_DRUGBANK": "data/drugbank_gene_set.txt",
        "DT_DTC": "data/dtc_gene_set.txt"
    }
elif experiment_to_run == 0:
    # Experiment 0
    paths_folder = {
        "lasso_features": "results/GDSC2/positive/lasso/features",
        "en_features": "results/GDSC2/positive/en/features",
        "ridge_features": "results/GDSC2/positive/ridge/features"
    }
    paths_file = {
     
    }
elif experiment_to_run == 2:
    # Experiment 2
    paths_folder = {
        "lasso_features": "results/pathway_enrichment/lasso_features",
        "drugbank": "results/pathway_enrichment/drugbank",
        "pagerank_from_lasso": "results/pathway_enrichment/pagerank_from_lasso",
        "pagerank_from_drugbank": "results/pathway_enrichment/pagerank_from_drugbank",
    }
    paths_file = {
        "pagerank": "results/pathway_enrichment/centrality_measures/pagerank_cancer.txt"
    }
elif experiment_to_run == 3:
    # Experiment 3
    paths_folder = {
        "lasso_features": "results/pathway_enrichment/lasso_features",
        "drugbank": "results/pathway_enrichment/drugbank",
        "dtc": "results/pathway_enrichment/dtc"
    }
    paths_file = {
        
    }
elif experiment_to_run == 4:
    # Experiment 4
    paths_folder = {

    }
    paths_file = {
        "pagerank": "results/pathway_enrichment/centrality_measures/pagerank_cancer.txt",
        "betwenness": "results/pathway_enrichment/centrality_measures/betweenness_cancer.txt",
        "degree": "results/pathway_enrichment/centrality_measures/degree_cancer.txt",
        "eigenvector": "results/pathway_enrichment/centrality_measures/eigenvector_cancer.txt"
    }
elif experiment_to_run == 5:
    # Experiment 5
    paths_folder = {

    }
    paths_file = {
        "pagerank": "results/centrality_measures/pagerank_top_50.txt",
        "betwenness": "results/centrality_measures/betweenness_top_50.txt",
        "degree": "results/centrality_measures/degree_top_50.txt",
        "eigenvector": "results/centrality_measures/eigenvector_top_50.txt"
    }

elif experiment_to_run == 6:
    # Experiment 6
    paths_folder = {
        "drugbank": "results/pathway_enrichment/drugbank",
        "pagerank_from_lasso": "results/pathway_enrichment/pagerank_from_lasso",
        "pagerank_from_drugbank": "results/pathway_enrichment/pagerank_from_drugbank",
    }
    paths_file = {
    }

elif experiment_to_run == 7:
    # Experiment 7
    paths_folder = {
        "pagerank_from_lasso": "results/pagerank_output/lasso/GDSC2/positive",
        "pagerank_from_drugbank": "results/pagerank_output/drugbank/drugwise",
        "pagerank_from_dtc": "results/pagerank_output/dtc/drugwise"
    }
    paths_file = {
        
    }
elif experiment_to_run == 8:
    # Experiment 8
    paths_folder = {
        
    }
    paths_file = {
        "betwenness": "results/centrality_measures/betweenness_top_all.txt",
        "degree": "results/centrality_measures/degree_top_all.txt",
        "eigenvector": "results/centrality_measures/eigenvector_top_all.txt",
        "pagerank": "results/centrality_measures/pagerank_top_all.txt"

    }
else:
    raise ValueError("Invalid experiment_to_run value. Choose 1 or 2.")

import matplotlib.pyplot as plt

def plot_venn_diagram(sets, set_labels):
    """
    Plot a Venn diagram for three sets.

    :param sets: List of three sets
    :param set_labels: List of labels for the sets
    """
    if len(sets) != 3 or len(set_labels) != 3:
        raise ValueError("This function only supports exactly three sets and three labels.")
    
    venn3([sets[0], sets[1], sets[2]], set_labels)
    plt.title("Overlap of feature sets")
    plt.show()

if experiment_to_run == 1 or experiment_to_run == 0:
    folder_sets = [load_set_from_folder(folder_path, suffix=".txt") for folder_path in paths_folder.values()]
elif experiment_to_run == 7:
    folder_sets = [load_set_from_folder(folder_path, suffix="50.txt") for folder_path in paths_folder.values()]
else:
    folder_sets = [load_set_from_folder(folder_path) for folder_path in paths_folder.values()]
file_sets = [read_lines_from_file(file_path) for file_path in paths_file.values()]

all_sets = folder_sets + file_sets

# Print the number of elements in each set
set_names = list(paths_folder.keys()) + list(paths_file.keys())
for name, s in zip(set_names, all_sets):
    print(f"Set '{name}' contains {len(s)} elements.")

jaccard_results = pairwise_jaccard_distances(all_sets)

# Print the results with set names instead of the sets themselves
for i, (set1, set2, distance) in enumerate(jaccard_results):
    name1 = set_names[all_sets.index(set1)]
    name2 = set_names[all_sets.index(set2)]
    print(f"Jaccard distance between {name1} and {name2}: {distance}")

if experiment_to_run == 0:
    # Compute Jaccard distance file-wise for matching names across folders
    folder_files_dicts = []
    folder_names = list(paths_folder.keys())
    for folder_path in paths_folder.values():
        folder_files = {filename: read_lines_from_file(os.path.join(folder_path, filename))
                        for filename in os.listdir(folder_path) if filename.endswith(".txt")}
        folder_files_dicts.append(folder_files)
    
    # Ensure there are at least two folders to compare
    if len(folder_files_dicts) < 2:
        raise ValueError("At least two folders are required to compute distances.")

    # Compute pairwise Jaccard distances for matching file names across folders
    for filename in folder_files_dicts[0].keys():
        if all(filename in folder_files for folder_files in folder_files_dicts[1:]):  # Check if file exists in all folders
            set1 = folder_files_dicts[0][filename]
            set2 = folder_files_dicts[1][filename]
            print(f"Jaccard distance between {filename} ({folder_names[0]}) and {filename} ({folder_names[1]}): {jaccard_distance(set1, set2)}")

# Example: Plot a Venn diagram for the first three sets
if len(all_sets) >= 3:
    plot_venn_diagram(all_sets[:3], set_names[:3])









from scipy.stats import wilcoxon
from statsmodels.stats.multitest import multipletests
import pandas as pd

# Input Pearson values
# Load the CSV file and focus on the PEARSON column

gdsc_pos_df_ridge = pd.read_csv("results/gCSI/positive/summary_results.csv")
gdsc_pos_ridge = gdsc_pos_df_ridge[gdsc_pos_df_ridge["Model"] == "ridge"]["PEARSON"].tolist()

gdsc_pos_df_lasso = pd.read_csv("results/gCSI/positive/summary_results.csv")
gdsc_pos_lasso = gdsc_pos_df_lasso[gdsc_pos_df_lasso["Model"] == "lasso"]["PEARSON"].tolist()

gdsc_pos_df_en = pd.read_csv("results/gCSI/positive/summary_results.csv")
gdsc_pos_en = gdsc_pos_df_en[gdsc_pos_df_en["Model"] == "en"]["PEARSON"].tolist()


gdsc_pos_df = pd.read_csv("results/GDSC2/positive/summary_results.csv")
gdsc_pos = gdsc_pos_df[gdsc_pos_df["Model"] == "en"]["PEARSON"].tolist()

gdsc_neg_df = pd.read_csv("results/GDSC2/negative/summary_results.csv")
gdsc_neg = gdsc_neg_df[gdsc_neg_df["Model"] == "en"]["PEARSON"].tolist()

gcsi_pos_df = pd.read_csv("results/gCSI/positive/summary_results.csv")
gcsi_pos = gcsi_pos_df[gcsi_pos_df["Model"] == "en"]["PEARSON"].tolist()

gcsi_neg_df = pd.read_csv("results/gCSI/negative/summary_results.csv")
gcsi_neg = gcsi_neg_df[gcsi_neg_df["Model"] == "en"]["PEARSON"].tolist()



# Wilcoxon tests
comparisons = {
    "gCSI_Pos_vs_Neg": wilcoxon(gcsi_pos, gcsi_neg),
    "GDSC_Pos_vs_Neg": wilcoxon(gdsc_pos, gdsc_neg),
    "gCSI_Pos_vs_GDSC_Pos": wilcoxon(gcsi_pos, gdsc_pos),
    "gCSI_Neg_vs_GDSC_Neg": wilcoxon(gcsi_neg, gdsc_neg),
    "lasso_en": wilcoxon(gdsc_pos_lasso, gdsc_pos_en),
    "lasso_ridge": wilcoxon(gdsc_pos_lasso, gdsc_pos_ridge),
    "ridge_en": wilcoxon(gdsc_pos_ridge, gdsc_pos_en),
}

# Extract p-values
p_values = [test.pvalue for test in comparisons.values()]
q_values = multipletests(p_values, method='fdr_bh')[1]
significance = q_values < 0.05

# Results table
results_df = pd.DataFrame({
    "Comparison": list(comparisons.keys()),
    "p-value": p_values,
    "q-value": q_values,
    "Significant (FDR < 0.05)": significance
})

print(results_df)
