import os
import pandas as pd

# Function to empty a file if it exists
def empty_file_if_exists(file_path):
    if os.path.exists(file_path):
        with open(file_path, 'w') as file:
            pass  # Truncate the file

# Define the folders and their corresponding column names
folders_mapping = {
    1: {
        "folders": {
            "results/pathway_enrichment/pagerank_from_drugbank": "COUNT_PR_DB",
            "results/pathway_enrichment/ppr_lasso_all": "COUNT_PR_LASSO",
            "results/pathway_enrichment/pagerank_from_dtc": "COUNT_PR_DTC",
        },
        "output_file": "results/pathway_enrichment/drug_pathway_counts_PPR.csv",
        "output_file_cancer": "results/pathway_enrichment/drug_pathway_counts_PPR_cancer.csv",
    },
    2: {
        "folders": {
            "results/pathway_enrichment/drugbank": "COUNT_DB",
            "results/pathway_enrichment/dtc": "COUNT_DTC",
            "results/pathway_enrichment/lasso_features": "COUNT_LASSO",
        },
        "output_file": "results/pathway_enrichment/drug_pathway_counts_DTC_DB_LASSO.csv",
        "output_file_cancer": "results/pathway_enrichment/drug_pathway_counts_DTC_DB_LASSO_cancer.csv",
    },
}

# Variable to determine which folders to use
experiment = 1  # Change this to 1 or 2
folders = folders_mapping[experiment]["folders"]
output_file = folders_mapping[experiment]["output_file"]
output_file_cancer = folders_mapping[experiment]["output_file_cancer"]

# Function to process pathway counts
# Function to process pathway counts with aggregation
def process_pathway_counts(folders, cancer_suffix=""):
    drug_pathway_counts = {}

    # Iterate through each folder
    for folder, column_name in folders.items():
        if os.path.exists(folder):
            for file_name in os.listdir(folder):
                if file_name.endswith(f"{cancer_suffix}.txt"):
                    drug_name = file_name.replace(f"{cancer_suffix}.txt", "")
                    file_path = os.path.join(folder, file_name)
                    with open(file_path, "r") as file:
                        pathways = file.readlines()
                        if drug_name not in drug_pathway_counts:
                            drug_pathway_counts[drug_name] = {column: 0 for column in folders.values()}
                        drug_pathway_counts[drug_name][column_name] += len(pathways)
        else:
            print(f"Folder {folder} does not exist.")

    # Ensure all drugs are accounted for, even if they have 0 pathways
    all_drugs = set(drug_pathway_counts.keys())
    for folder in folders.keys():
        if os.path.exists(folder):
            for file_name in os.listdir(folder):
                if file_name.endswith(f"{cancer_suffix}.txt"):
                    drug_name = file_name.replace(f"{cancer_suffix}.txt", "")
                    all_drugs.add(drug_name)

    # Add drugs with 0 pathways for missing entries
    for drug in all_drugs:
        if drug not in drug_pathway_counts:
            drug_pathway_counts[drug] = {column: 0 for column in folders.values()}

    # Create a DataFrame for the table
    df = pd.DataFrame.from_dict(drug_pathway_counts, orient="index").reset_index()
    df.rename(columns={"index": "DRUG"}, inplace=True)
    return df

# Process non-cancer and cancer-related pathways
df = process_pathway_counts(folders)
df = df[~df["DRUG"].str.endswith("_cancer")]
df_cancer = process_pathway_counts(folders, cancer_suffix="_cancer")




# Load the mapping file
mapping_file = "data/drugs_for_dtc_match.csv"
if os.path.exists(mapping_file):
    drug_mapping = pd.read_csv(mapping_file)
    dtc_to_drug_mapping = dict(zip(drug_mapping["name_dtc"], drug_mapping["name_drug"]))

    # Function to standardize drug names using the mapping
    def standardize_drug_name(drug_name):
        return dtc_to_drug_mapping.get(drug_name, drug_name)

    # Apply the mapping to standardize drug names in the DataFrames
    df["DRUG"] = df["DRUG"].apply(standardize_drug_name)
    df_cancer["DRUG"] = df_cancer["DRUG"].apply(standardize_drug_name)

# Save the DataFrames to CSV files
empty_file_if_exists(output_file)
df.to_csv(output_file, index=False)
print(f"Table saved to {output_file}")

empty_file_if_exists(output_file_cancer)
df_cancer.to_csv(output_file_cancer, index=False)
print(f"Table filtering cancer saved to {output_file_cancer}")

# Print the resulting DataFrames
print("Non-cancer pathway counts:")
print(df)

print("\nCancer pathway counts:")
print(df_cancer)

