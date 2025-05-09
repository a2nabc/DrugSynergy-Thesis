import os
import pandas as pd

# Function to empty a file if it exists
def empty_file_if_exists(file_path):
    if os.path.exists(file_path):
        with open(file_path, 'w') as file:
            pass  # Truncate the file

# Define the folders and their corresponding column names
folders2 = {
    "results/pathway_enrichment/drugbank": "COUNT_DB",
    "results/pathway_enrichment/dtc": "COUNT_DTC",
    "results/pathway_enrichment/lasso_features": "COUNT_LASSO"
}

folders = {
    "results/pathway_enrichment/pagerank_from_drugbank": "COUNT_PR_DB",
    "results/pathway_enrichment/pagerank_from_lasso": "COUNT_PR_LASSO",
}

# Initialize a dictionary to store drug pathway counts
drug_pathway_counts = {}

# Iterate through each folder
for folder, column_name in folders.items():
    if os.path.exists(folder):
        for file_name in os.listdir(folder):
            if file_name.endswith(".txt") and not file_name.endswith("_cancer.txt"):
                drug_name = file_name.replace(".txt", "")
                file_path = os.path.join(folder, file_name)
                with open(file_path, "r") as file:
                    pathways = file.readlines()
                    if drug_name not in drug_pathway_counts:
                        drug_pathway_counts[drug_name] = {column: 0 for column in folders.values()}
                    drug_pathway_counts[drug_name][column_name] = len(pathways)
    else:
        print(f"Folder {folder} does not exist.")

# Ensure all drugs are accounted for, even if they have 0 pathways
all_drugs = set(drug_pathway_counts.keys())
for folder in folders.keys():
    if os.path.exists(folder):
        for file_name in os.listdir(folder):
            if file_name.endswith(".txt") and not file_name.endswith("_cancer.txt"):
                drug_name = file_name.replace(".txt", "")
                all_drugs.add(drug_name)

# Add drugs with 0 pathways for missing entries
for drug in all_drugs:
    if drug not in drug_pathway_counts:
        drug_pathway_counts[drug] = {column: -1 for column in folders.values()}

# Create a DataFrame for the table
df = pd.DataFrame.from_dict(drug_pathway_counts, orient="index").reset_index()
df.rename(columns={"index": "DRUG"}, inplace=True)

# Save the table to a CSV file
output_file2 = "results/pathway_enrichment/drug_pathway_counts_DTC_DB_LASSO.csv"
output_file = "results/pathway_enrichment/drug_pathway_counts_PPR.csv"

# Empty the file if it exists
empty_file_if_exists(output_file)
df.to_csv(output_file, index=False)

print(f"Table saved to {output_file}")


################## CENCER RELATED PATHWAYS ##################
# Initialize a dictionary to store drug pathway counts
drug_pathway_counts_cancer = {}

# Iterate through each folder
for folder, column_name in folders.items():
    if os.path.exists(folder):
        for file_name in os.listdir(folder):
            if file_name.endswith("_cancer.txt"):
                drug_name = file_name.replace("_cancer.txt", "")
                file_path = os.path.join(folder, file_name)
                with open(file_path, "r") as file:
                    pathways_cancer = file.readlines()
                    if drug_name not in drug_pathway_counts_cancer:
                        drug_pathway_counts_cancer[drug_name] = {column: 0 for column in folders.values()}
                    drug_pathway_counts_cancer[drug_name][column_name] = len(pathways_cancer)
    else:
        print(f"Folder {folder} does not exist.")

# Ensure all drugs are accounted for, even if they have 0 pathways
all_drugs_cancer = set(drug_pathway_counts_cancer.keys())
for folder in folders.keys():
    if os.path.exists(folder):
        for file_name in os.listdir(folder):
            if file_name.endswith("_cancer.txt"):
                drug_name = file_name.replace("_cancer.txt", "")
                all_drugs_cancer.add(drug_name)

# Add drugs with 0 pathways for missing entries
for drug in all_drugs_cancer:
    if drug not in drug_pathway_counts_cancer:
        drug_pathway_counts_cancer[drug] = {column: -1 for column in folders.values()}

# Create a DataFrame for the table
df_cancer = pd.DataFrame.from_dict(drug_pathway_counts_cancer, orient="index").reset_index()
df_cancer.rename(columns={"index": "DRUG"}, inplace=True)

# Save the table to a CSV file
output_file2_cancer = "results/pathway_enrichment/drug_pathway_counts_DTC_DB_LASSO_cancer.csv"
output_file_cancer = "results/pathway_enrichment/drug_pathway_counts_PPR_cancer.csv"

# Empty the file if it exists
empty_file_if_exists(output_file_cancer)
df_cancer.to_csv(output_file_cancer, index=False)

print(f"Table filtering cancer saved to {output_file_cancer}")