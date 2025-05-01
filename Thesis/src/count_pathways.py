import os
import pandas as pd

# Define the folders and their corresponding column names
folders = {
    "results/pathway_enrichment/drugbank": "COUNT_DB",
    "results/pathway_enrichment/dtc": "COUNT_DTC",
    "results/pathway_enrichment/lasso_features": "COUNT_LASSO"
}

# Initialize a dictionary to store drug pathway counts
drug_pathway_counts = {}

# Iterate through each folder
for folder, column_name in folders.items():
    if os.path.exists(folder):
        for file_name in os.listdir(folder):
            if file_name.endswith(".txt"):
                drug_name = file_name.replace(".txt", "")
                file_path = os.path.join(folder, file_name)
                with open(file_path, "r") as file:
                    pathways = file.readlines()
                    if drug_name not in drug_pathway_counts:
                        drug_pathway_counts[drug_name] = {"COUNT_DTC": 0, "COUNT_DB": 0, "COUNT_LASSO": 0}
                    drug_pathway_counts[drug_name][column_name] = len(pathways)
    else:
        print(f"Folder {folder} does not exist.")

# Ensure all drugs are accounted for, even if they have 0 pathways
all_drugs = set(drug_pathway_counts.keys())
for folder in folders.keys():
    if os.path.exists(folder):
        for file_name in os.listdir(folder):
            if file_name.endswith(".txt"):
                drug_name = file_name.replace(".txt", "")
                all_drugs.add(drug_name)

# Add drugs with 0 pathways for missing entries
for drug in all_drugs:
    if drug not in drug_pathway_counts:
        drug_pathway_counts[drug] = {"COUNT_DTC": 0, "COUNT_DB": 0, "COUNT_LASSO": 0}

# Create a DataFrame for the table
df = pd.DataFrame.from_dict(drug_pathway_counts, orient="index").reset_index()
df.rename(columns={"index": "DRUG"}, inplace=True)

# Save the table to a CSV file
output_file = "results/pathway_enrichment/drug_pathway_counts_DTC_DB_LASSO.csv"
df.to_csv(output_file, index=False)

print(f"Table saved to {output_file}")