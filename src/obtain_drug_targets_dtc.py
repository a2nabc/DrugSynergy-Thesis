import requests
import pandas as pd

# Example drug name or CHEBI ID list
with open('./data/drugs_for_dtc.txt', 'r') as file:
    drugs = [line.strip() for line in file.readlines()]


# Load the drug-target interactions dataset
interactions_df = pd.read_csv('./data/DtcDrugTargetInteractions.csv')

# Check for rows where compound_name contains "DOCETAXEL"
docetaxel_rows = interactions_df[interactions_df['compound_name'].str.contains('DOCETAXEL', case=False, na=False)]

# Print the column names of the filtered rows
print(docetaxel_rows['compound_name'].unique())


# Filter the dataset for the drugs in the list
filtered_interactions = interactions_df[interactions_df['compound_name'].str.upper().isin([drug.upper() for drug in drugs])]

# Group by compound name and aggregate unique gene names
result = filtered_interactions.groupby('compound_name')['gene_names'].apply(
    lambda x: ', '.join(sorted(set(map(str, x.dropna()))))
).reset_index()

# Save the result to a CSV file or print it
result.to_csv('./results/dtc/AffectedGenesByDrug.csv', index=False)
print(result)

# Function to output a table with drug and count of unique gene_names
def save_drug_gene_count_table(interactions_df, output_path):
    # Group by compound name and count unique gene names
    gene_count = interactions_df.groupby('compound_name')['gene_names'].apply(
        lambda x: len(set(map(str, x.dropna())))
    ).reset_index()

    # Rename columns for clarity
    gene_count.columns = ['compound_name', 'unique_gene_count']

    # Save the result to a CSV file
    gene_count.to_csv(output_path, index=False)
    print(gene_count)
    print(f"Drug-gene count table saved to {output_path}")

# Call the function and save the table
save_drug_gene_count_table(filtered_interactions, './results/dtc/DrugGeneCount.csv')