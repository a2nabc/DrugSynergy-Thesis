import requests
import pandas as pd

# Example drug name or CHEBI ID list
with open('./data/drugs_for_dtc.txt', 'r') as file:
    drugs = [line.strip() for line in file.readlines()]


# Load the drug-target interactions dataset
interactions_df = pd.read_csv('./data/DtcDrugTargetInteractions.csv')

# Filter the dataset for the drugs in the list
filtered_interactions = interactions_df[interactions_df['compound_name'].str.upper().isin([drug.upper() for drug in drugs])]

# Group by compound name and aggregate unique gene names
result = filtered_interactions.groupby('compound_name')['gene_names'].apply(
    lambda x: ', '.join(sorted(set(map(str, x.dropna()))))
).reset_index()

# Save the result to a CSV file or print it
result.to_csv('./results/DTC/AffectedGenesByDrug.csv', index=False)
print(result)