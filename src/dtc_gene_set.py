import csv
import os

# Paths to the input files
file_paths = [
    'results/dtc/AffectedGenesByDrug.csv',
    'results/drugbank/AffectedGenesByDrug.txt'
]

# Separate unique genes into two sets for DTC and drugbank
dtc_genes = set()
drugbank_genes = set()

# Process each file again to separate genes
for file_path in file_paths:
    if file_path.endswith('.csv'):
        with open(file_path, 'r') as csvfile:
            reader = csv.DictReader(csvfile)
            for row in reader:
                genes = row['gene_names'].split(', ')
                dtc_genes.update(genes)
    elif file_path.endswith('.txt'):
        with open(file_path, 'r') as txtfile:
            next(txtfile)  # Skip the header line
            for line in txtfile:
                columns = line.strip().split('\t')
                if len(columns) >= 3:  # Ensure there are at least 3 columns
                    genes = columns[2].split(', ')  # 'GENES_AFFECTED' is the 3rd column
                    drugbank_genes.update(genes)

# Save the unique gene names for DTC
dtc_output_file_path = 'data/DTC_gene_set.txt'
with open(dtc_output_file_path, 'w') as dtc_outfile:
    dtc_outfile.write('\n'.join(sorted(dtc_genes)))

# Save the unique gene names for drugbank
drugbank_output_file_path = 'data/drugbank_gene_set.txt'
with open(drugbank_output_file_path, 'w') as drugbank_outfile:
    drugbank_outfile.write('\n'.join(sorted(drugbank_genes)))

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


# Load the drug name correspondences
drug_name_map = {}
with open('data/drugs_for_dtc_match.csv', 'r') as mapping_file:
    reader = csv.DictReader(mapping_file)
    for row in reader:
        # Ensure the correct column names are used
        if 'name_dtc' in row and 'name_drug' in row:
            drug_name_map[row['name_dtc']] = row.get('name_drug', '-')
        else:
            raise KeyError("Expected columns 'name_dtc' and 'name_drug' not found in the CSV file.")

# Load affected genes by drug for DTC
affected_genes_by_dtc_drug = {}
with open('results/dtc/AffectedGenesByDrug.csv', 'r') as dtc_file:
    reader = csv.DictReader(dtc_file)
    for row in reader:
        affected_genes_by_dtc_drug[row['compound_name']] = set(row['gene_names'].split(', '))

# Load affected genes by drug for DrugBank
affected_genes_by_drugbank_drug = {}
with open('results/drugbank/AffectedGenesByDrug.txt', 'r') as drugbank_file:
    next(drugbank_file)  # Skip the header line
    for line in drugbank_file:
        columns = line.strip().split('\t')
        if len(columns) >= 3:  # Ensure there are at least 3 columns
            drug_name = columns[1]  # 'DRUG_NAME' is the 2nd column
            genes = set(columns[2].split(', '))  # 'GENES_AFFECTED' is the 3rd column
            affected_genes_by_drugbank_drug[drug_name] = genes

# Compute Jaccard distances for each drug
jaccard_results = {}
for dtc_drug, drugbank_drug in drug_name_map.items():
    dtc_drug_genes = affected_genes_by_dtc_drug.get(dtc_drug, set())
    drugbank_drug_genes = affected_genes_by_drugbank_drug.get(drugbank_drug, set())
    distance = jaccard_distance(dtc_drug_genes, drugbank_drug_genes)
    jaccard_results[dtc_drug] = distance

# Save the Jaccard distances to a file
jaccard_output_file_path = 'results/jaccard/distances_dtc_db.txt'
with open(jaccard_output_file_path, 'w') as jaccard_outfile:
    jaccard_outfile.write('Drug,Jaccard_Distance\n')
    for drug, distance in sorted(jaccard_results.items()):
        jaccard_outfile.write(f'{drug},{distance:.4f}\n')

# Load affected genes for GDSC2 from multiple files
gdsc2_genes_by_drug = {}
gdsc2_folder_path = 'results/GDSC2/positive/lasso/features'
for file_name in os.listdir(gdsc2_folder_path):
    if file_name.endswith('.txt'):
        drug_name = os.path.splitext(file_name)[0]  # Extract drug name from file name
        file_path = os.path.join(gdsc2_folder_path, file_name)
        with open(file_path, 'r') as gdsc2_file:
            genes = set(line.strip() for line in gdsc2_file)
            gdsc2_genes_by_drug[drug_name] = genes

# Compute Jaccard distances between DTC and GDSC2 for each drug using the drug name mapping
dtc_gdsc2_jaccard_results = {}
for dtc_drug, dtc_drug_genes in affected_genes_by_dtc_drug.items():
    gdsc2_drug = drug_name_map.get(dtc_drug, None)  # Map DTC drug name to GDSC2 drug name
    if gdsc2_drug:
        gdsc2_drug_genes = gdsc2_genes_by_drug.get(gdsc2_drug, set())
        distance = jaccard_distance(dtc_drug_genes, gdsc2_drug_genes)
        dtc_gdsc2_jaccard_results[dtc_drug] = distance

# Save the Jaccard distances to a file
gdsc2_jaccard_output_file_path = 'results/jaccard/distances_dtc_gdsc2.txt'
with open(gdsc2_jaccard_output_file_path, 'w') as gdsc2_outfile:
    gdsc2_outfile.write('Drug,Jaccard_Distance\n')
    for drug, distance in sorted(dtc_gdsc2_jaccard_results.items()):
        gdsc2_outfile.write(f'{drug},{distance:.4f}\n')


# Create output directories if they don't exist
os.makedirs('results/dtc/features', exist_ok=True)
os.makedirs('results/drugbank/features', exist_ok=True)

# Generate files for DTC
for dtc_drug, genes in affected_genes_by_dtc_drug.items():
    standardized_name = drug_name_map.get(dtc_drug, dtc_drug)  # Use standardized name if available
    output_file_path = os.path.join('results/dtc/features', f'{standardized_name}.txt')
    with open(output_file_path, 'w') as outfile:
        outfile.write('\n'.join(sorted(genes)) + '\n')  # Add final newline

# Generate files for DrugBank
for drugbank_drug, genes in affected_genes_by_drugbank_drug.items():
    standardized_name = drugbank_drug  # DrugBank names are already standardized
    output_file_path = os.path.join('results/drugbank/features', f'{standardized_name}.txt')
    with open(output_file_path, 'w') as outfile:
        outfile.write('\n'.join(sorted(genes)) + '\n')  # Add final newline