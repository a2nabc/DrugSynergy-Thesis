import csv

# Paths to the input files
file_paths = [
    'results/DTC/AffectedGenesByDrug.csv',
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
                    genes = columns[2].split(', ')
                    drugbank_genes.update(genes)

# Save the unique gene names for DTC
dtc_output_file_path = 'data/DTC_gene_set.txt'
with open(dtc_output_file_path, 'w') as dtc_outfile:
    dtc_outfile.write('\n'.join(sorted(dtc_genes)))

# Save the unique gene names for drugbank
drugbank_output_file_path = 'data/drugbank_gene_set.txt'
with open(drugbank_output_file_path, 'w') as drugbank_outfile:
    drugbank_outfile.write('\n'.join(sorted(drugbank_genes)))