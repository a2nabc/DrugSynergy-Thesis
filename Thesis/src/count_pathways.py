import os
import csv

# Define the paths to the files
paths = [os.path.join("results/pathway_enrichment/centrality_measures", file) 
         for file in os.listdir("results/pathway_enrichment/centrality_measures") 
         if file.endswith(".txt")]

# Initialize a dictionary to store the counts
pathway_counts = {}

# Count the number of pathways in each file
for path in paths:
    if os.path.exists(path):
        with open(path, 'r') as file:
            lines = file.readlines()
            pathway_counts[path] = len(lines)
    else:
        pathway_counts[path] = "File not found"

# Print the table
print(f"{'File':<10} | {'Pathway Count':<15}")
print("-" * 30)
for path, count in pathway_counts.items():
    print(f"{path:<10} | {count:<15}")

    # Write the pathway counts to a CSV file
output_file = "results/pathway_enrichment/pathway_counts_CENTRALITY_MEASURES.csv"

with open(output_file, 'w', newline='') as csvfile:
    csvwriter = csv.writer(csvfile)
    csvwriter.writerow(["File", "Pathway Count"])
    for path, count in pathway_counts.items():
        csvwriter.writerow([path, count])