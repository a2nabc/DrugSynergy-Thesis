import os
import pandas as pd

# Define the directories and file patterns
base_dir = "results"
datasets = ["gCSI", "GDSC2"]
experiments = ["negative", "positive"]
methods = ["en", "lasso", "ridge"]
suffix = "_7.txt"
drug_names_file = './data/common_drugs.txt'

# Initialize a list to store the data
data = []

# Process each dataset and experiment
for dataset in datasets:
    for experiment in experiments:
        for method in methods:
            # Construct the directory path
            dir_path = os.path.join(base_dir, dataset, experiment, method, "correlations")
            
            dir_path_enrichment = os.path.join(base_dir, "pathway_enrichment", f"aac_correlations_{method}", dataset, experiment)


            # Collect files ending with the specified suffix
            if os.path.exists(dir_path):  # Ensure the directory exists
                for file_name in os.listdir(dir_path):
                    if file_name.endswith(suffix):
                        file_path = os.path.join(dir_path, file_name)
                        
                        # Read the file and concatenate all lines into a single string
                        with open(file_path, "r") as file:
                            contents = ", ".join(line.strip() for line in file.readlines())
                        
                        # Extract drug type from the file name
                        drug_type = file_name.split("_")[0]

                        # Collect pathway enrichment data, if there exists
                        if os.path.exists(dir_path_enrichment):
                            for file_name_enrichment in os.listdir(dir_path_enrichment):
                                if file_name_enrichment.startswith(drug_type):
                                    # Read the file and concatenate all lines into a single string
                                    file_path_enrichment = os.path.join(dir_path_enrichment, file_name_enrichment)
                                    with open(file_path_enrichment, "r") as file_enrichment:
                                        contents_pathway = ", ".join(line.strip() for line in file_enrichment.readlines())
                        
                        # Append the data to the list
                        data.append({
                            "DRUG": drug_type,
                            "DATASET": dataset,
                            "EXPERIMENT": experiment.capitalize(),
                            "METHOD": method.capitalize(),
                            "TOP AAC GENES": contents,
                            "PATHWAY ENRICHMENT": contents_pathway
                        })

# Convert the data to a DataFrame
df = pd.DataFrame(data)

# Save to CSV
output_file = "./results/highest_correlated_genes_table.csv"
df.to_csv(output_file, index=False)
print(f"Saved {output_file}")
