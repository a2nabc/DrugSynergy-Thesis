import os
import pandas as pd

# Define the folder paths
folder_paths = [
    "results/pathway_enrichment/lasso_features",
    "results/pathway_enrichment/en_features",
    "results/pathway_enrichment/ridge_features",
    "results/pathway_enrichment/pagerank_from_lasso",
    "results/pathway_enrichment/pagerank_from_drugbank"
]

for folder_path in folder_paths:
    # Initialize a dictionary to store the data
    data = {}

    # Loop through all files in the folder
    for file_name in os.listdir(folder_path):
        file_path = os.path.join(folder_path, file_name)
        if os.path.isfile(file_path) and file_name.endswith('.txt'):  # Ensure it's a .txt file
            # Extract the drug name from the file name
            name_drug = os.path.splitext(file_name)[0]

            # Read the pathways from the file
            with open(file_path, 'r') as file:
                pathways = file.read().strip().splitlines()

            # Combine pathways for the same drug
            if name_drug in data:
                data[name_drug].extend(pathways)
            else:
                data[name_drug] = pathways

    # Prepare the data for the DataFrame
    formatted_data = [{"NAME_DRUG": drug, "PATHWAY": ", ".join(set(pathways))} for drug, pathways in data.items()]

    # Create a DataFrame from the collected data
    df = pd.DataFrame(formatted_data)
    
    # Save the table to a CSV file inside the folder
    output_file = os.path.join(folder_path, "pathway_enrichment_table.csv")
    
    # Ensure the file is rewritten
    if os.path.exists(output_file):
        os.remove(output_file)
    
    df.to_csv(output_file, index=False)

    print(f"Table created and saved to {output_file}")
