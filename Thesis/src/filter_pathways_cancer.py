import os

# Load the list of cancer-related pathways
cancer_related_pathways_file = "data/cancer_related_pathways.txt"
if os.path.exists(cancer_related_pathways_file):
    with open(cancer_related_pathways_file, "r") as file:
        cancer_related_pathways = set(line.strip() for line in file.readlines())
else:
    print(f"File {cancer_related_pathways_file} does not exist.")
    cancer_related_pathways = set()

# Base directory containing the folders
base_directory = "results/pathway_enrichment"

# Dictionary to store the count of cancer-related pathways
drug_pathway_counts = {}

# Function to process files in a given folder
def process_folder(folder_path, folder_name):
    for file_name in os.listdir(folder_path):
        file_path = os.path.join(folder_path, file_name)
        if os.path.isfile(file_path) and file_name.endswith(".txt"):
            drug_name = file_name.replace(".txt", "")
            with open(file_path, "r") as file:
                pathways = [line.strip() for line in file.readlines()]
                cancer_pathways = [p for p in pathways if p in cancer_related_pathways]
            
            # Save the filtered pathways to a new file
            output_cancer_file = os.path.join(folder_path, f"{drug_name}_cancer.txt")
            with open(output_cancer_file, "w") as output_file:
                output_file.write("\n".join(cancer_pathways))
            
            # Update the count for cancer-related pathways
            if drug_name not in drug_pathway_counts:
                drug_pathway_counts[drug_name] = {}
            drug_pathway_counts[drug_name][folder_name] = len(cancer_pathways)

# Iterate through all subfolders in the base directory
for folder_name in os.listdir(base_directory):
    folder_path = os.path.join(base_directory, folder_name)
    if os.path.isdir(folder_path):  # Ensure it's a folder
        # Check if the folder contains subfolders like "Negative" and "Positive"
        subfolders = [f for f in os.listdir(folder_path) if os.path.isdir(os.path.join(folder_path, f))]
        if "Negative" in subfolders or "Positive" in subfolders:
            for subfolder_name in subfolders:
                subfolder_path = os.path.join(folder_path, subfolder_name)
                process_folder(subfolder_path, f"{folder_name}/{subfolder_name}")
        else:
            # Process files directly in the folder
            process_folder(folder_path, folder_name)
