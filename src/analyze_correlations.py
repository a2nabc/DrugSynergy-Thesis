import os
import pandas as pd

drug_names_file = './data/common_drugs.txt'

def extract_correlations(base_dir, condition, subtype, threshold=0.7):
    input_dir = os.path.join(base_dir, condition, subtype, "correlations")
    files = [f for f in os.listdir(input_dir) if f.endswith('.csv')]

    for file in files:
        file_path = os.path.join(input_dir, file)
        data = pd.read_csv(file_path)
        
        # Assuming the correlations are in a column named 'correlation'
        if 'x' not in data.columns:
            raise ValueError(f"'x' column not found in {file_path}")
        
        # Extract correlations > threshold
        filtered_data = data[data['x'] > threshold]
        output_file = os.path.join(input_dir, f"{file.split('.')[0]}_7.txt")
        filtered_data = filtered_data.drop(columns=['x'])
        filtered_data.iloc[:, 0].to_csv(output_file, index=False, header=False)

# Define directories, conditions, and subtypes
base_dir = './results'
conditions = ['gCSI/negative', 'gCSI/positive', 'GDSC2/negative', 'GDSC2/positive']
subtypes = ['en', 'lasso', 'ridge']

for condition in conditions:
    for subtype in subtypes:
        extract_correlations(base_dir, condition, subtype)
