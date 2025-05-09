import os
import matplotlib.pyplot as plt
from matplotlib_venn import venn3

def plot_venn_from_folder(folder_lasso, folder_en, folder_ridge, output_path, screen, experiment):
    def read_genes_from_folder(folder_path):
        gene_set = set()
        for filename in os.listdir(folder_path):
            file_path = os.path.join(folder_path, filename)
            if os.path.isfile(file_path):
                with open(file_path, 'r') as f:
                    gene_set.update(line.strip() for line in f)
        return gene_set

    # Read genes from all files in each folder
    genes_lasso = read_genes_from_folder(folder_lasso)
    genes_en = read_genes_from_folder(folder_en)
    genes_ridge = read_genes_from_folder(folder_ridge)
    
    # Create Venn diagram
    venn = venn3([genes_lasso, genes_en, genes_ridge], ('Lasso', 'Elastic Net', 'Ridge'))
    
    # Save the plot
    plt.title(f"Venn Diagram of Gene Sets ({screen}, {experiment})")
    plt.savefig(output_path)
    plt.close()

# Example usage:
for screen in ['gCSI', 'GDSC2']:
    for experiment in ['positive', 'negative']:
        folder_lasso = f'results/{screen}/{experiment}/lasso/features/'
        folder_en = f'results/{screen}/{experiment}/en/features/'
        folder_ridge = f'results/{screen}/{experiment}/ridge/features/'
        output_path = f'figs/LASSO_EN_RIDGE/venn_features_{screen}_{experiment}.png'
        os.makedirs(os.path.dirname(output_path), exist_ok=True)
        # Call the function to plot and save the Venn diagram
        plot_venn_from_folder(folder_lasso, folder_en, folder_ridge, output_path, screen, experiment)
