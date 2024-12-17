import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import os
import plotly.express as px # library for interactive plots

# Load and filter data

####################################### MONOTHERAPY DATA #######################################

data = pd.read_excel('Data/Monotherapy/GDSC/GDSC1_fitted_dose_response_27Oct23.xlsx', engine='openpyxl')

# Sort all drugs by increasing sensitivity
data_sorted = data.sort_values(by='LN_IC50')  

# Filter by one random drug name: Erlotinib
data_sorted_erlotinib = data_sorted[data_sorted['DRUG_NAME'] == 'Erlotinib']  

# Filter by the mode (drug with more rows of data) (one random criteria that can be used for selecting the 2 drugs)
def compute_mode_drug_id_monotherapy(data):
    mode = data['DRUG_ID'].mode().values[0]
    return mode
moda = compute_mode_drug_id_monotherapy(data_sorted)
moda_drug_name = data_sorted[data_sorted['DRUG_ID'] == moda]['DRUG_NAME'].iloc[0] # Get the drug name corresponding to the mode drug id
data_sorted_1drug = data_sorted[data_sorted['DRUG_ID'] == moda] #filter only the mode drug id



####################################### SYNERGY DATA #######################################

data_synergy = pd.read_csv('Data/Synergy/NCI/ComboDrugGrowth_Nov2017.csv')

# Choose a random drug pair for synergy plot
nsc1, nsc2 = 752, 3088  # Looked at NSC webpage (https://discover.nci.nih.gov/cellminerdata/html/NSC_LIST.html) and got the drugs names: 6-THIOGUANINE and Chlorambucil

# Get the most and the 2nd most common drug id for synergy data (one random criteria that can be used for selecting the 2 drugs)
def compute_2modes_nsc_synergy(data):
    mode1 = data['NSC1'].mode().values[0]
    mode2 = data[data['NSC1'] != mode1]['NSC1'].mode().values[0]
    return mode1, mode2

moda1, moda2 = compute_2modes_nsc_synergy(data_synergy)
# print(f"The mode of DRUG_ID is: {moda1}, {moda2}") # I used this print to look at the drug name in NSC webpage (https://discover.nci.nih.gov/cellminerdata/html/NSC_LIST.html)
data_sorted_2drugs = data_synergy[((data_synergy['NSC1'] == moda1) & (data_synergy['NSC2'] == moda2)) | ((data_synergy['NSC1'] == moda2) & (data_synergy['NSC2'] == moda1))]



####################################### PLOTTING FUNCTIONS #######################################

def plot_lnIC50_1drug(data_sorted, drug_name, output_dir):
    plt.figure(figsize=(12, 6))
    # Scatter plot for individual points
    sns.scatterplot(x=range(len(data_sorted)), y='LN_IC50', hue='TCGA_DESC', data=data_sorted, s=10, palette="tab10", legend=False)

    # Line plot for general trend
    sns.lineplot(x=range(len(data_sorted)), y='LN_IC50', data=data_sorted, color='red')

    # Styling the plot
    plt.xlabel('Ranked Cell Lines by Sensitivity')
    plt.ylabel('Sensitivity (ln IC50)')
    plt.title(f'Sensitivity of Cell Lines to {drug_name} (ln IC50)')
    plt.grid(True, linestyle='--', linewidth=0.5)
    plt.tight_layout()
    #plt.show()

    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    # Save the plot
    plt.savefig(os.path.join(output_dir, 'Monotherapy_plot_' + drug_name + '.png'))


def plot_synergy_2drugs(data, nsc1, nsc2, nameDrug1, nameDrug2, output_dir):
    # Filter the data for the selected drug pair
    drug_pair_data = data[((data['NSC1'] == nsc1) & (data['NSC2'] == nsc2)) | ((data['NSC1'] == nsc2) & (data['NSC2'] == nsc1))]

    # Sort by synergy score (SCORE column) TODO: Try to find what is the SCORE computed from
    drug_pair_data_sorted = drug_pair_data.sort_values(by='SCORE', ascending=True) 

    # Plot the synergy scores across cell lines
    plt.figure(figsize=(12, 6))
    sns.scatterplot(x=range(len(drug_pair_data_sorted)), y='SCORE', data=drug_pair_data_sorted, color='blue', s=10)
    #sns.lineplot(x=range(len(drug_pair_data_sorted)), y='SCORE', data=drug_pair_data_sorted, color='blue')

    # Add labels and title
    plt.xlabel('Cell Lines (ordered by synergy score)')
    plt.ylabel('Synergy Score')
    plt.title(f'Synergy Scores of Drug Pair NSC {nsc1} ({nameDrug1}) & NSC {nsc2} ({nameDrug2}) Across Cell Lines')
    plt.grid(True, linestyle='--', linewidth=0.5)
    plt.tight_layout()
    #plt.show()

    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    # Save the plot
    plt.savefig(os.path.join(output_dir, 'Monotherapy_plot_' + nameDrug1 + '_' + nameDrug2 + '.png'))


def plot_synergy_2drugs_interactive(data, nsc1, nsc2, output_dir):
    # Filter the data for the selected drugs only
    drug_pair_data = data[((data['NSC1'] == nsc1) & (data['NSC2'] == nsc2)) | ((data['NSC1'] == nsc2) & (data['NSC2'] == nsc1))]
    drug_pair_data_sorted = drug_pair_data.sort_values(by='SCORE', ascending=True)
    
    # Interactive plot creation 
    fig = px.scatter(
        drug_pair_data_sorted,
        x=range(len(drug_pair_data_sorted)),
        y='SCORE',
        labels={'x': 'Cell Lines (ordered by synergy score)', 'SCORE': 'Synergy Score'},
        title=f'Synergy Scores of Drug Pair NSC {nsc1} & NSC {nsc2} Across Cell Lines',
        hover_data={'CONCINDEX1': True, 'CONCINDEX2': True}
    )

    fig.update_traces(marker=dict(color='blue', size=8))
    # fig.show()

    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    
    # Save the plot
    fig.write_html(os.path.join(output_dir, f'Synergy_plot_{nsc1}_{nsc2}.html'))


####################################### PLOTTING #######################################

# Create directory if it doesn't exist
output_dir = 'Plots'

plot_lnIC50_1drug(data_sorted, "All Drugs", output_dir)
plot_lnIC50_1drug(data_sorted_erlotinib, "Erlotinib", output_dir)
plot_lnIC50_1drug(data_sorted_1drug, moda_drug_name, output_dir)

plot_synergy_2drugs(data_synergy,nsc1,nsc2,"6-THIOGUANINE","Chlorambucil", output_dir)
plot_synergy_2drugs_interactive(data_synergy,nsc1,nsc2, output_dir)

plot_synergy_2drugs(data_synergy,moda1,moda2,"Epirubicin","Idarubicin", output_dir)
plot_synergy_2drugs_interactive(data_sorted_2drugs,moda1,moda2, output_dir)
