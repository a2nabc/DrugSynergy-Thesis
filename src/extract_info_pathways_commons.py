import os

def create_new_tables(input_file, output_file_1, output_file_2):
    with open(input_file, 'r') as infile:
        lines = infile.readlines()

    # Creating interactions table
    with open(output_file_1, 'w') as outfile_1:
        outfile_1.write('\t'.join(lines[0].split('\t')[:3]) + '\n')  # Header
        for line in lines[1:10910]:  # First 10910 lines for interactions
            elements = line.split('\t')
            outfile_1.write('\t'.join(elements[:3]).strip() + '\n')

    # Creating participant details table
    with open(output_file_2, 'w') as outfile_2:
        outfile_2.write(lines[10910])  # Header
        for line in lines[10910:]:  # Remaining lines for participants
            elements = line.split('\t')
            if len(elements) > 1 and elements[1] == 'SmallMoleculeReference':
                outfile_2.write(line.strip() + '\n')


def filter_by_drug_names(input_file, drug_names_file, output_Drug_Chebi_file, output_missingDrugs_file):
    with open(drug_names_file, 'r') as dnfile:
        drug_names = {name.lower(): name for name in dnfile.read().splitlines()}

    with open(input_file, 'r') as infile:
        lines = infile.readlines()

    with open(output_Drug_Chebi_file, 'w') as outfile:
        header_elements = lines[0].split('\t')
        outfile.write('\t'.join([header_elements[0], header_elements[2]]) + '\n')  # Keep only CHEBI ID and drug name

        for line in lines[1:]:
            elements = line.split('\t')
            if len(elements) > 2 and elements[2].lower() in drug_names:
                original_name = drug_names[elements[2].lower()]
                outfile.write('\t'.join([elements[0], original_name]) + '\n')
            else:
                with open(output_missingDrugs_file, 'a') as not_found_file:
                    if not_found_file.tell() == 0:
                        not_found_file.write("Drugs in DetailsAboutParticipants file but not in the list\n")
                    not_found_file.write(f"{elements[0]}, {elements[2]}\n")


def count_chebi_interactions(chebi_file, interaction_file, output_file_number, output_file_genes):
    chebi_counts = {}
    chebi_drug_names = {}

    with open(chebi_file, 'r') as chebi_f:
        chebi_lines = chebi_f.readlines()[1:]  # Skip header line
        for line in chebi_lines:
            elements = line.split('\t')
            chebi_id = elements[0]
            drug_name = elements[1]
            chebi_drug_names[chebi_id] = drug_name.strip()

    with open(interaction_file, 'r') as interaction_f:
        interaction_lines = interaction_f.readlines()[1:]  # Skip header line
        for line in interaction_lines:
            elements = line.split('\t')
            chebi_id = elements[0]
            if chebi_id in chebi_drug_names:
                if chebi_id not in chebi_counts:
                    chebi_counts[chebi_id] = set()
                chebi_counts[chebi_id].add(elements[2].strip())

    with open(output_file_number, 'w') as output_num_f:
        output_num_f.write('CHEBI_ID\tDRUG_NAME\tCOUNT\n')
        for chebi_id, participants in chebi_counts.items():
            drug_name = chebi_drug_names[chebi_id]
            output_num_f.write(f'{chebi_id}\t{drug_name}\t{len(participants)}\n')

    with open(output_file_genes, 'w') as output_genes_f:
        output_genes_f.write('CHEBI_ID\tDRUG_NAME\tGENES_AFFECTED\n')
        for chebi_id, participants in chebi_counts.items():
            drug_name = chebi_drug_names[chebi_id]
            genes_affected = ', '.join(participants)
            output_genes_f.write(f'{chebi_id}\t{drug_name}\t{genes_affected}\n')


def count_unique_participant_a(interaction_file, output_file):
    participant_a_counts = {}

    with open(interaction_file, 'r') as infile:
        lines = infile.readlines()[1:]  # Skip header line
        for line in lines:
            elements = line.split('\t')
            participant_a = elements[0]
            if participant_a not in participant_a_counts:
                participant_a_counts[participant_a] = 0
            participant_a_counts[participant_a] += 1

    with open(output_file, 'w') as outfile:
        outfile.write('PARTICIPANT_A\tCOUNT\n')
        for participant_a, count in participant_a_counts.items():
            outfile.write(f'{participant_a}\t{count}\n')

        total_count = sum(participant_a_counts.values())
        outfile.write(f'TOTAL\t{total_count}\n')


# Define file paths
input_file = './data/raw/PathwayCommons12.drugbank.hgnc.txt'

os.makedirs('./results/drugbank/', exist_ok=True) 
output_file_1 = './results/drugbank/Interactions.txt'
output_file_2 = './results/drugbank/DetailsAboutParticipants.txt'

drug_names_file = './data/common_drugs.txt'
output_filtered_file = './results/drugbank/DrugName_vs_Chebi.txt'
missingDrugs_file = './results/drugbank/drugs_not_in_list.txt'

output_count_file = './results/drugbank/NumberGenesAffectedByDrug.txt'
output_drug_genes_affected = './results/drugbank/AffectedGenesByDrug.txt'
output_unique_participant_a_file = './results/drugbank/UniqueParticipantA_Count.txt'

# Run functions
create_new_tables(input_file, output_file_1, output_file_2)
filter_by_drug_names(output_file_2, drug_names_file, output_filtered_file, missingDrugs_file)
count_chebi_interactions(output_filtered_file, output_file_1, output_count_file, output_drug_genes_affected)
count_unique_participant_a(output_file_1, output_unique_participant_a_file)
