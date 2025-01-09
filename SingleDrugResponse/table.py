def create_new_table(input_file, output_file_1, output_file_2):
    with open(input_file, 'r') as infile:
        lines = infile.readlines()

    # creation of a table with interactions: PARTICIPANT_A (CHEBI ID), INTERACTION_TYPE, PARTICIPANT_B (Gene or protein affected)
    with open(output_file_1, 'w') as outfile_1:
        # Write the header line
        outfile_1.write('\t'.join(lines[0].split('\t')[:3]) + '\n')
        for line in lines[1:10910]:
            elements = line.split('\t')
            outfile_1.write('\t'.join(elements[:3]).strip() + '\n')

    
    #creation of a table with details about participants: PARTICIPANT_A (CHEBI ID), PARTICIPANT_B (Gene or protein affected), PARTICIPANT_C (Drug name)
    with open(output_file_2, 'w') as outfile_2:
        # Write the header line
        outfile_2.write(lines[10910])
        for line in lines[10910:]:
            elements = line.split('\t')
            if len(elements) > 1 and elements[1] == 'SmallMoleculeReference':
                outfile_2.write(line.strip() + '\n')


def filter_by_drug_names(input_file, drug_names_file, output_file):
    with open(drug_names_file, 'r') as dnfile:
        drug_names = set(name.lower() for name in dnfile.read().splitlines())

    with open(input_file, 'r') as infile:
        lines = infile.readlines()

    with open(output_file, 'w') as outfile:
        # Write the header line
        header_elements = lines[0].split('\t')
        outfile.write('\t'.join([header_elements[0], header_elements[2]]) + '\n')
        for line in lines[1:]:
            elements = line.split('\t')
            if len(elements) > 2 and elements[2].lower() in drug_names:
                elements[2] = elements[2].lower()
                outfile.write('\t'.join([elements[0], elements[2]]) + '\n')

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

        # Output another file with chebi, drug name, and the names of all genes affected
        with open(output_file_genes, 'w') as output_genes_f:
            output_genes_f.write('CHEBI_ID\tDRUG_NAME\tGENES_AFFECTED\n')
            for chebi_id, participants in chebi_counts.items():
                drug_name = chebi_drug_names[chebi_id]
                genes_affected = ', '.join(participants)
                #genes_affected_vector = genes_affected.split(', ')
                #print(genes_affected_vector)
                output_genes_f.write(f'{chebi_id}\t{drug_name}\t{genes_affected}\n')


input_file = './data/PathwayCommons12.drugbank.hgnc.txt'
output_file_1 = './results/Interactions.txt'
output_file_2 = './results/DetailsAboutParticipants.txt'

drug_names_file_GDSC = './data/DrugNamesGDSC.txt'
output_filtered_file_GDSC = './results/DrugName_vs_Chebi_GDSC.txt'

drug_names_file_NCI = './data/DrugNamesNCI.txt'
output_filtered_file_NCI = './results/DrugName_vs_Chebi_NCI.txt'

create_new_table(input_file, output_file_1, output_file_2)
filter_by_drug_names(output_file_2, drug_names_file_GDSC, output_filtered_file_GDSC)
filter_by_drug_names(output_file_2, drug_names_file_NCI, output_filtered_file_NCI)

output_count_file_GDSC = './results/NumberGenesAffectedByDrug_GDSC.txt'
output_drug_genes_affected_GDSC = './results/AffectedGenesByDrug_GDSC.txt'
count_chebi_interactions(output_filtered_file_GDSC, output_file_1, output_count_file_GDSC, output_drug_genes_affected_GDSC)
output_count_file_NCI = './results/NumberGenesAffectedByDrug_NCI.txt'
output_drug_genes_affected_NCI = './results/AffectedGenesByDrug_NCI.txt'
count_chebi_interactions(output_filtered_file_NCI, output_file_1, output_count_file_NCI, output_drug_genes_affected_NCI)