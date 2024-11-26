library(PharmacoGx)
library(SummarizedExperiment)

NCI_ALMANAC_2017 <- readRDS("NCI_ALMANAC_2017.rds")
NCI_gene_expression <- readRDS("PSet_NCI60.rds")
NCI_gene_expression <- updateObject(NCI_gene_expression) # I did this because an error line prompted


#visualitzo les dades q tinc:
NCI_ALMANAC_2017@treatmentResponse
NCI_ALMANAC_2017@treatmentResponse@rowData
NCI_ALMANAC_2017@treatmentResponse@colData
NCI_ALMANAC_2017@treatmentResponse@assays

NCI_gene_expression@treatmentResponse$info



#comprovacio de si puc relacionar els 2 datasets amb smapleid:
gene_sample_id <- NCI_gene_expression@treatmentResponse$info$sampleid
synergy_sample_id <- NCI_ALMANAC_2017@treatmentResponse@colData$sampleid
common_sample_ids <- gene_sample_id %in% synergy_sample_id
table(common_sample_ids)
intersection_sample_ids <- intersect(gene_sample_id, synergy_sample_id)


#intento exportar les dades de gene_expression de NCI_gene_expression i les de synergy de NCI_ALMANAC_2017

#filter and export gene expression data
rna_data <- NCI_gene_expression@molecularProfiles$rna
rna_df <- as.data.frame(assay(rna_data))  # Extract assay data

#extract sampleid and score columns for synergy data
synergy_data <- NCI_ALMANAC_2017@treatmentResponse$profiles
# Select the relevant columns: sampleid and SCORE
synergy_data_filtered <- synergy_data[, .(sampleid, treatment1id, treatment2id, SCORE)]
#synergy_data_filtered_common_sample_ids <- synergy_data_filtered[sampleid %in% intersection_sample_ids, ]



# Load required library
library(dplyr)
library(tidyr)
library(tibble)

# Reshape `rna_df` to a long format
rna_long <- rna_df %>%
  rownames_to_column(var = "gene") %>%   # Preserve gene names
  pivot_longer(cols = -gene,             # All other columns are samples
               names_to = "sampleid",    # New column for sample names
               values_to = "expression") # New column for expression values

head(rna_long)

write.csv(rna_long, "gene_expression.csv", row.names = FALSE)
write.csv(synergy_data_filtered, "synergy_data.csv", row.names = FALSE)

# THE FOLLOWING CODE EXPLODES BeCAUSE OF DIMENSIONALITY
# Join the reshaped RNA data with the second dataset

#merged_data <- synergy_data_filtered %>%
#  left_join(rna_long, by = "sampleid")  # Merge on sampleid

#head(merged_data)


# The idea was to output only 1 CSV with sampleid, gene, expression, treatment1id, treatment2id, and SCORE


