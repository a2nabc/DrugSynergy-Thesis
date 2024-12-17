library(PharmacoGx)
###################### Obtain gene expression NCI ########################

dataNCI60 <- readRDS("PSet_NCI60.rds")

# Extract gene expression data and convert to data frame
gene_expression_NCI <- dataNCI60@molecularProfiles$rna
gene_expression_matrix_NCI <- as.data.frame(assay(gene_expression_NCI))  # Extract assay data (gene expression)

#
gene_names <- rownames(gene_expression_matrix_NCI)  # Gene names
cell_line_names <- colnames(gene_expression_matrix_NCI)  # Cell line names

geneid_names <- gene_expression_NCI@rowRanges@elementMetadata@listData$gene_id

gene_expression_matrix_NCI_with_id <- cbind(gene_id = geneid_names, gene_expression_matrix_NCI)

# Transpose to have cell lines as rows and genes as columns
gene_expression_matrix_NCI_T <- as.data.frame(t(gene_expression_matrix_NCI_with_id))
gene_expression_matrix_NCI_T$cell_line <- rownames(gene_expression_matrix_NCI_T)  # Add cell line IDs


# Write to CSV
write.csv(gene_expression_matrix_NCI_T, "gene_expression_matrix_NCI_T.csv", row.names = FALSE)

# Additional operations to filter only relevant expression data for the drug selected: Camptothecin
columns_to_keep <- which(gene_expression_matrix_NCI_T["gene_id", ] == "ENSG00000198900")
columns_to_keep <- c("cell_line", colnames(gene_expression_matrix_NCI_T)[columns_to_keep])
filtered_matrix_NCI <- gene_expression_matrix_NCI_T[, columns_to_keep, drop = FALSE]

write.csv(filtered_matrix_NCI, "gene_expression_matrix_NCI_RELEVANT.csv", row.names = FALSE)

###################### Obtain ic50 and AAC ##################################
# Drug names are in the 'drugid' column within the @drug slot
drug_names <- unique(dataNCI60@drug$drugid)
print(drug_names)

drug_of_interest <- "Camptothecin"


# Filter rows corresponding to the selected drug
sensitivity_data_NCI <- dataNCI60@sensitivity$profiles
drug_ic50_NCI <- sensitivity_data_NCI$ic50_recomputed[dataNCI60@drug$drugid == drug_of_interest]
drug_aac_NCI <- sensitivity_data_NCI$aac_recomputed[dataNCI60@drug$drugid == drug_of_interest]

# Combine IC50 and AAC values into a data frame
ic50_aac_NCI_df <- data.frame(
  cell_line = dataNCI60@sensitivity$info$cellid[dataNCI60@drug$drugid == drug_of_interest],
  IC50 = drug_ic50_NCI,
  AAC = drug_aac_NCI
)

write.csv(ic50_aac_NCI_df, "ic50_aac_NCI.csv", row.names = FALSE)



###################### Merge gene expression and ic50/AAC ####################

# Read the saved gene expression and IC50/AAC files
gene_expression_NCI <- read.csv("gene_expression_matrix_NCI_T.csv")
ic50_aac_NCI <- read.csv("ic50_aac_NCI.csv")

# Merge by cell line to align the data
merged_data_NCI <- merge(ic50_aac_NCI, gene_expression_NCI, by = "cell_line")

# Save merged data to CSV for use in Python
write.csv(merged_data_NCI, "merged_NCI.csv", row.names = FALSE)

print("Data preprocessing complete. Merged file saved as 'merged_NCI.csv'")



# Merge relevant columns only
gene_expression_NCI_relevant <- read.csv("gene_expression_matrix_NCI_RELEVANT.csv")

ic50_aac_NCI <- read.csv("ic50_aac_NCI.csv")# Merge by cell line to align the data
merged_data_NCI_relevant <- merge(ic50_aac_NCI, gene_expression_NCI_relevant, by = "cell_line")

# Save merged data to CSV for use in Python
write.csv(merged_data_NCI_relevant, "merged_NCI_RELEVANT.csv", row.names = FALSE)


############################################################################## 

###################### Obtain gene expression GDSC ###########################

dataGDSC <- readRDS("GDSC2.rds")

# Extract gene expression data
gene_expression_GDSC <- dataGDSC@molecularProfiles$rna  # Assuming 'rna' contains gene expression data

# Convert to a matrix 
gene_expression_matrix_GDSC <- as.data.frame(assay(gene_expression_GDSC))

gene_names <- rownames(gene_expression_matrix_GDSC)  # Gene names
cell_line_names <- colnames(gene_expression_matrix_GDSC)  # Cell line names

cellid_values <- gene_expression_GDSC@colData@listData$cellid

colnames(gene_expression_matrix_GDSC) <- cellid_values

# Transpose to have cell lines as rows and genes as columns
gene_expression_matrix_GDSC_T <- as.data.frame(t(gene_expression_matrix_GDSC))
gene_expression_matrix_GDSC_T$cell_line <- rownames(gene_expression_matrix_GDSC_T)  # Add cell line IDs

write.csv(gene_expression_matrix_GDSC_T, "gene_expression_matrix_GDSC_T.csv", row.names = FALSE)

# Additional operations to filter only relevant expression data for the drug selected: Camptothecin
columns_to_keep2 <- c("cell_line", "ENSG00000198900")
filtered_matrix_GDSC <- gene_expression_matrix_GDSC_T[, columns_to_keep2, drop = FALSE]

write.csv(filtered_matrix_GDSC, "gene_expression_matrix_GDSC_RELEVANT.csv", row.names = FALSE)


###################### Obtain ic50 and AAC ##################################
drug_names_GDSC <- unique(dataGDSC@drug$drugid)
print(drug_names_GDSC)


# Sensitivity profiles for IC50 and AAC
sensitivity_data_GDSC <- dataGDSC@sensitivity$profiles

drug_of_interest2 <- "Camptothecin"

# Filter IC50 and AAC for the given drug
drug_ic50_GDSC <- sensitivity_data_GDSC$ic50_recomputed[dataGDSC@sensitivity$info$drugid == drug_of_interest2]
drug_aac_GDSC <- sensitivity_data_GDSC$aac_recomputed[dataGDSC@sensitivity$info$drugid == drug_of_interest2]

# Extract corresponding cell line names
cell_lines_GDSC <- dataGDSC@sensitivity$info$CELL_LINE_NAME[dataGDSC@sensitivity$info$drugid == drug_of_interest2]

# Combine into a data frame
ic50_aac_GDSC_df <- data.frame(
  cell_line = cell_lines_GDSC,
  IC50 = drug_ic50_GDSC,
  AAC = drug_aac_GDSC
)

# Write to CSV 
write.csv(ic50_aac_GDSC_df, "ic50_aac_GDSC.csv", row.names = FALSE)


###################### Merge gene expression and ic50/AAC ####################

# Read the saved gene expression and IC50/AAC files
gene_expression_GDSC <- read.csv("gene_expression_matrix_GDSC_T.csv")
ic50_aac_GDSC <- read.csv("ic50_aac_GDSC.csv")

# Merge by cell line to align the data
merged_data_GDSC <- merge(ic50_aac_GDSC, gene_expression_GDSC, by = "cell_line")

# Save merged data to CSV for use in Python
write.csv(merged_data_GDSC, "merged_GDSC.csv", row.names = FALSE)

print("Data preprocessing complete. Merged file saved as 'merged_GDSC.csv'")


gene_expression_GDSC_relevant <- read.csv("gene_expression_matrix_GDSC_RELEVANT.csv")
merged_data_GDSC_relevant <- merge(ic50_aac_GDSC, gene_expression_GDSC_relevant, by="cell_line")

write.csv(merged_data_GDSC_relevant, "merged_GDSC_RELEVANT.csv", row.names = FALSE)
print("Data preprocessing complete. Merged file saved as 'merged_GDSC_RELEVANT.csv'")
