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
write.csv(gene_expression_matrix_NCI_T, "gene_expression_matrix_NCI_allDrugs.csv", row.names = FALSE)

# Additional operations to filter only relevant expression data for the drug selected: 
columns_to_keep_Sorafenib <- which(gene_expression_matrix_NCI_T["gene_id", ] == "ENSG00000157764" |
                                     gene_expression_matrix_NCI_T["gene_id", ] == "ENSG00000077782" |
                                     gene_expression_matrix_NCI_T["gene_id", ] == "ENSG00000102755" |
                                     gene_expression_matrix_NCI_T["gene_id", ] == "ENSG00000122025" |
                                     gene_expression_matrix_NCI_T["gene_id", ] == "ENSG00000037280" |
                                     gene_expression_matrix_NCI_T["gene_id", ] == "ENSG00000128052" |
                                     gene_expression_matrix_NCI_T["gene_id", ] == "ENSG00000157404" |
                                     gene_expression_matrix_NCI_T["gene_id", ] == "ENSG00000113721" |
                                     gene_expression_matrix_NCI_T["gene_id", ] == "ENSG00000132155" |
                                     gene_expression_matrix_NCI_T["gene_id", ] == "ENSG00000165731")
columns_to_keep_Sorafenib <- c("cell_line", colnames(gene_expression_matrix_NCI_T)[columns_to_keep_Sorafenib])
filtered_matrix_NCI_Sorafenib <- gene_expression_matrix_NCI_T[, columns_to_keep_Sorafenib, drop = FALSE]

write.csv(filtered_matrix_NCI_Sorafenib, "gene_expression_matrix_NCI_Sorafenib.csv", row.names = FALSE)

###################### Obtain ic50 and AAC ##################################
# Drug names are in the 'drugid' column within the @drug slot
drug_names_NCI <- unique(dataNCI60@drug$drugid)
print(drug_names_NCI)
writeLines(as.character(drug_names_NCI), "DrugNamesNCI.txt")
which(drug_names_NCI == "Camptothecin")
which(drug_names_NCI == "Sorafenib")

drug_of_interest <- "Sorafenib"


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

# Merge by cell line to align the data
#gene_expression_NCI <- read.csv("gene_expression_matrix_NCI_T.csv")
#ic50_aac_NCI <- read.csv("ic50_aac_NCI.csv")
#merged_data_NCI <- merge(ic50_aac_NCI, gene_expression_NCI, by = "cell_line")
merged_data_NCI <- merge(ic50_aac_NCI_df, gene_expression_matrix_NCI_T, by = "cell_line")
write.csv(merged_data_NCI, "merged_NCI.csv", row.names = FALSE)

print("Data preprocessing complete. Merged file with all gene expression data + IC50 + AAC saved as 'merged_NCI.csv'")


# Prepare a csv for only relevant genes
#gene_expression_NCI_relevant <- read.csv("gene_expression_matrix_NCI_RELEVANT.csv")
#merged_data_NCI_relevant <- merge(ic50_aac_NCI, gene_expression_NCI_relevant, by = "cell_line")

############################### CAMPTOTHECIN #####################################
#merged_data_NCI_Camptothecin <- merge(ic50_aac_NCI_df, filtered_matrix_NCI_Camptothecin, by="cell_line")
#write.csv(merged_data_NCI_Camptothecin, "merged_NCI_Camptothecin.csv", row.names = FALSE)

############################### SORAFENIB #####################################
merged_data_NCI_Sorafenib <- merge(ic50_aac_NCI_df, filtered_matrix_NCI_Sorafenib, by="cell_line")
write.csv(merged_data_NCI_Sorafenib, "merged_NCI_Sorafenib.csv", row.names = FALSE)


print("Data preprocessing complete. Merged file with all gene expression data + IC50 + AAC saved as 'merged_NCI_Camptothecin.csv'")

############################################################################## 
###########################    GDSC    #######################################

dataGDSC <- readRDS("GDSC2.rds")


###################### Obtain ic50 and AAC GDSC ##################################

drug_names_GDSC <- unique(dataGDSC@drug$drugid)
print(drug_names_GDSC)
writeLines(as.character(drug_names_GDSC), "DrugNamesGDSC.txt")

# Sensitivity profiles for IC50 and AAC
sensitivity_data_GDSC <- dataGDSC@sensitivity$profiles

drug_of_interest2 <- "Camptothecin"

# Filter IC50 and AAC for the given drug
drug_ic50_GDSC <- sensitivity_data_GDSC$ic50_recomputed
drug_aac_GDSC <- sensitivity_data_GDSC$aac_recomputed
print(max(drug_ic50_GDSC, na.rm = TRUE))
# Extract corresponding cell line names
cell_lines_GDSC <- dataGDSC@sensitivity$info$CELL_LINE_NAME

# Combine into a data frame
ic50_aac_GDSC_df <- data.frame(
  cell_line = cell_lines_GDSC,
  IC50 = drug_ic50_GDSC,
  AAC = drug_aac_GDSC
)

# Group by cell line and calculate the average IC50 and AAC
library(dplyr)

# Define the threshold
threshold <- 1e6

# Filter out IC50 values above the threshold, then group by cell line and calculate averages
ic50_aac_GDSC_avg <- ic50_aac_GDSC_df %>%
  filter(IC50 < threshold) %>% 
  group_by(cell_line) %>%
  summarize(
    IC50 = mean(IC50, na.rm = TRUE),
    AAC = mean(AAC, na.rm = TRUE)
  )


# Write to CSV 
write.csv(ic50_aac_GDSC_avg, "ic50_aac_GDSC.csv", row.names = FALSE)



###################### Obtain gene expression GDSC ###########################

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

# Replace periods with dashes in row names (which are now the cell line names after transpose)
#rownames(gene_expression_matrix_GDSC_T) <- gsub("\\.", "-", rownames(gene_expression_matrix_GDSC_T))

gene_expression_matrix_GDSC_T$cell_line <- gsub("\\.", "-", rownames(gene_expression_matrix_GDSC_T))

write.csv(gene_expression_matrix_GDSC_T, "gene_expression_matrix_GDSC.csv", row.names = FALSE)

# Additional operations to filter only relevant expression data for the drug selected: Camptothecin
# Use read.delim for tab-separated files
ENSGS <- read.delim("../results/AffectedGenesByDrug_ENSG_GDSC.txt", 
                    stringsAsFactors = FALSE, 
                    sep = "\t", # Explicitly specify tab as the delimiter
                    header = TRUE) # Ensure the first row is treated as headers

# Check the structure of the dataframe
str(ENSGS)
# Loop through each drug in ENSGS
for (i in 1:nrow(ENSGS)) {
  # Extract drug information
  chebi_id <- ENSGS$CHEBI_ID[i]
  drug_name <- ENSGS$DRUG_NAME[i]
  gene_list <- unlist(strsplit(ENSGS$GENES_AFFECTED[i], ", "))
  
  # Filter columns to keep
  available_genes <- intersect(gene_list, colnames(gene_expression_matrix_GDSC_T))
  columns_to_keep <- c("cell_line", available_genes)
  
  # Check for missing genes and print a warning if any are missing
  missing_genes <- setdiff(gene_list, colnames(gene_expression_matrix_GDSC_T))
  if (length(missing_genes) > 0) {
    warning(paste("For drug", drug_name, "the following genes are missing and will be ignored:", 
                  paste(missing_genes, collapse = ", ")))
  }
  
  # Skip if no valid genes are found
  if (length(columns_to_keep) <= 1) {  # Only "cell_line" is present
    warning(paste("Skipping drug", drug_name, "because no valid genes are found in the matrix."))
    next
  }
  
  # Filter the gene expression matrix for this drug
  filtered_matrix <- gene_expression_matrix_GDSC_T[, columns_to_keep, drop = FALSE]
  
  # Merge with IC50 and AAC data
  merged_data <- merge(ic50_aac_GDSC_avg, filtered_matrix, by = "cell_line")
  
  # Create a dynamic filename
  output_filename <- paste0("merged_GDSC_", drug_name, ".csv")
  
  # Write the filtered matrix to a CSV file
  write.csv(merged_data, output_filename, row.names = FALSE)
}

###################### Merge gene expression and ic50/AAC ####################

# Merge by cell line to align the data
#gene_expression_GDSC <- read.csv("gene_expression_matrix_GDSC_T.csv")
#ic50_aac_GDSC <- read.csv("ic50_aac_GDSC.csv")
#merged_data_GDSC <- merge(ic50_aac_GDSC, gene_expression_GDSC, by = "cell_line")
merged_data_GDSC <- merge(ic50_aac_GDSC_avg, gene_expression_matrix_GDSC_T, by = "cell_line")
write.csv(merged_data_GDSC, "merged_GDSC.csv", row.names = FALSE)

print("Data preprocessing complete. Merged file with all gene expression data + IC50 + AAC saved as 'merged_GDSC.csv'")
