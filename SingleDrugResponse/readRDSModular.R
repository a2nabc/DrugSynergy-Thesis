library(PharmacoGx)
library(dplyr)

# Function to load dataset
load_dataset <- function(file_path) {
  readRDS(file_path)
}

# Function to extract gene expression data
extract_gene_expression <- function(data, is_GDSC, output_path) {
  gene_expression <- data@molecularProfiles$rna
  gene_expression_matrix <- as.data.frame(assay(gene_expression))
  gene_names <- rownames(gene_expression_matrix) #gene names

  
  if (is_GDSC==TRUE) {
    gene_names <- rownames(gene_expression_matrix)  # Gene names
    cell_line_names <- colnames(gene_expression_matrix)  # Cell line names
    cell_id_values <- gene_expression@colData@listData$cellid
    colnames(gene_expression_matrix) <- cell_id_values
    gene_expression_matrix_T <- as.data.frame(t(gene_expression_matrix))
    gene_expression_matrix_T$cell_line <- gsub("\\.", "-", rownames(gene_expression_matrix_T))
  }
  else {
    cell_line_names <- colnames(gene_expression_matrix)  # Cell line names
    geneid_names <- gene_expression@rowRanges@elementMetadata@listData$gene_id
    gene_expression_matrix_with_id <- cbind(gene_id = geneid_names, gene_expression_matrix)
    gene_expression_matrix_T <- as.data.frame(t(gene_expression_matrix_with_id))
    gene_expression_matrix_T$cell_line <- rownames(gene_expression_matrix_T)
  }
  write.csv(gene_expression_matrix_T, output_path, row.names = FALSE)
  return(gene_expression_matrix_T)
}

# Function to filter gene expression data for specific genes
filter_gene_expression <- function(gene_expression_matrix_T, genes_of_interest, output_path) {
  columns_to_keep <- c("cell_line", colnames(gene_expression_matrix_T)[
    gene_expression_matrix_T["gene_id", ] %in% genes_of_interest
  ])
  filtered_matrix <- gene_expression_matrix_T[, columns_to_keep, drop = FALSE]
  write.csv(filtered_matrix, output_path, row.names = FALSE)
  return(filtered_matrix)
}

# Function to extract IC50 and AAC data
extract_ic50_aac <- function(data, drug_name, output_path) {
  sensitivity_data <- data@sensitivity$profiles
  drug_ic50 <- sensitivity_data$ic50_recomputed[data@drug$drugid == drug_name]
  drug_aac <- sensitivity_data$aac_recomputed[data@drug$drugid == drug_name]
  ic50_aac_df <- data.frame(
    cell_line = data@sensitivity$info$cellid[data@drug$drugid == drug_name],
    IC50 = drug_ic50,
    AAC = drug_aac
  )
  write.csv(ic50_aac_df, output_path, row.names = FALSE)
  return(ic50_aac_df)
}

# Function to merge gene expression data with IC50/AAC
merge_data <- function(ic50_aac_df, gene_expression_matrix, output_path) {
  merged_data <- merge(ic50_aac_df, gene_expression_matrix, by = "cell_line")
  write.csv(merged_data, output_path, row.names = FALSE)
  return(merged_data)
}

# Additional operations to filter only relevant expression data
process_GDSC <- function(dataGDSC, matrix, output_path_prefix, ENSGS_path, threshold = 1e6) {
  ENSGS <- read.delim(ENSGS_path, 
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
    available_genes <- intersect(gene_list, colnames(matrix))
    columns_to_keep <- c("cell_line", available_genes)
    
    print("COLS TO KEEP")
    print(columns_to_keep)
    
    # Check for missing genes and print a warning if any are missing
    missing_genes <- setdiff(gene_list, colnames(matrix))
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
    filtered_matrix <- matrix[, columns_to_keep, drop = FALSE]
    
    
    # Filter IC50 and AAC for the given drug
    drug_id <- dataGDSC@drug[dataGDSC@drug$drugid == drug_name, "drugid"]
    
    sensitivity_data_GDSC <- dataGDSC@sensitivity$profiles
    drug_ic50_GDSC <- sensitivity_data_GDSC$ic50_recomputed[dataGDSC@drug$drugid == drug_id]
    drug_aac_GDSC <- sensitivity_data_GDSC$aac_recomputed[dataGDSC@drug$drugid == drug_id]
    
    # Extract corresponding cell line names
    cell_lines_GDSC <- dataGDSC@sensitivity$info$CELL_LINE_NAME[dataGDSC@drug$drugid == drug_id]
    
    # Combine into a data frame
    ic50_aac_GDSC_df <- data.frame(
      cell_line = cell_lines_GDSC,
      IC50 = drug_ic50_GDSC,
      AAC = drug_aac_GDSC
    )
    
    # Group by cell line and calculate the average IC50 and AAC

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
    
    # ONLY RELEVANT GENES
    # Merge with IC50 and AAC data
    merged_data_RelevantGenes <- merge(ic50_aac_GDSC_avg, filtered_matrix, by = "cell_line")
    
    # Create a dynamic filename
    output_filename <- paste0("merged_GDSC_ONLY_RELEVANT_", drug_name, ".csv")
    
    # Write the filtered matrix to a CSV file
    write.csv(merged_data_RelevantGenes, output_filename, row.names = FALSE)
    
    
    # ALL GENE EXPRESSION DATA
    # Merge with IC50 and AAC data
    merged_data <- merge(ic50_aac_GDSC_avg, gene_expression_matrix_GDSC_T, by = "cell_line")
    
    # Create a dynamic filename
    output_filename <- paste0("merged_GDSC_all_gene_expression_", drug_name, ".csv")
    
    # Write the filtered matrix to a CSV file
    write.csv(merged_data, output_filename, row.names = FALSE)
  }
}

# Main Function
main <- function() {
  # Load datasets
  dataNCI60 <- load_dataset("PSet_NCI60.rds")
  dataGDSC <- load_dataset("GDSC2.rds")
  
  # Process NCI60 data
  isGDSC <- FALSE
  gene_expression_matrix_NCI <- extract_gene_expression(dataNCI60, isGDSC, "gene_expression_matrix_NCI_allDrugs.csv")
  
  # Process GDSC data
  isGDSC <- TRUE
  gene_expression_matrix_GDSC <- extract_gene_expression(dataGDSC, isGDSC, "gene_expression_matrix_GDSC.csv")
  
  # Filter only relevant data
  process_GDSC(dataGDSC, gene_expression_matrix_GDSC, "GDSC_output", "../results/AffectedGenesByDrug_ENSG_GDSC.txt")
}

# Run the main function
main()
