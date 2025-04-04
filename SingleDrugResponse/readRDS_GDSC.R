library(PharmacoGx)
library(dplyr)

# Function to load dataset
load_dataset <- function(file_path) {
  readRDS(file_path)
}

# Function to extract all different Drug Names in the PSets
extract_drug_names <- function(data, output_path){
  drug_names <- drugNames(data)
  writeLines(as.character(drug_names), output_path)
}


extract_gene_expression <- function(data, output_path, profile_name = "Kallisto_0.46.1.rnaseq") {
  experiment <- molecularProfiles(data)[[profile_name]]
  
  # Extract gene expression matrix
  gene_expression_matrix <- as.data.frame(assay(experiment))
  
  # Filter only protein-coding genes
  relevant_gene_info <- as.data.frame(rowData(experiment)) %>%
    filter(gene_type == "protein_coding") %>%
    select(gene_name)
  
  gene_expression_matrix <- merge(relevant_gene_info, gene_expression_matrix, by=0) 
  
  # Aggregate duplicate gene names using median before transposing
  gene_expression_matrix <- gene_expression_matrix %>%
    group_by(gene_name) %>%
    summarize(across(where(is.numeric), median, na.rm = TRUE)) %>%
    as.data.frame()
  
  # Store gene names before transposing
  gene_names <- gene_expression_matrix$gene_name
  gene_expression_matrix <- gene_expression_matrix %>% select(-gene_name)  # Remove gene_name column before transposing
  
  # Transpose matrix (genes as columns, cell lines as rows)
  gene_expression_matrix_T <- as.data.frame(t(gene_expression_matrix))
  
  # Restore gene names as column names
  colnames(gene_expression_matrix_T) <- gene_names
  
  # Ensure rownames of gene_expression_matrix_GDSC_T are stored
  gene_expression_matrix_T$EGAR_ID <- rownames(gene_expression_matrix_T)
  
  # Extract mapping of EGAR_ID to sampleid from colData
  egar_to_sampleid <- data.frame(
    EGAR_ID = rownames(colData(experiment)), #EGAR_NAMES
    cell_line = colData(experiment)$sampleid
  )
  
  sample_info <- as.data.frame(sampleInfo(data))
  
  keep_sample_ids <- sample_info %>%
    filter(!tissueid %in% c("Myeloid", "Lymphoid")) %>%
    select(sampleid) %>%
    pull(sampleid)
  
  # Merge with sample ID mapping
  gene_expression_matrix_T <- merge(egar_to_sampleid, gene_expression_matrix_T, by = "EGAR_ID")
  gene_expression_matrix_T$EGAR_ID <- NULL  # Drop redundant ID column
  
  # Filter to keep only cell lines in keep_sample_ids
  gene_expression_matrix_T <- gene_expression_matrix_T %>%
    filter(cell_line %in% keep_sample_ids)
  
  # Move sampleid (cell_line) column to the first position
  gene_expression_matrix_T <- gene_expression_matrix_T[, c("cell_line", setdiff(names(gene_expression_matrix_T), "cell_line"))]
  
  # Clean column names (remove version numbers from gene IDs)
  colnames(gene_expression_matrix_T) <- gsub("\\.\\d+$", "", colnames(gene_expression_matrix_T))
 
  # Save to CSV
  #write.csv(gene_expression_matrix_T, output_path, row.names = FALSE)
  
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
                      sep = "\t", # tab is the delimiter
                      header = TRUE) # First row is treated as header
  
  # Check structure of dataframe
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
    
    print(paste("GENE COLUMNS TO KEEP for drug", drug_name))
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
    drug_ic50_GDSC <- summarizeSensitivityProfiles(dataGDSC, sensitivity.measure = "ic50_recomputed")[drug_name, ]
    drug_aac_GDSC <- summarizeSensitivityProfiles(dataGDSC, sensitivity.measure = "aac_recomputed")[drug_name, ]
    
    # Extract corresponding cell line names
    cell_lines_GDSC <- colnames(summarizeSensitivityProfiles(dataGDSC, sensitivity.measure = "ic50_recomputed"))
    
    # Combine into a data frame
    ic50_aac_GDSC_df <- data.frame(
      cell_line = cell_lines_GDSC,
      IC50 = as.numeric(drug_ic50_GDSC),
      AAC = as.numeric(drug_aac_GDSC)
    )
    
    # Group by cell line and calculate the average IC50 and AAC
    
    # Filter out IC50 values above the threshold, then group by cell line and calculate averages
    ic50_aac_GDSC_avg <- ic50_aac_GDSC_df %>%
      #filter(IC50 < threshold) %>% 
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
    merged_data <- merge(ic50_aac_GDSC_avg, matrix, by = "cell_line")
    
    # Create a dynamic filename
    output_filename <- paste0("merged_GDSC_all_gene_expression_", drug_name, ".csv")
    
    # Write the filtered matrix to a CSV file
    write.csv(merged_data, output_filename, row.names = FALSE)
  }
}

# Main Function
main <- function() {
  # Load datasets
  dataGDSC <- load_dataset("GDSC2.rds")
  dataGDSC <- updateObject(dataGDSC)
  
  # Extract drug names of each PSet
  extract_drug_names(dataGDSC, "DrugNamesGDSC.txt")
  
  # Process GDSC data
  gene_expression_matrix_GDSC <- extract_gene_expression(dataGDSC, "gene_expression_matrix_GDSC.csv")
  
  # Filter only relevant data
  process_GDSC(dataGDSC, gene_expression_matrix_GDSC, "GDSC_output", "../results/AffectedGenesByDrug_GDSC.txt")
}

# Run the main function
main()
