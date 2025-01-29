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

# Function to process GDSC data
process_GDSC <- function(dataGDSC, output_path_prefix, ENSGS_path, threshold = 1e6) {
  drug_names_GDSC <- unique(dataGDSC@drug$drugid)
  sensitivity_data <- dataGDSC@sensitivity$profiles
  
  # Load ENSGS data
  ENSGS <- read.delim(ENSGS_path, stringsAsFactors = FALSE, sep = "\t", header = TRUE)
  
  for (i in 1:nrow(ENSGS)) {
    drug_name <- ENSGS$DRUG_NAME[i]
    gene_list <- unlist(strsplit(ENSGS$GENES_AFFECTED[i], ", "))
    
    # Filter gene expression matrix
    filtered_matrix <- filter_gene_expression(
      gene_expression_matrix_GDSC_T, 
      gene_list, 
      paste0(output_path_prefix, "_filtered_", drug_name, ".csv")
    )
    
    # Extract IC50 and AAC
    ic50_aac_df <- extract_ic50_aac(dataGDSC, drug_name, paste0(output_path_prefix, "_ic50_aac_", drug_name, ".csv"))
    ic50_aac_avg <- ic50_aac_df %>%
      filter(IC50 < threshold) %>%
      group_by(cell_line) %>%
      summarize(IC50 = mean(IC50, na.rm = TRUE), AAC = mean(AAC, na.rm = TRUE))
    
    # Merge data
    merged_data <- merge(ic50_aac_avg, filtered_matrix, by = "cell_line")
    write.csv(merged_data, paste0(output_path_prefix, "_merged_", drug_name, ".csv"), row.names = FALSE)
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
  gene_expression_matrix_GDSC_T <- extract_gene_expression(dataGDSC, isGDSC, "gene_expression_matrix_GDSC.csv")
  
  # Filter only relevant data
  process_GDSC(dataGDSC, "GDSC_output", "../results/AffectedGenesByDrug_ENSG_GDSC.txt")
}

# Run the main function
main()
