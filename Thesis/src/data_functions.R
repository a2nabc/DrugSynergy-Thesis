library(dplyr)
library(tidyr)

# Function to get 
get_gene_union <- function(gene_subsets, path) {
  all_genes <- character(0)
  for (subset in gene_subsets) {
    genes <- readLines(paste0(path, subset, ".txt"))
    all_genes <- union(all_genes, genes)
  }
  return(all_genes)
}

load_data <- function(expression_path, response_path) {
  expression_data <- read.csv(expression_path)
  response_data <- read.csv(response_path)
  return(list(expression = expression_data, response = response_data))
}

create_output_folders <- function(target_set, model_types, experiment_positive, results_dir, features_dir, eval_dir, correlations_dir, cv_dir) {
  results_subdir <- if (experiment_positive == "true") paste0(target_set, "/positive/") else paste0(target_set, "/negative/")
  dir.create(results_dir, showWarnings = FALSE)
  dir.create(paste0(results_dir, target_set, "/"), showWarnings = FALSE)
  dir.create(paste0(results_dir, results_subdir), showWarnings = FALSE)
  dir.create(paste0(results_dir, results_subdir, "metadata/"), showWarnings = FALSE)
  for (model in model_types) {
    dir.create(paste0(results_dir, results_subdir, model), showWarnings = FALSE)
    dir.create(paste0(results_dir, results_subdir, model, "/",  features_dir), showWarnings = FALSE)
    dir.create(paste0(results_dir, results_subdir, model, "/", eval_dir), showWarnings = FALSE)
    dir.create(paste0(results_dir, results_subdir, model, "/", correlations_dir), showWarnings = FALSE)
    dir.create(paste0(results_dir, results_subdir, model, "/", cv_dir), showWarnings = FALSE)
  }
  return(results_subdir)
}

create_df_metadata <- function(results_dir, results_metadata_dir){
  data_dimensions <- data.frame(
    Drug = character(),
    Model_Type = character(),
    Training_Samples = integer(),
    Training_Features = integer(),
    Testing_Samples = integer(),
    Testing_Features = integer(),
    stringsAsFactors = FALSE
  )
  return(data_dimensions)
}

filter_cell_lines <- function(is_consistency_experiment, source_data, target_data) {
  source_cells <- source_data$expression$Cell_line
  target_cells <- target_data$expression$Cell_line
  
  if(is_consistency_experiment){ # consistency experiment --> we want to test on common cell lines
    common_cells <- intersect(source_cells, target_cells)
    target_data$expression <- target_data$expression %>%
      filter(Cell_line %in% common_cells)
    target_data$response <- target_data$response %>%
      filter(Cell_line %in% common_cells)
    print("COMMON cell lines")
  } else{ # generalization -> we want to test only in unseen cell_lines
    target_data$expression <- target_data$expression %>%
      filter(!Cell_line %in% source_cells)
    target_data$response <- target_data$response %>%
      filter(!Cell_line %in% source_cells)
    print("DIFFERENT cell lines")
  }
  
  return(target_data)
}

load_pagerank_feature_list <- function(gene_list_path) {
  if (file.exists(gene_list_path)) {
    return(readLines(gene_list_path))
  } else{
    print(paste("Pagerank feature list not found:", gene_list_path))
    return(NULL)
  }
}