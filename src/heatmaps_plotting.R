library(reshape2)
library(dplyr)
library(tidyr)
library(ggplot2)
library(reshape2)
library(pheatmap)
library(purrr)


generate_model_heatmaps <- function(screen_path, correlations_subdir, drugs, plot_dir) {
  #We need to generate tables for both experiments: positive and negative_
  subfolders <- list.dirs(screen_path, recursive = FALSE)
  for (subfolder in subfolders) {
    lasso_path <- paste0(subfolder, "/lasso/", correlations_subdir)
    en_path <- paste0(subfolder, "/en/", correlations_subdir)
    ridge_path <- paste0(subfolder, "/ridge/", correlations_subdir)
    output_folder_path <- paste0(plot_dir, "/", basename(subfolder))
    print(lasso_path)
    
    generate_experiment_heatmaps(lasso_path, en_path, ridge_path, drugs, output_folder_path)
    
    if (grepl("positive", subfolder, ignore.case = TRUE)) {
      positive_lasso_path <- paste0(subfolder, "/lasso/", correlations_subdir)
      positive_en_path <- paste0(subfolder, "/en/", correlations_subdir)
      positive_ridge_path <- paste0(subfolder, "/ridge/", correlations_subdir)
    }
    else if (grepl("negative", subfolder, ignore.case = TRUE)) {
      # If negative, create paths for negative correlations
      negative_lasso_path <- paste0(subfolder, "/lasso/", correlations_subdir)
      negative_en_path <- paste0(subfolder, "/en/", correlations_subdir)
      negative_ridge_path <- paste0(subfolder, "/ridge/", correlations_subdir)
    }
  }
  generate_negative_positive_heatmaps(positive_lasso_path, negative_lasso_path, drugs, plot_dir)
  generate_negative_positive_heatmaps(positive_en_path, negative_en_path, drugs, plot_dir)
  generate_negative_positive_heatmaps(positive_ridge_path, negative_ridge_path, drugs, plot_dir)
  
}

generate_experiment_heatmaps <- function(lasso_path, en_path, ridge_path, drugs, plot_dir){
  for (drug in drugs){
    # Define output directory for plots
    plot_directory <- paste0(plot_dir,"/", drug)
    if (!dir.exists(plot_directory)) {
      dir.create(plot_directory, recursive = TRUE)
    }
    
    lasso_corrs <- read.csv(paste0(lasso_path, "/", drug, ".csv"), row.names = 1)
    en_corrs <- read.csv(paste0(en_path, "/", drug, ".csv"), row.names = 1)
    ridge_corrs <- read.csv(paste0(ridge_path, "/", drug, ".csv"), row.names = 1)
    
    # Convert row names back into a column for merging
    lasso_corrs <- lasso_corrs %>% mutate(genes = rownames(lasso_corrs))
    en_corrs <- en_corrs %>% mutate(genes = rownames(en_corrs))
    ridge_corrs <- ridge_corrs %>% mutate(genes = rownames(ridge_corrs))
    
    lasso_corrs <- rename(lasso_corrs, LASSO = x)
    en_corrs <- rename(en_corrs, EN = x)
    ridge_corrs <- rename(ridge_corrs, RIDGE = x)
    
    merged_data <-reduce(list(ridge_corrs, lasso_corrs, en_corrs), full_join, by="genes")
    
    # Replace NA values with 0 (or another value if needed)
    merged_data[is.na(merged_data)] <- 0  # Fill missing values
    
    # Set 'genes' as row names for heatmap compatibility
    rownames(merged_data) <- merged_data$genes
    
    merged_data <- merged_data %>% select(-"genes")  # Remove genes column
    
    merged_data_2 <-reduce(list(lasso_corrs, en_corrs), full_join, by="genes")
    
    # Replace NA values with 0 (or another value if needed)
    merged_data_2[is.na(merged_data_2)] <- 0  # Fill missing values
    
    # Set 'genes' as row names for heatmap compatibility
    rownames(merged_data_2) <- merged_data_2$genes
    merged_data_2 <- merged_data_2 %>% select(-"genes")  # Remove genes column
    
    # Convert to matrix for heatmap
    heatmap_matrix <- as.matrix(merged_data)
    heatmap_matrix_2 <- as.matrix(merged_data_2)
    
    # Generate heatmap
    h1 <- pheatmap(
      heatmap_matrix,
      cluster_rows = TRUE,
      cluster_cols = TRUE,
      scale = "column",
      color = colorRampPalette(c("blue", "white", "red"))(50),
      main = paste("Heatmap of", drug, "Model Correlations"),
      show_rownames = FALSE
    )
    
    h2 <- pheatmap(
      heatmap_matrix_2,
      cluster_rows = TRUE,
      cluster_cols = TRUE,
      scale = "column",
      color = colorRampPalette(c("blue", "white", "red"))(50),
      main = paste("Heatmap of", drug, "Model Correlations"),
      show_rownames = FALSE
    )
    ggsave(filename = paste0(plot_dir, "/", drug, "/heatmap_3_models.png"), plot = h1, width = 6, height = 4, dpi = 300)
    ggsave(filename = paste0(plot_dir, "/", drug, "/heatmap_2_models.png"), plot = h2, width = 6, height = 4, dpi = 300)
    message("heatmaps for comparing models saved at: ", plot_dir, "/", drug)
  }
}

generate_negative_positive_heatmaps <- function(path_positive_model, path_negative_model, drugs, plot_dir) {
  for (drug in drugs){
    
    lasso_corrs <- read.csv(paste0("results/GDSC2/negative/lasso/correlations/", drug, ".csv"), row.names = 1)
    en_corrs <- read.csv(paste0("results/GDSC2/positive/lasso/correlations/", drug, ".csv"), row.names = 1)
    
    # Convert row names back into a column for merging
    lasso_corrs <- lasso_corrs %>% mutate(genes = rownames(lasso_corrs))
    en_corrs <- en_corrs %>% mutate(genes = rownames(en_corrs))
    
    
    colnames(en_corrs)[colnames(en_corrs) == "x"] <- "POSITIVE"
    colnames(lasso_corrs)[colnames(lasso_corrs) == "x"] <- "NEGATIVE"
    
    
    merged_data <-reduce(list(lasso_corrs, en_corrs), full_join, by="genes")
    
    # Replace NA values with 0 (or another value if needed)
    merged_data[is.na(merged_data)] <- 0  # Fill missing values
    
    # Set 'genes' as row names for heatmap compatibility
    rownames(merged_data) <- merged_data$genes
    
    merged_data <- merged_data %>% select(-"genes")  # Remove genes column
    
    # Convert to matrix for heatmap
    heatmap_matrix <- as.matrix(merged_data)
    
    # Generate heatmap
    h1 <- pheatmap(
      heatmap_matrix,
      cluster_rows = TRUE,
      cluster_cols = TRUE,
      scale = "column",
      color = colorRampPalette(c("blue", "white"))(50),
      main = paste("Heatmap of", drug, "Model Correlations"),
      show_rownames = FALSE
    )
    
    ggsave(filename = paste0(plot_dir, "/", drug, "_heatmap_pos_neg.png"), plot = h1, width = 6, height = 4, dpi = 300)
    message("heatmaps for comparing positive/negative experiments saved at: ", plot_dir)
  }
}


