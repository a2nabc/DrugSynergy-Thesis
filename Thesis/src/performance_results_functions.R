library(VennDiagram)
library(ggplot2)
library(dplyr)
library(pheatmap)

# Function to plot a venn diagram with cell_lines across the datasets 
# ccl_data_list is a list with an entry for each screen
venn_plot <- function(ccl_data_list, output_file) {
  
  data_cell_lines <- list()
  for (file in ccl_data_list){
    dataset_name <- basename(dirname(file))
    sample_data <- read.csv(file)
    data_cell_lines[[dataset_name]] <- sample_data$sampleid
  }
  
  # Generate the Venn diagram
  venn.plot <- venn.diagram(
    x = data_cell_lines,
    category.names = names(data_cell_lines),
    filename = output_file,
    output = TRUE,
    fill = c("red", "blue", "green"), 
    alpha = 0.5,
    cat.col = c("red", "blue", "green"),
    cat.cex = 1.2,
    margin = 0.1
  )
  
  message("Venn diagram saved at: ", output_file)
}

generate_tables_one_experiment <- function(file_path, output_folder_path) {
  final_results <- read.csv(file_path)
  final_results <- final_results %>%
    arrange(Drug)
  
  drug_summary <- final_results %>%
    group_by(Drug) %>%
    summarise(across(c(MSE, RMSE, MAE, R2, PEARSON), mean, na.rm = TRUE))
  
  drug_summary <- drug_summary %>% arrange(desc(PEARSON))
  
  model_summary <- final_results %>%
    group_by(Model) %>%
    summarise(across(c(MSE, RMSE, MAE, R2, PEARSON), mean, na.rm = TRUE))
  
  model_summary <- model_summary %>% arrange(desc(PEARSON))  
  
  write.csv(final_results, paste0(output_folder_path, "summary_results_sorted.csv"))
  write.csv(drug_summary,  paste0(output_folder_path, "summary_per_drug.csv"))
  write.csv(model_summary, paste0(output_folder_path, "summary_per_model.csv"))
  
  message("Tables with summary of performance saved at: ", output_folder_path)
}

generate_tables_summary_performance <- function(csv_path, csv_name) {
  #We need to generate tables for both experiments: positive and negative_
  subfolders <- list.dirs(csv_path, recursive = FALSE)
  for (subfolder in subfolders) {
    file_path <- paste0(subfolder, "/", csv_name)
    output_folder_path <- paste0(subfolder, "/")
    generate_tables_one_experiment(file_path, output_folder_path)
  }
}

generate_one_violin <- function(csv_path, plot_dir) {
  final_results <- read.csv(csv_path)
  final_results$Model <- as.factor(final_results$Model)
  
  # Define output directory for plots
  if (!dir.exists(plot_dir)) {
    dir.create(plot_dir, recursive = TRUE)
  }
  
  # List of performance metrics to plot
  metrics <- c("MSE", "RMSE", "R2", "PEARSON")
  
  # Loop through metrics and create violin plots
  for (metric in metrics) {
    p <- ggplot(final_results, aes(x = Model, y = final_results[[metric]], fill = Model)) +
      geom_violin(trim = TRUE, alpha = 0.6) +  # Violin plot
      geom_boxplot(width = 0.1, fill = "white", outlier.shape = NA) +  # Boxplot overlay
      theme_minimal() +
      labs(title = paste("Model Performance:", metric), x = "Model Type", y = metric) +
      theme(legend.position = "none")  # Hide legend
    
    # Save plots as PNG and PDF
    ggsave(filename = paste0(plot_dir, "/", metric, "_violin_plot.png"), plot = p, width = 6, height = 4, dpi = 300)
    ggsave(filename = paste0(plot_dir, "/", metric, "_violin_plot.pdf"), plot = p, width = 6, height = 4)
    
    message("Saved", metric, "violin plot to", plot_dir, "\n")
  }
}

generate_violin_plots <- function(csv_path, csv_name, plot_dir) {
  #We need to generate tables for both experiments: positive and negative_
  subfolders <- list.dirs(csv_path, recursive = FALSE)
  for (subfolder in subfolders) {
    file_path <- paste0(subfolder, "/", csv_name)
    output_folder_path <- paste0(plot_dir, "/", basename(subfolder))
    
    generate_one_violin(file_path, output_folder_path)
  }
}

compute_jaccard_index <- function(set1, set2) {
  
  intersection_size <- length(intersect(set1, set2))
  union_size <- length(union(set1, set2))
  jaccard_index <- intersection_size / union_size
  
  return(jaccard_index)
}

#computes the jaccard index for positive and negative experiments given a screen
jaccard_index <- function(screen_path, features_subdir, drugs) {
  results <- data.frame(
    DRUG = character(),
    TYPE_OF_EXPERIMENT = character(),
    J_lasso_en = numeric(),
    J_lasso_ridge = numeric(),
    J_en_ridge = numeric(),
    stringsAsFactors = FALSE
  )
  
  subfolders <- list.dirs(screen_path, recursive = FALSE) #subfolders are positive and negative
  for (subfolder in subfolders) {
    #lasso
    lasso_features_path <- paste0(subfolder, "/lasso/", features_subdir)
    en_features_path <- paste0(subfolder, "/en/", features_subdir)
    ridge_features_path <- paste0(subfolder, "/ridge/", features_subdir)
    
    for (drug in drugs) {
      lasso_features <- readLines(paste0(lasso_features_path, drug, ".txt"))
      en_features <- readLines(paste0(en_features_path, drug, ".txt"))
      ridge_features <- readLines(paste0(ridge_features_path, drug, ".txt"))
      
      # lasso vs en
      lasso_en_index <- compute_jaccard_index(lasso_features, en_features)
      # lasso vs ridge
      lasso_ridge_index <- compute_jaccard_index(lasso_features, ridge_features)
      # ridge vs en
      ridge_en_index <- compute_jaccard_index(ridge_features, en_features)
      
      # Append results to the data frame
      results <- rbind(results, data.frame(
        DRUG = drug,
        TYPE_OF_EXPERIMENT = basename(subfolder),
        J_lasso_en = lasso_en_index,
        J_lasso_ridge = lasso_ridge_index,
        J_en_ridge = ridge_en_index
      ))
    }
  }
  
  results <- results %>% arrange(DRUG)
  
  # Generate the output file name based on the basename of the screen path
  output_filename <- paste0(screen_path, "/", basename(screen_path), "_jaccard.csv")
  
  # Save the results to a CSV file
  write.csv(results, file = output_filename, row.names = FALSE)
  
}


