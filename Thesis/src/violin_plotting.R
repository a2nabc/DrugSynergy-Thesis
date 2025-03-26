library(ggplot2)
library(dplyr)
library(pheatmap)


# Load the final results table
final_results <- read.csv("results/gCSI/negative/summary_results.csv")

# Convert Model to a factor for proper ordering
final_results$Model <- as.factor(final_results$Model)

# Define output directory for plots
plot_dir <- "figs/gCSI/negative"
if (!dir.exists(plot_dir)) {
  dir.create(plot_dir, recursive = TRUE)
}

# List of performance metrics to plot
metrics <- c("MSE", "RMSE", "R2", "PEARSON")
summary(final_results$PEARSON)
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
  
  cat("Saved", metric, "violin plot to", plot_dir, "\n")
}



######################################################################

library(reshape2)
library(dplyr)
library(tidyr)
library(ggplot2)
library(reshape2)
library(pheatmap)
library(purrr)
drugs <- readLines("./data/common_drugs.txt")
models <- c("lasso", "en", "ridge")

for (drug in drugs){
  # Define output directory for plots
  plot_directory <- paste0(plot_dir,"/", drug)
  if (!dir.exists(plot_directory)) {
    dir.create(plot_directory, recursive = TRUE)
  }

  lasso_corrs <- read.csv(paste0("results/gCSI/negative/lasso/correlations/", drug, ".csv"), row.names = 1)
  en_corrs <- read.csv(paste0("results/gCSI/negative/en/correlations/", drug, ".csv"), row.names = 1)
  ridge_corrs <- read.csv(paste0("results/gCSI/negative/ridge/correlations/", drug, ".csv"), row.names = 1)
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
  ggsave(filename = paste0(plot_dir, "/", drug, "/heatmap_3models.png"), plot = h1, width = 6, height = 4, dpi = 300)
  ggsave(filename = paste0(plot_dir, "/", drug, "/heatmap_2models.png"), plot = h2, width = 6, height = 4, dpi = 300)
  
}


