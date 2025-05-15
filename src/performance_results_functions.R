library(VennDiagram)
library(ggplot2)
library(dplyr)
library(pheatmap)

compute_jaccard_index <- function(set1, set2) {
  
  intersection_size <- length(intersect(set1, set2))
  union_size <- length(union(set1, set2))
  jaccard_index <- intersection_size / union_size
  
  return(jaccard_index)
}

jaccard_db_dtc_lasso <- function(drugbank_genes, dtc_genes, lasso_genes, output_file) {
  # Initialize results data frame
  jaccard_results <- data.frame(
    DRUG = character(),
    J_DTC_DB = numeric(),
    J_DTC_Lasso = numeric(),
    J_DB_Lasso = numeric(),
    stringsAsFactors = FALSE
  )
  
  # Get the unique list of all drugs
  all_drugs <- unique(c(names(drugbank_genes), names(dtc_genes), names(lasso_genes)))
  print(all_drugs)
  for (drug in all_drugs) {
    drug_sets <- list(
      DrugBank = drugbank_genes[[drug]] %||% character(0),
      DTC = dtc_genes[[drug]] %||% character(0),
      Lasso = lasso_genes[[drug]] %||% character(0)
    )
    # Compute Jaccard indices for each pair
    j_dtc_db <- if (length(drug_sets$DTC) > 0 && length(drug_sets$DrugBank) > 0) {
      compute_jaccard_index(drug_sets$DTC, drug_sets$DrugBank)
    } else {
      NA
    }
    
    j_dtc_lasso <- if (length(drug_sets$DTC) > 0 && length(drug_sets$Lasso) > 0) {
      compute_jaccard_index(drug_sets$DTC, drug_sets$Lasso[[1]])
    } else {
      NA
    }
    
    j_db_lasso <- if (length(drug_sets$DrugBank) > 0 && length(drug_sets$Lasso) > 0) {
      compute_jaccard_index(drug_sets$DrugBank, drug_sets$Lasso[[1]])
    } else {
      NA
    }
    
    # Add results to the data frame
    jaccard_results <- rbind(jaccard_results, data.frame(
      DRUG = drug,
      J_DTC_DB = j_dtc_db,
      J_DTC_Lasso = j_dtc_lasso,
      J_DB_Lasso = j_db_lasso,
      stringsAsFactors = FALSE
    ))
  }
  
  # Save results to a CSV file
  write.csv(jaccard_results, output_file, row.names = FALSE)
  
  return(jaccard_results)
}


save_venn_plots_db_dtc_lasso <- function(drugbank_genes, dtc_genes, lasso_genes, name_output) {
  
  # Create a list of gene sets
  gene_sets <- list(
    Lasso = unique(lasso_genes),
    DTC = unique(dtc_genes),
    DrugBank = unique(drugbank_genes)
  )
  # Create the Venn diagram with custom colors
  venn.plot <- venn.diagram(
    x = gene_sets,
    category.names = c("Lasso", "DTC", "DrugBank"),
    filename = NULL,  # Don't specify a file to prevent immediate saving
    col = "black",    # Color of the border
    fill = c("red", "blue", "green"),  # Colors for each category (Lasso, DTC, DrugBank)
    alpha = 0.5,      # Transparency (0 is fully transparent, 1 is fully opaque)
    label.col = "black",  # Color of the labels
    cex = 2,             # Font size of the labels
    cat.col = c("red", "blue", "green"),  # Colors of the category labels
    cat.cex = 2,      # Font size of the category labels
    cat.dist = c(0.05, 0.05, 0.05),  # Category distance from the circle
    margin = 0.1       # Margin between the diagram and the edges of the plot
  )
  
  # Draw the plot
  png(name_output, width = 800, height = 800)
  grid.draw(venn.plot)
  dev.off()  # Close the device
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
  #We need to generate tables for both experiments: positive and negative
  subfolders <- list.dirs(csv_path, recursive = FALSE)
  for (subfolder in subfolders) {
    file_path <- paste0(subfolder, "/", csv_name)
    output_folder_path <- paste0(plot_dir, "/", basename(subfolder))
    
    generate_one_violin(file_path, output_folder_path)
  }
  print("DONE VIOLIN")
}
generate_violin_plots_all_in_one <- function(csv_path_screen1, csv_path_screen2, csv_name, plot_dir) {
  # Initialize an empty data frame to store results from both screens
  all_results <- data.frame()
  
  # Process the first screen (GDSC2)
  subfolders_screen1 <- list.dirs(csv_path_screen1, recursive = FALSE)
  for (subfolder in subfolders_screen1) {
    file_path <- paste0(subfolder, "/", csv_name)
    experiment_results <- read.csv(file_path)
    experiment_results$Screen <- "GDSC2"
    experiment_results$Experiment <- basename(subfolder)
    experiment_results$Condition <- ifelse(grepl("positive", subfolder, ignore.case = TRUE), "Positive", "Negative")
    all_results <- rbind(all_results, experiment_results)
  }
  
  # Process the second screen (gCSI)
  subfolders_screen2 <- list.dirs(csv_path_screen2, recursive = FALSE)
  for (subfolder in subfolders_screen2) {
    file_path <- paste0(subfolder, "/", csv_name)
    experiment_results <- read.csv(file_path)
    experiment_results$Screen <- "gCSI"
    experiment_results$Experiment <- basename(subfolder)
    experiment_results$Condition <- ifelse(grepl("positive", subfolder, ignore.case = TRUE), "Positive", "Negative")
    all_results <- rbind(all_results, experiment_results)
  }
  
  # Add grouping information
  all_results$Model <- as.factor(all_results$Model)
  all_results$Group <- paste(all_results$Model, all_results$Screen, all_results$Condition, sep = "_")
  
  # Define custom colors for models
  model_colors <- c(
    "lasso" = "blue",
    "en" = "green",
    "ridge" = "red"
  )
  
  # Add a color column based on Model
  all_results$ColorGroup <- all_results$Model
  
  # Ensure consistent ordering of groups for consecutive colors
  all_results$Group <- factor(all_results$Group, levels = unique(all_results$Group))
  
  # Add separators for clear grouping
  all_results$Group <- factor(all_results$Group, levels = c(
    "lasso_GDSC2_Negative", "en_GDSC2_Negative", "ridge_GDSC2_Negative",
    "---",  # Unique separator
    "lasso_GDSC2_Positive", "en_GDSC2_Positive", "ridge_GDSC2_Positive",
    "----",  # Unique separator
    "lasso_gCSI_Negative", "en_gCSI_Negative", "ridge_gCSI_Negative",
    "-----",  # Unique separator
    "lasso_gCSI_Positive", "en_gCSI_Positive", "ridge_gCSI_Positive"
  ))
  
  # List of performance metrics to plot
  metrics <- c("PEARSON", "RMSE")
  
  for (metric in metrics) {
    p <- ggplot(all_results, aes(x = Group, y = all_results[[metric]], fill = ColorGroup)) +
      geom_violin(trim = TRUE, alpha = 0.6) +
      geom_boxplot(width = 0.1, fill = "white", outlier.shape = NA) +
      scale_fill_manual(values = model_colors) +
      theme_minimal() +
      labs(title = paste("Model Performance:", metric), x = "Group", y = metric) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none") +
      geom_vline(xintercept = c(3.5, 6.5, 9.5), linetype = "dashed", color = "black")  # Add vertical separators
    
    # Save combined plot as PNG and PDF
    ggsave(filename = paste0(plot_dir, "/", metric, "_all_in_one_violin_plot.png"), plot = p, width = 12, height = 6, dpi = 300)
    ggsave(filename = paste0(plot_dir, "/", metric, "_all_in_one_violin_plot.pdf"), plot = p, width = 12, height = 6)
    
    message("Saved combined", metric, "violin plot to", plot_dir, "\n")
  }
  print("DONE VIOLIN")
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

# Generalized function to compute the Jaccard index for different methods given a screen
new_jaccard_index <- function(screen_path, features_subdir, drugs, methods) {
  results <- data.frame(
    DRUG = character(),
    TYPE_OF_EXPERIMENT = character(),
    stringsAsFactors = FALSE
  )
  
  # Dynamically create column names for pairwise comparisons
  method_combinations <- combn(methods, 2, simplify = FALSE)
  for (comb in method_combinations) {
    col_name <- paste0("J_", comb[1], "_", comb[2])
    results[[col_name]] <- numeric()
  }
  
  subfolders <- list.dirs(screen_path, recursive = FALSE) # subfolders are positive and negative
  for (subfolder in subfolders) {
    method_features <- list()
    
    # Read features for each method
    for (method in methods) {
      method_features[[method]] <- list()
      for (drug in drugs) {
        method_features[[method]][[drug]] <- readLines(paste0(subfolder, "/", method, "/", features_subdir, drug, ".txt"))
      }
    }
    
    for (drug in drugs) {
      row <- list(DRUG = drug, TYPE_OF_EXPERIMENT = basename(subfolder))
      
      # Compute Jaccard index for each pair of methods
      for (comb in method_combinations) {
        jaccard_index_value <- compute_jaccard_index(
          method_features[[comb[1]]][[drug]],
          method_features[[comb[2]]][[drug]]
        )
        col_name <- paste0("J_", comb[1], "_", comb[2])
        row[[col_name]] <- jaccard_index_value
      }
      
      # Append the row to the results data frame
      results <- rbind(results, as.data.frame(row, stringsAsFactors = FALSE))
    }
  }
  
  results <- results %>% arrange(DRUG)
  
  # Generate the output file name based on the basename of the screen path
  output_filename <- paste0(screen_path, "/", basename(screen_path), "_jaccard.csv")
  
  # Save the results to a CSV file
  write.csv(results, file = output_filename, row.names = FALSE)
}
