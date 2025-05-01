library(VennDiagram)
library(ggplot2)
source("src/heatmaps_plotting.R")

# Function to plot a venn diagram with cell_lines across the datasets 
# ccl_data_list is a list with an entry for each screen
venn_plot <- function(ccl_data_list, output_file) {
  
  data_cell_lines <- ccl_data_list
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

# ============================== pagerank_analysis_utils.R ==============================

library(dplyr)
library(tidyr)
library(ggplot2)
library(RColorBrewer)

load_and_label <- function(size, dataset, method, subdir, prefix = "lasso_performance_top_") {
  file <- paste0(subdir, ifelse(dataset != "", paste0(dataset, "/"), ""), prefix, size, "_features.csv")
  return(load_cv_performance_and_annotate(file, method, size))
}

load_all_performance_data <- function(feature_sizes, config) {
  all_data <- list()
  for (size in feature_sizes) {
    all_data[[paste0("lasso_gdsc_", size)]] <- load_and_label(size, "GDSC2", "Lasso_gdsc_Pagerank", config$results.cv.lasso.subdir)
    all_data[[paste0("lasso_gcsi_", size)]] <- load_and_label(size, "gCSI", "Lasso_gcsi_Pagerank", config$results.cv.lasso.subdir)
    all_data[[paste0("drugbank_", size)]]   <- load_and_label(size, "", "Drugbank_Pagerank", config$results.cv.drugbank.subdir, "drugbank_performance_top_")
    all_data[[paste0("random_", size)]]     <- load_and_label(size, "", "Random", config$results.cv.random.subdir, prefix = "random_performance_")
  }
  return(bind_rows(all_data))
}

filter_valid_drugs <- function(data) {
  valid_drugs <- data %>%
    filter(Method == "Drugbank_Pagerank") %>%
    pull(Drug) %>%
    unique()
  return(data %>% filter(Drug %in% valid_drugs))
}

plot_pearson_barplot <- function(data, size, output_path, palette = "Set2", width = 10, height = 6) {
  p <- data %>%
    filter(FEATURE_SIZE == size) %>%
    ggplot(aes(x = Drug, y = PEARSON_Mean, fill = Method)) +
    geom_bar(stat = "identity", position = "dodge") +
    theme_minimal() +
    labs(title = paste("Pearson Correlation per Drug (Top", size, "Features)"),
         y = "Pearson Correlation", x = "Drug") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    scale_fill_brewer(palette = palette)
  ggsave(output_path, plot = p, width = width, height = height)
  return(p)
}

plot_pearson_boxplot <- function(data, size, output_path, width = 10, height = 6) {
  p <- data %>%
    filter(FEATURE_SIZE == size) %>%
    ggplot(aes(x = Method, y = PEARSON_Mean, fill = Method)) +
    geom_boxplot() +
    theme_minimal() +
    labs(title = paste("Distribution of Pearson Correlation (Top", size, "Features)"),
         y = "Pearson Correlation", x = "Method") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  ggsave(output_path, plot = p, width = width, height = height)
  return(p)
}

run_wilcoxon_tests <- function(data, size) {
  data_filtered <- data %>% filter(FEATURE_SIZE == size)
  pearson_wide <- data_filtered %>%
    select(Drug, Method, PEARSON_Mean) %>%
    pivot_wider(names_from = Method, values_from = PEARSON_Mean)
  
  methods <- c("Lasso_gdsc_Pagerank", "Lasso_gcsi_Pagerank", "Drugbank_Pagerank", "Random")
  
  comparisons <- combn(methods, 2, simplify = FALSE)
  results <- lapply(comparisons, function(pair) {
    test <- wilcox.test(pearson_wide[[pair[1]]], pearson_wide[[pair[2]]], paired = TRUE)
    data.frame(Comparison = paste(pair, collapse = " vs "), P_Value = test$p.value)
  })
  
  # Combine the results
  results_df <- bind_rows(results)
  
  # Apply FDR correction using the Benjamini-Hochberg method
  results_df$FDR_P_Value <- p.adjust(results_df$P_Value, method = "BH")
  
  return(results_df)
}
