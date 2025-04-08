library(VennDiagram)
library(ggplot2)
source("src/heatmaps_plotting.R")

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




library(ggplot2)
library(dplyr)
library(tidyr)
library(RColorBrewer)
library(purrr)

### --- Helper Function: Load and Annotate Performance Data --- ###
load_performance_data <- function(size, config) {
  list(
    lasso_gdsc = load_cv_performance_and_annotate(
      file = paste0(config$results.cv.lasso.subdir, "GDSC2/lasso_performance_top_", size, "_features.csv"),
      method_name = "Lasso_gdsc_Pagerank", size = size
    ),
    lasso_gcsi = load_cv_performance_and_annotate(
      file = paste0(config$results.cv.lasso.subdir, "gCSI/lasso_performance_top_", size, "_features.csv"),
      method_name = "Lasso_gcsi_Pagerank", size = size
    ),
    drugbank = load_cv_performance_and_annotate(
      file = paste0(config$results.cv.drugbank.subdir, "drugbank_performance_top_", size, "_features.csv"),
      method_name = "Drugbank_Pagerank", size = size
    ),
    random = load_cv_performance_and_annotate(
      file = paste0(config$results.cv.random.subdir, "random_performance_", size, "_features.csv"),
      method_name = "Random", size = size
    )
  )
}

### --- Combine All Performance Data --- ###
combine_all_performances <- function(feature_sizes, config) {
  all_data <- feature_sizes %>%
    map(~ load_performance_data(.x, config)) %>%
    flatten() %>%
    bind_rows()
  return(all_data)
}

### --- Filter Valid Drugs --- ###
get_valid_drugs <- function(data) {
  data %>%
    filter(Method == "Drugbank_Pagerank") %>%
    pull(Drug) %>%
    unique()
}

### --- Plot Pearson Barplot per Feature Size --- ###
plot_pearson_bars <- function(data, size) {
  data %>%
    filter(FEATURE_SIZE == size) %>%
    ggplot(aes(x = Drug, y = PEARSON_Mean, fill = Method)) +
    geom_bar(stat = "identity", position = "dodge") +
    theme_minimal() +
    labs(
      title = paste("Pearson Correlation per Drug (Top", size, "Features)"),
      y = "Pearson Correlation", x = "Drug"
    ) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    scale_fill_brewer(palette = "Set2")
}

### --- Plot Pearson Boxplot per Feature Size --- ###
plot_pearson_boxplot <- function(data, size) {
  data %>%
    filter(FEATURE_SIZE == size) %>%
    ggplot(aes(x = Method, y = PEARSON_Mean, fill = Method)) +
    geom_boxplot() +
    theme_minimal() +
    labs(
      title = paste("Distribution of Pearson Correlation (Top", size, "Features)"),
      y = "Pearson Correlation", x = "Method"
    ) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
}

### --- Wilcoxon Test Wrapper --- ###
run_wilcoxon_tests <- function(data, size) {
  wide_data <- data %>%
    filter(FEATURE_SIZE == size) %>%
    select(Drug, Method, PEARSON_Mean) %>%
    pivot_wider(names_from = Method, values_from = PEARSON_Mean)
  
  cat("Wilcoxon test results (Top", size, "features):\n")
  list(
    Lasso_GDSC_vs_Drugbank = wilcox.test(wide_data$Lasso_gdsc_Pagerank, wide_data$Drugbank_Pagerank, paired = TRUE),
    Lasso_gCSI_vs_Drugbank = wilcox.test(wide_data$Lasso_gcsi_Pagerank, wide_data$Drugbank_Pagerank, paired = TRUE),
    Lasso_GDSC_vs_gCSI     = wilcox.test(wide_data$Lasso_gdsc_Pagerank, wide_data$Lasso_gcsi_Pagerank, paired = TRUE),
    Random_vs_Drugbank     = wilcox.test(wide_data$Random, wide_data$Drugbank_Pagerank, paired = TRUE),
    Lasso_GDSC_vs_Random   = wilcox.test(wide_data$Lasso_gdsc_Pagerank, wide_data$Random, paired = TRUE),
    Lasso_gCSI_vs_Random   = wilcox.test(wide_data$Lasso_gcsi_Pagerank, wide_data$Random, paired = TRUE)
  )
}
