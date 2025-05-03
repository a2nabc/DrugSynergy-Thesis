library(glmnet)
library(PharmacoGx)
library(dplyr)
library(caret)
source("src/model_functions.R")
source("src/data_functions.R")
source("src/plots.R")
source("src/performance_results_functions.R")
#source("src/heatmaps_plotting.R")

############# 
# READ THIS:
# 		it is a VERY bad idea to use library(config)
# 		this will overwrite important baseline R functions
#		it is best practice to use config::get()
###########

args = commandArgs(trailingOnly = TRUE)
config.file <- args[1]

config<-config::get(file=config.file)
print(config)

## I suggest you to the 21 drugs that we have screened across all the datasets 
# (gCSI, GDSC, etc.)


# The code could look something like
#	Read in the data for the source/training screen
# 	Read in the data for the target/testing screen
# 	Either:
#		(a) Remove the common cell lines from the testing screen
#		(b) keep only the commmon cell lines from the testing screen
# 	Depending on if you want to test generalization or consistency. It would make sense
# 	to put a boolean for what kind of experiment in your config file.

# Then:
#	for each model:
#		for each drug:
#			train the model on the training screen
#			record the genes with non-zero coefficients
#			evaluate the model on the testing screen
#				predict(X_test, model)
#			for each of the genes with non-zero coefficient
#				record the correlation between that gene's expression
#				and the measure of response.
#			do the following 10 times:
#				split the test screen into a training subset (80%)
# 				and a testing subset (20%). 
#				train a linear regression (not lasso, simple regression)
#				on the training subset and evaluate on the testing subset/
#				record the performance
#					What we're trying to capture here is if the biomarkers you found
#					with lasso/enet/ridge actually generalize

config <- config::get(file = "src/configs/Broad_to_gCSI.yml")

# load source screen and target screen data
source_data <- load_data(config$source.expression, config$source.response)
target_data <- load_data(config$target.expression, config$target.response)

source_data <- load_data(config$source.expression, config$source.response)
target_data_GDSC2 <- load_data(config$target.expression, config$target.response)
target_data_gCSI <-  load_data(config$target.expression, config$target.response)

# Generate Venn plot for common cell lines
ccl_list <- list(
  GDSC2 = unique(target_data_GDSC2$expression$Cell_line),
  gCSI = unique(target_data_gCSI$expression$Cell_line),
  CTRPv2_CCLE = unique(source_data$expression$Cell_line)
)
venn_plot(ccl_list, paste0(config$figs.dir, "venn_plot_cell_lines.png"))

# Generate Venn plot for common drugs
drug_list <- list(
  GDSC2 = unique(target_data_GDSC2$response$Drug),
  gCSI = unique(target_data_gCSI$response$Drug),
  CTRPv2_CCLE = unique(source_data$response$Drug)
)
venn_plot(drug_list, paste0(config$figs.dir, "venn_plot_drugs.png"))

# filter source and target just with genes from KEGG and LINCS
keep_genes <- get_gene_union(config$gene.sets, config$gene.sets.path)
keep_genes <- keep_genes[keep_genes %in% colnames(source_data$expression)]

writeLines(keep_genes, "data/background_genes.txt")

source_data$expression <- source_data$expression %>% select(c("Cell_line", keep_genes))
target_data$expression <- target_data$expression %>% select(c("Cell_line", keep_genes))

# filter cell lines across the two screens depending on the type of experiment
target_data <- filter_cell_lines(config$experiment.type.is.positive, source_data, target_data)

# Print dataset sizes after filtering
#cat("Training Data:", nrow(source_data$expression), "samples\n")
#cat("Testing Data:", nrow(target_data$expression), "samples\n")


# Read drug universe and initialize list and folders to store results
common_drugs <- readLines(config$drugs) 
results.subdir <- create_output_folders(config$target.screen, config$models, config$experiment.type.is.positive, config$results.dir, config$results.features.dir, config$results.eval.dir, config$results.correlations.dir, config$results.cv.random.subdir, config$results.cv.lasso.subdir, config$results.cv.lasso.gdsc, config$results.cv.lasso.gcsi, config$results.cv.drugbank.subdir,  config$results.cv.dir)
data_dimensions_metadata <- create_df_metadata(config$results.dir, config$results.metadata.dir)
model_results <- list()
debug <- list()
summary_results <- data.frame(Drug = character(), Model = character(), MSE = numeric(), RMSE = numeric(), R2 = numeric(), Biomarkers = character(), stringsAsFactors = FALSE)


for (model.type in config$models) {
  for (drug in common_drugs) {
    output_features <- paste0(config$results.dir, results.subdir, model.type, "/", config$results.features.dir, drug, ".txt")
    output_eval <- paste0(config$results.dir, results.subdir,  model.type, "/", config$results.eval.dir, drug, ".csv")
    output_corrs <- paste0(config$results.dir, results.subdir, model.type, "/", config$results.correlations.dir, drug, ".csv")
    
    results <- train.model.and.get.results(source_data, drug, model.type, output_features, target_data, output_eval, output_corrs)
    
    summary_results <- rbind(summary_results, results$new_row_summary_results)
    model_results[[drug]][[model.type]] <- results$new_row_model_results
    data_dimensions_metadata <- rbind(data_dimensions_metadata, results$new_row_dimensions)
  }
}

write.csv(data_dimensions_metadata, paste0(config$results.dir, results.subdir, config$results.metadata.dir, "data_dimensions.csv"))

write.csv(summary_results, paste0(config$results.dir, results.subdir, "summary_results.csv"))

############################# RESULTS AND PLOTS TO COMPARE LASSO, EN AND RIDGE ###############################
ccl_list <- list()
for (screen in config$screens) {
  ccl_list <- rbind(ccl_list, paste0(config$path.to.processed.screens, screen, "/", config$name.ccl.file))
}
venn_plot(ccl_list, paste0(config$figs.dir, "venn_plot.png"))

for (screen in config$target.screens){
  generate_tables_summary_performance(paste0(config$results.dir, screen), "summary_results.csv") 
  generate_violin_plots(paste0(config$results.dir, screen), "summary_results.csv", paste0(config$figs.dir, screen))
  
  generate_model_heatmaps(paste0(config$results.dir, screen), config$results.correlations.dir, common_drugs, paste0(config$figs.dir, screen))
 # jaccard_index(paste0(config$results.dir, screen), config$results.features.dir, common_drugs)
  screen_path <- paste0(config$results.dir, screen)
  features_subdir <- config$results.features.dir
  methods <- c('lasso', 'en', 'ridge')
  new_jaccard_index(screen_path, features_subdir, common_drugs, methods)
}

generate_violin_plots_all_in_one(paste0(config$results.dir, "GDSC2"), paste0(config$results.dir, "gCSI"), "summary_results.csv", paste0(config$figs.dir, screen))

############################################################# K-FOLD CROSS VALIDATION FOR ALL THE GENES IN THE DATA ########################
# Perform 10-fold CV for both screens (GDSC2 and gCSI) using all genes in the original data
results_all_genes <- list()

for (screen in config$target.screens) {
  data_expression <- paste0("./data/processed/", screen, "/expression.csv")
  data_response <- paste0("./data/processed/", screen, "/responses.csv")
  data <- load_data(data_expression, data_response)
  
  num_folds <- 10
  for (drug in common_drugs) {
    result_drug <- k.fold.linear.model(data, drug, num_folds)
    results_all_genes[[paste0(screen, "_", drug)]] <- result_drug
  }
}

# Combine results into a single data frame
results_all_genes_df <- do.call(rbind, results_all_genes)

# Save results to a CSV file
write.csv(results_all_genes_df, file = paste0(config$results.dir, "/cv_performance/PC_genes/", "cv_results_all_genes.csv"), row.names = FALSE)

############################################################# K-FOLD CROSS VALIDATION FOR GIVEN GENE SETS ######################## --> OK

results_biological_priors_drugbank <- list()
results_biological_priors_dtc <- list()


# Training on Drugbank/DTC drug targets
for (screen in config$target.screens) {
  data_expression <- paste0("./data/processed/", screen, "/expression.csv")
  data_response <- paste0("./data/processed/", screen, "/responses.csv")
  data <- load_data(data_expression, data_response)
  results_biological_priors_drugbank <- run.cv.gene.set.per.drug(results_biological_priors_drugbank, screen, common_drugs, "drugbank", config$drugbank.gene.set, data, config, 10)
  results_biological_priors_dtc <- run.cv.gene.set.per.drug(results_biological_priors_dtc, screen, common_drugs, "dtc", config$dtc.gene.set, data, config, 10)
}
write.results.drug.targets(results_biological_priors_drugbank, config$target.screens, "drugbank", config)
write.results.drug.targets(results_biological_priors_dtc, config$target.screens, "dtc", config)
# Read the data
df1 <- read.csv("./results/cv_performance/drugbank/gCSI/performance_drugbank_drug_targets.csv")
df2 <- read.csv("./results/cv_performance/drugbank/GDSC2/performance_drugbank_drug_targets.csv")
df3 <- read.csv("./results/cv_performance/dtc/gCSI/performance_dtc_drug_targets.csv")
df4 <- read.csv("./results/cv_performance/dtc/GDSC2/performance_dtc_drug_targets.csv")

# Add a column to identify the dataset
df1$Dataset <- "Drugbank_gCSI"
df2$Dataset <- "Drugbank_GDSC2"
df3$Dataset <- "DTC_gCSI"
df4$Dataset <- "DTC_GDSC2"

# Combine the data frames
combined_df <- bind_rows(df1, df2, df3, df4)

combined_df12 <- bind_rows(df1, df2)
combined_df34 <- bind_rows(df3, df4)

# Violin plot for PEARSON_Mean
ggplot(combined_df, aes(x = Dataset, y = PEARSON_Mean, fill = Dataset)) +
  geom_violin(trim = TRUE, alpha = 0.6) +    # Violin plot
  geom_boxplot(width = 0.1, fill = "white", outlier.shape = NA) + # Boxplot overlay
  theme_minimal() +
  labs(title = "Distribution of PEARSON across all drugs",
       x = "Gene set and training dataset",
       y = "PEARSON mean values after 10-fold cv") +
  theme(legend.position = "none")   # Hide legend if you don't need it
ggsave(filename = "figs/DRUGBANK_DTC/pearson_mean.png", width = 8, height = 6)


# Violin plot for RMSE_Mean
ggplot(combined_df12, aes(x = Dataset, y = RMSE_Mean, fill = Dataset)) +
  geom_violin(trim = TRUE, alpha = 0.6) +    # Violin plot
  geom_boxplot(width = 0.1, fill = "white", outlier.shape = NA) + # Boxplot overlay
  theme_minimal() +
  labs(title = "Distribution of RMSE across all drugs",
       x = "Gene set and training dataset",
       y = "RMSE mean values after 10-fold cv") +
  theme(legend.position = "none")   # Hide legend if you don't need it
ggsave(filename = "figs/DRUGBANK_DTC/rmse_mean_drugbank.png", width = 8, height = 6)

# Violin plot for RMSE_Mean
ggplot(combined_df34, aes(x = Dataset, y = RMSE_Mean, fill = Dataset)) +
  geom_violin(trim = TRUE, alpha = 0.6) +    # Violin plot
  geom_boxplot(width = 0.1, fill = "white", outlier.shape = NA) + # Boxplot overlay
  theme_minimal() +
  labs(title = "Distribution of RMSE across all drugs",
       x = "Gene set and training dataset",
       y = "RMSE mean values after 10-fold cv") +
  theme(legend.position = "none")   # Hide legend if you don't need it
ggsave(filename = "figs/DRUGBANK_DTC/rmse_mean_dtc.png", width = 8, height = 6)


# Load DrugBank data
drugbank_data <- read.delim("results/drugbank/AffectedGenesByDrug.txt", stringsAsFactors = FALSE)
drugbank_genes <- lapply(split(drugbank_data$GENES_AFFECTED, drugbank_data$DRUG_NAME), function(x) unlist(strsplit(x, ", ")))

# Load DTC data
dtc_data <- read.csv("results/dtc/AffectedGenesByDrug.csv", stringsAsFactors = FALSE)
dtc_genes <- lapply(split(dtc_data$gene_names, dtc_data$compound_name), function(x) unlist(strsplit(x, ", ")))

# Load Lasso data
lasso_files <- list.files("results/GDSC2/positive/lasso/features", full.names = TRUE)
lasso_genes <- lapply(lasso_files, function(file) list(readLines(file)))
names(lasso_genes) <- gsub(".txt", "", basename(lasso_files))

# Function to extract all unique genes from lasso_genes and save to a file
save_unique_lasso_genes <- function(lasso_genes, output_file) {
  # Flatten the list and extract unique genes
  all_genes <- unique(unlist(lasso_genes))
  
  # Write the unique genes to the specified file
  writeLines(all_genes, output_file)
}

# Call the function
save_unique_lasso_genes(lasso_genes, "data/lasso_gene_set.txt")


jaccard_db_dtc_lasso(drugbank_genes, dtc_genes, lasso_genes, "results/jaccard/dtc_drugbank_jaccard_indices.csv")

# Load gene sets
lasso_genes <- readLines("data/lasso_gene_set.txt")
dtc_genes <- readLines("data/dtc_gene_set.txt")
drugbank_genes <- readLines("data/drugbank_gene_set.txt")
save_venn_plots_db_dtc_lasso(drugbank_genes, dtc_genes, lasso_genes, "./figs/DRUGBANK_DTC/venn_plot_gene_sets_colored.png")


############################################################### 10-FOLD CV FOR CENTRALIZED MEASURES ##########################################################
results_central_genes <- list()
central_gene_sets <- list(
  betweenness = config$betweenness.gene.set,
  eigenvector = config$eigenvector.gene.set,
  degree = config$degree.gene.set,
  uniform_pagerank = config$uniform.pagerank.gene.set
)

for (screen in config$target.screens) {
  data_expression <- paste0("./data/processed/", screen, "/expression.csv")
  data_response <- paste0("./data/processed/", screen, "/responses.csv")
  data <- load_data(data_expression, data_response)
  
  for (gene_set_name in names(central_gene_sets)) {
    gene_set <- central_gene_sets[[gene_set_name]]
    results_central_genes <- run.cv.gene.set.per.drug(
      results_central_genes, 
      screen, 
      common_drugs, 
      gene_set_name, 
      gene_set, 
      data, 
      config, 
      10
    )
  }
}
write.results.central.genes(results_central_genes, config$target.screens, central_gene_sets, config)

# Combine all centrality data into a single data frame
combined_centrality <- bind_rows(
  betweenness_gCSI, degree_gCSI, eigenvector_gCSI, pagerank_gCSI,
  betweenness_GDSC, degree_GDSC, eigenvector_GDSC, pagerank_GDSC
)

# Add a column to identify the centrality measure and dataset
combined_centrality <- combined_centrality %>%
  mutate(Method = factor(Method, levels = c("Betweenness", "Degree", "Eigenvector", "Pagerank")),
         Screen = factor(Screen, levels = c("gCSI", "GDSC2")))

# Filter data for valid PEARSON_Mean and RMSE_Mean values
filtered_centrality <- combined_centrality %>%
  filter(is.finite(PEARSON_Mean), is.finite(RMSE_Mean))

# Generate a single-axis violin plot for PEARSON_Mean
pearson_plot <- ggplot(filtered_centrality, aes(x = interaction(Screen, Method), y = PEARSON_Mean, fill = Method)) +
  geom_violin(trim = TRUE, alpha = 0.6) +
  geom_boxplot(width = 0.1, fill = "white", outlier.shape = NA) +
  theme_minimal() +
  labs(title = "PEARSON_Mean for Centrality Measures Across Screens",
       x = "Screen and Centrality Measure",
       y = "PEARSON Mean") +
  theme(legend.position = "bottom") +
  scale_fill_brewer(palette = "Dark2")
  scale_x_discrete(labels = function(x) gsub("\\.", "\n", x))  # Format x-axis labels

# Save the PEARSON_Mean plot
ggsave(filename = "figs/CENTRALITY_MEASURES/Combined_Violin_PEARSON_Mean_SingleAxis_centrality.png",
       plot = pearson_plot, width = 12, height = 8)

# Generate a single-axis violin plot for RMSE_Mean
rmse_plot <- ggplot(filtered_centrality, aes(x = interaction(Screen, Method), y = RMSE_Mean, fill = Method)) +
  geom_violin(trim = TRUE, alpha = 0.6) +
  geom_boxplot(width = 0.1, fill = "white", outlier.shape = NA) +
  theme_minimal() +
  labs(title = "RMSE_Mean for Centrality Measures Across Screens",
       x = "Screen and Centrality Measure",
       y = "RMSE Mean") +
  theme(legend.position = "bottom") +
  scale_fill_brewer(palette = "Dark2") +
  scale_x_discrete(labels = function(x) gsub("\\.", "\n", x))  # Format x-axis labels

# Save the RMSE_Mean plot
ggsave(filename = "figs/CENTRALITY_MEASURES/Combined_Violin_RMSE_Mean_SingleAxis_centrality.png",
       plot = rmse_plot, width = 12, height = 8)


############ now we do the same for pathway analysis ###########

# Load DrugBank data
drugbank_files <- list.files("results/pathway_enrichment/drugbank", full.names = TRUE, pattern = "\\.txt$")
drugbank_pathways <- lapply(drugbank_files, function(file) list(readLines(file)))
names(drugbank_pathways) <- gsub(".txt", "", basename(drugbank_files))

# Load DTC data
dtc_files <- list.files("results/pathway_enrichment/dtc", full.names = TRUE, pattern = "\\.txt$")
dtc_pathways <- lapply(dtc_files, function(file) list(readLines(file)))
names(dtc_pathways) <- gsub(".txt", "", basename(dtc_files))

# Load Lasso data
lasso_files <- list.files("results/pathway_enrichment/lasso_features", full.names = TRUE, pattern = "\\.txt$")
lasso_pathways <- lapply(lasso_files, function(file) list(readLines(file)))
names(lasso_pathways) <- gsub(".txt", "", basename(lasso_files))

jaccard_db_dtc_lasso(drugbank_pathways, dtc_pathways, lasso_pathways, "results/jaccard/dtc_drugbank_pathways_jaccard_indices.csv")

# for the venn plot we have to merge all pathways (eliminate separation per drug)
merged_drugbank_pathways <- unique(unlist(drugbank_pathways))
merged_dtc_pathways <- unique(unlist(dtc_pathways))
merged_lasso_pathways <- unique(unlist(lasso_pathways))

save_venn_plots_db_dtc_lasso(merged_drugbank_pathways, merged_dtc_pathways, merged_lasso_pathways, "./figs/DRUGBANK_DTC/venn_plot_pathways_colored.png")




##################################################################### DRUGBANK AND DTC AND RANDOM ##################################################

# Initialize
results <- init.results(results, c("lasso", "random", "drugbank", "dtc"), config$target.screens, feature_sizes)

# Lasso
for (screen in config$target.screens) {
  data_expression <- paste0("./data/processed/", screen, "/expression.csv")
  data_response <- paste0("./data/processed/", screen, "/responses.csv")
  data <- load_data(data_expression, data_response)
  results <- run.lasso.cv(results, screen, feature_sizes, common_drugs, source_data, data, config, 10)
}
write.lasso.results(results, config$target.screens, feature_sizes, config)

# Random + Drugbank + DTC
for (screen in config$target.screens){
  data_expression <- paste0("./data/processed/", screen, "/expression.csv")
  data_response <- paste0("./data/processed/", screen, "/responses.csv")
  data <- load_data(data_expression, data_response)
  results <- run.other.cv(results, feature_sizes, common_drugs, screen, data, config, 10)
  write.other.results(results, feature_sizes, config, screen)
}

#ultima prova. pagerank all genes:
for (screen in config$target.screens){
  data_expression <- paste0("./data/processed/", screen, "/expression.csv")
  data_response <- paste0("./data/processed/", screen, "/responses.csv")
  data <- load_data(data_expression, data_response)
  
  genes_to_keep_dtc <- load_pagerank_feature_list("./results/pagerank_output/dtc/pagerank_dtc_50.txt")
  genes_to_keep_dtc <- genes_to_keep_dtc[genes_to_keep_dtc %in% colnames(data$expression)]
  genes_to_keep_drugbank <- load_pagerank_feature_list("./results/pagerank_output/drugbank/pagerank_db_50.txt")
  genes_to_keep_drugbank <- genes_to_keep_drugbank[genes_to_keep_drugbank %in% colnames(data$expression)]
  
  # Ensure filtered_1 and filtered_2 are initialized as lists
  filtered_1 <- list()
  filtered_2 <- list()
  
  # Correctly subset the data
  filtered_1$expression <- data$expression[, c("Cell_line", genes_to_keep_dtc), drop = FALSE]
  filtered_1$response <- data$response
  filtered_2$expression <- data$expression[, c("Cell_line", genes_to_keep_drugbank), drop = FALSE]
  filtered_2$response <- data$response
  
  num_folds <- 10
  for (drug in common_drugs) {
    
    result_drug_dtc <- k.fold.linear.model(filtered_1, drug, num_folds)
    results[[paste0("dtc_", screen, "_all")]] <- rbind(results[[paste0("dtc_", screen, "_all")]], result_drug_dtc)
    
    result_drug_drugbank <- k.fold.linear.model(filtered_2, drug, num_folds)
    results[[paste0("drugbank_", screen, "_all")]] <- rbind(results[[paste0("drugbank_", screen, "_all")]], result_drug_drugbank)
  }
  
  
  
  write.other.results(results, feature_sizes, config, screen)
}


############################# PLOTS FOR PAGERANK VS RANDOM PERFORMANCES #############################
methods_to_plot <- c("random", "lasso", "dtc", "drugbank")
library(readr)
library(patchwork)

# Initialize lists to store data
all_pearson_data <- list()
all_rmse_data <- list()

# Loop through each method and collect data
for (method in methods_to_plot) {
  # Safely load the CSV file and handle potential errors
  csv_path <- paste0("../../resum_top50_", method, ".csv")
  if (file.exists(csv_path)) {
    method_data <- tryCatch({
      read_csv(csv_path) %>%
        pivot_longer(cols = starts_with("RMSE_Mean"), names_to = "Screen", values_to = "RMSE_Mean") %>%
        pivot_longer(cols = starts_with("PEARSON_Mean"), names_to = "Screen_Pearson", values_to = "PEARSON_Mean") %>%
        filter(if ("PEARSON_Mean" %in% colnames(.)) is.finite(PEARSON_Mean) else FALSE,
               if ("RMSE_Mean" %in% colnames(.)) is.finite(RMSE_Mean) else FALSE) %>%
        mutate(Screen = ifelse(grepl("GDSC", Screen_Pearson), "GDSC", "gCSI"),
               Method = method)
    }, error = function(e) {
      message(paste("Error reading file:", csv_path, " - ", e$message))
      NULL
    })
    
    if (!is.null(method_data)) {
      all_pearson_data[[method]] <- method_data %>% select(Screen, PEARSON_Mean, Method)
      all_rmse_data[[method]] <- method_data %>% select(Screen, RMSE_Mean, Method)
    }
  } else {
    message(paste("File does not exist:", csv_path))
  }
}

# Combine all data into single data frames
combined_pearson_data <- bind_rows(all_pearson_data)
combined_rmse_data <- bind_rows(all_rmse_data)

# PEARSON plot
pearson_plot <- ggplot(combined_pearson_data, aes(x = interaction(Screen, Method), y = PEARSON_Mean, fill = Method)) +
  geom_violin(trim = TRUE, alpha = 0.6) +
  geom_boxplot(width = 0.1, fill = "white", outlier.shape = NA) +
  theme_minimal() +
  labs(title = "PEARSON_Mean for All Methods",
       x = "Screen and Method",
       y = "PEARSON Mean") +
  theme(legend.position = "bottom") +
  scale_fill_brewer(palette = "Set3") +
  scale_x_discrete(labels = function(x) gsub("\\.", "\n", x))  # Format x-axis labels

ggsave(filename = "figs/LASSO_PAGERANK1_PAGERANK2/Combined_Violin_PEARSON_Mean_SingleAxis.png",
       plot = pearson_plot, width = 12, height = 8)

# RMSE plot
rmse_plot <- ggplot(combined_rmse_data, aes(x = interaction(Screen, Method), y = RMSE_Mean, fill = Method)) +
  geom_violin(trim = TRUE, alpha = 0.6) +
  geom_boxplot(width = 0.1, fill = "white", outlier.shape = NA) +
  theme_minimal() +
  labs(title = "RMSE_Mean for All Methods",
       x = "Screen and Method",
       y = "RMSE Mean") +
  theme(legend.position = "bottom") +
  scale_fill_brewer(palette = "Set3") +
  scale_x_discrete(labels = function(x) gsub("\\.", "\n", x))  # Format x-axis labels

ggsave(filename = "figs/LASSO_PAGERANK1_PAGERANK2/Combined_Violin_RMSE_Mean_SingleAxis.png",
       plot = rmse_plot, width = 12, height = 8)




#################################################### COMPARISONNNNNNNNNNNNNNNNN 
# Load required libraries
library(ggplot2)
library(dplyr)
library(readr)

# Define file paths
file_paths <- list(
  Lasso_derived_features = "results/GDSC2/positive/summary_results.csv",
  Drugbank_features = "results/cv_performance/drugbank/GDSC2/performance_drugbank_drug_targets.csv",
  Pagerank_features = "results/cv_performance/central_genes/GDSC2/performance_uniform_pagerank_central_genes.csv",
  Lasso_Pagerank_features = "../../resum_top50_lasso.csv",
  Drugbank_Pagerank_features = "../../resum_top50_drugbank.csv",
  Random_features = "../../resum_top50_random.csv"
)
# Load and process data
data_list <- lapply(names(file_paths), function(name) {
  file <- file_paths[[name]]
  if (file.exists(file)) {
    data <- read_csv(file, show_col_types = FALSE)
    if (name == "Lasso_derived_features") {
      if ("PEARSON" %in% colnames(data)) {
        data <- data %>%
          select(Drug, PEARSON_Mean = PEARSON) %>%
          filter(is.finite(PEARSON_Mean)) %>%
          mutate(Method = name)
      } else {
        data <- tibble()  # Return an empty tibble if PEARSON is missing
      }
    } else if (name %in% c("Drugbank_features", "Pagerank_features")) {
      if ("PEARSON_Mean" %in% colnames(data)) {
        data <- data %>%
          select(Drug, PEARSON_Mean) %>%
          filter(is.finite(PEARSON_Mean)) %>%
          mutate(Method = name)
      } else {
        data <- tibble()  # Return an empty tibble if PEARSON_Mean is missing
      }
    } else if (name %in% c("Lasso_Pagerank_features", "Drugbank_Pagerank_features", "Random_features")) {
      if ("PEARSON_Mean_GDSC" %in% colnames(data)) {
        data <- data %>%
          select(Drug, PEARSON_Mean = PEARSON_Mean_GDSC) %>%
          filter(is.finite(PEARSON_Mean)) %>%
          mutate(Method = name)
      } else {
        data <- tibble()  # Return an empty tibble if PEARSON_Mean_GDSC is missing
      }
    } else {
      data <- tibble()  # Return an empty tibble for unsupported formats
    }
  } else {
    message(paste("File does not exist:", file))
    data <- tibble()  # Return an empty tibble if the file does not exist
  }
  return(data)
})

# Combine all data into a single data frame
combined_data <- bind_rows(data_list)

# Ensure combined_data has the expected structure
if (nrow(combined_data) == 0) {
  stop("No valid data was loaded. Please check the file paths and data structure.")
}

# Generate violin and boxplot for PEARSON_Mean
plot <- ggplot(combined_data, aes(x = Method, y = PEARSON_Mean, fill = Method)) +
  geom_violin(trim = TRUE, alpha = 0.6) +
  geom_boxplot(width = 0.1, fill = "white", outlier.shape = NA) +
  theme_minimal() +
  labs(title = "Distribution of PEARSON_Mean Across Methods",
       x = "Method",
       y = "PEARSON Mean") +
  theme(legend.position = "none") +
  scale_fill_brewer(palette = "Set3")

# Save the plot
ggsave(filename = "figs/Violin_Boxplot_PEARSON_Mean_Distributions.png", plot = plot, width = 10, height = 6)




################33 END COMPARISON

# Plots for drug bank / DTC data
generate_violin_plots(paste0(config$results.dir, screen), "summary_results.csv", paste0(config$figs.dir, screen))


# Jaccard distances to compare the feature sets in results/pagerank_output/
pagerank_output_path <- "results/pagerank_output/"
methods <- c("lasso", "en", "ridge")
subfolders <- c("GDSC2", "gCSI")
subsubfolders <- c("positive", "negative")
feature_suffix <- "_50.txt"

# Get the list of drugs from the drugbank folder
drug_files <- list.files(paste0(pagerank_output_path, "drugbank/"), pattern = feature_suffix, full.names = TRUE)
drug_names <- gsub(feature_suffix, "", basename(drug_files))

# Initialize a data frame to store pairwise Jaccard distances
pairwise_jaccard_results <- data.frame(Drug = character(), Method1 = character(), Subfolder1 = character(), Subsubfolder1 = character(),
                     Method2 = character(), Subfolder2 = character(), Subsubfolder2 = character(), Jaccard = numeric(), stringsAsFactors = FALSE)
method2 <- "drugbank"
# Compare pairwise Jaccard distances for each drug
for (drug in drug_names) {
  for (method1 in methods) {
  for (subfolder1 in subfolders) {
    for (subsubfolder1 in subsubfolders) {  # We got the ML loop, we just compare it with drugbank data
       
        file1 <- paste0(pagerank_output_path, method1, "/", subfolder1, "/", subsubfolder1, "/", drug, feature_suffix)
        file2 <- paste0(pagerank_output_path, method2, "/", drug, feature_suffix)
        
        if (file.exists(file1) && file.exists(file2)) {
          features1 <- readLines(file1)
          features2 <- readLines(file2)
          
          jaccard <- compute_jaccard_index(features1, features2)
          
          pairwise_jaccard_results <- rbind(pairwise_jaccard_results, data.frame(
          Drug = drug,
          Method1 = method1,
          Subfolder1 = subfolder1,
          Subsubfolder1 = subsubfolder1,
          Method2 = method2,
          Jaccard = jaccard,
          stringsAsFactors = FALSE
          ))
        }
    }
  }
  }
}

# Save pairwise Jaccard results to a CSV file
write.csv(pairwise_jaccard_results, file = paste0(pagerank_output_path, "pairwise_jaccard_distances.csv"), row.names = FALSE)

# Generate a heatmap for pairwise Jaccard distances
plot_pairwise_jaccard_heatmap(pairwise_jaccard_results, paste0(pagerank_output_path, "pairwise_jaccard_heatmap.png"))

# First we group the lasso_GDSC, lasso_GCSI, drugbank and random performances in a single table:
# Load all data
combined2 <- load_all_performance_data(feature_sizes, config)
filtered_combined2 <- filter_valid_drugs(combined2)

# Prepare combined_data for Wilcoxon input
prepare_combined_data <- function(combined_data, file_paths) {
  # Initialize an empty data frame to store the results
  prepared_data <- data.frame(
    Drug = character(),
    MSE_Mean = numeric(),
    MSE_SD = numeric(),
    RMSE_Mean = numeric(),
    RMSE_SD = numeric(),
    MAE_Mean = numeric(),
    MAE_SD = numeric(),
    R2_Mean = numeric(),
    R2_SD = numeric(),
    PEARSON_Mean = numeric(),
    PEARSON_SD = numeric(),
    Method = character(),
    stringsAsFactors = FALSE
  )
  
  # Iterate through each method in file_paths
  for (method_name in names(file_paths)) {
    file_path <- file_paths[[method_name]]
    
    # Safely load the data
    if (file.exists(file_path)) {
      method_data <- read.csv(file_path)
      
      # Check if required columns exist
      if (method_name == "Lasso_derived_features") {
        if ("PEARSON" %in% colnames(method_data)) {
          method_data <- method_data %>%
            select(Drug, PEARSON_Mean = PEARSON) %>%
            mutate(
              MSE_Mean = NA,
              MSE_SD = NA,
              RMSE_Mean = NA,
              RMSE_SD = NA,
              MAE_Mean = NA,
              MAE_SD = NA,
              R2_Mean = NA,
              R2_SD = NA,
              PEARSON_SD = NA,
              Method = method_name
            )
          
          # Append to the prepared_data
          prepared_data <- rbind(prepared_data, method_data)
        }
      } else if (method_name %in% c("Drugbank_features", "Pagerank_features")) {
        if (all(c("Drug", "PEARSON_Mean") %in% colnames(method_data))) {
          method_data <- method_data %>%
            mutate(
              MSE_Mean = ifelse("MSE_Mean" %in% colnames(.), MSE_Mean, NA),
              MSE_SD = ifelse("MSE_SD" %in% colnames(.), MSE_SD, NA),
              RMSE_Mean = ifelse("RMSE_Mean" %in% colnames(.), RMSE_Mean, NA),
              RMSE_SD = ifelse("RMSE_SD" %in% colnames(.), RMSE_SD, NA),
              MAE_Mean = ifelse("MAE_Mean" %in% colnames(.), MAE_Mean, NA),
              MAE_SD = ifelse("MAE_SD" %in% colnames(.), MAE_SD, NA),
              R2_Mean = ifelse("R2_Mean" %in% colnames(.), R2_Mean, NA),
              R2_SD = ifelse("R2_SD" %in% colnames(.), R2_SD, NA),
              PEARSON_SD = ifelse("PEARSON_SD" %in% colnames(.), PEARSON_SD, NA)
            ) %>%
            mutate(Method = method_name) %>%
            select(
              Drug,
              MSE_Mean,
              MSE_SD,
              RMSE_Mean,
              RMSE_SD,
              MAE_Mean,
              MAE_SD,
              R2_Mean,
              R2_SD,
              PEARSON_Mean,
              PEARSON_SD,
              Method
            )
          
          # Append to the prepared_data
          prepared_data <- rbind(prepared_data, method_data)
        }
      } else if (method_name %in% c("Lasso_Pagerank_features", "Drugbank_Pagerank_features", "Random_features")) {
        if ("PEARSON_Mean_GDSC" %in% colnames(method_data)) {
          method_data <- method_data %>%
            select(Drug, PEARSON_Mean = PEARSON_Mean_GDSC) %>%
            mutate(
              MSE_Mean = NA,
              MSE_SD = NA,
              RMSE_Mean = NA,
              RMSE_SD = NA,
              MAE_Mean = NA,
              MAE_SD = NA,
              R2_Mean = NA,
              R2_SD = NA,
              PEARSON_SD = NA,
              Method = method_name
            )
          
          # Append to the prepared_data
          prepared_data <- rbind(prepared_data, method_data)
        }
      }
    }
  }
  
  return(prepared_data)
}

# Call the function to prepare the data
prepared_combined_data <- prepare_combined_data(combined_data, file_paths)

# Save the prepared data to a CSV file
write.csv(prepared_combined_data, "results/prepared_combined_data.csv", row.names = FALSE)

wilcoxon <- run_wilcoxon_tests_no_size(prepared_combined_data)

# Barplots
######   plot_pearson_barplot(filtered_combined, 10, "figs/Barplot_10_features.png", "RdYlGn")
######   plot_pearson_barplot(filtered_combined, 20,  "figs/Barplot_20_features.png", "RdBu")
######   plot_pearson_barplot(filtered_combined, 50,  "figs/Barplot_50_features.png", "PuOr")
######   
######   # Boxplots
######   plot_pearson_boxplot(filtered_combined, 10, "figs/Boxplot_10_features.png")
######   plot_pearson_boxplot(filtered_combined, 20, "figs/Boxplot_20_features.png")
######   plot_pearson_boxplot(filtered_combined, 50, "figs/Boxplot_50_features.png")

# Wilcoxon test results
wilcoxon_results <- run_wilcoxon_tests(prepared_combined_data, 50)
write.csv(wilcoxon_results, file = "results/cv_performance/wilcoxon_fdr_results.csv", row.names = FALSE)
print(wilcoxon_results)
