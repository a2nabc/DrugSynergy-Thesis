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


############################################################# K-FOLD CROSS VALIDATION FOR GIVEN GENE SETS ######################## --> OK

results_lasso <- list()
results_en <- list()
results_ridge <- list()

load_set_from_folder <- function(folder_path, suffix = ".txt") {
  # Initialize an empty set to store the union of all elements
  result_set <- character()
  
  # List all files in the folder with the specified suffix
  files <- list.files(folder_path, pattern = paste0(suffix, "$"), full.names = TRUE)
  
  # Iterate through each file and add its contents to the result set
  for (file in files) {
    file_contents <- readLines(file)
    result_set <- union(result_set, file_contents)
  }
  
  return(result_set)
}

path_lasso <- "results/GDSC2/positive/lasso/features"
path_en <- "results/GDSC2/positive/en/features"
path_ridge <- "results/GDSC2/positive/ridge/features"

# Load the union of all features from the lasso folder
lasso_feature_set <- load_set_from_folder(path_lasso)

# Save the result to a .txt file in the same folder
output_file_lasso <- file.path(path_lasso, "lasso_feature_set.txt")
writeLines(lasso_feature_set, output_file_lasso)

# Load the union of all features from the ridge folder
ridge_feature_set <- load_set_from_folder(path_ridge)

# Save the result to a .txt file in the same folder
output_file_ridge <- file.path(path_ridge, "ridge_feature_set.txt")
writeLines(ridge_feature_set, output_file_ridge)

# Load the union of all features from the en folder
en_feature_set <- load_set_from_folder(path_en)

# Save the result to a .txt file in the same folder
output_file_en <- file.path(path_en, "en_feature_set.txt")
writeLines(en_feature_set, output_file_en)

# Training on lasso/en/ridge drug targets using the saved feature sets
for (screen in config$target.screens) {
  data_expression <- paste0("./data/processed/", screen, "/expression.csv")
  data_response <- paste0("./data/processed/", screen, "/responses.csv")
  data <- load_data(data_expression, data_response)
  
  # Load the saved feature sets
  lasso_gene_set <- "results/GDSC2/positive/lasso/features/"
  ridge_gene_set <- "results/GDSC2/positive/ridge/features/"
  en_gene_set <- "results/GDSC2/positive/en/features/"
  
  # Run cross-validation for each feature set
  results_lasso <- run.cv.gene.set.per.drug(results_lasso, screen, common_drugs, "lasso", lasso_gene_set, data, config, 10)
  results_ridge <- run.cv.gene.set.per.drug(results_ridge, screen, common_drugs, "ridge", ridge_gene_set, data, config, 10)
  results_en <- run.cv.gene.set.per.drug(results_en, screen, common_drugs, "en", en_gene_set, data, config, 10)
}

# Write results for each feature set

write.results.lasso.en.ridge(results_lasso, config$target.screens, "lasso", config)
write.results.lasso.en.ridge(results_ridge, config$target.screens, "ridge", config)
write.results.lasso.en.ridge(results_en, config$target.screens, "en", config)

################################################# PLOT DRUG-WISE PERFORMANCE ##################################
# Define a function to generate the drug-wise PEARSON barplot
generate_pearson_barplot <- function(folder_paths, output_file) {
  # Load required libraries
  library(ggplot2)
  library(dplyr)
  
  # Initialize an empty data frame to store results
  pearson_data <- data.frame(Drug = character(), PEARSON_Mean = numeric(), Model = character(), stringsAsFactors = FALSE)
  
  # Loop through each folder and calculate the mean of eval_metrics.PEARSON
  for (model in names(folder_paths)) {
    folder_path <- folder_paths[[model]]
    csv_files <- list.files(folder_path, pattern = "\\.csv$", full.names = TRUE)
    
    for (file in csv_files) {
      # Extract the drug name from the file name
      drug_name <- gsub("\\.csv$", "", basename(file))
      
      # Read the CSV file
      data <- read.csv(file)
      
      # Check if the column exists and calculate the mean
      if ("eval_metrics.PEARSON" %in% colnames(data)) {
        pearson_mean <- mean(data$eval_metrics.PEARSON, na.rm = TRUE)
        pearson_data <- rbind(pearson_data, data.frame(
          Drug = drug_name,
          PEARSON_Mean = pearson_mean,
          Model = model
        ))
      }
    }
  }
  
  # Create the barplot without vertical lines
  plot <- ggplot(pearson_data, aes(x = Drug, y = PEARSON_Mean, fill = Model)) +
    geom_bar(stat = "identity", position = "dodge") +
    theme_minimal() +
    labs(title = "Mean PEARSON Values for Each Drug and Model",
         x = "Drug",
         y = "Mean PEARSON",
         fill = "Model") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))  # Rotate x-axis labels for better readability
  
  # Save the plot
  ggsave(filename = output_file, plot = plot, width = 12, height = 6)
}

# Call the function with the specified folder paths and output file
folder_paths <- list(
  lasso = "/home/anna/Desktop/BHK/Thesis/results/GDSC2/negative/lasso/evaluation",
  en = "/home/anna/Desktop/BHK/Thesis/results/GDSC2/negative/en/evaluation",
  ridge = "/home/anna/Desktop/BHK/Thesis/results/GDSC2/negative/ridge/evaluation"
)
output_file <- "/home/anna/Desktop/BHK/Thesis/results/gCSI/positive/evaluation/pearson_barplot_combined.png"

generate_pearson_barplot(folder_paths, output_file)

#################################################




results_biological_priors_drugbank <- list()
results_biological_priors_dtc <- list()


# Training on Drugbank/DTC drug targets
for (screen in config$target.screens) {
  data_expression <- paste0("./data/processed/", screen, "/expression.csv")
  data_response <- paste0("./data/processed/", screen, "/responses.csv")
  data <- load_data(data_expression, data_response)
  folder_drugbank <- "results/drugbank/features/"
  folder_dtc <- "results/dtc/features/"
  results_biological_priors_drugbank <- run.cv.gene.set.per.drug(results_biological_priors_drugbank, screen, common_drugs, "drugbank", folder_drugbank, data, config, 10)
  results_biological_priors_dtc <- run.cv.gene.set.per.drug(results_biological_priors_dtc, screen, common_drugs, "dtc", folder_dtc, data, config, 10)
}
write.results.drug.targets(results_biological_priors_drugbank, config$target.screens, "drugbank", config)
write.results.drug.targets(results_biological_priors_dtc, config$target.screens, "dtc", config)


################################################# PLOT DRUG-WISE PERFORMANCE ##################################

plot_pearson_distribution_filtering(results_biological_priors_drugbank, "GDSC2", "drugbank", "figs/DRUGBANK_DTC/pearson_distribution_drugbank_GDSC_drugwise.png")
plot_pearson_distribution_filtering(results_biological_priors_drugbank, "gCSI", "drugbank", "figs/DRUGBANK_DTC/pearson_distribution_drugbank_gCSI_drugwise.png")
plot_pearson_distribution_filtering(results_biological_priors_dtc, "GDSC2", "dtc", "figs/DRUGBANK_DTC/pearson_distribution_dtc_GDSC_drugwise.png")
plot_pearson_distribution_filtering(results_biological_priors_dtc, "gCSI", "dtc", "figs/DRUGBANK_DTC/pearson_distribution_dtc_gCSI_drugwise.png")


# Call the function with the specified folder paths and output file

###############################################################################################################

# Read the data
df1 <- read.csv("./results/cv_performance/drugbank/gCSI/performance_drugbank_drug_targets.csv")
df2 <- read.csv("./results/cv_performance/drugbank/GDSC2/performance_drugbank_drug_targets.csv")
df3 <- read.csv("./results/cv_performance/dtc/gCSI/performance_dtc_drug_targets.csv")
df4 <- read.csv("./results/cv_performance/dtc/GDSC2/performance_dtc_drug_targets.csv")
#df5 <- read.csv("./results/cv_performance/random_pagerank_input/gCSI/random_performance_50_features.csv")
#df6 <- read.csv("./results/cv_performance/random_pagerank_input/GDSC2/random_performance_50_features.csv")
df5 <- read.csv("results/cv_performance/ridge/gCSI/performance_ridge_lasso_en_ridge.csv")
df6 <- read.csv("results/cv_performance/ridge/GDSC2/performance_ridge_lasso_en_ridge.csv")

# Add a column to identify the dataset
df1$Dataset <- "Drugbank_gCSI"
df2$Dataset <- "Drugbank_GDSC2"
df3$Dataset <- "DTC_gCSI"
df4$Dataset <- "DTC_GDSC2"
df5$Dataset <- "All_features_gCSI"
df6$Dataset <- "All_features_GDSC2"

# Combine the data frames
combined_df <- bind_rows(df1, df2, df3, df4, df5, df6)

# Define colors for gCSI and GDSC2
dataset_colors <- c("gCSI" = "#1f77b4", "GDSC2" = "#ff7f0e")

# Violin plot for PEARSON_Mean
ggplot(combined_df, aes(x = Dataset, y = PEARSON_Mean, fill = ifelse(grepl("gCSI", Dataset), "gCSI", "GDSC2"))) +
  geom_violin(trim = TRUE, alpha = 0.6) +    # Violin plot
  geom_boxplot(width = 0.1, fill = "white", outlier.shape = NA) + # Boxplot overlay
  theme_minimal() +
  scale_fill_manual(values = dataset_colors) +
  labs(title = "Distribution of PEARSON across all drugs",
       x = "Gene set and training dataset",
       y = "PEARSON mean values after 10-fold cv",
       fill = "Dataset") +
  theme(legend.position = "bottom")   # Show legend
ggsave(filename = "figs/DRUGBANK_DTC/pearson_mean_colored_all.png", width = 8, height = 6)

# Violin plot for RMSE_Mean
ggplot(combined_df, aes(x = Dataset, y = RMSE_Mean, fill = ifelse(grepl("gCSI", Dataset), "gCSI", "GDSC2"))) +
  geom_violin(trim = TRUE, alpha = 0.6) +    # Violin plot
  geom_boxplot(width = 0.1, fill = "white", outlier.shape = NA) + # Boxplot overlay
  theme_minimal() +
  scale_fill_manual(values = dataset_colors) +
  labs(title = "Distribution of RMSE across all drugs",
       x = "Gene set and training dataset",
       y = "RMSE mean values after 10-fold cv",
       fill = "Dataset") +
  theme(legend.position = "bottom")   # Show legend
ggsave(filename = "figs/DRUGBANK_DTC/rmse_mean_colored_all.png", width = 8, height = 6)



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
#write.results.central.genes(results_central_genes, config$target.screens, central_gene_sets, config)
write.results.centrality(results_central_genes, config$screens, c("betweenness", "degree", "uniform_pagerank", "eigenvector"), config) 


################################################# PLOT DRUG-WISE PERFORMANCE ##################################

plot_pearson_distribution_filtering(results_central_genes, "GDSC2", "betweenness", "figs/CENTRALITY_MEASURES/pearson_distribution_betweenness_GDSC_drugwise.png")
plot_pearson_distribution_filtering(results_central_genes, "gCSI", "betweenness", "figs/CENTRALITY_MEASURES/pearson_distribution_betweenness_gCSI_drugwise.png")
plot_pearson_distribution_filtering(results_central_genes, "GDSC2", "degree", "figs/CENTRALITY_MEASURES/pearson_distribution_degree_GDSC_drugwise.png")
plot_pearson_distribution_filtering(results_central_genes, "gCSI", "degree", "figs/CENTRALITY_MEASURES/pearson_distribution_degree_gCSI_drugwise.png")
plot_pearson_distribution_filtering(results_central_genes, "GDSC2", "uniform_pagerank", "figs/CENTRALITY_MEASURES/pearson_distribution_uniform_pagerank_GDSC_drugwise.png")
plot_pearson_distribution_filtering(results_central_genes, "gCSI", "uniform_pagerank", "figs/CENTRALITY_MEASURES/pearson_distribution_uniform_pagerank_gCSI_drugwise.png")
plot_pearson_distribution_filtering(results_central_genes, "GDSC2", "eigenvector", "figs/CENTRALITY_MEASURES/pearson_distribution_eigenvector_GDSC_drugwise.png")
plot_pearson_distribution_filtering(results_central_genes, "gCSI", "eigenvector", "figs/CENTRALITY_MEASURES/pearson_distribution_eigenvector_gCSI_drugwise.png")


################################################### 10-FOLD CV FOR PPR ##########################################################
results_ppr_lasso <- list()
results_ppr_db <- list()
results_ppr_dtc <- list()

#training on feature sets drug-wise
for (screen in config$target.screens) {
  data_expression <- paste0("./data/processed/", screen, "/expression.csv")
  data_response <- paste0("./data/processed/", screen, "/responses.csv")
  data <- load_data(data_expression, data_response)
  
  # Load the saved feature sets
  lasso_ppr_gene_set <- "results/pagerank_output/lasso/GDSC2/positive"
  drugbank_ppr_gene_set <- "results/pagerank_output/drugbank/drugwise"
  dtc_ppr_gene_set <- "results/pagerank_output/dtc/drugwise"
  
  # Run cross-validation for each feature set
  results_ppr_lasso <- run.cv.gene.set.per.drug.ppr(results_ppr_lasso, screen, common_drugs, "ppr_lasso", lasso_ppr_gene_set, data, config, 10)
  results_ppr_db <- run.cv.gene.set.per.drug.ppr(results_ppr_db, screen, common_drugs, "ppr_db", drugbank_ppr_gene_set, data, config, 10)
  results_ppr_dtc <- run.cv.gene.set.per.drug.ppr(results_ppr_dtc, screen, common_drugs, "ppr_dtc", dtc_ppr_gene_set, data, config, 10)
}

write.results.ppr(results_ppr_lasso, config$target.screens, "ppr_lasso", config)
write.results.ppr(results_ppr_db, config$target.screens, "ppr_db", config)
write.results.ppr(results_ppr_dtc, config$target.screens, "ppr_dtc", config)

################################################# PLOT DRUG-WISE PERFORMANCE ##################################
plot_pearson_distribution_filtering(results_ppr_lasso, "GDSC2", "ppr_lasso", "figs/PPR/pearson_distribution_ppr_lasso_GDSC_drugwise.png")
plot_pearson_distribution_filtering(results_ppr_lasso, "gCSI", "ppr_lasso", "figs/PPR/pearson_distribution_ppr_lasso_gCSI_drugwise.png")
plot_pearson_distribution_filtering(results_ppr_db, "GDSC2", "ppr_db", "figs/PPR/pearson_distribution_ppr_db_GDSC_drugwise.png")
plot_pearson_distribution_filtering(results_ppr_db, "gCSI", "ppr_db", "figs/PPR/pearson_distribution_ppr_db_gCSI_drugwise.png")
plot_pearson_distribution_filtering(results_ppr_dtc, "GDSC2", "ppr_dtc", "figs/PPR/pearson_distribution_ppr_dtc_GDSC_drugwise.png")
plot_pearson_distribution_filtering(results_ppr_dtc, "gCSI", "ppr_dtc", "figs/PPR/pearson_distribution_ppr_dtc_gCSI_drugwise.png")


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

# Load additional data from ../../resums/resum_ridge_LM.csv
# Load additional data from ../../resums/resum_random50.csv
ridge_data <- tryCatch({
  read.csv("../../resums/resum_random50.csv")
}, error = function(e) {
  stop("Error loading ridge_data: ", e$message)
})

# Check if required columns exist
if (!all(c("PEARSON_Mean_GDSC", "PEARSON_Mean_gCSI") %in% colnames(ridge_data))) {
  stop("The required columns 'PEARSON_Mean_GDSC' and 'PEARSON_Mean_gCSI' are missing in the ridge_data.")
}

ridge_data_gdsc <- ridge_data %>%
  select(PEARSON_Mean = PEARSON_Mean_GDSC) %>%
  mutate(Method = "Random 50", Screen = "GDSC2")
ridge_data_gcsi <- ridge_data %>%
  select(PEARSON_Mean = PEARSON_Mean_gCSI) %>%
  mutate(Method = "Random 50", Screen = "gCSI")

# Combine ridge data with centrality data
filtered_centrality <- bind_rows(filtered_centrality, ridge_data_gdsc, ridge_data_gcsi)

# Generate a single-axis violin plot for PEARSON_Mean
pearson_plot <- ggplot(filtered_centrality, aes(x = interaction(Screen, Method), y = PEARSON_Mean, fill = Method)) +
  geom_violin(trim = TRUE, alpha = 0.6) +
  geom_boxplot(width = 0.1, fill = "white", outlier.shape = NA) +
  theme_minimal() +
  labs(title = "PEARSON_Mean for Centrality Measures and Random sets across screens",
       x = "Screen and Method",
       y = "PEARSON Mean") +
  theme(legend.position = "bottom") +
  scale_fill_brewer(palette = "Dark2") +
  scale_x_discrete(labels = function(x) gsub("\\.", "\n", x))  # Format x-axis labels

# Save the PEARSON_Mean plot
ggsave(filename = "figs/CENTRALITY_MEASURES/Combined_Violin_PEARSON_Mean_SingleAxis_centrality_rnd.png",
       plot = pearson_plot, width = 12, height = 8)

# Generate a single-axis violin plot for RMSE_Mean
rmse_plot <- ggplot(filtered_centrality, aes(x = interaction(Screen, Method), y = RMSE_Mean, fill = Method)) +
  geom_violin(trim = TRUE, alpha = 0.6) +
  geom_boxplot(width = 0.1, fill = "white", outlier.shape = NA) +
  theme_minimal() +
  labs(title = "RMSE_Mean for Centrality Measures and Random sets across screens",
       x = "Screen and Method",
       y = "RMSE Mean") +
  theme(legend.position = "bottom") +
  scale_fill_brewer(palette = "Dark2") +
  scale_x_discrete(labels = function(x) gsub("\\.", "\n", x))  # Format x-axis labels

# Save the RMSE_Mean plot
ggsave(filename = "figs/CENTRALITY_MEASURES/Combined_Violin_RMSE_Mean_SingleAxis_centrality_rnd.png",
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

################### now a venn plot only for relevant pathways

# Load DrugBank data
drugbank_files_cancer <- list.files("results/pathway_enrichment/drugbank", full.names = TRUE, pattern = "\\_cancer.txt$")
drugbank_pathways_cancer <- lapply(drugbank_files_cancer, function(file) list(readLines(file)))
names(drugbank_pathways_cancer) <- gsub("_Cancer.txt", "", basename(drugbank_files_cancer))

# Load DTC data
dtc_files_cancer <- list.files("results/pathway_enrichment/dtc", full.names = TRUE, pattern = "\\_cancer.txt$")
dtc_pathways_cancer <- lapply(dtc_files_cancer, function(file) list(readLines(file)))
names(dtc_pathways_cancer) <- gsub("_cancer.txt", "", basename(dtc_files_cancer))

# Load Lasso data
lasso_files_cancer <- list.files("results/pathway_enrichment/lasso_features", full.names = TRUE, pattern = "\\_cancer.txt$")
lasso_pathways_cancer <- lapply(lasso_files_cancer, function(file) list(readLines(file)))
names(lasso_pathways_cancer) <- gsub("_cancer.txt", "", basename(lasso_files_cancer))

merged_drugbank_pathways_cancer <- unique(unlist(drugbank_pathways_cancer))
merged_dtc_pathways_cancer <- unique(unlist(dtc_pathways_cancer))
merged_lasso_pathways_cancer <- unique(unlist(lasso_pathways_cancer))

save_venn_plots_db_dtc_lasso(merged_drugbank_pathways_cancer, merged_dtc_pathways_cancer, merged_lasso_pathways_cancer, "./figs/DRUGBANK_DTC/venn_plot_pathways_colored_cancer.png")

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


resultats_lasso_derived <- init.results.wo.methods.and.sizes("lasso", resultats_lasso_derived, config$target.screens)
resultats_ridge_derived <- init.results.wo.methods.and.sizes("ridge", resultats_ridge_derived, config$target.screens)
resultats_en_derived <- init.results.wo.methods.and.sizes("en", resultats_en_derived, config$target.screens)

for (screen in config$target.screens) {
  data_expression <- paste0("./data/processed/", screen, "/expression.csv")
  data_response <- paste0("./data/processed/", screen, "/responses.csv")
  data <- load_data(data_expression, data_response)
  
  resultats_lasso_derived <- run.lasso.derived.cv("lasso", resultats_lasso_derived, screen, common_drugs, source_data, data, 10)
  resultats_ridge_derived <- run.lasso.derived.cv("ridge", resultats_ridge_derived, screen, common_drugs, source_data, data, 10)
  resultats_en_derived <- run.lasso.derived.cv("en", resultats_en_derived, screen, common_drugs, source_data, data, 10)
}

write.lasso.derived.results("lasso", resultats_lasso_derived, config$target.screens, config)
write.lasso.derived.results("ridge", resultats_ridge_derived, config$target.screens, config)
write.lasso.derived.results("en", resultats_en_derived, config$target.screens, config)


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



##############################3 VIOLKINS PPR
# Load required libraries
library(ggplot2)
library(dplyr)
library(readr)
# Define file paths
file_paths <- list(
  Lasso_Pagerank_features = "../../resums/resum_lasso_PPR.csv",
  Drugbank_Pagerank_features = "../../resums/resum_drugbank_PPR.csv",
  DTC_Pagerank_features = "../../resums/resum_dtc_PPR.csv",
  all_features = "../../resums/resum_ridge_LM.csv",
  Random_features = "../../resums/resum_random50.csv"
)

# Load and process data for GDSC
data_list_gdsc <- lapply(names(file_paths), function(name) {
  file <- file_paths[[name]]
  if (file.exists(file)) {
    data <- read_csv(file, show_col_types = FALSE)
    if ("PEARSON_Mean_GDSC" %in% colnames(data)) {
      data <- data %>%
        select(Drug, PEARSON_Mean = PEARSON_Mean_GDSC) %>%
        filter(is.finite(PEARSON_Mean)) %>%
        mutate(Method = name)
    } else {
      data <- tibble()  # Return an empty tibble if PEARSON_Mean_GDSC is missing
    }
  } else {
    message(paste("File does not exist:", file))
    data <- tibble()  # Return an empty tibble if the file does not exist
  }
  return(data)
})

combined_data_gdsc <- bind_rows(data_list_gdsc)

# Generate violin and boxplot for GDSC
plot_gdsc <- ggplot(combined_data_gdsc, aes(x = Method, y = PEARSON_Mean, fill = Method)) +
  geom_violin(trim = TRUE, alpha = 0.6) +
  geom_boxplot(width = 0.1, fill = "white", outlier.shape = NA) +
  theme_minimal() +
  labs(title = "Distribution of PEARSON_Mean Across Methods (GDSC)",
       x = "Method",
       y = "PEARSON Mean") +
  theme(legend.position = "none") +
  scale_fill_brewer(palette = "Set3")

# Save the GDSC plot
ggsave(filename = "figs/PPR/Violin_Boxplot_PEARSON_Mean_Distributions_GDSC.png", plot = plot_gdsc, width = 10, height = 6)

# Load and process data for gCSI
data_list_gcsi <- lapply(names(file_paths), function(name) {
  file <- file_paths[[name]]
  if (file.exists(file)) {
    data <- read_csv(file, show_col_types = FALSE)
    if ("PEARSON_Mean_gCSI" %in% colnames(data)) {
      data <- data %>%
        select(Drug, PEARSON_Mean = PEARSON_Mean_gCSI) %>%
        filter(is.finite(PEARSON_Mean)) %>%
        mutate(Method = name)
    } else {
      data <- tibble()  # Return an empty tibble if PEARSON_Mean_gCSI is missing
    }
  } else {
    message(paste("File does not exist:", file))
    data <- tibble()  # Return an empty tibble if the file does not exist
  }
  return(data)
})

combined_data_gcsi <- bind_rows(data_list_gcsi)
# Generate violin and boxplot for gCSI
plot_gcsi <- ggplot(combined_data_gcsi, aes(x = Method, y = PEARSON_Mean, fill = Method)) +
  geom_violin(trim = TRUE, alpha = 0.6) +
  geom_boxplot(width = 0.1, fill = "white", outlier.shape = NA) +
  theme_minimal() +
  labs(title = "Distribution of PEARSON_Mean Across Methods (gCSI)",
       x = "Method",
       y = "PEARSON Mean") +
  theme(legend.position = "none") +
  scale_fill_brewer(palette = "Set3")

# Save the gCSI plot
ggsave(filename = "figs/PPR/Violin_Boxplot_PEARSON_Mean_Distributions_gCSI.png", plot = plot_gcsi, width = 10, height = 6)


# Load required libraries
library(ggplot2)
library(dplyr)
library(readr)
# Define file paths
file_paths <- list(
  Lasso_Pagerank_features = "../../resums/resum_lasso_PPR.csv",
  Drugbank_Pagerank_features = "../../resums/resum_drugbank_PPR.csv",
  DTC_Pagerank_features = "../../resums/resum_dtc_PPR.csv",
  all_features = "../../resums/resum_ridge_LM.csv",
  Random_features = "../../resums/resum_random50.csv"
)

# Load and process data for GDSC
data_list_gdsc <- lapply(names(file_paths), function(name) {
  file <- file_paths[[name]]
  if (file.exists(file)) {
    data <- read_csv(file, show_col_types = FALSE)
    if ("PEARSON_Mean_GDSC" %in% colnames(data)) {
      data <- data %>%
        select(Drug, PEARSON_Mean = PEARSON_Mean_GDSC) %>%
        filter(is.finite(PEARSON_Mean)) %>%
        mutate(Method = name)
    } else {
      data <- tibble()  # Return an empty tibble if PEARSON_Mean_GDSC is missing
    }
  } else {
    message(paste("File does not exist:", file))
    data <- tibble()  # Return an empty tibble if the file does not exist
  }
  return(data)
})

combined_data_gdsc <- bind_rows(data_list_gdsc)

# Generate violin and boxplot for GDSC
plot_gdsc <- ggplot(combined_data_gdsc, aes(x = Method, y = PEARSON_Mean, fill = Method)) +
  geom_violin(trim = TRUE, alpha = 0.6) +
  geom_boxplot(width = 0.1, fill = "white", outlier.shape = NA) +
  theme_minimal() +
  labs(title = "Distribution of PEARSON_Mean Across Methods (GDSC)",
       x = "Method",
       y = "PEARSON Mean") +
  theme(legend.position = "none") +
  scale_fill_brewer(palette = "Set3")

# Save the GDSC plot
ggsave(filename = "figs/PPR/Violin_Boxplot_PEARSON_Mean_Distributions_GDSC.png", plot = plot_gdsc, width = 10, height = 6)

# Load and process data for gCSI
data_list_gcsi <- lapply(names(file_paths), function(name) {
  file <- file_paths[[name]]
  if (file.exists(file)) {
    data <- read_csv(file, show_col_types = FALSE)
    if ("PEARSON_Mean_gCSI" %in% colnames(data)) {
      data <- data %>%
        select(Drug, PEARSON_Mean = PEARSON_Mean_gCSI) %>%
        filter(is.finite(PEARSON_Mean)) %>%
        mutate(Method = name)
    } else {
      data <- tibble()  # Return an empty tibble if PEARSON_Mean_gCSI is missing
    }
  } else {
    message(paste("File does not exist:", file))
    data <- tibble()  # Return an empty tibble if the file does not exist
  }
  return(data)
})

combined_data_gcsi <- bind_rows(data_list_gcsi)
# Generate violin and boxplot for gCSI
plot_gcsi <- ggplot(combined_data_gcsi, aes(x = Method, y = PEARSON_Mean, fill = Method)) +
  geom_violin(trim = TRUE, alpha = 0.6) +
  geom_boxplot(width = 0.1, fill = "white", outlier.shape = NA) +
  theme_minimal() +
  labs(title = "Distribution of PEARSON_Mean Across Methods (gCSI)",
       x = "Method",
       y = "PEARSON Mean") +
  theme(legend.position = "none") +
  scale_fill_brewer(palette = "Set3")

# Save the gCSI plot
ggsave(filename = "figs/PPR/Violin_Boxplot_PEARSON_Mean_Distributions_gCSI.png", plot = plot_gcsi, width = 10, height = 6)

############################3violins for other files

# Load required libraries
library(ggplot2)
library(dplyr)
library(readr)

# Define file paths for ridge, lasso, and en
file_paths <- list(
  EN = "../../resums/resum_en_derived.csv",
  Lasso = "../../resums/resum_lasso_derived.csv",
  Ridge = "../../resums/resum_ridge_derived.csv",
  All_features = "../../resums/resum_ridge_LM.csv"
)

# Initialize lists to store data
data_list <- list()

# Loop through each method and collect data
for (method in names(file_paths)) {
  file <- file_paths[[method]]
  if (file.exists(file)) {
    method_data <- read_csv(file, show_col_types = FALSE) %>%
      select(Drug, PEARSON_Mean_GDSC, PEARSON_Mean_gCSI) %>%
      pivot_longer(cols = starts_with("PEARSON_Mean"), names_to = "Screen", values_to = "PEARSON_Mean") %>%
      filter(is.finite(PEARSON_Mean)) %>%
      mutate(Screen = ifelse(grepl("GDSC", Screen), "GDSC", "gCSI"),
             Method = method)
    data_list[[method]] <- method_data
  } else {
    message(paste("File does not exist:", file))
  }
}

# Combine all data into a single data frame
combined_data <- bind_rows(data_list)

# Generate violin plot for PEARSON_Mean
violin_plot <- ggplot(combined_data, aes(x = interaction(Screen, Method), y = PEARSON_Mean, fill = Method)) +
  geom_violin(trim = TRUE, alpha = 0.6) +
  geom_boxplot(width = 0.1, fill = "white", outlier.shape = NA) +
  theme_minimal() +
  labs(title = "PEARSON_Mean for Ridge, Lasso, and ElasticNet",
       x = "Screen and Method",
       y = "PEARSON Mean") +
  theme(legend.position = "bottom") +
  scale_fill_brewer(palette = "Set3") +
  scale_x_discrete(labels = function(x) gsub("\\.", "\n", x))  # Format x-axis labels

# Save the violin plot
ggsave(filename = "figs/LASSO_EN_RIDGE/Violin_PEARSON_Mean.png", plot = violin_plot, width = 12, height = 8)
#################################################### COMPARISONNNNNNNNNNNNNNNNN 
# Load required libraries
library(ggplot2)
library(dplyr)
library(readr)
# Define file paths
file_paths <- list(
  Lasso_derived_features = "../../resums/resum_lasso_LM.csv",
  Drugbank_features = "../../resums/resum_drugbank_drug_targets.csv",
  DTC_features = "../../resums/resum_dtc_drug_targets.csv",
  Pagerank_features = "../../resums/resum_uniform_pagerank_central_genes.csv",
  Lasso_Pagerank_features = "../../resums/resum_lasso_PPR.csv",
  Drugbank_Pagerank_features = "../../resums/resum_drugbank_PPR.csv",
  Random_features = "../../resums/resum_random50.csv",
  All_genes = "../../resums/resum_ridge_LM.csv"
)

# Load and process data for GDSC
data_list_gdsc <- lapply(names(file_paths), function(name) {
  file <- file_paths[[name]]
  if (file.exists(file)) {
    data <- read_csv(file, show_col_types = FALSE)
    if ("PEARSON_Mean_GDSC" %in% colnames(data)) {
      data <- data %>%
        select(Drug, PEARSON_Mean = PEARSON_Mean_GDSC) %>%
        filter(is.finite(PEARSON_Mean)) %>%
        mutate(Method = name)
    } else {
      data <- tibble()  # Return an empty tibble if PEARSON_Mean_GDSC is missing
    }
  } else {
    message(paste("File does not exist:", file))
    data <- tibble()  # Return an empty tibble if the file does not exist
  }
  return(data)
})

# Load and process data for gCSI
data_list_gcsi <- lapply(names(file_paths), function(name) {
  file <- file_paths[[name]]
  if (file.exists(file)) {
    data <- read_csv(file, show_col_types = FALSE)
    if ("PEARSON_Mean_gCSI" %in% colnames(data)) {
      data <- data %>%
        select(Drug, PEARSON_Mean = PEARSON_Mean_gCSI) %>%
        filter(is.finite(PEARSON_Mean)) %>%
        mutate(Method = name)
    } else {
      data <- tibble()  # Return an empty tibble if PEARSON_Mean_gCSI is missing
    }
  } else {
    message(paste("File does not exist:", file))
    data <- tibble()  # Return an empty tibble if the file does not exist
  }
  return(data)
})

# Combine all data into single data frames
combined_data_gdsc <- bind_rows(data_list_gdsc)
combined_data_gcsi <- bind_rows(data_list_gcsi)

# Ensure combined_data_gdsc and combined_data_gcsi have the expected structure
if (nrow(combined_data_gdsc) == 0) {
  stop("No valid GDSC data was loaded. Please check the file paths and data structure.")
}
if (nrow(combined_data_gcsi) == 0) {
  stop("No valid gCSI data was loaded. Please check the file paths and data structure.")
}

# Generate violin and boxplot for GDSC
plot_gdsc <- ggplot(combined_data_gdsc, aes(x = Method, y = PEARSON_Mean, fill = Method)) +
  geom_violin(trim = TRUE, alpha = 0.6) +
  geom_boxplot(width = 0.1, fill = "white", outlier.shape = NA) +
  theme_minimal() +
  labs(title = "Distribution of PEARSON_Mean Across Methods (GDSC)",
       x = "Method",
       y = "PEARSON Mean") +
  theme(legend.position = "none") +
  scale_fill_brewer(palette = "Set3")

# Save the GDSC plot
ggsave(filename = "figs/Overall/Violin_Boxplot_PEARSON_Mean_Distributions_GDSC.png", plot = plot_gdsc, width = 10, height = 6)

# Generate violin and boxplot for gCSI
plot_gcsi <- ggplot(combined_data_gcsi, aes(x = Method, y = PEARSON_Mean, fill = Method)) +
  geom_violin(trim = TRUE, alpha = 0.6) +
  geom_boxplot(width = 0.1, fill = "white", outlier.shape = NA) +
  theme_minimal() +
  labs(title = "Distribution of PEARSON_Mean Across Methods (gCSI)",
       x = "Method",
       y = "PEARSON Mean") +
  theme(legend.position = "none") +
  scale_fill_brewer(palette = "Set3")

# Save the gCSI plot
ggsave(filename = "figs/Overall/Violin_Boxplot_PEARSON_Mean_Distributions_gCSI.png", plot = plot_gcsi, width = 10, height = 6)


# Perform Wilcoxon tests for GDSC and gCSI
wilcoxon_results_GDSC <- run_wilcoxon_tests_no_size(filtered_data_GDSC)
wilcoxon_results_gCSI <- run_wilcoxon_tests_no_size(filtered_data_gCSI)

# Save the results to CSV files
write.csv(wilcoxon_results_GDSC, "results/wilcoxon_results_GDSC.csv", row.names = FALSE)
write.csv(wilcoxon_results_gCSI, "results/wilcoxon_results_gCSI.csv", row.names = FALSE)

# Print the results
print(wilcoxon_results_GDSC)
print(wilcoxon_results_gCSI)

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


##############################################################################################
##############################################################################################
        ############################## WILCOXON #######################################
##############################################################################################
##############################################################################################
##############################################################################################

# Load required libraries
library(dplyr)
library(tidyr)

# Define the folder path
folder_path <- "../../resums/"

# Get the list of CSV files in the folder
# Prompt the user to choose the pattern for filtering files

choice <- 2

# Set the pattern based on the user's choice
if (choice == 1) {
  pattern <- "\\_(drug_targets|ridge_LM)\\.csv$"
} else if (choice == 2) {
  pattern <- "\\_(PPR|ridge_LM)\\.csv$"
} else if (choice == 3){
  pattern <- "\\_(central_genes|ridge_LM)\\.csv$"
} else if (choice == 4) {
  pattern <- "\\_(PPR|ridge_LM)\\.csv$"
} else {
  stop("Invalid choice. Please run the script again and choose 1 or 2.")
}

# List files based on the chosen pattern
csv_files <- list.files(folder_path, pattern = pattern, full.names = TRUE)

csv_files <- list(
  Lasso_derived_features = "../../resums/resum_lasso_derived.csv",
  Drugbank_features = "../../resums/resum_drugbank_drug_targets.csv",
  Pagerank_features = "../../resums/resum_uniform_pagerank_central_genes.csv",
  Lasso_Pagerank_features = "../../resums/resum_lasso_PPR.csv",
  Drugbank_Pagerank_features = "../../resums/resum_drugbank_PPR.csv",
  Random_features = "../../resums/resum_random50.csv",
  All_features = "../../resums/resum_ridge_LM.csv",
)

# Initialize empty data frames to store combined data for GDSC and gCSI
combined_data_GDSC <- data.frame(
  Drug = character(),
  PEARSON_Mean = numeric(),
  Method = character(),
  stringsAsFactors = FALSE
)

combined_data_gCSI <- data.frame(
  Drug = character(),
  PEARSON_Mean = numeric(),
  Method = character(),
  stringsAsFactors = FALSE
)

# Loop through each file and extract PEARSON_Mean data for GDSC and gCSI
for (file in csv_files) {
  # Extract the method name from the file name
  method_name <- gsub("\\.csv$", "", basename(file))
  
  # Read the CSV file
  data <- read.csv(file)
  
  # Process GDSC data
  if ("PEARSON_Mean_GDSC" %in% colnames(data)) {
    gdsc_data <- data %>%
      select(Drug, PEARSON_Mean = PEARSON_Mean_GDSC) %>%
      mutate(Method = method_name)
    
    # Append to the combined GDSC data
    combined_data_GDSC <- bind_rows(combined_data_GDSC, gdsc_data)
  }
  
  # Process gCSI data
  if ("PEARSON_Mean_gCSI" %in% colnames(data)) {
    gcsi_data <- data %>%
      select(Drug, PEARSON_Mean = PEARSON_Mean_gCSI) %>%
      mutate(Method = method_name)
    
    # Append to the combined gCSI data
    combined_data_gCSI <- bind_rows(combined_data_gCSI, gcsi_data)
  }
}

# Ensure combined_data_GDSC and combined_data_gCSI have valid PEARSON_Mean values
filtered_data_GDSC <- combined_data_GDSC %>%
  filter(is.finite(PEARSON_Mean))

filtered_data_gCSI <- combined_data_gCSI %>%
  filter(is.finite(PEARSON_Mean))
# Perform Wilcoxon tests for GDSC and gCSI
wilcoxon_results_GDSC <- run_wilcoxon_tests_no_size(filtered_data_GDSC)
wilcoxon_results_gCSI <- run_wilcoxon_tests_no_size(filtered_data_gCSI)

# Save the results to CSV files
write.csv(wilcoxon_results_GDSC, "results/wilcoxon_results_GDSC.csv", row.names = FALSE)
write.csv(wilcoxon_results_gCSI, "results/wilcoxon_results_gCSI.csv", row.names = FALSE)

print(wilcoxon_results_GDSC)
print(wilcoxon_results_gCSI)

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
