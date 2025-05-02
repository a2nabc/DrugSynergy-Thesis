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

# Load centrality measure performance data for gCSI
betweenness_gCSI <- read.csv("results/cv_performance/central_genes/gCSI/performance_betweenness_central_genes.csv")
degree_gCSI <- read.csv("results/cv_performance/central_genes/gCSI/performance_degree_central_genes.csv")
eigenvector_gCSI <- read.csv("results/cv_performance/central_genes/gCSI/performance_eigenvector_central_genes.csv")
pagerank_gCSI <- read.csv("results/cv_performance/central_genes/gCSI/performance_uniform_pagerank_central_genes.csv")

# Load centrality measure performance data for GDSC
betweenness_GDSC <- read.csv("results/cv_performance/central_genes/GDSC2/performance_betweenness_central_genes.csv")
degree_GDSC <- read.csv("results/cv_performance/central_genes/GDSC2/performance_degree_central_genes.csv")
eigenvector_GDSC <- read.csv("results/cv_performance/central_genes/GDSC2/performance_eigenvector_central_genes.csv")
pagerank_GDSC <- read.csv("results/cv_performance/central_genes/GDSC2/performance_uniform_pagerank_central_genes.csv")

# Add a column to identify the centrality measure and dataset
betweenness_gCSI$Method <- "Betweenness"
degree_gCSI$Method <- "Degree"
eigenvector_gCSI$Method <- "Eigenvector"
pagerank_gCSI$Method <- "Pagerank"
betweenness_GDSC$Method <- "Betweenness"
degree_GDSC$Method <- "Degree"
eigenvector_GDSC$Method <- "Eigenvector"
pagerank_GDSC$Method <- "Pagerank"

betweenness_gCSI$Screen <- "gCSI"
degree_gCSI$Screen <- "gCSI"
eigenvector_gCSI$Screen <- "gCSI"
pagerank_gCSI$Screen <- "gCSI"
betweenness_GDSC$Screen <- "GDSC2"
degree_GDSC$Screen <- "GDSC2"
eigenvector_GDSC$Screen <- "GDSC2"
pagerank_GDSC$Screen <- "GDSC2"

# Combine all data into a single data frame
combined_centrality <- bind_rows(
  betweenness_gCSI, degree_gCSI, eigenvector_gCSI, pagerank_gCSI,
  betweenness_GDSC, degree_GDSC, eigenvector_GDSC, pagerank_GDSC
)

# Save the combined data for further analysis
write.csv(combined_centrality, "results/cv_performance/central_genes/combined_centrality_performance.csv", row.names = FALSE)

# Generate violin plots for PEARSON_Mean by centrality measures and screens
# Filter data for Betweenness and Degree methods
methods_to_plot <- c("Betweenness", "Degree", "Eigenvector", "Pagerank")

# Loop through each method and generate plots
for (method in methods_to_plot) {
  method_data <- combined_centrality %>% 
    filter(Method == method) %>% 
    filter(is.finite(PEARSON_Mean), is.finite(RMSE_Mean))  # Remove non-finite values
  
  # PEARSON plot for gCSI and GDSC2
  pearson_plot <- ggplot(method_data, aes(x = Screen, y = PEARSON_Mean, fill = Screen)) +
    geom_violin(trim = TRUE, alpha = 0.6) +
    geom_boxplot(width = 0.1, fill = "white", outlier.shape = NA) +
    theme_minimal() +
    labs(title = paste("PEARSON_Mean for", method, "Method"),
         x = "Screen",
         y = "PEARSON Mean") +
    theme(legend.position = "none") +
    scale_fill_brewer(palette = "Set3")
  ggsave(filename = paste0("figs/CENTRALITY_MEASURES/Violin_PEARSON_Mean_", method, ".png"), plot = pearson_plot, width = 8, height = 6)
  
  # RMSE plot for gCSI and GDSC2
  rmse_plot <- ggplot(method_data, aes(x = Screen, y = RMSE_Mean, fill = Screen)) +
    geom_violin(trim = TRUE, alpha = 0.6) +
    geom_boxplot(width = 0.1, fill = "white", outlier.shape = NA) +
    theme_minimal() +
    labs(title = paste("RMSE_Mean for", method, "Method"),
         x = "Screen",
         y = "RMSE Mean") +
    theme(legend.position = "none") +
    scale_fill_brewer(palette = "Set3")
  ggsave(filename = paste0("figs/CENTRALITY_MEASURES/Violin_RMSE_Mean_", method, ".png"), plot = rmse_plot, width = 8, height = 6)
}


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




#####################################################################

# Initialize
results <- init.results(results, c("lasso", "random", "drugbank"), config$target.screens, feature_sizes)

# Lasso
for (screen in config$target.screens) {
  results <- run.lasso.cv(results, screen, feature_sizes, common_drugs, source_data, target_data, config, 10)
}
write.lasso.results(results, config$target.screens, feature_sizes, config)

# Random + Drugbank
results <- run.other.cv(results, feature_sizes, common_drugs, source_data, target_data, config, 10)
write.other.results(results, feature_sizes, config)

############################# PLOTS FOR PAGERANK VS RANDOM PERFORMANCES #############################
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
combined <- load_all_performance_data(feature_sizes, config)
filtered_combined <- filter_valid_drugs(combined)

# Barplots
plot_pearson_barplot(filtered_combined, 10, "figs/Barplot_10_features.png", "RdYlGn")
plot_pearson_barplot(filtered_combined, 20,  "figs/Barplot_20_features.png", "RdBu")
plot_pearson_barplot(filtered_combined, 50,  "figs/Barplot_50_features.png", "PuOr")

# Boxplots
plot_pearson_boxplot(filtered_combined, 10, "figs/Boxplot_10_features.png")
plot_pearson_boxplot(filtered_combined, 20, "figs/Boxplot_20_features.png")
plot_pearson_boxplot(filtered_combined, 50, "figs/Boxplot_50_features.png")

# Wilcoxon test results
wilcoxon_results <- run_wilcoxon_tests(filtered_combined, 50)
write.csv(wilcoxon_results, file = "results/cv_performance/wilcoxon_fdr_results.csv", row.names = FALSE)
print(wilcoxon_results)
