library(glmnet)
library(PharmacoGx)
library(dplyr)
library(caret)
source("src/model_functions.R")
source("src/data_functions.R")
source("src/performance_results_functions.R")
source("src/heatmaps_plotting.R")

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

# filter source and target just with genes from KEGG and LINCS
keep_genes <- get_gene_union(config$gene.sets, config$gene.sets.path)
keep_genes <- keep_genes[keep_genes %in% colnames(source_data$expression)]

source_data$expression <- source_data$expression %>% select(c("Cell_line", keep_genes))
target_data$expression <- target_data$expression %>% select(c("Cell_line", keep_genes))

# filter cell lines across the two screens depending on the type of experiment
target_data <- filter_cell_lines(config$experiment.type.is.positive, source_data, target_data)

# Print dataset sizes after filtering
#cat("Training Data:", nrow(source_data$expression), "samples\n")
#cat("Testing Data:", nrow(target_data$expression), "samples\n")


# Read drug universe and initialize list and folders to store results
common_drugs <- readLines(config$drugs) 
results.subdir <- create_output_folders(config$target.screen, config$models, config$experiment.type.is.positive, config$results.dir, config$results.features.dir, config$results.eval.dir, config$results.correlations.dir, config$results.cv.dir)
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

############################# RESULTS AND PLOTS ###############################
ccl_list <- list()
for (screen in config$screens) {
  ccl_list <- rbind(ccl_list, paste0(config$path.to.processed.screens, screen, "/", config$name.ccl.file))
}
venn_plot(ccl_list, paste0(config$figs.dir, "venn_plot.png"))

for (screen in config$target.screens){
  generate_tables_summary_performance(paste0(config$results.dir, screen), "summary_results.csv") 
  generate_violin_plots(paste0(config$results.dir, screen), "summary_results.csv", paste0(config$figs.dir, screen))
  generate_model_heatmaps(paste0(config$results.dir, screen), config$results.correlations.dir, common_drugs, paste0(config$figs.dir, screen))
  jaccard_index(paste0(config$results.dir, screen), config$results.features.dir, common_drugs)
}

####################### K-FOLD CROSS VALIDATION FOR PAGERANK GENES ########################

for (screen in config$target.screens) {
  results_pagerank_10 <- list()
  results_random_10 <- list()
  results_pagerank_20 <- list()
  results_random_20 <- list()
  results_pagerank_50 <- list()
  results_random_50 <- list()
  
  for (drug in common_drugs) {
    path_10 <- paste0(config$results.dir, screen, "/positive/lasso/", config$results.pagerank.features.dir, drug, "_10.txt")
    results_drug_10 <- compute.cv.for.pagerank.input(path_10, source_data, target_data, drug, 10)
    results_pagerank_10 <- rbind(results_pagerank_10, results_drug_10)
    
    results_random_drug_10 <- compute.cv.for.random.input(source_data, target_data, drug, 10, 10)
    results_random_10 <- rbind(results_random_10, results_random_drug_10)
    
    
    path_20 <- paste0(config$results.dir, screen, "/positive/lasso/", config$results.pagerank.features.dir, drug, "_20.txt")
    results_drug_20 <- compute.cv.for.pagerank.input(path_20, source_data, target_data, drug, 10)
    results_pagerank_20 <- rbind(results_pagerank_20, results_drug_20)
    
    results_random_drug_20 <- compute.cv.for.random.input(source_data, target_data, drug, 10, 20)
    results_random_20 <- rbind(results_random_20, results_random_drug_20)
    
    
    path_50 <- paste0(config$results.dir, screen, "/positive/lasso/", config$results.pagerank.features.dir, drug, "_50.txt")  
    results_drug_50 <- compute.cv.for.pagerank.input(path_50, source_data, target_data, drug, 10)
    results_pagerank_50 <- rbind(results_pagerank_50, results_drug_50)
    
    results_random_drug_50 <- compute.cv.for.random.input(source_data, target_data, drug, 10, 50)
    results_random_50 <- rbind(results_random_50, results_random_drug_50)
  }
  
  write.csv(results_pagerank_10, paste0(config$results.dir, screen, "/positive/lasso", "/pagerank_performance_top_10_features.csv"))
  write.csv(results_pagerank_20, paste0(config$results.dir, screen , "/positive/lasso", "/pagerank_performance_top_20_features.csv"))
  write.csv(results_pagerank_50, paste0(config$results.dir, screen, "/positive/lasso", "/pagerank_performance_top_50_features.csv"))
  write.csv(results_random_10, paste0(config$results.dir, screen, "/positive/lasso", "/random_performance_10_features.csv"))
  write.csv(results_random_20, paste0(config$results.dir, screen , "/positive/lasso", "/random_performance_20_features.csv"))
  write.csv(results_random_50, paste0(config$results.dir, screen, "/positive/lasso", "/random_performance_50_features.csv"))
  
}

