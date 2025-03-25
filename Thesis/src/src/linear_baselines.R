library(glmnet)
library(PharmacoGx)
library(dplyr)
source("src/model_functions.R")
source("src/data_functions.R")

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

#config <- config::get(file = "src/configs/Broad_to_gCSI.yml")

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
    # Train model
    train_data <- prepare.data.for.ml(source_data, drug)
    print(paste("Training data is", nrow(train_data), "samples and", ncol(train_data), "features"))
    model_info <- train.model(model.type, train_data, drug)
    
    # Save features
    features <- model_info$features
    features <- features[features != "(Intercept)"]  # Exclude "(Intercept)"
    writeLines(features, paste0(config$results.dir, results.subdir, model.type, "/", config$results.features.dir, drug, ".txt"))
    
    # Evaluate model
    test_data <- prepare.data.for.ml(target_data, drug)
    print(paste("Testing data is", nrow(test_data), "samples and", ncol(test_data), "features"))

    eval_info <- evaluate.model(model_info$model, test_data, drug)
    write.csv(eval_info, paste0(config$results.dir, results.subdir,  model.type, "/", config$results.eval.dir, drug, ".csv"))
    
    # Compute gene-expression correlation
    gene_corrs <- compute.gene.correlations(features, test_data, eval_info$y_test, drug)
    write.csv(gene_corrs, paste0(config$results.dir, results.subdir, model.type, "/", config$results.correlations.dir, drug, ".csv"))
    
    # Perform cross-validation
    cross_val_avg <- perform.cv(test_data)
    write.csv(cross_val_avg$mean, paste0(config$results.dir, results.subdir, model.type, "/", config$results.cv.dir, drug, "_mean.csv"))
    write.csv(cross_val_avg$std, paste0(config$results.dir, results.subdir, model.type, "/", config$results.cv.dir, drug, "_std.csv"))
    
    # Store results
    model_results[[drug]][[model.type]] <- list(
      model = model_info$model,
      eval = eval_info$eval_metrics,
      features_selected = features,
      gene_correlations = gene_corrs,
      cross_validation_avg = cross_val_avg
    )
    
    # Store dimensions of training / testing data
    new_row_dimensions <- data.frame(
      Drug = drug,
      Model_Type = model.type,
      Training_Samples = nrow(train_data),
      Training_Features = ncol(train_data),
      Testing_Samples = nrow(test_data),
      Testing_Features = ncol(test_data)
    )
    data_dimensions_metadata <- rbind(data_dimensions_metadata, new_row_dimensions)
    
    biomarkers_str <- paste(features, collapse = ", ")
    
    new_row_summary_results <- data.frame(
      Drug= drug,
      Model = model.type,
      MSE = eval_info$eval_metrics$MSE,
      RMSE = eval_info$eval_metrics$RMSE,
      MAE = eval_info$eval_metrics$MAE,
      R2 = eval_info$eval_metrics$R2,
      PEARSON = eval_info$eval_metrics$PEARSON,
      Biomarkers = biomarkers_str,
      stringsAsFactors = FALSE
    )
    summary_results <- rbind(summary_results, new_row_summary_results)
    
  }
}
write.csv(data_dimensions_metadata, paste0(config$results.dir, results.subdir, config$results.metadata.dir, "data_dimensions.csv"))
write.csv(summary_results, paste0(config$results.dir, results.subdir, "summary_results.csv"))
# if(config$experiment.type.is.positive) {
#   saveRDS(model_results, "results_positive_2.rds")
# } else{
#   saveRDS(model_results, "results_negative_2.rds")
# }

# results_positive <- readRDS("results_positive.rds")
# results_negative <- readRDS("results_negative.rds")
