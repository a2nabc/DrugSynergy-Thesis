library(glmnet)
library(PharmacoGx)
library(dplyr)
source("src/model_functions.R")
#install.packages("config")
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

load_data <- function(expression_path, response_path) {
  expression_data <- read.csv(expression_path)
  response_data <- read.csv(response_path)
  return(list(expression = expression_data, response = response_data))
}

source_data <- load_data(config$source.expression, config$source.response)

target_data <- load_data(config$target.expression, config$target.response)

# filter source and target just with genes from KEGG and LINCS
keep_genes <- get_gene_union(config$gene.sets, config$gene.sets.path)
keep_genes <- keep_genes[keep_genes %in% colnames(source_data$expression)]
source_data$expression <- source_data$expression %>%
  select(c("Cell_line", keep_genes))
target_data$expression <- target_data$expression %>%
  select(c("Cell_line", keep_genes))

# Extract cell line identifiers
source_cells <- source_data$expression$Cell_line
target_cells <- target_data$expression$Cell_line

if(config$experiment.type.positive){ #consistency experiment --> we want to test on common cell lines
  common_cells <- intersect(source_cells, target_cells)
  target_data$expression <- target_data$expression %>%
    filter(Cell_line %in% common_cells)
  target_data$response <- target_data$response %>%
    filter(Cell_line %in% common_cells)
} else{ # generalization -> we want to test only in unseen cell_lines
  target_data$expression <- target_data$expression %>%
    filter(!Cell_line %in% source_cells)
  target_data$response <- target_data$response %>%
    filter(!Cell_line %in% source_cells)
}

# Print dataset sizes after filtering
#cat("Training Data:", nrow(source_data$expression), "samples\n")
#cat("Testing Data:", nrow(target_data$expression), "samples\n")

common_drugs <- readLines(config$drugs) 


# Initialize list to store results

model_results <- list()

for (model.type in config$models) {
  for (drug in common_drugs) {
    cat("\nTraining", model.type, "model for", drug, "...\n")
    
    # Prepare training data
    train_data <- prepare.data.for.ml(source_data, drug)
    cat("Training data dimensions for", drug, ",", model.type, "model:", dim(train_data), "\n")
    X_train <- train_data %>% select(-AAC)  # Gene expression
    y_train <- train_data$AAC  # Drug response
    
    # Train model
    trained_model <- fit.linear.model(X_train, y_train, model.type)
    
    # Extract non-zero coefficient genes
    features <- coef(trained_model$model)
    features <- as.data.frame(as.matrix(features))
    features <- rownames(features)[features[, 1] != 0]
    
    # ho guardem a results/MODEL/Drug
    ifelse(!dir.exists(config$results.features),dir.create(config$results.features, recursive=TRUE), FALSE)
    writeLines(features, paste0(config$results.features, drug, ".txt"))
    
    # Prepare testing data
    test_data <- prepare.data.for.ml(target_data, drug)
    cat("Testing data dimensions for", drug, ":", dim(test_data), "\n")
    X_test <- as.matrix(test_data %>% select(-AAC))  # Gene expression
    y_test <- test_data$AAC  # Drug response
    
    # Predict on test set
    y_pred <- predict(trained_model$model, newx = X_test)
    
    # Evaluate model performance
    eval_metrics <- evaluate_regression(y_test, y_pred)
    write.csv(eval_metrics, paste0(config$results.features, drug, "_metrics.csv"))
    
    # Compute gene-expression correlation with response
    gene_corrs <- sapply(features, function(gene) {
      if (gene %in% colnames(test_data)) {
        cor(test_data[[gene]], y_test, use = "complete.obs")
      } else {
        NA
      }
    })
    write.csv(gene_corrs, paste0(config$results.features, drug, "_genes_corr.csv"))
    
    # Store results
    model_results[[drug]][[model.type]] <- list(
      model = trained_model$model,
      eval = eval_metrics,
      features_selected = features,
      gene_correlations = gene_corrs
    )
   
    # Perform 10x cross-validation (80/20 splits)
    cross_val_results <- perform.cross.validation(10, test_data)
    
    # Save cross-validation results
    model_results[[drug]][[model.type]]$cross_validation <- cross_val_results
  }
}
