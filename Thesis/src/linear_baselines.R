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

config <- config::get(file = "configs/Broad_to_gCSI.yml")

load_data <- function(expression_path, response_path) {
  expression_data <- read.csv(expression_path)
  response_data <- read.csv(response_path)
  return(list(expression = expression_data, response = response_data))
}

source_data <- load_data(config$source.expression, config$source.response)

target_data <- load_data(config$target.expression, config$target.response)

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

# Iterate over models and drugs
for (model.type in config$models) {
  for (drug in config$drugs) {
    
  }
}
