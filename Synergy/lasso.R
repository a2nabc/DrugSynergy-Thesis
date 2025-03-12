library(dplyr)
library(tidyr)
library(glmnet)
library(Metrics)  
library(caret)  
library(pROC)  

# METRIC PERFORMANCE EVALUATION for Lasso, Elastic Net, Ridge, Kernel Ridge Regression
evaluate_regression <- function(y_true, y_pred) { 
  mse <- mean((y_true - y_pred)^2)
  rmse <- sqrt(mse)
  mae <- mean(abs(y_true - y_pred))
  r2 <- 1 - (sum((y_true - y_pred)^2) / sum((y_true - mean(y_true))^2))
  
  return(list(MSE = mse, RMSE = rmse, MAE = mae, R2 = r2))
}

# METRIC PERFORMANCE EVALUATION for Lasso Classification, Support Vector Machine (SVM)
evaluate_classification <- function(y_true, y_pred) {
  cm <- confusionMatrix(as.factor(y_pred), as.factor(y_true))
  accuracy <- cm$overall["Accuracy"]
  precision <- cm$byClass["Pos Pred Value"]
  recall <- cm$byClass["Sensitivity"]
  f1 <- 2 * (precision * recall) / (precision + recall)
  
  roc_obj <- roc(y_true, as.numeric(y_pred))
  auc <- auc(roc_obj)
  
  return(list(Accuracy = accuracy, Precision = precision, Recall = recall, F1_Score = f1, AUC = auc))
}

run_lasso <- function(X, y) {
  X <- as.matrix(X)  # Ensure X is a matrix
  lasso_cv <- cv.glmnet(X, y, alpha = 1)  # Perform cross-validation
  best_lambda <- lasso_cv$lambda.min  # Best lambda found by CV
  print(paste("Optimal lambda:", best_lambda))
  
  # Fit Lasso with optimal lambda
  lasso_model <- glmnet(X, y, alpha = 1, lambda = best_lambda)
  
  y_pred <- predict(lasso_model, X, s = best_lambda)
  
  # Evaluate model performance
  metrics <- evaluate_regression(y, y_pred)
  print(metrics)
  

  
  # Extract coefficients
  coefficients <- coef(lasso_model)
  coefs_df <- data.frame(
    Feature = rownames(coefficients),
    Coefficient = as.numeric(coefficients)
  )
  
  # Filter non-zero coefficients
  non_zero_coefs <- coefs_df[coefs_df$Coefficient != 0, ]
  
  print(non_zero_coefs)
  return(list(metrics = metrics, non_zero_coefs = non_zero_coefs))
}

run_elastic_net <- function(X, y) {
  alpha_value = 0.5
  X <- as.matrix(X)  # Ensure X is a matrix
  
  # Perform cross-validation for Elastic Net
  elastic_cv <- cv.glmnet(X, y, alpha = alpha_value)  
  best_lambda <- elastic_cv$lambda.min  # Best lambda found by CV
  print(paste("Optimal lambda:", best_lambda))
  
  # Fit Elastic Net with optimal lambda
  elastic_model <- glmnet(X, y, alpha = alpha_value, lambda = best_lambda)
  
  y_pred <- predict(elastic_model, X, s = best_lambda)
  
  # Evaluate model performance
  metrics <- evaluate_regression(y, y_pred)
  print(metrics)
  
  # Extract coefficients and filter non-zero ones
  coefficients <- coef(elastic_model)
  
  # Convert to data frame
  coefs_df <- data.frame(
    Feature = rownames(coefficients),
    Coefficient = as.numeric(coefficients)
  )
  
  non_zero_coefs <- coefs_df[coefs_df$Coefficient != 0, ]
  
  print(non_zero_coefs)
  return(list(metrics = metrics, non_zero_coefs = non_zero_coefs))
}

run_ridge <- function(X, y) {
  X <- as.matrix(X)  # Ensure X is a matrix
  
  # Perform cross-validation for Ridge Regression (alpha = 0)
  ridge_cv <- cv.glmnet(X, y, alpha = 0)  
  best_lambda <- ridge_cv$lambda.min  # Best lambda found by CV
  print(paste("Optimal lambda:", best_lambda))
  
  # Fit Ridge Regression with optimal lambda
  ridge_model <- glmnet(X, y, alpha = 0, lambda = best_lambda)
  
  y_pred <- predict(ridge_model, X, s = best_lambda)
  
  # Evaluate model performance
  metrics <- evaluate_regression(y, y_pred)
  print(metrics)
  
  # Extract coefficients
  coefficients <- coef(ridge_model)
  
  # Convert to data frame
  coefs_df <- data.frame(
    Feature = rownames(coefficients),
    Coefficient = as.numeric(coefficients)
  )
  
  non_zero_coefs <- coefs_df[coefs_df$Coefficient != 0, ]
  
  print(non_zero_coefs)
  return(list(metrics = metrics, non_zero_coefs = non_zero_coefs))
}



run_lasso_classification <- function(X, y) {
  # Ensure X is a matrix
  X <- as.matrix(X)
  
  # Perform cross-validation for logistic regression with Lasso regularization
  lasso_cv <- cv.glmnet(X, y, family = "binomial", alpha = 1)  # Cross-validation for logistic regression
  best_lambda <- lasso_cv$lambda.min  # Best lambda found by CV
  
  print(paste("Optimal lambda:", best_lambda))
  
  # Fit Lasso with optimal lambda
  lasso_model <- glmnet(X, y, family = "binomial", alpha = 1, lambda = best_lambda)
  
  y_pred_prob <- predict(lasso_model, X, s = best_lambda, type = "response")
  y_pred <- ifelse(y_pred_prob > 0.5, 1, 0)  # Convert probabilities to class labels
  
  metrics <- evaluate_classification(y, y_pred)
  
  
  
  # Extract coefficients
  coefficients <- coef(lasso_model)
  coefs_df <- data.frame(
    Feature = rownames(coefficients),
    Coefficient = as.numeric(coefficients)
  )
  
  # Filter non-zero coefficients
  non_zero_coefs <- coefs_df[coefs_df$Coefficient != 0, ]
  
  print(non_zero_coefs)
  
  return(list(metrics = metrics, non_zero_coefs = non_zero_coefs))
}

library(e1071)  # For SVM
run_svm <- function(X, y) {
  X <- as.matrix(X)  # Ensure X is a matrix
  y <- factor(y, levels=c(0,1))
  # Train an SVM model with RBF kernel
  svm_model <- svm(X, y, kernel = "radial", cost = 1, gamma = 1/ncol(X), class.weights = c("0" = 9, "1" = 1))
  
  # Check if the model was successfully trained
  if (is.null(svm_model)) {
    return(NULL)
  }
  
  # Extract support vectors
  support_vectors <- svm_model$SV
  
  # Predict on training data
  y_pred <- predict(svm_model, X)
  
  metrics <- evaluate_classification(y, y_pred)
  return(list(metrics = metrics, non_zero_coefs = NULL))
}

library(kernlab)  # For Kernel Ridge Regression

run_krr <- function(X, y) {
  X <- as.matrix(X)  # Ensure X is a matrix
  
  # Train KRR model with RBF kernel
  krr_model <- ksvm(X, y, kernel = "rbfdot", kpar = "automatic", C = 1)
  
  # Extract learned parameters
  alpha_values <- krr_model@alpha
  
  # Predict on training data
  y_pred <- predict(krr_model, X)
  
  metrics <- evaluate_regression(y, y_pred)
  
  return(list(metrics = metrics, non_zero_coefs = NULL))
}

# drugs to analyze individually
drug_names <- c("Bortezomib", "Axitinib", "Cisplatin", "Erlotinib", "Lapatinib", "Ruxolitinib", "Sirolimus", "Vinorelbine", "Vorinostat")

gene_subset_names <- c("ProteinCoding", "COSMIC", "KEGG", "LINCS", "MDSig_Hallmarks")


##################### Read single drug files #############################

# One file for each drug and each gene subset

drug_data_lists <- list()

for (subset in gene_subset_names) {
  drug_data_lists[[subset]] <- list()
  for (drug in drug_names) {
    file_suffix <- ifelse(subset == "ProteinCoding", "", paste0("_", subset))
    file_path <- paste0("data/single_drug/", drug, file_suffix, ".csv")
    drug_data_lists[[subset]][[drug]] <- read.csv(file_path)
  }
}


######################### Read synergy files #############################

# in synergy files I have the following combinations

#synergy1 -> Bortezomib, Erlotinib
#synergy2 -> Cisplatin, Lapatinib
#synergy3 -> Sirolimus, Axitinib
#synergy4 -> Vorinostat, Ruxolitinib
#synergy5 -> Vorinostat, Vinorelbine

# List of synergy file names
synergy_files <- paste0("data/synergy/synergy", 1:5, ".csv")
synergy_names <- paste0("synergy", 1:5)

# Read synergy data dynamically
synergy_data_list <- list()
for (i in seq_along(synergy_names)) {
  synergy_data_list[[synergy_names[i]]] <- read.csv(synergy_files[i])
}


############ Find common genes in synergy / single drug files ################

# Find common genes 
common_cols_list <- list()
for (subset in gene_subset_names) {
  common_cols_list[[subset]] <- intersect(names(synergy_data_list$synergy1), names(drug_data_lists[[subset]]$Bortezomib))
}

# Update synergy data in filtered_synergy_list
filtered_synergy_list <- list()
for (subset in gene_subset_names) {
  filtered_synergy_list[[subset]] <- list()
  for (name in synergy_names) {
    filtered_synergy_list[[subset]][[name]] <- synergy_data_list[[name]] %>%
      select(sampleid, treatment1id, treatment2id, median_SCORE, all_of(common_cols_list[[subset]]))
  }
}

# Update single drug data in filtered_drug_data_list
filtered_drug_data_list <- list()
for (subset in gene_subset_names) {
  filtered_drug_data_list[[subset]] <- list()
  for (drug in drug_names) {
    filtered_drug_data_list[[subset]][[drug]] <- drug_data_lists[[subset]][[drug]] %>%
      select(sampleid, aac, label, all_of(common_cols_list[[subset]]))
  }
}


############ Find common cell lines in synergy / single drug files ################

common_cell_lines <- list()
filtered_common_cell_lines_synergy_data_list <- list()

for (subset in gene_subset_names) {
  common_cell_lines[[subset]] <- list()
  filtered_common_cell_lines_synergy_data_list[[subset]] <- list()
  for (i in seq_along(synergy_names)) {
    cell_lines_synergy <- synergy_data_list[[synergy_names[i]]]$sampleid
    drug1 <- synergy_data_list[[synergy_names[i]]]$treatment1id[1]
    drug2 <- synergy_data_list[[synergy_names[i]]]$treatment2id[1]
    cell_lines_drug1 <- drug_data_lists[[subset]][[drug1]]$sampleid
    cell_lines_drug2 <- drug_data_lists[[subset]][[drug2]]$sampleid
    
    common_cell_lines[[subset]][[synergy_names[i]]] <- intersect(cell_lines_synergy, cell_lines_drug1)
    common_cell_lines[[subset]][[synergy_names[i]]] <- intersect(common_cell_lines[[subset]][[synergy_names[i]]], cell_lines_drug2)
    
    filtered_common_cell_lines_synergy_data_list[[subset]][[synergy_names[i]]] <- filtered_synergy_list[[subset]][[synergy_names[i]]] %>%
      filter(sampleid %in% common_cell_lines[[subset]][[synergy_names[i]]])
  }
}

############################################## LASSO MODEL ######################################

######################### Synergy (all data available in NCI ALMANAC) ###########################

results_synergy_list <- list()
results_synergy_list_classification <- list()

for (i in seq_along(synergy_names)) {
  synergy_data <- synergy_data_list[[synergy_names[i]]]  
  y <- synergy_data$median_SCORE

  y <- (y - min(y)) / (max(y) - min(y)) 
  X <- synergy_data[, 5:ncol(synergy_data)]
  X <- log2(X+1)
  print(synergy_names[[i]])
  result <- tryCatch({
    run_lasso(X, y)
  }, error = function(e) {
    message(paste("Skipping drug:", drug, "due to error:", e$message))
    return(NULL)
  })
  results_synergy_list[[synergy_names[i]]] <- result

  z <- synergy_data$median_label
  X <- synergy_data[, 5:ncol(synergy_data)]
  X <- log2(X+1)
  print(synergy_names[[i]])
  result <- tryCatch({
    run_lasso_classification(X, z)
  }, error = function(e) {
    message(paste("Skipping drug:", drug, "due to error:", e$message))
    return(NULL)
  })
  results_synergy_list_classification[[synergy_names[i]]] <- result
}

for (entry_name in names(results_synergy_list)) {
  feature_names <- results_synergy_list[[entry_name]]$Feature  # Extract feature names
  feature_names <- feature_names[feature_names != "(Intercept)"]  # Exclude "(Intercept)"
  feature_list <- paste(feature_names, collapse = " ")  # Join feature names into a string
  cat(paste(entry_name, ":", feature_list, "\n\n"))  # Print the result
}


##################### Filtered synergy (common cell lines and common genes) ########################

results_synergy_list_filtered <- list()
for (subset in gene_subset_names) {
  results_synergy_list_filtered[[subset]] <- list()
  for (i in seq_along(synergy_names)) {
    filtered_synergy_data <- filtered_synergy_list[[subset]][[synergy_names[i]]]
    y <- filtered_synergy_data$median_SCORE
    y <- (y - min(y)) / (max(y) - min(y))
    X <- filtered_synergy_data[, 5:ncol(filtered_synergy_data)]
    X <- log2(X + 1)
    print(synergy_names[[i]])
    
    result <- tryCatch({
      run_lasso(X, y)
    }, error = function(e) {
      message(paste("Skipping drug:", drug, "due to error:", e$message))
      return(NULL)
    })
    
    results_synergy_list_filtered[[subset]][[synergy_names[i]]] <- result
  }
  
  for (entry_name in names(results_synergy_list_filtered[[subset]])) {
    feature_names <- results_synergy_list_filtered[[subset]][[entry_name]]$Feature  # Extract feature names
    feature_names <- feature_names[feature_names != "(Intercept)"]  # Exclude "(Intercept)"
    feature_list <- paste(feature_names, collapse = " ")  # Join feature names into a string
    cat(paste(entry_name, ":", feature_list, "\n\n"))  # Print the result
  }
}

################################# Filtered single drug ######################################

results_list_lasso <- list()
results_list_elastic_net <- list()
results_list_ridge <- list()
results_list_classification <- list()
results_list_svm <- list()
results_list_krr <- list()

for (subset in gene_subset_names) {
  subset_results_lasso <- list()
  subset_results_elastic_net <- list()
  subset_results_ridge <- list()
  subset_results_classification <- list()
  subset_results_svm <- list()
  subset_results_krr <- list()
  
  data_regression <- filtered_drug_data_list[[subset]]
  
  for (drug in drug_names) {
    drug_data <- data_regression[[drug]]
    y <- drug_data$aac
    z <- drug_data$label
    X <- drug_data[, 4:ncol(drug_data)]
    X <- log2(X + 1)
    
    print(paste(subset, "gene subset\n\n"))
    print(drug)
    print("LASSO: ")
    result <- tryCatch({
      run_lasso(X, y)
    }, error = function(e) {
      message(paste("Skipping drug:", drug, "due to error:", e$message))
      return(NULL)
    })
    subset_results_lasso[[drug]] <- result
    
    print("ELASTIC NET: ")
    result <- tryCatch({
      run_elastic_net(X, y)
    }, error = function(e) {
      message(paste("Skipping drug:", drug, "due to error:", e$message))
      return(NULL)
    })
    subset_results_elastic_net[[drug]] <- result
    
    print("RIDGE REGRESSION: ")
    result <- tryCatch({
      run_ridge(X, y)
    }, error = function(e) {
      message(paste("Skipping drug:", drug, "due to error:", e$message))
      return(NULL)
    })
    subset_results_ridge[[drug]] <- result
    
    print("CLASSIFICATION (LASSO): ")
    result_classification <- tryCatch({
      run_lasso_classification(X, z)
    }, error = function(e) {
      message(paste("Skipping classification for drug:", drug, "due to error:", e$message))
      return(NULL)
    })
    subset_results_classification[[drug]] <- result_classification
    
    print("SVM: ")
    print(table(z))
    result_svm <- tryCatch({
      run_svm(X, z)
    }, error = function(e) {
      message(paste("Skipping drug:", drug, "due to error:", e$message))
      return(NULL)
    })
    subset_results_svm[[drug]] <- result_svm
    
    print("KRR: ")
    result_krr <- tryCatch({
      run_krr(X, y)
    }, error = function(e) {
      message(paste("Skipping drug:", drug, "due to error:", e$message))
      return(NULL)
    })
    subset_results_krr[[drug]] <- result_krr
    
  }
  
  results_list_lasso[[subset]] <- subset_results_lasso
  results_list_elastic_net[[subset]] <- subset_results_elastic_net
  results_list_ridge[[subset]] <- subset_results_ridge
  results_list_classification[[subset]] <- subset_results_classification
  results_list_svm[[subset]] <- subset_results_svm
  results_list_krr[[subset]] <- subset_results_krr
}

model_list <- c("lasso", "elastic_net", "ridge", "classification", "svm", "krr")
  
# Print and save files with features
save_results <- function(gene_subset) {
  for (drug_name in drug_names) {
    for (model in model_list){
      results <- switch(model,
                        "lasso" = results_list_lasso,
                        "elastic_net" = results_list_elastic_net,
                        "ridge" = results_list_ridge,
                        "classification" = results_list_classification,
                        "svm" = results_list_svm,
                        "krr" = results_list_krr)
      feature_names <- results[[gene_subset]][[drug_name]]$non_zero_coefs$Feature  # Extract feature names <=> gene names
      feature_names <- feature_names[feature_names != "(Intercept)"]  # Exclude "(Intercept)"
      
      if (length(feature_names) == 0) {
        feature_names <- "No features selected"
      }
      
      # Define directory and file paths
      dir_path <- file.path("Results2", "Single_Drug", gene_subset, toupper(model))  # Convert model to uppercase for consistency
      file_path <- file.path(dir_path, paste0(drug_name, ".txt"))
      
      # Create directory if it doesn't exist
      if (!dir.exists(dir_path)) {
        dir.create(dir_path, recursive = TRUE)
      }
      
      # Save feature names to a text file
      writeLines(feature_names, file_path)
      
      # Print status message
      cat(paste("Saved:", file_path, "\n"))
      
      # Save metrics
      metrics_path <- file.path(dir_path, paste0(drug_name, "_metrics.txt"))
      writeLines(capture.output(print(results[[gene_subset]][[drug_name]]$metrics)), metrics_path)
      cat(paste("Saved metrics:", metrics_path, "\n"))
      
      # Print results on console
    #  feature_list <- paste(feature_names, collapse = " ")  # Join feature names into a string
    #  cat(paste(entry_name, model, ":", feature_list, "\n\n")) 
    }
  }
}

# Extract and print feature names for each gene subset
for (subset in gene_subset_names) {
  print(subset)
  save_results(subset)
}
