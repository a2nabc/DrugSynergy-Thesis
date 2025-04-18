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
  adj_r2 <- 1 - ((1 - r2) * (length(y_true) - 1) / (length(y_true) - ncol(y_pred) - 1))
  mape <- mean(abs((y_true - y_pred) / y_true)) * 100
  
  return(list(MSE = mse, RMSE = rmse, MAE = mae, R2 = r2, Adjusted_R2 = adj_r2, MAPE = mape))
}

# METRIC PERFORMANCE EVALUATION for Lasso Classification, Support Vector Machine (SVM)
evaluate_classification <- function(y_true, y_pred) {
  cm <- confusionMatrix(as.factor(y_pred), as.factor(y_true))
  accuracy <- cm$overall["Accuracy"]
  precision <- cm$byClass["Pos Pred Value"]
  recall <- cm$byClass["Sensitivity"]
  specificity <- cm$byClass["Specificity"]
  balanced_accuracy <- (recall + specificity) / 2
  f1 <- 2 * (precision * recall) / (precision + recall)
  
  # Compute AUC only if y_pred contains probabilities, otherwise NA
  auc <- tryCatch({
    roc_obj <- roc(y_true, as.numeric(y_pred))
    auc(roc_obj)
  }, error = function(e) NA)
  
  cm_table <- as.table(cm$table)
  
  return(list(Accuracy = accuracy, Precision = precision, Recall = recall, 
              Specificity = specificity, Balanced_Accuracy = balanced_accuracy,
              F1_Score = f1, AUC = auc, Confusion_Matrix = cm_table))
}

# External function to run a model on K different splits and return averaged metrics
run_model_k_splits <- function(model_function, X, y, k = 5, classification = FALSE) {
  if (classification) {
    y <- as.factor(y)  # Ensure y is a factor for classification
  }
  
  folds <- createFolds(y, k = k, list = TRUE)  # Create K random splits
  metrics_list <- list()
  features_list <- list()
  
  for (fold_idx in seq_along(folds)) {
    print(paste("Processing fold", fold_idx, "of", k))
    
    test_idx <- folds[[fold_idx]]
    train_idx <- setdiff(seq_len(nrow(X)), test_idx)
    
    X_train <- X[train_idx, , drop = FALSE]
    y_train <- y[train_idx]
    X_test <- X[test_idx, , drop = FALSE]
    y_test <- y[test_idx]
    
    # Handle classification tasks separately
    if (classification) {
      print("It's classification")
      unique_classes <- length(unique(y_train))
      print(paste("Unique classes in training data:", unique_classes))
      
      if (length(unique(y_train)) < 2) {  
        message(paste("Skipping fold", fold_idx, "- Not enough classes in training data"))
        next  # Skip this fold if there's only one class
      }
    }
    
    # Try running the model
    result <- tryCatch({
      model_function(X_train, y_train)
    }, error = function(e) {
      message(paste("Skipping fold", fold_idx, "due to error:", e$message))
      return(NULL)
    })
    
    # Ensure result is valid before adding it to metrics_list
    if (!is.null(result) && "metrics" %in% names(result)) {
      print(paste("Fold", fold_idx, "- Metrics collected successfully"))
      metrics_list[[fold_idx]] <- result$metrics
    } else {
      print(paste("Fold", fold_idx, "- No valid metrics returned"))
    }
    
    if (!is.null(result) && "non_zero_coefs" %in% names(result)){
      print(paste("Fold", fold_idx, "- Feature selection collected successfully"))
      features_list[[fold_idx]] <- result$non_zero_coefs
    } else {
      print(paste("Fold", fold_idx, "- No valid features returned"))
    }
  }
  # If no valid results, return NULL
  if (length(metrics_list) == 0) {
    message("No valid metrics obtained from any fold. Returning NULL.")
    return(NULL)
  }
  
  # Compute mean & std of metrics across splits
  print("Aggregating metrics across all folds... aka computing mean and std")
  #print(metrics_list)
  final_metrics <- aggregate_metrics(metrics_list)
  print("Aggregating metrics complete!")
  #print(final_metrics)
  
  return(list(features_list = features_list, average_metrics = final_metrics, metrics_list = metrics_list))
}

extract_numeric <- function(x) {
  if (inherits(x, "auc")) {
    return(as.numeric(x))  # Convert AUC object to numeric
  } else if (is.character(x) && grepl("Area under the curve", x)) {
    return(as.numeric(sub(".*: ", "", x)))  # Extract numeric value after the colon
  } else if (is.list(x) && length(x) == 1) {
    return(as.numeric(x[[1]]))  # Extract numeric from named lists
  } else if (is.numeric(x)) {
    return(x)  # Already numeric
  } else {
    return(NA)  # If non-numeric, set to NA
  }
}

aggregate_metrics <- function(metrics_list) {
  if (length(metrics_list) == 0) {
    message("No metrics provided.")
    return(NULL)
  }
  
  # Convert metrics list into a proper dataframe
  metrics_df <- lapply(metrics_list, function(metrics) {
    metrics <- lapply(metrics, extract_numeric)  # Extract only numeric values
    metrics$Confusion_Matrix <- NULL  # Remove confusion matrices
    as.data.frame(metrics, stringsAsFactors = FALSE)
  }) %>% bind_rows()
  
  # Remove columns that are entirely NA
  metrics_df <- metrics_df[, colSums(!is.na(metrics_df)) > 0, drop = FALSE]
  
  if (ncol(metrics_df) == 0) {
    message("All metric columns are non-numeric or missing.")
    return(NULL)
  }
  
  # Compute mean and standard deviation
  mean_metrics <- colMeans(metrics_df, na.rm = TRUE)
  std_metrics <- apply(metrics_df, 2, sd, na.rm = TRUE)
  
  print("Aggregation complete!")
  return(list(mean = mean_metrics, std = std_metrics))
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
  #print(metrics)
  

  
  # Extract coefficients
  coefficients <- coef(lasso_model)
  coefs_df <- data.frame(
    Feature = rownames(coefficients),
    Coefficient = as.numeric(coefficients)
  )
  
  # Filter non-zero coefficients
  non_zero_coefs <- coefs_df[coefs_df$Coefficient != 0, ]
  
  #print(non_zero_coefs)
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
  #print(metrics)
  
  # Extract coefficients and filter non-zero ones
  coefficients <- coef(elastic_model)
  
  # Convert to data frame
  coefs_df <- data.frame(
    Feature = rownames(coefficients),
    Coefficient = as.numeric(coefficients)
  )
  
  non_zero_coefs <- coefs_df[coefs_df$Coefficient != 0, ]
  
  #print(non_zero_coefs)
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
  #print(metrics)
  
  # Extract coefficients
  coefficients <- coef(ridge_model)
  
  # Convert to data frame
  coefs_df <- data.frame(
    Feature = rownames(coefficients),
    Coefficient = as.numeric(coefficients)
  )
  
  non_zero_coefs <- coefs_df[coefs_df$Coefficient != 0, ]
  
  #print(non_zero_coefs)
  return(list(metrics = metrics, non_zero_coefs = non_zero_coefs))
}



run_lasso_classification <- function(X, y) {
  # Ensure X is a matrix
  X <- as.matrix(X)
  y <- as.factor(y)
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
  
  #print(non_zero_coefs)
  
  return(list(metrics = metrics, non_zero_coefs = non_zero_coefs))
}

library(e1071)  # For SVM
run_svm <- function(X, y) {
  X <- as.matrix(X)  # Ensure X is a matrix
  y <- factor(y, levels=c(0,1))
  
  class_counts <- table(y)
  class_weights <- 1 / class_counts  # Inverse weighting for imbalance
  
  # Train an SVM model with RBF kernel
  svm_model <- svm(X, y, kernel = "radial", cost = 1, gamma = 1/ncol(X), class.weights = class_weights)
  
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

# results_synergy_list <- list()
# results_synergy_list_classification <- list()
# 
# for (i in seq_along(synergy_names)) {
#   synergy_data <- synergy_data_list[[synergy_names[i]]]  
#   y <- synergy_data$median_SCORE
# 
#   y <- (y - min(y)) / (max(y) - min(y)) 
#   X <- synergy_data[, 5:ncol(synergy_data)]
#   X <- log2(X+1)
#   print(synergy_names[[i]])
#   result <- tryCatch({
#     run_lasso(X, y)
#   }, error = function(e) {
#     message(paste("Skipping drug:", drug, "due to error:", e$message))
#     return(NULL)
#   })
#   results_synergy_list[[synergy_names[i]]] <- result
# 
#   z <- synergy_data$median_label
#   X <- synergy_data[, 5:ncol(synergy_data)]
#   X <- log2(X+1)
#   print(synergy_names[[i]])
#   result <- tryCatch({
#     run_lasso_classification(X, z)
#   }, error = function(e) {
#     message(paste("Skipping drug:", drug, "due to error:", e$message))
#     return(NULL)
#   })
#   results_synergy_list_classification[[synergy_names[i]]] <- result
# }
# 
# for (entry_name in names(results_synergy_list)) {
#   feature_names <- results_synergy_list[[entry_name]]$Feature  # Extract feature names
#   feature_names <- feature_names[feature_names != "(Intercept)"]  # Exclude "(Intercept)"
#   feature_list <- paste(feature_names, collapse = " ")  # Join feature names into a string
#   cat(paste(entry_name, ":", feature_list, "\n\n"))  # Print the result
# }


##################### Filtered synergy (common cell lines and common genes) ########################

# results_synergy_list_filtered <- list()
# for (subset in gene_subset_names) {
#   results_synergy_list_filtered[[subset]] <- list()
#   for (i in seq_along(synergy_names)) {
#     filtered_synergy_data <- filtered_synergy_list[[subset]][[synergy_names[i]]]
#     y <- filtered_synergy_data$median_SCORE
#     y <- (y - min(y)) / (max(y) - min(y))
#     X <- filtered_synergy_data[, 5:ncol(filtered_synergy_data)]
#     X <- log2(X + 1)
#     print(synergy_names[[i]])
#     
#     result <- tryCatch({
#       run_lasso(X, y)
#     }, error = function(e) {
#       message(paste("Skipping drug:", drug, "due to error:", e$message))
#       return(NULL)
#     })
#     
#     results_synergy_list_filtered[[subset]][[synergy_names[i]]] <- result
#   }
#   
#   for (entry_name in names(results_synergy_list_filtered[[subset]])) {
#     feature_names <- results_synergy_list_filtered[[subset]][[entry_name]]$Feature  # Extract feature names
#     feature_names <- feature_names[feature_names != "(Intercept)"]  # Exclude "(Intercept)"
#     feature_list <- paste(feature_names, collapse = " ")  # Join feature names into a string
#     cat(paste(entry_name, ":", feature_list, "\n\n"))  # Print the result
#   }
# }



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
    subset_results_lasso[[drug]] <- run_model_k_splits(run_lasso, X, y, k = 5)
    
    print("ELASTIC NET: ")
    subset_results_elastic_net[[drug]] <- run_model_k_splits(run_elastic_net, X, y, k = 5)
    
    print("RIDGE REGRESSION: ")
    subset_results_ridge[[drug]] <- run_model_k_splits(run_ridge, X, y, k = 5)
    
    print("CLASSIFICATION (LASSO): ")
    subset_results_classification[[drug]] <- run_model_k_splits(run_lasso_classification, X, z, k = 5, classification = TRUE)
    
    print("SVM: ")
    print(table(z))
    subset_results_svm[[drug]] <- run_model_k_splits(run_svm, X, z, k = 5,classification = TRUE)
    
    print("KRR: ")
    subset_results_krr[[drug]] <- run_model_k_splits(run_krr, X, y, k = 5)
  }
  
  results_list_lasso[[subset]] <- subset_results_lasso
  results_list_elastic_net[[subset]] <- subset_results_elastic_net
  results_list_ridge[[subset]] <- subset_results_ridge
  results_list_classification[[subset]] <- subset_results_classification
  results_list_svm[[subset]] <- subset_results_svm
  results_list_krr[[subset]] <- subset_results_krr
}

model_list <- c("lasso", "elastic_net", "ridge", "classification", "svm", "krr")

save_results_auxiliar <- function(gene_subset){
  for (drug_name in drug_names) {
    for (model in model_list){
      results <- switch(model,
                        "lasso" = results_list_lasso,
                        "elastic_net" = results_list_elastic_net,
                        "ridge" = results_list_ridge,
                        "classification" = results_list_classification,
                        "svm" = results_list_svm,
                        "krr" = results_list_krr)
      # Define directory and file paths
      
      # Extract mean and std for every drug
      mean_values <- results[[gene_subset]][[drug_name]]$average_metrics$mean
      std_values <- results[[gene_subset]][[drug_name]]$average_metrics$std
      
      # Combine into a data frame
      df <- data.frame(
        Metric = names(mean_values),
        Mean = mean_values,
        Std = std_values
      )
      
      dir_path <- file.path("Results3", "Single_Drug", gene_subset, toupper(model))  # Convert model to uppercase for consistency
      file_path <- file.path(dir_path, paste0(drug_name, "_metrics", ".csv"))
      
      # Create directory if it doesn't exist
      if (!dir.exists(dir_path)) {
        dir.create(dir_path, recursive = TRUE)
      }
      
      # Save feature names to a text file
      write.csv(df, file_path, row.names = FALSE)
      
      # Print status message
      cat(paste("Saved:", file_path, "\n"))
    }
  }
}
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
      feature_names <- results[[gene_subset]][[drug_name]]$featrues_list[[1]]$Feature  # Extract feature names <=> gene names
      feature_names <- feature_names[feature_names != "(Intercept)"]  # Exclude "(Intercept)"

      if (length(feature_names) == 0) {
        feature_names <- "No features selected"
      }
      
      # Define directory and file paths
      dir_path <- file.path("Results6", "Single_Drug", gene_subset, toupper(model))  # Convert model to uppercase for consistency
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
      metrics_path <- file.path(dir_path, paste0(drug_name, "_metrics.csv"))
      
      mean_values <- results[[gene_subset]][[drug_name]]$average_metrics$mean
      std_values <- results[[gene_subset]][[drug_name]]$average_metrics$std
      
      # Combine into a data frame
      df <- data.frame(
        Metric = names(mean_values),
        Mean = mean_values,
        Std = std_values
      )
      
      # Save feature names to a text file
      write.csv(df, metrics_path, row.names = FALSE)
      
      # Print status message
      cat(paste("Saved metrics:", metrics_path, "\n"))
      
      #writeLines(capture.output(print(results[[gene_subset]][[drug_name]]$average_metrics)), metrics_path) #COMPROVAR!!!!!!!!!!!!!!
      #cat(paste("Saved metrics:", metrics_path, "\n"))
      
    }
  }
}

# Extract and print feature names for each gene subset
for (subset in gene_subset_names) {
  print(subset)
  save_results(subset)
}
