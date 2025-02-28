library(dplyr)
library(tidyr)
library(glmnet)

run_lasso <- function(X, y) {
  X <- as.matrix(X)  # Ensure X is a matrix
  lasso_cv <- cv.glmnet(X, y, alpha = 1)  # Perform cross-validation
  best_lambda <- lasso_cv$lambda.min  # Best lambda found by CV
  print(paste("Optimal lambda:", best_lambda))
  
  # Fit Lasso with optimal lambda
  lasso_model <- glmnet(X, y, alpha = 1, lambda = best_lambda)
  
  # Extract coefficients
  coefficients <- coef(lasso_model)
  coefs_df <- data.frame(
    Feature = rownames(coefficients),
    Coefficient = as.numeric(coefficients)
  )
  
  # Filter non-zero coefficients
  non_zero_coefs <- coefs_df[coefs_df$Coefficient != 0, ]
  
  print(non_zero_coefs)
  return(non_zero_coefs)
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
  
  # Extract coefficients and filter non-zero ones
  coefficients <- coef(elastic_model)
  
  # Convert to data frame
  coefs_df <- data.frame(
    Feature = rownames(coefficients),
    Coefficient = as.numeric(coefficients)
  )
  
  non_zero_coefs <- coefs_df[coefs_df$Coefficient != 0, ]
  
  print(non_zero_coefs)
  return(non_zero_coefs)
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
  
  # Extract coefficients
  coefficients <- coef(lasso_model)
  coefs_df <- data.frame(
    Feature = rownames(coefficients),
    Coefficient = as.numeric(coefficients)
  )
  
  # Filter non-zero coefficients
  non_zero_coefs <- coefs_df[coefs_df$Coefficient != 0, ]
  
  print(non_zero_coefs)
  
  return(non_zero_coefs)
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
results_list_classification <- list()

for (subset in gene_subset_names) {
  subset_results_lasso <- list()
  subset_results_elastic_net <- list()
  subset_results_classification <- list()
  
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
    
    print("CLASSIFICATION: ")
    result_classification <- tryCatch({
      run_lasso_classification(X, z)
    }, error = function(e) {
      message(paste("Skipping classification for drug:", drug, "due to error:", e$message))
      return(NULL)
    })
    subset_results_classification[[drug]] <- result_classification
  }
  
  results_list_lasso[[subset]] <- subset_results_lasso
  results_list_elastic_net[[subset]] <- subset_results_elastic_net
  results_list_classification[[subset]] <- subset_results_classification
}

# Extract and print feature names for each gene subset
for (subset in names(results_list_lasso)) {
  print(subset)
  for (entry_name in names(results_list_lasso[[subset]])) {
    feature_names <- results_list_lasso[[subset]][[entry_name]]$Feature  # Extract feature names <=> gene names
    feature_names <- feature_names[feature_names != "(Intercept)"]  # Exclude "(Intercept)"
    feature_list <- paste(feature_names, collapse = " ")  # Join feature names into a string
    cat(paste(entry_name, "LASSO", ":", feature_list, "\n\n"))  
    
    feature_names <- results_list_elastic_net[[subset]][[entry_name]]$Feature  # Extract feature names <=> gene names
    feature_names <- feature_names[feature_names != "(Intercept)"]  # Exclude "(Intercept)"
    feature_list <- paste(feature_names, collapse = " ")  # Join feature names into a string
    cat(paste(entry_name, "ELASTIC_NET", ":", feature_list, "\n\n"))  
    
    feature_names <- results_list_classification[[subset]][[entry_name]]$Feature  # Extract feature names <=> gene names
    feature_names <- feature_names[feature_names != "(Intercept)"]  # Exclude "(Intercept)"
    feature_list <- paste(feature_names, collapse = " ")  # Join feature names into a string
    cat(paste(entry_name, "CLASSIFICATION", ":", feature_list, "\n\n"))  
  }
}
