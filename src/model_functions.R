library(dplyr)
library(tidyr)
library(glmnet)
library(biglm)
library(caret)
# library(Metrics)  
# library(caret)  
# library(pROC)  


evaluate_regression <- function(y_true, y_pred) { 
  mse <- mean((y_true - y_pred)^2)
  rmse <- sqrt(mse)
  mae <- mean(abs(y_true - y_pred))
  r2 <- 1 - (sum((y_true - y_pred)^2) / sum((y_true - mean(y_true))^2))
  pearson <-  cor(y_true, y_pred)
  
  return(list(MSE = mse, RMSE = rmse, MAE = mae, R2 = r2, PEARSON = pearson))
}


aggregate_cv_metrics <- function(cv_results) {
  if (length(cv_results) == 0) {
    message("No cross-validation results provided.")
    return(NULL)
  }
  
  aggregated_results <- list()
  
  # Extract metrics and keep CV structure
  for (cv in names(cv_results)) {
    metrics <- cv_results[[cv]]
    metrics <- lapply(metrics, as.numeric)  # Ensure numeric values
    aggregated_results[[cv]] <- metrics
  }
  
  # Convert list to a data frame for statistical computations
  metrics_df <- do.call(rbind, lapply(aggregated_results, function(metrics) {
    as.data.frame(metrics, stringsAsFactors = FALSE)
  }))
  
  # Compute mean and standard deviation
  mean_metrics <- colMeans(metrics_df, na.rm = TRUE)
  std_metrics <- apply(metrics_df, 2, sd, na.rm = TRUE)
  
  result <- list(
    mean = mean_metrics,
    std = std_metrics,
    iterations_list = aggregated_results
  )
  
  print("Aggregation of CV metrics complete!")
  return(result)
}


prepare.data.for.ml <- function(data, drug) {
  data <- data$expression %>%
    inner_join(data$response %>% filter(Drug == drug), by = "Cell_line") %>%
    select(-Cell_line, -Drug)
  return(data)
}

get_gene_union <- function(gene_subsets, path) {
  all_genes <- character(0)
  for (subset in gene_subsets) {
    genes <- readLines(paste0(path, subset, ".txt"))
    all_genes <- union(all_genes, genes)
  }
  return(all_genes)
}


fit.linear.model <- function(X, y, model.type = 'lasso', alpha = 0.5){
  ## Wrapper function to fit linear models on 
  
  if(model.type =='lasso'){
    alpha <- 1
  } else if (model.type == 'en'){
    alpha <- alpha
  } else if (model.type == 'ridge'){
    alpha <- 0
  } else {
    stop("Invalid model type. Choose 'lasso', 'en' (elastic net), or 'ridge'.")
  }
  
  X <- as.matrix(X)  # Ensure X is a matrix
  
  model <- cv.glmnet(X, y, alpha = alpha)  # Perform cross-validation
  
  best_lambda <- model$lambda.min  # Best lambda found by CV
  model <- glmnet(X, y, alpha = alpha, lambda = best_lambda)
  
  return(list(model=model))
}



train.model <- function(model_type, train_data, drug) {
  cat("\nTraining", model_type, "model for", drug, "...\n")
  
  # Prepare training data
  X_train <- train_data %>% select(-AAC)  # Gene expression
  y_train <- train_data$AAC  # Drug response
  
  # Train model
  trained_model <- fit.linear.model(X_train, y_train, model_type)
  
  # Extract non-zero coefficient genes
  features <- coef(trained_model$model)
  features <- as.data.frame(as.matrix(features))
  features <- rownames(features)[features[, 1] != 0]
  
  return(list(model = trained_model$model, features = features))
}


evaluate.model <- function(model, test_data, drug) {
  cat("Evaluating model for", drug, "...\n")
  
  X_test <- as.matrix(test_data %>% select(-AAC))  # Gene expression
  y_test <- test_data$AAC  # Drug response
  
  # Predict on test set
  y_pred <- predict(model, newx = X_test)
  y_pred <- as.vector(y_pred)
  
  # Evaluate model performance
  eval_metrics <- evaluate_regression(y_test, y_pred)
  
  return(list(y_test = y_test, y_pred = y_pred, eval_metrics = eval_metrics))
}



compute.gene.correlations <- function(features, test_data, y_test, drug) {
  gene_corrs <- sapply(features, function(gene) {
    if (gene %in% colnames(test_data)) {
      cor(test_data[[gene]], y_test, use = "complete.obs")
    } else {
      NA
    }
  })
  return(gene_corrs)
}


perform.cv <- function(test_data, drug) {
  cross_val_results <- perform.cross.validation(10, test_data, drug)
  cross_val_avg <- aggregate_cv_metrics(cross_val_results$metrics)
  return(list(cv=cross_val_avg, features=cross_val_results$features))
}


perform.cross.validation.protein.coding <- function(n, test_data){
  # should I do PCA? Or how can I address problems of applying lm to a 700 x 19900 matrix??
  # I have already tried biglm library, as well as lm.fit... none worked...
}

perform.cross.validation <- function(n, test_data, drug){
  cross_val_results <- list()
  print("Starting cross validation")
  for (i in 1:n) {
    train_idx <- sample(1:nrow(test_data), size = 0.8 * nrow(test_data))
    
    train_subset <- test_data[train_idx, ]
    test_subset <- test_data[-train_idx, ]
    
    X_cv_train <- as.data.frame(train_subset %>% select(-AAC))
    y_cv_train <- (train_subset$AAC)
    X_cv_test <- as.data.frame(test_subset %>% select(-AAC))
    y_cv_test <- (test_subset$AAC)
    
    # Train simple linear regression
    lm_model <- lm(y_cv_train ~ ., data=X_cv_train)
    
    # Predict and evaluate
    y_cv_pred <- predict(lm_model, newdata= X_cv_test)
    
    cv_eval <- evaluate_regression(y_cv_test, y_cv_pred)
    cross_val_results[[paste0("CV_", i)]] <- cv_eval
    features <- names(coef(lm_model))
    features <- features[features != "(Intercept)"]  # Exclude "(Intercept)"
  #  gene_corrs <- compute.gene.correlations(features, train_subset, cv_eval$y_test, drug)
    
  }
  return(list(metrics = cross_val_results, features = features))
}




train.model.and.get.results <- function(source_data, drug, model_type, output_features, target_data, output_eval, output_corrs){
  # Train model
  train_data <- prepare.data.for.ml(source_data, drug)
  print(paste("Training data is", nrow(train_data), "samples and", ncol(train_data), "features"))
  model_info <- train.model(model_type, train_data, drug)
  
  # Save features
  features <- model_info$features
  features <- features[features != "(Intercept)"]  # Exclude "(Intercept)"
  writeLines(features, output_features)
  
  # Evaluate model
  test_data <- prepare.data.for.ml(target_data, drug)
  print(paste("Testing data is", nrow(test_data), "samples and", ncol(test_data), "features"))
  
  eval_info <- evaluate.model(model_info$model, test_data, drug)
  write.csv(eval_info, output_eval)
  
  # Compute gene-expression correlation
  gene_corrs <- compute.gene.correlations(features, test_data, eval_info$y_test, drug)
  write.csv(gene_corrs, output_corrs)
  
  # Extract all Pearson correlation observations
  all_pearson_observations <- cor(eval_info$y_test, eval_info$y_pred, method = "pearson")
  
  # Perform cross-validation
  #    cv <- perform.cv(test_data, drug)
  
  # WE IGNORE THIS PART OF THE CODE :)
  #    write.csv(cv$cv$mean, paste0(config$results.dir, results.subdir, model.type, "/", config$results.cv.dir, drug, "_mean.csv"))
  #    write.csv(cv$cv$std, paste0(config$results.dir, results.subdir, model.type, "/", config$results.cv.dir, drug, "_std.csv"))
  #    write.csv(cv$features, paste0(config$results.dir, results.subdir, model.type, "/", config$results.cv.dir, drug, "_features.csv"))
  # write.csv(cv$corrs, paste0(config$results.dir, results.subdir, model.type, "/", config$results.cv.dir, drug, "_correlations.csv"))
  
  # Store results
  
  new_row_model_results <- list(
    model = model_info$model,
    eval = eval_info$eval_metrics,
    features_selected = features,
    gene_correlations = gene_corrs,
    all_pearson_observations = all_pearson_observations
  )
  
  # Store dimensions of training / testing data
  new_row_dimensions <- data.frame(
    Drug = drug,
    Model_Type = model_type,
    Training_Samples = nrow(train_data),
    Training_Features = ncol(train_data),
    Testing_Samples = nrow(test_data),
    Testing_Features = ncol(test_data)
  )
  
  biomarkers_str <- paste(features, collapse = ", ")
  
  new_row_summary_results <- data.frame(
    Drug= drug,
    Model = model_type,
    MSE = eval_info$eval_metrics$MSE,
    RMSE = eval_info$eval_metrics$RMSE,
    MAE = eval_info$eval_metrics$MAE,
    R2 = eval_info$eval_metrics$R2,
    PEARSON = eval_info$eval_metrics$PEARSON,
    Biomarkers = biomarkers_str,
    stringsAsFactors = FALSE
  )
  
  return(list(new_row_model_results = new_row_model_results, new_row_summary_results = new_row_summary_results, new_row_dimensions = new_row_dimensions))
}

one.fold.linear.model <- function(train_subset, test_subset, drug){
  
  x_train <- as.data.frame(train_subset %>% select(-AAC))

  y_train <- train_subset$AAC
  
  x_test <- as.data.frame(test_subset %>% select(-AAC))
  y_test <- test_subset$AAC
  
  model <- lm(y_train ~ ., data=x_train)
  
  y_pred <- predict(model, newdata = x_test)
  
  return(evaluate_regression(y_test, y_pred))
}

k.fold.linear.model <- function(data_to_split, drug, num_folds) {
  data <- prepare.data.for.ml(data_to_split, drug)
  responses <- data$AAC
  cv_folds <- createFolds(responses, k = num_folds, list = TRUE)
  
  results <- data.frame(Drug = character(), MSE = numeric(), RMSE = numeric(), MAE = numeric(), R2 = numeric(), PEARSON = numeric())
  for (fold_idx in seq_along(cv_folds)) {
    
    # Split data for this fold
    train_indices <- unlist(cv_folds[-fold_idx])  # All but one fold for training
    test_indices <- unlist(cv_folds[fold_idx])  # The held-out fold for testing
    
    train_data <- data[train_indices, ]
    test_data <- data[test_indices, ]
    
    fold_results <- one.fold.linear.model(train_data, test_data, drug)
    results <- rbind(results, fold_results)
  }
  
  numeric_results <- results[, sapply(results, is.numeric)]
  
  # Compute column-wise mean and standard deviation
  mean_results <- colMeans(numeric_results, na.rm = TRUE)
  sd_results <- apply(numeric_results, 2, sd, na.rm = TRUE)
  
  # Create summary dataframe
  summary_results <- data.frame(
    Drug = drug,
    MSE_Mean = mean_results["MSE"],
    MSE_SD = sd_results["MSE"],
    RMSE_Mean = mean_results["RMSE"],
    RMSE_SD = sd_results["RMSE"],
    MAE_Mean = mean_results["MAE"],
    MAE_SD = sd_results["MAE"],
    R2_Mean = mean_results["R2"],
    R2_SD = sd_results["R2"],
    PEARSON_Mean = mean_results["PEARSON"],
    PEARSON_SD = sd_results["PEARSON"]
  )
  
  return(summary_results)
}
evaluate.feature.set <- function(path_to_feature_set, data_to_split, drug, num_folds) {
  if (!file.exists(path_to_feature_set)) {
    cat("Error: The file at path", path_to_feature_set, "does not exist.\n")
    return(list(
      summary = data.frame(
        Drug = drug,
        MSE_Mean = -1000,
        MSE_SD = -1000,
        RMSE_Mean = -1000,
        RMSE_SD = -1000,
        MAE_Mean = -1000,
        MAE_SD = -1000,
        R2_Mean = -1000,
        R2_SD = -1000,
        PEARSON_Mean = -1000,
        PEARSON_SD = -1000
      ),
      fold_results = list()
    ))
  }

  feature_set <- readLines(path_to_feature_set)
  feature_set <- unlist(feature_set)  # Ensure feature_set is a character vector
  filtered_data <- list()
  filtered_data$expression <- data_to_split$expression %>%
    select(any_of(c("Cell_line", feature_set)))
  filtered_data$response <- data_to_split$response

  if (is.null(filtered_data) || !("expression" %in% names(filtered_data)) || !("response" %in% names(filtered_data))) {
    stop("filtered_data is NULL or does not contain the required 'expression' and 'response' components for drug:", drug) 
  }
  
  data <- prepare.data.for.ml(filtered_data, drug)
  responses <- data$AAC
  cv_folds <- createFolds(responses, k = num_folds, list = TRUE)
  
  results <- data.frame(Drug = character(), MSE = numeric(), RMSE = numeric(), MAE = numeric(), R2 = numeric(), PEARSON = numeric())
  fold_results_list <- list()
  
  for (fold_idx in seq_along(cv_folds)) {
    tryCatch({
      # Split data for this fold
      train_indices <- unlist(cv_folds[-fold_idx])  # All but one fold for training
      test_indices <- unlist(cv_folds[fold_idx])  # The held-out fold for testing
      
      train_data <- data[train_indices, ]
      test_data <- data[test_indices, ]
      
      # Prepare training and testing data
      X_train <- as.matrix(train_data %>% select(-AAC))
      y_train <- train_data$AAC
      X_test <- as.matrix(test_data %>% select(-AAC))
      y_test <- test_data$AAC
      
      # Ensure X_train has at least two columns
      if (ncol(X_train) < 2) {
        stop("X_train must have at least two columns for glmnet.")
      }
      
      # Fit Lasso model
      lasso_model <- cv.glmnet(X_train, y_train, alpha = 1)  # Lasso regression
      best_lambda <- lasso_model$lambda.min
      final_model <- glmnet(X_train, y_train, alpha = 1, lambda = best_lambda)
      
      # Predict on test data
      y_pred <- predict(final_model, newx = X_test)
      y_pred <- as.vector(y_pred)
      
      # Evaluate model performance
      eval_metrics <- evaluate_regression(y_test, y_pred)
      
      # Append results for this fold
      fold_results <- data.frame(
        Drug = drug,
        MSE = eval_metrics$MSE,
        RMSE = eval_metrics$RMSE,
        MAE = eval_metrics$MAE,
        R2 = eval_metrics$R2,
        PEARSON = eval_metrics$PEARSON
      )
      results <- rbind(results, fold_results)
      fold_results_list[[paste0("Fold_", fold_idx)]] <- eval_metrics
    }, error = function(e) {
      cat("Training stopped for drug:", drug, "at fold:", fold_idx, "\n")
      # Append -1 for all metrics in case of failure
      fold_results <- data.frame(
        Drug = drug,
        MSE = -1000,
        RMSE = -1000,
        MAE = -1000,
        R2 = -1000,
        PEARSON = -1000
      )
      results <- rbind(results, fold_results)
      fold_results_list[[paste0("Fold_", fold_idx)]] <- list(MSE = -1000, RMSE = -1000, MAE = -1000, R2 = -1000, PEARSON = -1000)
    })
  }
  
  # Compute mean and standard deviation of metrics
  numeric_results <- results[, sapply(results, is.numeric)]
  mean_results <- colMeans(numeric_results, na.rm = TRUE)
  sd_results <- apply(numeric_results, 2, sd, na.rm = TRUE)
  
  # Create summary dataframe
  summary_results <- data.frame(
    Drug = drug,
    MSE_Mean = mean_results["MSE"],
    MSE_SD = sd_results["MSE"],
    RMSE_Mean = mean_results["RMSE"],
    RMSE_SD = sd_results["RMSE"],
    MAE_Mean = mean_results["MAE"],
    MAE_SD = sd_results["MAE"],
    R2_Mean = mean_results["R2"],
    R2_SD = sd_results["R2"],
    PEARSON_Mean = mean_results["PEARSON"],
    PEARSON_SD = sd_results["PEARSON"]
  )
  
  return(list(summary = summary_results, fold_results = fold_results_list))
}

compute.cv.for.pagerank.input <- function(path_pagerank_genes, data, drug, num_folds) {
  pagerank_genes <- load_pagerank_feature_list(path_pagerank_genes)
  pagerank_genes <- pagerank_genes[pagerank_genes %in% colnames(source_data$expression)]
  if (length(pagerank_genes) == 0 || is.null(pagerank_genes)) {
    # Handle the case when gene_list is NULL or empty
    cat("gene_list is empty or NULL\n")
    return(NULL)
  } else {
    filtered_data <- list()
    filtered_data$expression <- data$expression %>% select(c("Cell_line", pagerank_genes))
    filtered_data$response <- data$response
    #results_drug <- k.fold.linear.model(filtered_data, drug, num_folds)
    results_drug <- evaluate.feature.set(filtered_data, drug, num_folds)
    results_drug <- results_drug$summary
    return(results_drug)
  }

}

compute.cv.for.random.input <- function(source_data, target_data, drug, num_folds, num_genes) {
  random_genes_to_keep <- sample(2:ncol(target_data$expression), num_genes)

  filtered_source$expression <- source_data$expression[, c("Cell_line", names(source_data$expression)[random_genes_to_keep])]
  filtered_target$expression <- target_data$expression[, c("Cell_line", names(target_data$expression)[random_genes_to_keep])]
  filtered_source$response <- source_data$response
  filtered_target$response <- target_data$response
  
  results_drug <- k.fold.linear.model(filtered_target, drug, num_folds)
  
  return(results_drug)
}


process.results.pagerank <- function(path_pagerank_output, drug, feature_size, target_data, num_folds) {
  pagerank_output <- paste0(path_pagerank_output, "_", feature_size, ".txt")
  results_cv <- compute.cv.for.pagerank.input(pagerank_output, target_data, drug, num_folds) 
  return(results_cv)
}


# pagerank_cv_helpers.R

init.results <- function(results, methods, screens, sizes) {
  for (method in methods) {
    for (size in sizes) {
      if (method == "lasso") {
        for (screen in screens) {
          key <- paste0("lasso_", screen, "_", size)
          results[[key]] <- list()
        }
      } else {
        key <- paste0(method, "_", screen, "_", size)
        results[[key]] <- list()
      }
    }
  }
  return(results)
}

init.results.wo.methods.and.sizes <- function(method, results, screens) {
  results <-  list()
  for (screen in screens) {
    key <- paste0(method, "_", screen, "_derived")
    results[[key]] <- list()
  }
  return(results)
}



#given a path gene_set_file_path with a subset of genes, performs 10-fold cv per drug
run.cv.gene.set.per.drug.ppr <- function(results, screen, drugs, gene_set_name, gene_set_file_path, target_data, config, num_folds) {
  
  for (drug in drugs) {
  #cat("Processing drug:", drug, "with size:", size, "\n")  #  ADD THIS
  if (grepl("\\.txt$", gene_set_file_path)) {
    gene_set_path <- gene_set_file_path
  } else {
    gene_set_path <- file.path(gene_set_file_path, paste0(drug, "_50.txt"))
  }
  result <- evaluate.feature.set(gene_set_path, target_data, drug, num_folds)
  
  # Extract summary and per-fold results
  summary_result <- result$summary
  fold_results <- result$fold_results
  
  # Store summary results
  key <- paste0(gene_set_name, "_", screen)
  results[[key]] <- rbind(results[[key]], summary_result)
  
  # Store per-fold results
  fold_key <- paste0(gene_set_name, "_", screen, "_folds")
  if (is.null(results[[fold_key]])) {
    results[[fold_key]] <- list()
  }
  results[[fold_key]][[drug]] <- fold_results
  }
  
  return(results)
}

#given a path gene_set_file_path with a subset of genes, performs 10-fold cv per drug
run.cv.gene.set.per.drug <- function(results, screen, drugs, gene_set_name, gene_set_file_path, target_data, config, num_folds) {
  
  for (drug in drugs) {
  #cat("Processing drug:", drug, "with size:", size, "\n")  #  ADD THIS
  if (grepl("\\.txt$", gene_set_file_path)) {
    gene_set_path <- gene_set_file_path
  } else {
    gene_set_path <- file.path(gene_set_file_path, paste0(drug, ".txt"))
  }
  result <- evaluate.feature.set(gene_set_path, target_data, drug, num_folds)
  
  # Extract summary and per-fold results
  summary_result <- result$summary
  fold_results <- result$fold_results
  
  # Store summary results
  key <- paste0(gene_set_name, "_", screen)
  results[[key]] <- rbind(results[[key]], summary_result)
  
  # Store per-fold results
  fold_key <- paste0(gene_set_name, "_", screen, "_folds")
  if (is.null(results[[fold_key]])) {
    results[[fold_key]] <- list()
  }
  results[[fold_key]][[drug]] <- fold_results
  }
  
  return(results)
}
plot_pearson_distribution_filtering <- function(results, screen, gene_set_name, output_file) {
  # Extract fold results for the specified screen and gene set
  fold_key <- paste0(gene_set_name, "_", screen, "_folds")
  fold_results <- results[[fold_key]]
  
  if (is.null(fold_results)) {
    stop("No fold results found for the specified screen and gene set.")
  }
  
  # Prepare data for plotting
  plot_data <- do.call(rbind, lapply(names(fold_results), function(drug) {
    pearson_values <- sapply(fold_results[[drug]], function(fold) fold$PEARSON)
    if (length(pearson_values) > 0) {
      data.frame(
        Drug = drug,
        Pearson = pearson_values
      )
    } else {
      NULL
    }
  }))
  
  if (is.null(plot_data) || nrow(plot_data) == 0) {
    stop("No valid Pearson values found for plotting.")
  }
  
  # Create the plot
  library(ggplot2)
  plot <- ggplot(plot_data, aes(x = Drug, y = Pearson, fill = Drug)) +
    geom_boxplot(alpha = 0.7) +
    theme_minimal() +
    labs(
      title = paste("Distribution of Pearson values across folds for", gene_set_name, "feature set in", screen),
      x = "Drug",
      y = "Pearson Correlation"
    ) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  # Save the plot
  ggsave(filename = output_file, plot = plot, width = 12, height = 6)
}


plot_pearson_distribution <- function(results, screen, gene_set_name) {
  # Extract fold results for the specified screen and gene set
  fold_key <- paste0(gene_set_name, "_", screen, "_folds")
  fold_results <- results[[fold_key]]
  
  if (is.null(fold_results)) {
    stop("No fold results found for the specified screen and gene set.")
  }
  
  # Prepare data for plotting
  plot_data <- do.call(rbind, lapply(names(fold_results), function(drug) {
    data.frame(
      Drug = drug,
      Pearson = sapply(fold_results[[drug]], function(fold) fold$PEARSON)
    )
  }))
  
  # Create the plot
  library(ggplot2)
  ggplot(plot_data, aes(x = Drug, y = Pearson, fill = Drug)) +
    geom_boxplot(alpha = 0.7) +
    theme_minimal() +
    labs(
      title = "Distribution of Pearson Values Across Folds",
      x = "Drug",
      y = "Pearson Correlation"
    ) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
}

run.lasso.cv <- function(results, screen, sizes, drugs, source_data, target_data, config, num_folds) {
  for (drug in drugs) {
    for (size in sizes) {
      #cat("Processing drug:", drug, "with size:", size, "\n")  
      path <- file.path(paste0(config$results.pagerank.output, "lasso"), screen, "positive", drug)
      result <- process.results.pagerank(path, drug, size, target_data, num_folds)
      key <- paste0("lasso_", screen, "_", size)
      results[[key]] <- rbind(results[[key]], result)
    }
  }
  return(results)
}

#method is wither lasso, ridge or en
run.lasso.derived.cv <- function(method, results, screen, drugs, source_data, target_data, num_folds) {
  for (drug in drugs) {
    #cat("Processing drug:", drug, "\n")  
    path <- file.path("results", screen, "positive", method, "features", paste0(drug, ".txt"))
    result <- compute.cv.for.pagerank.input(path, target_data, drug, num_folds)
    key <- paste0(method, "_", screen, "_derived")
    results[[key]] <- rbind(results[[key]], result)
  }
  return(results)
}
write.results.lasso.en.ridge <- function(results, screens, method, config) {
  for (screen in screens) {
    key <- paste0(method, "_", screen)
    out_path <- paste0("results/cv_performance/", method, "/", screen, "/performance_", method, "_lasso_en_ridge.csv")
    
    # Select only the required columns
    filtered_results <- results[[key]][, c("Drug", "PEARSON_Mean", "PEARSON_SD", "RMSE_Mean", "RMSE_SD")]
    
    write.csv(filtered_results, out_path, row.names = FALSE)
    print(paste0("Written results of 10-fold CV for ", method, " in screen ", screen, " at ", out_path))
  }
}
#methos is ppr-lasso, ppr-db, ppr-dtc
#method is drugbank or dtc
write.results.ppr <- function(results, screens, method, config) {
  for (screen in screens) {
    key <- paste0(method, "_", screen)
    out_path <- paste0("results/cv_performance/PPR/", method, "/", screen, "/performance_", method, ".csv")
    
    # Select only the required columns
    filtered_results <- results[[key]][, c("Drug", "PEARSON_Mean", "PEARSON_SD", "RMSE_Mean", "RMSE_SD")]
    
    write.csv(filtered_results, out_path, row.names = FALSE)
    print(paste0("Written results of 10-fold CV for ", method, " in screen ", screen, " at ", out_path))
  }
}

#method is drugbank or dtc
write.results.drug.targets <- function(results, screens, method, config) {
  for (screen in screens) {
    key <- paste0(method, "_", screen)
    out_path <- paste0("results/cv_performance/", method, "/", screen, "/performance_", method, "_drug_targets.csv")
    
    # Select only the required columns
    filtered_results <- results[[key]][, c("Drug", "PEARSON_Mean", "PEARSON_SD", "RMSE_Mean", "RMSE_SD")]
    
    write.csv(filtered_results, out_path, row.names = FALSE)
    print(paste0("Written results of 10-fold CV for ", method, " in screen ", screen, " at ", out_path))
  }
}
write.results.centrality <- function(results, screens, centrality_measures, config) {
  for (screen in screens) {
    for (measure in centrality_measures) {
      key <- paste0(measure, "_", screen)
      out_path <- file.path("results", "cv_performance", "centrality", screen, paste0("performance_", measure, "_centrality.csv"))
      dir.create(dirname(out_path), recursive = TRUE, showWarnings = FALSE)
      write.csv(results[[key]], out_path, row.names = FALSE)
      print(paste0("Written results of 10-fold CV for ", measure, " centrality in screen ", screen, " at ", out_path))
    }
  }
}


write.results.drug.targets2 <- function(results, screens, drug_target_source, config) {
  for (screen in screens) {
    key <- paste0(drug_target_source, "_", screen)
    out_path <- paste0("results/cv_performance/", drug_target_source, "/", screen, "/performance_", drug_target_source, "_drug_targets.csv")
    write.csv(results[[key]], out_path, row.names = FALSE)
    print(paste0("Written results of 10-fold CV for ", drug_target_source, " drug targets in screen ", screen," at ", out_path))
  }
}

write.results.central.genes <- function(results, screens, central_gene_sets, config) {
  for (screen in screens) {
    for (gene_set_name in names(central_gene_sets)) {
      key <- paste0(gene_set_name, "_", screen)
      out_path <- paste0("results/cv_performance/central_genes/", screen, "/performance_", gene_set_name, "_central_genes.csv")
      write.csv(results[[key]], out_path, row.names = FALSE)
      print(paste0("Written results of 10-fold CV for ", gene_set_name, " central genes in screen ", screen, " at ", out_path))
    }
  }
}

write.lasso.results <- function(results, screens, sizes, config) {
  for (screen in screens) {
    for (size in sizes) {
      key <- paste0("lasso_", screen, "_", size)
      out_path <- file.path(config$results.cv.lasso.subdir, screen, paste0("lasso_performance_top_", size, "_features.csv"))
      print(out_path)
      write.csv(results[[key]], out_path, row.names = FALSE)
    }
  }
}

write.lasso.derived.results <- function(method, results, screens, config) {
  for (screen in screens) {
    key <- paste0(method, "_", screen, "_derived")
    out_path <- file.path("results", "cv_performance", paste0(method, "_derived"), screen, paste0(method, "_performance_derived_features.csv"))
    print(out_path)
    write.csv(results[[key]], out_path, row.names = FALSE)
  }
}

run.other.cv <- function(results, sizes, drugs, screen, screen_data, config, num_folds) {
  for (drug in drugs) {
    for (size in sizes) {
      # Drugbank
      db_path <- paste0(config$results.pagerank.output, config$drugbank.subfolder, drug)
      result_db <- process.results.pagerank(db_path, drug, size, screen_data, num_folds) #dbpath = "results/pagerank_output/drugbank/drugwise"
      results[[paste0("drugbank_", screen, "_", size)]] <- rbind(results[[paste0("drugbank_", screen, "_", size)]], result_db)
      
      # DTC
      dtc_path <- paste0(config$results.pagerank.output, config$dtc.subfolder, drug)
      result_dtc <- process.results.pagerank(dtc_path, drug, size, screen_data, num_folds) 
      results[[paste0("dtc_", screen, "_", size)]] <- rbind(results[[paste0("dtc_", screen, "_", size)]], result_dtc)
      
      # Random
      result_rand <- compute.cv.for.random.input(source_data, screen_data, drug, num_folds, size)
      results[[paste0("random_", screen, "_", size)]] <- rbind(results[[paste0("random_", screen, "_", size)]], result_rand)
    }
  }
  return(results)
}

write.other.results <- function(results, sizes, config, screen) {
  for (size in sizes) {
    write.csv(results[[paste0("drugbank_", screen, "_", size)]],
              file.path(config$results.cv.drugbank.subdir, screen, paste0("drugbank_performance_top_", size, "_features.csv")),
              row.names = FALSE)
    
    write.csv(results[[paste0("dtc_", screen, "_", size)]],
              file.path(config$results.cv.dtc.subdir, screen, paste0("dtc_performance_top_", size, "_features.csv")),
              row.names = FALSE)
    
    write.csv(results[[paste0("random_", screen, "_", size)]],
              file.path(config$results.cv.random.subdir, screen, paste0("random_performance_", size, "_features.csv")),
              row.names = FALSE)
  }
  write.csv(results[[paste0("dtc_", screen, "_all")]],
            file.path(config$results.cv.dtc.subdir, screen, "dtc_performance_all_gene_set.csv"),
            row.names = FALSE)
  write.csv(results[[paste0("drugbank_", screen, "_all")]],
            file.path(config$results.cv.drugbank.subdir, screen, "drugbank_performance_all_gene_set.csv"),
            row.names = FALSE)
}
