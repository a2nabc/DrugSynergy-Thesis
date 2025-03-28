library(dplyr)
library(tidyr)
library(glmnet)
library(biglm)
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
