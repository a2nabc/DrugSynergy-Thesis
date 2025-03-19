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
  
  return(list(MSE = mse, RMSE = rmse, MAE = mae, R2 = r2))
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


prepare.data.for.ml <- function(data, drug) {
  data <- data$expression %>%
    inner_join(data$response %>% filter(Drug == drug), by = "Cell_line") %>%
    select(-Cell_line, -Drug)
}

get_gene_union <- function(gene_subsets, path) {
  all_genes <- character(0)
  for (subset in gene_subsets) {
    genes <- readLines(paste0(path, subset, ".txt"))
    all_genes <- union(all_genes, genes)
  }
  return(all_genes)
}


perform.cross.validation.protein.coding <- function(n, test_data){
  # should I do PCA? Or how can I address problems of applying lm to a 700 x 19900 matrix??
  # I have already tried biglm library, as well as lm.fit... none worked...
}

perform.cross.validation <- function(n, test_data){
  cross_val_results <- list()
  
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
  }
  
  return(cross_val_results)
}
