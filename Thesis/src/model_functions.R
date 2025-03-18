library(dplyr)
library(tidyr)
library(glmnet)
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




fit.linear.model <- function(X,y, model.type = 'lasso',alpha = 0.5){
	## Wrapper function to fit linear models on 
	# ensure lower case
	model.type <- toLower(type)
	
	if(type =='lasso'){
		alpha <- 1
	} else if (model.type == 'en'){
		alpha <- alpha
	} else if (model.type == 'ridge'){
		alpha <- 0
	}
	
	X <- as.matrix(X)  # Ensure X is a matrix
  
  	model <- cv.glmnet(X, y, alpha = alpha)  # Perform cross-validation

 	best_lambda <- model$lambda.min  # Best lambda found by CV
 	model <- glmnet(X, y, alpha = alpha, lambda = best_lambda)
 	return(list(model=model))
}
