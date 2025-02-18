run_lasso <- function(X, Y, lambda_value = 10) {
  # Load necessary library
  library(glmnet)
  
  # Convert X to matrix (required by glmnet)
  X <- as.matrix(X)
  
  # Fit Lasso model
  lasso_model <- glmnet(X, Y, alpha = 1, lambda = lambda_value)
  
  # Extract coefficients
  coefficients <- coef(lasso_model)
  
  # Convert to data frame
  coefs_df <- data.frame(
    Feature = rownames(coefficients),
    Coefficient = as.numeric(coefficients)
  )
  
  # Filter non-zero coefficients
  non_zero_coefs <- coefs_df[coefs_df$Coefficient != 0, ]
  
  # Print and return non-zero coefficients
  print(non_zero_coefs)
  return(non_zero_coefs)
}

# drugs to analyze individually
drug_names <- c("Bortezomib", "Axitinib", "Cisplatin", "Erlotinib", "Lapatinib", "Ruxolitinib", "Sirolimus", "Vinorelbine", "Vorinostat")

##################### LASSO MODEL FOR SINGLE DRUG #############################
drug_data_list <- list()
for (drug in drug_names) {
  drug_data_list[[drug]] <- read.csv(paste0("data/single_drug/", drug, ".csv"))
}

# read synergy
#synergy1 # Bortezomib, Erlotinib
#synergy2 # Cisplatin, Lapatinib
#synergy3 # Sirolimus, Axitinib
#synergy4 # Vorinostat, Ruxolitinib
#synergy5 # Vorinostat, Vinorelbine

# List of synergy file names
synergy_files <- paste0("data/synergy/synergy", 1:5, ".csv")
synergy_names <- paste0("synergy", 1:5)

# Read synergy data dynamically
synergy_data_list <- list()
for (i in seq_along(synergy_names)) {
  synergy_data_list[[synergy_names[i]]] <- read.csv(synergy_files[i])
}

library(dplyr)

# Find common genes 
common_cols <- intersect(names(synergy_data_list$synergy1), names(drug_data_list$Bortezomib))

# Filter synergy data by common genes with single drug data
filtered_synergy_list <- list()
for (name in synergy_names) {
  filtered_synergy_list[[name]] <- synergy_data_list[[name]] %>%
    select(sampleid, treatment1id, treatment2id, median_SCORE, all_of(common_cols))
}

# Filter single drug data by common genes with synergy data
filtered_drug_data_list <- list()
for (drug in drug_names) {
  filtered_drug_data_list[[drug]] <- drug_data_list[[drug]] %>%
    select(sampleid, aac, all_of(common_cols))
}

# Only cell lines in common
for (i in seq_along(synergy_names)) {
  synergy_key <- synergy_names[i]
  drug1 <- synergy_data_list[[synergy_key]]$treatment1id[1]  # First drug in synergy pair
  drug2 <- synergy_data_list[[synergy_key]]$treatment2id[1]  # Second drug in synergy pair
  
  # Keep only common cell lines
  filtered_synergy_list[[synergy_key]] <- filtered_synergy_list[[synergy_key]] %>%
    filter(sampleid %in% filtered_drug_data_list[[drug1]]$cell_line)
  
  filtered_drug_data_list[[drug1]] <- filtered_drug_data_list[[drug1]] %>%
    filter(cell_line %in% filtered_synergy_list[[synergy_key]]$sampleid)
  
  filtered_drug_data_list[[drug2]] <- filtered_drug_data_list[[drug2]] %>%
    filter(cell_line %in% filtered_synergy_list[[synergy_key]]$sampleid)
}

############################################## LASSO MODEL ######################################

############################################## Synergy ##########################################

results_synergy_list <- list()

for (i in seq_along(synergy_names)) {
  synergy_data <- filtered_synergy_list[[synergy_names[i]]]  
  y <- synergy_data$median_SCORE
  X <- synergy_data[, 5:ncol(synergy_data)]
  results_synergy_list[[synergy_names[i]]] <- run_lasso(X, y)
}

############################################## Single drug ##########################################

results_list <- list()

for (drug in drug_names) {
  drug_data <- filtered_drug_data_list[[drug]] 
  y <- drug_data$AAC
  X <- drug_data[, 4:ncol(drug_data)]
  X <- log2(X + 1)
  result <- run_lasso(X, y)
  results_list[[drug]] <- result
}


