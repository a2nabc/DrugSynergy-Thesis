# Load necessary libraries
library(dplyr)
library(readr)

# Load drug response data
broad_resp <- read.csv("../data/processed/Broad/responses.csv")
gcsi_resp <- read.csv("../data/processed/gCSI/responses.csv")
gdsc_resp <- read.csv("../data/processed/GDSC2/responses.csv")

# Extract unique drug lists
broad_drugs <- unique(broad_resp$Drug)
gcsi_drugs <- unique(gcsi_resp$Drug)
gdsc_drugs <- unique(gdsc_resp$Drug)

# Find drugs common across all datasets
common_drugs <- Reduce(intersect, list(broad_drugs, gcsi_drugs, gdsc_drugs))

# Print common drugs
cat("Common drugs across datasets:\n", paste(common_drugs, collapse=", "), "\n")

# Save to file for reference
writeLines(common_drugs, "../data/common_drugs.txt")
write.csv(data.frame(Drug = common_drugs), "../data/processed/common_drugs.csv", row.names = FALSE)