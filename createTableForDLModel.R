library(PharmacoGx)
library(SummarizedExperiment)
library(dplyr)     # For the pipe operator and data manipulation
library(tidyr)     # For the `pivot_longer` function
library(tibble)   

NCI_ALMANAC_2017 <- readRDS("NCI_ALMANAC_2017.rds")
NCI_gene_expression <- readRDS("PSet_NCI60.rds")
NCI_gene_expression <- updateObject(NCI_gene_expression) # I did this because an error line prompted


######################### 1 obtain sampleid, drug1, drug2, synergy score from NCI_ALMANAC###########################

# Extract synergy data from NCI_ALMANAC_2017
synergy_data <- NCI_ALMANAC_2017@treatmentResponse$profiles  # Contains synergy scores

# Rename and select relevant columns
synergy_table <- synergy_data %>%
  dplyr::rename(
    sampleid = sampleid,      # Keep sampleid as is
    drugid1 = treatment1id,   # Rename treatment1id to drugid1
    drugid2 = treatment2id,   # Rename treatment2id to drugid2
    synergy_score = SCORE     # Rename SCORE to synergy_score
  ) %>%
  dplyr::select(sampleid, drugid1, drugid2, synergy_score)

# Clean incorrect rows --> no drugid2 or no drugid1 or NA synergy SCORE
cleaned_synergy_table <- synergy_table %>%
  dplyr::filter(!is.na(drugid1) & !is.na(drugid2) & !is.na(synergy_score) & !is.na(sampleid))


###################### 2 Obtain IC50 from monotherapy dataset on the corresponding samples #######################

# Extract treatment response data from NCI_gene_expression
treatment_info <- NCI_gene_expression@treatmentResponse$info  # Contains sampleID and NSC
treatment_profiles <- NCI_gene_expression@treatmentResponse$profiles  # Contains ic50_recomputed

# Step 1: Extract IC50 data from treatment_profiles
ic50_data <- treatment_profiles %>%
  dplyr::select(ic50_recomputed) %>%
  dplyr::rename(ic50 = ic50_recomputed)

# Step 2: Add sampleid and treatmentid to ic50_data from treatment_info
ic50_data <- cbind(treatment_info[, c("sampleid", "treatmentid")], ic50_data)

# Clean the rows with NA values or with Inf values
cleaned_ic50_data <- ic50_data %>%
  dplyr::filter(!is.na(sampleid) & !is.na(treatmentid) & !is.na(ic50) & !is.infinite(ic50))


############################### CHECK DUPLICATES IN BOTH TABLES BEFORE JOIN ###############################
# Check for duplicates in ic50_data
ic50_duplicates <- cleaned_ic50_data %>% 
  dplyr::group_by(sampleid, treatmentid) %>% 
  dplyr::filter(n() > 1)

# Check for duplicates in cleaned_synergy_table
synergy_duplicates <- cleaned_synergy_table %>%
  dplyr::group_by(sampleid, drugid1, drugid2) %>%
  dplyr::filter(n() > 1)

# I confirmed that they are not unique values so I decided to take the mean:

# Resolve duplicates in ic50_data by taking the mean IC50
ic50_data_mean <- cleaned_ic50_data %>%
  dplyr::group_by(sampleid, treatmentid) %>%
  dplyr::summarize(ic50 = mean(ic50, na.rm = TRUE), .groups = "drop")

# Resolve duplicates in cleaned_synergy_table by taking the mean synergy_score
synergy_table_mean <- cleaned_synergy_table %>%
  dplyr::group_by(sampleid, drugid1, drugid2) %>%
  dplyr::summarize(synergy_score = mean(synergy_score, na.rm = TRUE), .groups = "drop")

 ############################### CHECK DUPLICATES MEAN ###############################
# Check for duplicates in ic50_data
ic50_duplicates_mean <- ic50_data_mean %>% 
  dplyr::group_by(sampleid, treatmentid) %>% 
  dplyr::filter(n() > 1)

# Check for duplicates in cleaned_synergy_table
synergy_duplicates_mean <- synergy_table_mean %>%
  dplyr::group_by(sampleid, drugid1, drugid2) %>%
  dplyr::filter(n() > 1)

# Now it seems ok

# Step 3: Filter the data for drugid1 (drug 1) and drugid2 (drug 2) from synergy table
synergy_with_ic50 <- synergy_table_mean %>%
  dplyr::filter(drugid1 != drugid2) %>% 
  dplyr::inner_join(
    ic50_data_mean %>% dplyr::rename(ic50_drugid1 = ic50), 
    by = c("drugid1" = "treatmentid", "sampleid" = "sampleid")  # Match drugid1 and sampleid
  ) %>%
  dplyr::inner_join(
    ic50_data_mean %>% dplyr::rename(ic50_drugid2 = ic50), 
    by = c("drugid2" = "treatmentid", "sampleid" = "sampleid")  # Match drugid2 and sampleid
  )

write.csv(synergy_with_ic50, "synergy_with_ic50.csv", row.names = FALSE)

###################### 3 Obtain gene expression ########################

# Extract RNA data and convert to data frame
rna_data <- NCI_gene_expression@molecularProfiles$rna
rna_df <- as.data.frame(assay(rna_data))  # Extract assay data (gene expression)

# Gene expression data: sample IDs as columns, genes as rows
rna_long <- rna_df %>%
  rownames_to_column(var = "gene") %>%  # Add gene names as a column
  pivot_longer(cols = -gene,           # All columns except 'gene' are sample IDs
               names_to = "sampleID", 
               values_to = "expression") %>%
  group_by(sampleID) %>%
  summarise(gene_expression_vector = list(expression))

# Flatten gene expression vector for CSV export
rna_long_flat <- rna_long %>%
  mutate(gene_expression_vector = sapply(gene_expression_vector, paste, collapse = ","))

# Write to CSV
write.csv(rna_long_flat, "rna_expression.csv", row.names = FALSE)
