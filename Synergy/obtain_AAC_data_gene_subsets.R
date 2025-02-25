library(PharmacoGx)
library(dplyr)
library(tidyr)

data <- readRDS("../Single_Drug_Response/data/GDSC2.rds")
data <- updateObject(data)

drug_names <- c("Bortezomib", "Axitinib", "Cisplatin", "Erlotinib", 
                "Lapatinib", "Ruxolitinib", "Sirolimus", "Vinorelbine", "Vorinostat")

########################### OBTAIN RESPONSE DATA ##############################

sample_info <- as.data.frame(sampleInfo(data))
keep_tissues <- names(which(table(sample_info$tissueid) > 5))

keep_sample_ids <- sample_info %>%
  filter(tissueid %in% keep_tissues) %>%
  select(sampleid) %>%
  pull(sampleid)

experiment_info <- as.data.frame(treatmentResponse(data)$info) %>%
  filter(treatmentid %in% drug_names & sampleid %in% keep_sample_ids)

response_vars <- treatmentResponse(data)$profiles %>%
  select(aac_recomputed) %>%
  drop_na()

sensitivity <- merge(experiment_info, response_vars, by=0) %>%
  select(sampleid, treatmentid, aac_recomputed) %>%
  group_by(sampleid, treatmentid) %>%
  summarize(aac = mean(aac_recomputed), .groups = "drop")

# Binarize AAC: Top 10% performing drugs per cell line = Sensitive (1), others = Resistant (0)
sensitivity <- sensitivity %>%
  group_by(sampleid) %>%
  mutate(threshold = quantile(aac, 0.9, na.rm = TRUE),  # Compute 90th percentile per cell line
         label = ifelse(aac >= threshold, 1, 0)) %>%    # Assign labels based on threshold
  select(-threshold) %>%
  ungroup()

######################## OBTAIN GENE EXPRESSION MATRICES ########################

# Read gene lists
gene_sets <- list(
  COSMIC = read.table("data/single_drug/gene_subsets/COSMIC.txt", header = FALSE, stringsAsFactors = FALSE)$V1,
  KEGG = read.table("data/single_drug/gene_subsets/kegg.txt", header = FALSE, stringsAsFactors = FALSE)$V1,
  LINCS = read.table("data/single_drug/gene_subsets/LINCS.txt", header = FALSE, stringsAsFactors = FALSE)$V1,
  MDSig_Hallmarks = read.table("data/single_drug/gene_subsets/mdsig_hallmarks.txt", header = FALSE, stringsAsFactors = FALSE)$V1
)

expression <- molecularProfiles(data)$Kallisto_0.46.1.rnaseq.counts

# Store expression matrices for each gene set
expression_matrices <- list()

for (set_name in names(gene_sets)) {
  gene_info <- as.data.frame(rowData(expression)) %>%
    filter(gene_name %in% gene_sets[[set_name]]) %>%
    select(gene_name)
  
  filtered_expression <- merge(gene_info, as.data.frame(assay(expression)), by=0) %>%
    select(-Row.names) %>%
    group_by(gene_name) %>%
    summarise_at(vars(matches("EGAR")), median, .groups = "drop") %>%
    as.data.frame() %>%
    tibble::column_to_rownames("gene_name") %>%
    t()
  
  expression_matrices[[set_name]] <- filtered_expression
}

# Sample metadata
keep.ccls <- as.data.frame(colData(expression)) %>%
  filter(sampleid %in% keep_sample_ids) %>%
  select(sampleid, Primary_Tissue = Factor.Value.organism.part.)

# Process each gene set expression matrix
for (set_name in names(expression_matrices)) {
  expression_data <- keep.ccls %>%
    select(sampleid) %>%
    merge(expression_matrices[[set_name]], by=0) %>%
    select(-Row.names) %>%
    group_by(sampleid) %>%
    summarize_all(median) %>%
    as.data.frame()
  
  expression_matrices[[set_name]] <- expression_data  # Store processed matrix
}

############### CREATE COMPLETE MATRICES FOR EACH DRUG ####################

for (drug in drug_names) {
  drug_data <- sensitivity %>% 
    filter(treatmentid == drug) %>%
    arrange(sampleid) 
  
  assign(paste0("drug_data_", drug), drug_data)
  
  for (set_name in names(expression_matrices)) {
    complete_data <- merge(drug_data, expression_matrices[[set_name]], by = "sampleid")
    assign(paste0("complete_data_", drug, "_", set_name), complete_data)
    write.csv(complete_data, paste0("data/single_drug/", drug, "_", set_name, ".csv"), row.names = FALSE)
  }
}
