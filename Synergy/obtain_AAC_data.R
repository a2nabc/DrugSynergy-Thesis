library(PharmacoGx)
library(dplyr)
library(tidyr)


data <- readRDS("../Single_Drug_Response/data/GDSC2.rds")
data <- updateObject(data)

drug_names <- c("Bortezomib", "Axitinib", "Cisplatin", "Erlotinib", 
                "Lapatinib", "Ruxolitinib", "Sirolimus", "Vinorelbine", "Vorinostat")

########################### OBTAIN RESPONSE DATA ##############################

sample_info <- as.data.frame(sampleInfo(data))
keep_tissues <- names(which(table(sample_info$tissueid)>5))

keep_sample_ids <- sample_info %>%
  filter(tissueid %in% keep_tissues & !tissueid %in% c("Myeloid, Lymphoid")) %>%
  select(sampleid) %>%
  pull(sampleid)

experiment_info <- as.data.frame(treatmentResponse(data)$info)
experiment_info <- experiment_info %>% filter(treatmentid %in% drug_names & sampleid %in% keep_sample_ids)
response_vars <- treatmentResponse(data)$profiles %>% select(aac_recomputed) %>% drop_na()

sensitivity <- merge(experiment_info,response_vars, by=0) %>%
  select(sampleid, treatmentid, aac_recomputed) %>%
  group_by(sampleid,treatmentid) %>%
  summarize(aac = mean(aac_recomputed))

# Binarize AAC: Top 10% performing drugs per cell line = Sensitive (1), others = Resistant (0)
sensitivity <- sensitivity %>%
  group_by(sampleid) %>%
  mutate(threshold = quantile(aac, 0.9, na.rm = TRUE),  # Compute 90th percentile per cell line
         label = ifelse(aac >= threshold, 1, 0)) %>%    # Assign labels based on threshold
  select(-threshold) %>%
  ungroup()

######################## OBTAIN GENE EXPRESSION MATRIX ########################

keep_types = c("protein_coding")
expression <- molecularProfiles(data)$Kallisto_0.46.1.rnaseq.counts
gene_info <- as.data.frame(rowData(expression)) %>% filter(gene_type %in% keep_types) %>% select(gene_name)
expression.data <- merge(gene_info, as.data.frame(assay(expression)), by=0)%>%
  select(-Row.names) %>%
  group_by(gene_name) %>%
  summarise_at(vars(matches("EGAR")),median)

keep.ccls <- as.data.frame(colData(expression)) %>%
  filter(sampleid %in% keep_sample_ids) %>%
  select(sampleid, Primary_Tissue = Factor.Value.organism.part.)

str(keep.ccls$sampleid)
expression_data <- expression.data %>% as.data.frame() %>% tibble::column_to_rownames("gene_name") %>% t()

expression_data <- keep.ccls %>%
  select(sampleid) %>%
  merge(expression_data, by=0) %>%
  select(-Row.names) %>%
  group_by(sampleid) %>%
  summarize_all(median) %>%
  as.data.frame()

############### create a complete matrix for each drug ####################

for (drug in drug_names) {
  drug_data <- sensitivity %>% 
    filter(treatmentid == drug) %>%
    arrange(sampleid) 
  # Dynamically create a variable name for each drug
  assign(paste0("drug_data_", drug, "_bin"), drug_data)
  complete_data <- merge(drug_data, expression_data, by = "sampleid")
  assign(paste0("complete_data_", drug, "_bin"), complete_data)
  write.csv(complete_data, paste0("data/single_drug/", drug, "_bin.csv"), row.names = FALSE)
}


