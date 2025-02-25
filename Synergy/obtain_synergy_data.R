library(PharmacoGx)
library(dplyr)

###################### Obtain Synergy for All Pairs of Drugs ########################

data <- readRDS("NCI_ALMANAC_2017.rds")
data <- updateObject(data)

x <- treatmentResponse(data)
synergy <- x@assays$profiles

assayindex <- x@.intern$assayIndex
assayindex <- assayindex[, c("rowKey", "colKey", ".profiles")]

synergy_data <- merge(synergy, assayindex, by=".profiles")

cell_line_data <- x@colData[, c("sampleid", "colKey")]
merged_synergy_cell_lines <- merge(synergy_data, cell_line_data, by="colKey")

drug_pairs <- x@rowData %>%
  select(treatment1id, treatment2id, rowKey) %>%
  na.omit()

# Merge with synergy scores
synergy_data <- merge(merged_synergy_cell_lines, drug_pairs, by = "rowKey") %>%
  filter(!is.na(SCORE))

synergy_data <- synergy_data %>%
  select(sampleid, treatment1id, treatment2id, SCORE)

######################## Binarizing Synergy Scores Per Cell Line ########################
synergy_data <- synergy_data %>%
  group_by(sampleid) %>%
  mutate(threshold = quantile(SCORE, 0.9, na.rm = TRUE),  # Top 10% per cell line
         label = ifelse(SCORE >= threshold, 1, 0)) %>%  
  select(-threshold) %>%
  ungroup()

######################### GENE EXPRESSION ##################################
y <- molecularProfiles(data, "rnaseq.comp")

y <- t(y)
y <- as.data.frame(y)
y$sampleid <- rownames(y)

########################### FILTER ONLY 5 GIVEN DRUG PAIRS ###########################

selected_pairs <- list(
  c("Bortezomib", "5-Fluorouracil"),
  c("Cisplatin", "Lapatinib"),
  c("Sirolimus", "Axitinib"),
  c("Vorinostat", "Ruxolitinib"),
  c("Vorinostat", "Vinorelbine")
)

for (i in seq_along(selected_pairs)) {
  pair <- selected_pairs[[i]]
  
  # Filter for the specific drug pair
  selected_drug_pairs <- synergy_data %>%
    filter(treatment1id == pair[1] & treatment2id == pair[2]) %>%
    group_by(sampleid, treatment1id, treatment2id) %>%
    summarise(median_SCORE = median(SCORE, na.rm = TRUE),
              median_label = median(label, na.rm = TRUE),  # Keep binarized value
              .groups = "drop")
  
  # Merge with gene expression data
  merged_data <- merge(selected_drug_pairs, y, by = "sampleid")
  
  # Save dataset
  write.csv(merged_data, paste0("data/synergy/synergy_bin", i, ".csv"), row.names=FALSE)
}
