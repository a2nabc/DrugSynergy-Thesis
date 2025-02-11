library(PharmacoGx)
library(dplyr)

###################### Obtain synergy for all pairs of drugs ########################

data <- readRDS("NCI_ALMANAC_2017.rds")
data <- updateObject(data)

x <- treatmentResponse(data)
synergy_data <- x@assays$profiles


assayindex <- x@.intern$assayIndex
assayindex <- assayindex[, c("rowKey", "colKey", ".profiles")]

merged_synergy_rowData <- merge(synergy_data, assayindex, by=".profiles")

cell_line_data <- x@colData[, c("sampleid", "colKey")]
merged_synergy_cell_lines <- merge(merged_synergy_rowData, cell_line_data, by="colKey")

drug_pairs <- x@rowData %>%
  select(treatment1id, treatment2id, rowKey) %>%
  na.omit()

# Merge with synergy scores
merged_synergy_data <- merge(merged_synergy_cell_lines, drug_pairs, by = "rowKey")

# Count the number of observations for each drug pair
drug_pair_counts <- merged_synergy_data %>%
  group_by(treatment1id, treatment2id) %>%
  summarise(count = n(), .groups = "drop") %>%
  arrange(desc(count))



#######################3 GENE EXPRESSION ##################################
y<-molecularProfiles(data, "rnaseq.comp")

y <- t(y)
y <- as.data.frame(y)
y$sampleid <- rownames(y)
merged_data_with_y <- merge(merged_synergy_data, y, by = "sampleid")
