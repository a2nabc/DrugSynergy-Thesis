library(PharmacoGx)
library(dplyr)

###################### Obtain synergy for all pairs of drugs ########################

data <- readRDS("NCI_ALMANAC_2017.rds")
data <- updateObject(data)

x <- treatmentResponse(data)
synergy_data <- x@assays$profiles


assayindex <- x@.intern$assayIndex
assayindex <- assayindex[, c("rowKey", ".profiles")]

merged_synergy_rowData <- merge(synergy_data, assayindex, by=".profiles")

drug_pairs <- x@rowData %>%
  select(treatment1id, treatment2id, rowKey) %>%
  na.omit()

# Merge with synergy scores
merged_synergy_data <- merge(merged_synergy_rowData, drug_pairs, by = "rowKey")

# Count the number of observations for each drug pair
drug_pair_counts <- merged_synergy_data %>%
  group_by(treatment1id, treatment2id) %>%
  summarise(count = n(), .groups = "drop") %>%
  arrange(desc(count))


