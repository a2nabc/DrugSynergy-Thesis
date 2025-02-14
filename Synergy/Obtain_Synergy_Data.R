library(PharmacoGx)
library(dplyr)

###################### Obtain synergy for all pairs of drugs ########################

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

synergy_data$rowKey <- NULL
synergy_data$colKey <- NULL
synergy_data$.profiles <- NULL

# Count the number of observations for each drug pair
drug_pair_counts_sampleid <- synergy_data %>%
  group_by(treatment1id, treatment2id, sampleid) %>%
  summarise(count = n(), .groups = "drop") %>%
  arrange(desc(count))

drug_pair_counts <- synergy_data %>%
  group_by(treatment1id, treatment2id) %>%
  summarise(count = n(), .groups = "drop") %>%
  arrange(desc(count))

######################### GENE EXPRESSION ##################################
y<-molecularProfiles(data, "rnaseq.comp")

y <- t(y)
y <- as.data.frame(y)
y$sampleid <- rownames(y)




##################### RUNNING THE MERGE LOCALLY CRASHES...################
#merged_data_with_y <- merge(merged_synergy_data, y, by = "sampleid")

#join is more efficient than merge:
#merged_data_with_y <- synergy_data %>%
#  inner_join(y, by = "sampleid")


################ FILTER ONLY 5 GIVEN DRUG PAIRS ###########################

selected_drug_pairs <- synergy_data %>%
  filter((treatment1id == "Bortezomib" & treatment2id == "5-Fluorouracil")) %>%
  group_by(treatment1id, treatment2id, sampleid) %>%
  summarise(median_SCORE = median(SCORE, na.rm = TRUE))

merged_data_1 <- merge(selected_drug_pairs, y, by = "sampleid")

selected_drug_pairs <- synergy_data %>%
  filter((treatment1id == "Bortezomib" & treatment2id == "Gemcitabine")) %>%
  group_by(treatment1id, treatment2id, sampleid) %>%
  summarise(median_SCORE = median(SCORE, na.rm = TRUE))

merged_data_2 <- merge(selected_drug_pairs, y, by = "sampleid")

selected_drug_pairs <- synergy_data %>%
  filter((treatment1id == "Bortezomib" & treatment2id == "Axitinib")) %>%
  group_by(treatment1id, treatment2id, sampleid) %>%
  summarise(median_SCORE = median(SCORE, na.rm = TRUE))

merged_data_3 <- merge(selected_drug_pairs, y, by = "sampleid")

selected_drug_pairs <- synergy_data %>%
  filter((treatment1id == "Bortezomib" & treatment2id == "Erlotinib")) %>%
  group_by(treatment1id, treatment2id, sampleid) %>%
  summarise(median_SCORE = median(SCORE, na.rm = TRUE))

merged_data_4 <- merge(selected_drug_pairs, y, by = "sampleid")

selected_drug_pairs <- synergy_data %>%
  filter((treatment1id == "Cisplatin" & treatment2id == "Ruxolitinib")) %>%
  group_by(treatment1id, treatment2id, sampleid) %>%
  summarise(median_SCORE = median(SCORE, na.rm = TRUE))

merged_data_5 <- merge(selected_drug_pairs, y, by = "sampleid")

selected_drug_pairs <- synergy_data %>%
  filter((treatment1id == "Cisplatin" & treatment2id == "Nelarabine")) %>%
  group_by(treatment1id, treatment2id, sampleid) %>%
  summarise(median_SCORE = median(SCORE, na.rm = TRUE))

merged_data_6 <- merge(selected_drug_pairs, y, by = "sampleid")

selected_drug_pairs <- synergy_data %>%
  filter((treatment1id == "Cisplatin" & treatment2id == "Lapatinib")) %>%
  group_by(treatment1id, treatment2id, sampleid) %>%
  summarise(median_SCORE = median(SCORE, na.rm = TRUE))

merged_data_7 <- merge(selected_drug_pairs, y, by = "sampleid")

selected_drug_pairs <- synergy_data %>%
  filter((treatment1id == "Erlotinib" & treatment2id == "Crizotinib")) %>%
  group_by(treatment1id, treatment2id, sampleid) %>%
  summarise(median_SCORE = median(SCORE, na.rm = TRUE))

merged_data_8 <- merge(selected_drug_pairs, y, by = "sampleid")

selected_drug_pairs <- synergy_data %>%
  filter((treatment1id == "Erlotinib" & treatment2id == "Docetaxel")) %>%
  group_by(treatment1id, treatment2id, sampleid) %>%
  summarise(median_SCORE = median(SCORE, na.rm = TRUE))

merged_data_9 <- merge(selected_drug_pairs, y, by = "sampleid")

selected_drug_pairs <- synergy_data %>%
  filter((treatment1id == "Sirolimus" & treatment2id == "Gefitinib")) %>%
  group_by(treatment1id, treatment2id, sampleid) %>%
  summarise(median_SCORE = median(SCORE, na.rm = TRUE))

merged_data_10 <- merge(selected_drug_pairs, y, by = "sampleid")

selected_drug_pairs <- synergy_data %>%
  filter((treatment1id == "Sirolimus" & treatment2id == "Axitinib")) %>%
  group_by(treatment1id, treatment2id, sampleid) %>%
  summarise(median_SCORE = median(SCORE, na.rm = TRUE))

merged_data_11 <- merge(selected_drug_pairs, y, by = "sampleid")

selected_drug_pairs <- synergy_data %>%
  filter((treatment1id == "Sirolimus" & treatment2id == "Mitoxantrone")) %>%
  group_by(treatment1id, treatment2id, sampleid) %>%
  summarise(median_SCORE = median(SCORE, na.rm = TRUE))

merged_data_12 <- merge(selected_drug_pairs, y, by = "sampleid")

selected_drug_pairs <- synergy_data %>%
  filter((treatment1id == "Vorinostat" & treatment2id == "Ruxolitinib")) %>%
  group_by(treatment1id, treatment2id, sampleid) %>%
  summarise(median_SCORE = median(SCORE, na.rm = TRUE))

merged_data_13 <- merge(selected_drug_pairs, y, by = "sampleid")

selected_drug_pairs <- synergy_data %>%
  filter((treatment1id == "Vorinostat" & treatment2id == "Vinorelbine")) %>%
  group_by(treatment1id, treatment2id, sampleid) %>%
  summarise(median_SCORE = median(SCORE, na.rm = TRUE))

merged_data_14 <- merge(selected_drug_pairs, y, by = "sampleid")

selected_drug_pairs <- synergy_data %>%
  filter((treatment1id == "Vorinostat" & treatment2id == "Gefitinib")) %>%
  group_by(treatment1id, treatment2id, sampleid) %>%
  summarise(median_SCORE = median(SCORE, na.rm = TRUE))

merged_data_15 <- merge(selected_drug_pairs, y, by = "sampleid")

write.csv(merged_data_4, "Single_drug_Lasso/synergy1.csv", row.names=FALSE)
write.csv(merged_data_7, "Single_drug_Lasso/synergy2.csv", row.names=FALSE)
write.csv(merged_data_11, "Single_drug_Lasso/synergy3.csv", row.names=FALSE)
write.csv(merged_data_13, "Single_drug_Lasso/synergy4.csv", row.names=FALSE)
write.csv(merged_data_14, "Single_drug_Lasso/synergy5.csv", row.names=FALSE)

