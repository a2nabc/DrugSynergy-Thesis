#install.packages(c("httr", "readr", "dplyr", "igraph", "tidyr", "ggplot2"))
library(httr)
library(readr)
library(dplyr)
library(igraph)
library(tidyr)
library(ggplot2)

query_string_api <- function(target, species = 9606, score_threshold = 750) { #species id for humans
  base_url <- "https://string-db.org/api/tsv/network"
  query_url <- paste0(
    base_url, 
    "?identifiers=", URLencode(target),  # Note: "identifiers" instead of "identifier"
    "&species=", species, 
    "&required_score=", score_threshold
  )
  
  # send request to the API
  response <- GET(query_url)
  
  # Check if the request was successful
  if (status_code(response) == 200) {
    content <- content(response, "text", encoding = "UTF-8")
    
    # Check if content is empty
    if (nchar(content) > 0) {
      data <- read_tsv(content, col_types = cols())
      return(data %>% mutate(query_target = target))
    } else {
      cat(sprintf("No interactions found for %s.\n", target))
      return(NULL)
    }
  } else {
    cat(sprintf("Error querying STRING for %s: HTTP %s\n", target, status_code(response)))
    return(NULL)
  }
}

# Load the CSV file
data_matrix <- read.csv("../data/single_drug/Bortezomib.csv", header = TRUE, check.names = FALSE)

# Extract gene names (column names)
gene_list <- colnames(data_matrix)[-(1:4)]
writeLines(gene_list, "./python/gene_list.txt")

# Query STRING for all genes in the dataset
all_interactions <- lapply(unique(gene_list), function(gene) {
  cat(sprintf("Querying STRING for %s...\n", gene))
  query_string_api(gene)
})

# Combine non-null results into a single dataframe
interactions_df <- do.call(rbind, Filter(Negate(is.null), all_interactions)) %>%
  filter(score > 0.75)  # Filter for high-confidence interactions

# Select relevant columns for network (source and target genes)
edges <- interactions_df %>% select(preferredName_A, preferredName_B)

# Create a graph object using igraph
g_protein_coding <- graph_from_data_frame(edges, directed = FALSE)

saveRDS(g_protein_coding, file = "g_protein_coding.rds")
# Plot the network

node_degree <- degree(g_protein_coding)
threshold <- quantile(node_degree, 0.90)

V(g_protein_coding)$label <- ifelse(node_degree >= threshold, V(g_protein_coding)$name, NA)

plot(g_protein_coding, vertex.size = 8, vertex.label = V(g_protein_coding)$label, vertex.label.cex = 0.8, edge.color = "gray",
     main = "Gene-Gene Interaction Network")



################################## GENE SUBSETS ################################

gene_subsets <- c("COSMIC", "kegg", "LINCS", "mdsig_hallmarks")


# Loop through each subset to construct graphs
for (subset_name in gene_subsets) {
  
  # Construct file path
  file_path <- paste0("../data/single_drug/gene_subsets/", subset_name, ".txt")
  
  # Load subset of genes
  if (file.exists(file_path)) {
    subset_genes <- readLines(file_path)
    
    # Filter interactions for only subset genes
    subset_interactions <- interactions_df %>%
      filter(preferredName_A %in% subset_genes & preferredName_B %in% subset_genes)
    
    # Create a graph object
    assign(paste0("g_", subset_name), graph_from_data_frame(subset_interactions, directed = FALSE))
    #saveRDS(paste0("g_", subset_name), file = paste0("g_", subset_name, ".rds"))
  } else {
    cat(sprintf("File %s not found. Skipping...\n", file_path))
  }
}

# Save an RDS file for each graph
for (subset_name in gene_subsets){
  saveRDS(paste0("g_", subset_name), file = paste0("g_", subset_name, ".rds"))
  #write_graph(paste0("g_", subset_name), file = paste0("g_", subset_name, ".graphml"), format = "graphml") # --------> TO USE OUTSIDE OF R
}

# Plot every graph
for (subset_name in gene_subsets) {
  file_path <- paste0("../data/single_drug/gene_subsets/", subset_name, ".txt")
  
  if (file.exists(file_path)) {
    subset_genes <- readLines(file_path)
    # Plot the subset graph
    graph <- get(paste0("g_", subset_name))
    node_degree <- degree(graph)
    threshold <- quantile(node_degree, 0.90)
    
    V(graph)$label <- ifelse(node_degree >= threshold, V(graph)$name, NA)
    
    plot(graph, vertex.size = 8, vertex.label = V(graph)$label, vertex.label.cex = 0.8, edge.color = "gray",
         main = paste(subset_name, "Gene-Gene Interaction Network"))
  }    
}



library(biomaRt)

# Connect to Ensembl BioMart
mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

# Extract Ensembl Protein IDs (removing "9606.")
ensembl_protein_ids <- sub("9606.", "", V(g_kegg)$name)

# Retrieve gene symbols
mapping <- getBM(
  attributes = c("ensembl_peptide_id", "external_gene_name"),
  filters = "ensembl_peptide_id",
  values = ensembl_protein_ids,
  mart = mart
)

# Replace vertex names with gene symbols
V(g_kegg)$name <- mapping$external_gene_name[match(sub("9606.", "", V(g_kegg)$name), mapping$ensembl_peptide_id)]

# Check updated names
V(g_kegg)$name

