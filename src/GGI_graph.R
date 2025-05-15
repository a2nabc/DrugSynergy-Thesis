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


# Construct file path
KEGG <- readLines("./data/gene_subsets/KEGG.txt")
LINCS <- readLines("./data/gene_subsets/LINCS.txt")
subset_genes <- union(KEGG, LINCS)


# Query STRING for all genes in the dataset
all_interactions <- lapply(unique(subset_genes), function(gene) {
  cat(sprintf("Querying STRING for %s...\n", gene))
  query_string_api(gene)
})

# Combine non-null results into a single dataframe
interactions_df <- do.call(rbind, Filter(Negate(is.null), all_interactions)) %>%
  filter(score > 0.75)  # Filter for high-confidence interactions


# Select relevant columns for network (source and target genes)
edges <- interactions_df %>% select(preferredName_A, preferredName_B)

# Create a graph object using igraph
g <- graph_from_data_frame(edges, directed = FALSE)

saveRDS(g, file = "data/g_KEGG_union_LINCS.rds")
write_graph(g, file = "data/g_KEGG_UNION_LINCS.graphml", format = "graphml")

  
       