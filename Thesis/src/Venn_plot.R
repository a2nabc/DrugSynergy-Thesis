# library
#install.packages("VennDiagram")
library(VennDiagram)

# Read sample IDs from processed data
broad_samples <- read.csv("../data/processed/Broad/ccl_data.csv")$sampleid
gdsc_samples <- read.csv("../data/processed/GDSC2/ccl_data.csv")$sampleid
gCSI_samples <- read.csv("../data/processed/gCSI/ccl_data.csv")$sampleid

# Create a list of sample IDs for the Venn diagram
sample_sets <- list(
  "Broad" = broad_samples,
  "GDSC2" = gdsc_samples,
  "gCSI" = gCSI_samples
)

# Generate and save the Venn diagram
venn.plot <- venn.diagram(
  x = sample_sets,
  category.names = c("Broad", "GDSC2", "gCSI"),
  filename = "../data/sample_venn_diagram.png",
  output = TRUE,
  fill = c("red", "blue", "green"),
  alpha = 0.5,
  cat.col = c("red", "blue", "green"),
  cat.cex = 1.2,
  margin = 0.1
)

# Display the Venn Diagram
grid::grid.draw(venn.plot)
