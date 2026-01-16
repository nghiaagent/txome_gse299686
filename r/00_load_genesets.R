here::i_am("r/00_load_genesets.R")

# Attach necessary packages
library("GSEABase")
library("purrr")

# Load GMT files
## Populate list of files
## Construct file path
paths_gmt <- list(
  "mh" = "mh.all.v2025.1.Mm.entrez.gmt",
  "m2_cgp" = "m2.cgp.v2025.1.Mm.entrez.gmt",
  "m2_cp" = "m2.cp.v2025.1.Mm.entrez.gmt",
  "m5_go" = "m2.cp.v2025.1.Mm.entrez.gmt"
) %>%
  map(\(filename) {
    here::here(
      "input",
      "gene_sets",
      filename
    )
  })

## Load gene sets
gsva_gmt <- paths_gmt %>%
  map(\(x) getGmt(con = x))

# Clean up
rm(paths_gmt)

message("Sourced ", here::here("r/00_load_genesets.R"))
