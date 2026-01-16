here::i_am("r/05_convert_dds_to_eset.R")

# Attach packages
library("Biobase")
library("DESeq2")

# Load data
dge_quant <- readRDS(
  file = here::here("output/dge/dge_quant.rds")
)

# Change identifiers from ENSEMBL to ENTREZID
# Sort dataset by decreasing expression
# Remove genes with duplicate IDs and lower expression

## Sort
order <- order(
  mcols(dge_quant)$baseMean,
  decreasing = TRUE
)

dge_quant <- dge_quant[order, ]

## Remove genes with no ENTREZ ID
has_entrezid <- dge_quant %>%
  mcols() %>%
  as.data.frame() %>%
  .$ENTREZID %>%
  is.na(.) %>%
  !.

dge_quant <- dge_quant[has_entrezid, ]

## Remove genes with dupllicated ENTREZ ID
entrezid <- dge_quant %>%
  mcols() %>%
  as.data.frame() %>%
  .$ENTREZID

dge_quant <- dge_quant[!duplicated(entrezid), ]
entrezid <- entrezid[!duplicated(entrezid)]

## Rename rownames to entrezid
rownames(dge_quant) <- entrezid

# Convert data to ExpressionSet for GSVA
dge_eset <- dge_quant %$%
  ExpressionSet(
    assayData = counts(.),
    phenoData = colData(.) %>%
      as.data.frame() %>%
      AnnotatedDataFrame(),
    featureData = mcols(.) %>%
      as.data.frame() %>%
      AnnotatedDataFrame()
  )

# Save data
saveRDS(
  dge_eset,
  file = here::here("output/gsva/dge_eset.rds")
)
