# Declare location
here::i_am("r/01_load_expression.R")

# Attach packages
library("DESeq2")
library("edgeR")
library("magrittr")
library("Mus.musculus")
library("sva")
library("tidyverse")

# Load expression data
dge_counts <- read_tsv(here::here(
  "input/txome/GSE299686_RNAseq_Raw_Count.txt"
)) %>%
  column_to_rownames(var = "gene_id") %>%
  dplyr::mutate(gene_id = NULL) %>%
  # Modify row names to be in expected MSigDB format
  set_rownames(
    rownames(.) %>%
      str_sub(start = 1, end = 18)
  ) %>%
  # Turn into matrix, change into integers
  as.matrix() %>%
  round()

# Load and transform metadata to be sample table format
dge_coldata <- read_tsv(
  here::here(
    "input",
    "txome",
    "GSE299686_series_matrix.txt"
  ),
  # Skip header rows
  skip = 40
) %>%
  # Transpose to long table format
  t() %>%
  as.data.frame() %>%
  # Adjust column names
  set_colnames(.[1, ]) %>%
  repair_names() %>%
  dplyr::rename_with(\(colname) {
    colname %>%
      str_replace("!", "") %>%
      tolower()
  }) %>%
  # Remove redundant first row
  extract(2:nrow(.), ) %>%
  # Select columns
  dplyr::select(
    sample_characteristics_ch11,
    sample_geo_accession,
    sample_description,
    sample_description1,
  ) %>%
  # Rename treatment response status
  dplyr::rename(treatment_response = sample_description1) %>%
  # Adjust treatment and sample description columns
  dplyr::mutate(
    treatment = sample_characteristics_ch11 %>%
      str_replace("treatment: ", ""),
    sample_characteristics_ch11 = NULL,
    sample_id = sample_description %>%
      str_replace("Library name: ", ""),
    sample_description = NULL
  ) %>%
  # Split sample id column to yield mice ID and timepoint variables
  separate_wider_delim(
    cols = sample_id,
    delim = "_",
    names = c("sample_number", "day", "biological_replicate"),
    cols_remove = FALSE
  ) %>%
  as.data.frame() %>%
  set_rownames(.$sample_id) %>%
  # Filter and arrange columns
  dplyr::select(
    sample_id,
    sample_geo_accession,
    biological_replicate,
    day,
    treatment,
    treatment_response
  ) %>%
  # Arrange rows following counts table
  dplyr::arrange(match(.$sample_id, colnames(dge_counts))) %>%
  # Turn variables into factors
  dplyr::mutate(
    treatment_response = treatment_response %>%
      factor(levels = c("NonResponder", "Responder")),
    biological_replicate = biological_replicate %>%
      factor() %>%
      fct_inseq(),
    day = day %>%
      factor(levels = c("D4", "D9", "D17", "D24"))
  ) %>%
  dplyr::mutate(
    treatment_day = str_c(day, treatment_response, sep = "_") %>%
      factor(
        levels = c(
          "D4_NonResponder",
          "D4_Responder",
          "D9_NonResponder",
          "D9_Responder",
          "D17_NonResponder",
          "D17_Responder",
          "D24_NonResponder",
          "D24_Responder"
        )
      )
  )

# Get rowData
dge_rowdata <- AnnotationDbi::select(
  Mus.musculus,
  rownames(dge_counts),
  columns = c("SYMBOL", "GENENAME", "GENETYPE", "ENTREZID"),
  keytype = "ENSEMBL"
) %>%
  magrittr::extract(!duplicated(.$ENSEMBL), )

# Create DESeq dataset, filter for low expression genes
dge_quant <- DESeqDataSetFromMatrix(
  countData = dge_counts,
  colData = dge_coldata,
  rowData = dge_rowdata,
  design = ~ treatment_day + biological_replicate
) %>%
  magrittr::extract(filterByExpr(., group = .$treatment_day), )

# Run sva to calculate technical covariates
dge_sva_design_full <- model.matrix(
  ~ treatment_day + biological_replicate,
  data = colData(dge_quant)
)

dge_sva_design_null <- model.matrix(
  ~biological_replicate,
  data = colData(dge_quant)
)

dge_sva_n_surrogate_vars <- num.sv(
  counts(dge_quant),
  dge_sva_design_full,
  method = "be"
)

dge_sva <- sva(
  dat = counts(dge_quant),
  mod = dge_sva_design_full,
  mod0 = dge_sva_design_null,
  n.sv = dge_sva_n_surrogate_vars
)

# Update DESeq object to include surrogate variable
colData(dge_quant) <- cbind(
  colData(dge_quant),
  var_surrogate = dge_sva$sv[, 1]
)

design(dge_quant) <- ~ treatment_day + biological_replicate

# Create DESeqDataset object, run DESeq2
dge_quant_lrt <- dge_quant %>%
  DESeq(
    test = "LRT",
    reduced = ~biological_replicate
  )

dge_quant <- dge_quant %>%
  DESeq()

# Save DESeqDataset object
saveRDS(
  dge_quant,
  file = here::here("output/dge/dge_quant.rds")
)

saveRDS(
  dge_quant_lrt,
  file = here::here("output/dge/dge_quant_lrt.rds")
)
