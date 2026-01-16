here::i_am("r/06_gsva_04_volcano3d.R")

# Attach packages
library("GSVA")
library("limma")
library("magrittr")
library("tidyverse")
library("volcano3D")

# Define labels
volcano3d_labels <- list(
  mh = c(
    "HALLMARK_INFLAMMATORY_RESPONSE",
    "HALLMARK_KRAS_SIGNALING_UP",
    "HALLMARK_COMPLEMENT",
    "HALLMARK_TNFA_SIGNALING_VIA_NFKB",
    "HALLMARK_INTERFERON_ALPHA_RESPONSE"
  ),
  m2_cgp = NULL,
  m2_cp = NULL,
  m5_go = NULL
)

# Load data
## GSVA ES scores
gsva_quant <- readRDS(
  here::here("output/gsva/gsva_quant.rds")
)

## GSVA limma model fit
gsva_fit_contrasts <- readRDS(
  here::here("output/gsva/gsva_fit_contrasts.rds")
)

# Define outcome variable
outcome <- phenoData(gsva_quant$mh)$day %>%
  factor(levels = c("D4", "D9", "D17"))

# Define axis specs
breaks <- seq(
  from = 0,
  to = 3.5,
  by = 0.5
)

# Create matrix with appropriate samples
## Analyse on mh only (for now)

## Create GSVA quant object with appropriate samples
volcano3d_gsva_quant <- gsva_quant$mh

### Extract D4, D9, D17 samples
volcano3d_gsva_quant <- volcano3d_gsva_quant %>%
  .[, .$day %in% c("D4", "D9", "D17")]

### Adjust metadata
volcano3d_gsva_quant@phenoData$day <- volcano3d_gsva_quant@phenoData$day %>%
  factor(levels = c("D4", "D9", "D17"))

# Get list of appropriate p values
volcano3d_pvals <- map(
  c("D4" = 1, "D9" = 2, "D17" = 3),
  \(coef) {
    pvals <- topTable(
      gsva_fit_contrasts$mh,
      coef = coef,
      number = Inf,
      sort.by = "none"
    ) %>%
      .$P.Value

    return(pvals)
  }
) %>%
  as.data.frame() %>%
  as.matrix()

volcano3d_padj <- map(
  c("D4" = 1, "D9" = 2, "D17" = 3),
  \(coef) {
    padj <- topTable(
      gsva_fit_contrasts$mh,
      coef = coef,
      number = Inf,
      sort.by = "none"
    ) %>%
      .$adj.P.Val

    return(padj)
  }
) %>%
  as.data.frame() %>%
  as.matrix()

# Create volcano3d object with appropriate samples
volcano3d_polar <- polar_coords_2x3(
  data = exprs(volcano3d_gsva_quant) %>%
    t(),
  metadata = volcano3d_gsva_quant@phenoData %>%
    as("data.frame"),
  outcome = "treatment_response",
  group = "day",
  pvals = volcano3d_pvals,
  padj = volcano3d_padj
)
