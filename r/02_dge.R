here::i_am("r/02_dge.R")

# Attach necessary packages
library("ashr")
library("DESeq2")
library("tidyverse")

# Load data
dge_quant <- readRDS(
  file = here::here("output/dge/dge_quant.rds")
)

dge_quant_lrt <- readRDS(
  file = here::here("output/dge/dge_quant_lrt.rds")
)

# Obtain results (Wald test)
## No LFC shrinking
## Filter for protein-coding genes
dge_results <- dge_deseq2_contrasts %>%
  map(
    \(contrast) {
      # Get DGE results
      results <- results(
        dge_quant,
        contrast = contrast,
        alpha = 0.05
      )

      # Filter for protein-coding genes
      results <- results[
        !is.na(mcols(dge_quant)$GENETYPE) &
          mcols(dge_quant)$GENETYPE == "protein-coding",
      ]

      # Readjust p values
      results$padj <- results$pvalue %>%
        p.adjust(method = "BH")

      # Return data
      return(results)
    },
    .progress = TRUE
  )

## With LFC shrinking
dge_results_lfcshrink <- dge_deseq2_contrasts %>%
  map2(
    .x = .,
    .y = dge_results,
    \(contrast, results) {
      lfcShrink(
        dge_quant,
        contrast = contrast,
        res = results,
        type = "ashr"
      )
    },
    .progress = TRUE
  )

# Obtain results (LR test)
## Get p values, subset to protein-coding genes
dge_results_lrt <- dge_quant_lrt %>%
  results(alpha = 0.05)

dge_results_lrt <- dge_results_lrt[
  !is.na(mcols(dge_quant_lrt)$GENETYPE) &
    mcols(dge_quant_lrt)$GENETYPE == "protein-coding",
]

## Readjust p values
dge_results_lrt$padj <- dge_results_lrt$pvalue %>%
  p.adjust(method = "BH")

# Save data
## rlog expression

## Results
### LFC
saveRDS(dge_results, file = here::here("output/dge/dge_results.rds"))

### Shrunken LFC
saveRDS(
  dge_results,
  file = here::here("output/dge/dge_results_lfcshrink.rds")
)

### LRT
saveRDS(
  dge_results,
  file = here::here("output/dge/dge_results_lrt.rds")
)
