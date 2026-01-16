here::i_am("r/06_gsva_01_run.R")

# Attach package
library("Biobase")
library("GSVA")
library("limma")
library("magrittr")
library("tidyverse")

# Load data
dge_eset <- readRDS(
  file = here::here("output/gsva/dge_eset.rds")
)

# Define design and contrasts
## Design
gsva_design <- model.matrix(
  ~ treatment_response +
    day +
    treatment_response:day +
    biological_replicate +
    var_surrogate,
  data = dge_eset %>%
    phenoData() %>%
    pData()
) %>%
  set_colnames(
    colnames(.) %>%
      make.names
  )

## Contrasts
gsva_contrasts <- makeContrasts(
  # Treatment at each day
  treat_d4 = treatment_responseResponder,
  treat_d9 = treatment_responseResponder + treatment_responseResponder.dayD9,
  treat_d17 = treatment_responseResponder + treatment_responseResponder.dayD17,
  treat_d24 = treatment_responseResponder + treatment_responseResponder.dayD24,

  # Between days in non-responders
  d9_vs_d4_non_responders = dayD9 - 0,
  d17_vs_d4_non_responders = dayD17 - 0,
  d24_vs_d4_non_responders = dayD24 - 0,
  d17_vs_d9_non_responders = dayD17 - dayD9,
  d24_vs_d9_non_responders = dayD24 - dayD9,
  d24_vs_d17_non_responders = dayD24 - dayD17,

  # Between days in responders
  d9_vs_d4_responders = dayD9 + treatment_responseResponder.dayD9,
  d17_vs_d4_responders = dayD17 + treatment_responseResponder.dayD17,
  d24_vs_d4_responders = dayD24 + treatment_responseResponder.dayD24,
  d17_vs_d9_responders = (dayD17 + treatment_responseResponder.dayD17) -
    (dayD9 + treatment_responseResponder.dayD9),
  d24_vs_d9_responders = (dayD24 + treatment_responseResponder.dayD24) -
    (dayD9 + treatment_responseResponder.dayD9),
  d24_vs_d17_responders = (dayD24 + treatment_responseResponder.dayD24) -
    (dayD17 + treatment_responseResponder.dayD17),

  # Interaction term
  d9_vs_d4_treat = treatment_responseResponder.dayD9,
  d17_vs_d4_treat = treatment_responseResponder.dayD17,
  d24_vs_d4_treat = treatment_responseResponder.dayD24,
  d17_vs_d9_treat = treatment_responseResponder.dayD17 -
    treatment_responseResponder.dayD9,
  d24_vs_d9_treat = treatment_responseResponder.dayD24 -
    treatment_responseResponder.dayD9,
  d24_vs_d17_treat = treatment_responseResponder.dayD24 -
    treatment_responseResponder.dayD17,

  levels = gsva_design
)

# Perform GSVA
## Calculate GSVA enrichment scores
gsva_quant <- gsva_gmt %>%
  map(
    \(collection) {
      dge_eset %>%
        gsvaParam(
          geneSets = collection,
          minSize = 5,
          maxSize = 500,
          kcdf = "auto"
        ) %>%
        gsva()
    },
    .progress = TRUE
  )

# Perform limma model fit
## Use limma-trend method with robust = TRUE
gsva_fit_contrasts <- gsva_quant %>%
  map(\(quant) {
    quant %>%
      lmFit(gsva_design) %>%
      eBayes(robust = TRUE, trend = geneSetSizes(quant)) %>%
      contrasts.fit(gsva_contrasts) %>%
      eBayes(robust = TRUE, trend = geneSetSizes(quant))
  })

# Save data
## GSVA ES scores
saveRDS(
  gsva_quant,
  here::here("output/gsva/gsva_quant.rds")
)

## GSVA limma model fit
saveRDS(
  gsva_fit_contrasts,
  here::here("output/gsva/gsva_fit_contrasts.rds")
)
