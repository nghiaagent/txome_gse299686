here::i_am("r/06_gsva_04_volcano3d.R")

# Attach package
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
  to = 2.5,
  by = 0.5
)

# Define cooefficients of interest
volcano3d_coefficients <- list(
  "d4_v_d9_v_d17_non_responders" = c(
    "d9_vs_d4_non_responders",
    "d17_vs_d4_non_responders",
    "d17_vs_d9_non_responders"
  ),
  "d4_v_d9_v_d17_responders" = c(
    "d9_vs_d4_responders",
    "d17_vs_d4_responders",
    "d17_vs_d9_responders"
  )
)

# Derive GSVA results
## Expression matrices
gsva_exprs <- gsva_quant %>%
  map(\(eset) exprs(eset))

## Supply pvalues and padj
gsva_polar_pvals <- map(
  gsva_fit_contrasts,
  \(fit_contrasts) {
    map(
      volcano3d_coefficients,
      \(coef) {
        cbind(
          topTable(
            fit_contrasts,
            number = Inf,
            sort.by = "none"
          )$P.Value,
          topTable(
            fit_contrasts,
            number = Inf,
            sort.by = "none",
            coef = coef[[1]]
          )$P.Value,
          topTable(
            fit_contrasts,
            number = Inf,
            sort.by = "none",
            coef = coef[[2]]
          )$P.Value,
          topTable(
            fit_contrasts,
            number = Inf,
            sort.by = "none",
            coef = coef[[3]]
          )$P.Value
        )
      }
    )
  }
)

gsva_polar_padj <- map(
  gsva_fit_contrasts,
  \(fit_contrasts) {
    map(
      volcano3d_coefficients,
      \(coef) {
        cbind(
          topTable(
            fit_contrasts,
            number = Inf,
            sort.by = "none"
          )$adj.P.Val,
          topTable(
            fit_contrasts,
            number = Inf,
            sort.by = "none",
            coef = coef[[1]]
          )$adj.P.Val,
          topTable(
            fit_contrasts,
            number = Inf,
            sort.by = "none",
            coef = coef[[2]]
          )$adj.P.Val,
          topTable(
            fit_contrasts,
            number = Inf,
            sort.by = "none",
            coef = coef[[3]]
          )$adj.P.Val
        )
      }
    )
  }
)

# Get rownames
gsva_rownames <- map(gsva_fit_contrasts, \(fit_contrasts) {
  rownames(
    topTable(
      fit_contrasts,
      number = Inf,
      sort.by = "none"
    )
  )
})

# Construct volcano3D object
volcano3d_polar <- pmap(
  list(
    gsva_exprs,
    gsva_polar_pvals,
    gsva_polar_padj,
    gsva_rownames
  ),
  \(
    gsva_exprs,
    gsva_polar_pvals,
    gsva_polar_padj,
    gsva_rownames
  ) {
    map2(
      gsva_polar_pvals,
      gsva_polar_padj,
      \(gsva_polar_pvals, gsva_polar_padj) {
        polar_manual <- polar_coords(
          outcome = outcome,
          data = t(gsva_exprs),
          pvals = gsva_polar_pvals,
          padj = gsva_polar_padj
        )

        # Add rownames
        rownames(polar_manual@pvals) <- gsva_rownames
        rownames(polar_manual@padj) <- gsva_rownames

        # Return
        return(polar_manual)
      }
    )
  }
)

# Plot
volcano3d_allplots <- map2(
  volcano3d_polar,
  volcano3d_labels,
  \(polar, labels) {
    map(polar, \(polar) {
      # Build 3D plot
      volcano3d <- volcano3D(
        polar,
        type = 1,
        label_rows = labels,
        label_size = 24,
        grid_width = 2.5,
        axis_width = 2.5,
        z_axis_title_size = 20,
        radial_axis_title_size = 20,
        r_axis_ticks = breaks
      )

      # Build radial plots
      radial_plotly <- radial_plotly(
        polar,
        type = 1,
        label_rows = labels,
        label_size = 12,
        grid_width = 2,
        axis_width = 2,
        z_axis_title_size = 20,
        radial_axis_title_size = 20,
        r_axis_ticks = breaks
      )

      radial_ggplot <- radial_ggplot(
        polar,
        type = 1,
        r_axis_ticks = breaks
      )

      # Return list
      return(list(
        "volcano3d" = volcano3d,
        "radial_plotly" = radial_plotly,
        "radial_ggplot" = radial_ggplot
      ))
    })
  },
  .progress = TRUE
)

# Save data
saveRDS(
  volcano3d_polar,
  file = here::here(
    "output",
    "gsva",
    "volcano3d_polar.RDS"
  )
)

saveRDS(
  volcano3d_allplots,
  file = here::here(
    "output",
    "gsva",
    "volcano3d_allplots_passages.RDS"
  )
)
