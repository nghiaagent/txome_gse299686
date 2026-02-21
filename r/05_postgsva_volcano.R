here::i_am("r/05_postgsva_volcano.R")

######################
# Build volcano plots of GSVA results for selected gene sets
######################

# Attach packages
library("Biobase")
library("DESeq2")
library("ggrepel")
library("magrittr")
library("Mus.musculus")
library("tidyverse")

# Load data
dge_quant <- readRDS(
  file = here::here("output/dge/dge_quant.rds")
)

dge_results_lfcshrink <- readRDS(
  file = here::here("output/dge/dge_results_lfcshrink.rds")
)

# Load list of genes for the gene sets of interest
## Hallmark MYC targets V2
## Hallmark interferon alpha

### Define gene sets
hallmark_sel <- c(
  "HALLMARK_MYC_TARGETS_V1",
  "HALLMARK_MYC_TARGETS_V2",
  "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION"
)

### Get list of ENSEMBL IDs within gene set
genes_sel <- hallmark_sel %>%
  set_names(., .) %>%
  map(\(hallmark) {
    gsva_gmt %>%
      use_series(mh) %>%
      geneIds() %>%
      extract2(hallmark) %>%
      mapIds(
        x = Mus.musculus,
        keys = .,
        keytype = "ENTREZID",
        column = "ENSEMBL"
      )
  })

# Build volcano plots
## Highlight genes in the selected hallmark gene sets
## Compare between days in non-responders

plots_volcano_postgsva <- map(dge_results_lfcshrink, \(results) {
  map(genes_sel, \(genes_sel) {
    ## Split into 2 groups: Not highlighted and highlighted
    toptable_nothighlight <- results %>%
      extract(!rownames(.) %in% genes_sel, )

    toptable_highlight <- results %>%
      extract(rownames(.) %in% genes_sel, )

    ## Extract values for plotting
    nothighlight_xvals <- toptable_nothighlight$log2FoldChange
    nothighlight_yvals <- -log10(toptable_nothighlight$padj)

    highlight_xvals <- toptable_highlight$log2FoldChange
    highlight_yvals <- -log10(toptable_highlight$padj)
    highlight_colourmapping <- case_when(
      toptable_highlight$padj < 0.05 &
        abs(toptable_highlight$log2FoldChange) > 1 ~ "highlight_both",
      toptable_highlight$padj < 0.05 ~ "highlight_padj",
      .default = "highlight_ns"
    )

    ## Clip values: pvals to -30 and lfc to +-5
    ### Highlight type of clip via shape
    nothighlight_shape <- toptable_nothighlight %$%
      case_when(
        .$log2FoldChange >= 5 ~ "clip_right",
        .$log2FoldChange <= -5 ~ "clip_left",
        .$padj <= 1e-30 ~ "clip_up",
        .default = "unclipped"
      ) %>%
      set_names(case_when(
        . == "clip_right" ~ str_c("logFC > 5"),
        . == "clip_left" ~ str_c("logFC < -5"),
        . == "clip_up" ~ str_c("padj < 1e-30"),
        . == "unclipped" ~ "Unclipped"
      ))

    highlight_shape <- toptable_highlight %$%
      case_when(
        .$log2FoldChange >= 5 ~ "clip_right",
        .$log2FoldChange <= -5 ~ "clip_left",
        .$padj <= 1e-30 ~ "clip_up",
        .default = "unclipped"
      ) %>%
      set_names(case_when(
        . == "clip_right" ~ str_c("logFC > 5"),
        . == "clip_left" ~ str_c("logFC < -5"),
        . == "clip_up" ~ str_c("padj < 1e-30"),
        . == "unclipped" ~ "Unclipped"
      ))

    ### Clip x and y values
    nothighlight_xvals_clipped <- nothighlight_xvals %>%
      case_when(
        . >= 5 ~ 5,
        . <= -5 ~ -5,
        .default = .
      )
    nothighlight_yvals_clipped <- nothighlight_yvals %>%
      case_when(
        . <= 1e-30 ~ 1e-30,
        .default = .
      )

    highlight_xvals_clipped <- highlight_xvals %>%
      case_when(
        . >= 5 ~ 5,
        . <= -5 ~ -5,
        .default = .
      )
    highlight_yvals_clipped <- highlight_yvals %>%
      case_when(
        . <= 1e-30 ~ 1e-30,
        .default = .
      )

    ## Define labels of highlighted genes
    highlight_labels <- rownames(toptable_highlight) %>%
      recode_values(
        from = mcols(dge_quant)$ENSEMBL,
        to = mcols(dge_quant)$SYMBOL
      )

    # Construct plot
    plot <- ggplot() +
      geom_hline(
        yintercept = -log10(0.05),
        linetype = 2
      ) +
      geom_vline(
        xintercept = c(-1, 1),
        linetype = 2
      ) +
      geom_point(
        aes(
          x = nothighlight_xvals_clipped,
          y = nothighlight_yvals_clipped,
          shape = nothighlight_shape
        ),
        colour = "grey80",
        alpha = 0.4
      ) +
      geom_point(
        aes(
          x = highlight_xvals_clipped,
          y = highlight_yvals_clipped,
          shape = highlight_shape,
          colour = highlight_colourmapping
        ),
        alpha = 0.5
      ) +
      geom_label_repel(
        aes(
          x = highlight_xvals_clipped,
          y = highlight_yvals_clipped,
          label = highlight_labels
        ),
        size = 3,
        force = 1.5,
        max.overlaps = 5,
        label.padding = 0.1
      ) +
      scale_shape_manual(
        values = c(
          "clip_right" = -9658,
          "clip_left" = -9668,
          "clip_up" = 17,
          "unclipped" = 19
        )
      ) +
      scale_colour_manual(
        values = c(
          "highlight_both" = "red",
          "highlight_padj" = "blue",
          "highlight_ns" = "grey30"
        )
      ) +
      xlim(-5, 5) +
      ylim(NA, 30) +
      labs(
        x = bquote(~ Log[2] ~ "fold change"),
        y = bquote(~ -Log[10] ~ "adjusted P value")
      ) +
      theme_classic() +
      theme(legend.position = "none")

    return(plot)
  })
}) %>%
  unlist(recursive = FALSE)

# Save data
saveRDS(
  plots_volcano_postgsva,
  file = here::here(
    "output/dge/plots_volcano_postgsva/plots_volcano_postgsva.rds"
  )
)

# Export plots
plots_volcano_postgsva %>%
  extract(c(
    "d9_vs_d4_nonresponder.HALLMARK_MYC_TARGETS_V1",
    "d9_vs_d4_nonresponder.HALLMARK_MYC_TARGETS_V2",
    "d9_vs_d4_nonresponder.HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION",
    "d17_vs_d4_nonresponder.HALLMARK_MYC_TARGETS_V1",
    "d17_vs_d4_nonresponder.HALLMARK_MYC_TARGETS_V2",
    "d17_vs_d4_nonresponder.HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION",
    "d17_vs_d9_nonresponder.HALLMARK_MYC_TARGETS_V1",
    "d17_vs_d9_nonresponder.HALLMARK_MYC_TARGETS_V2",
    "d17_vs_d9_nonresponder.HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION"
  )) %>%
  iwalk(
    \(plot, name) {
      ggsave(
        filename = here::here(
          "output/dge/plots_volcano_postgsva",
          str_c(name, ".png")
        ),
        plot = plot,
        width = 6,
        height = 8,
        dpi = 600,
        scale = 0.6
      )
    }
  )
