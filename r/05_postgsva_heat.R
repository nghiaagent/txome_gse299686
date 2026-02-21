here::i_am("r/05_postgsva_heat.R")

######################
# Build heatmap of GSVA results for selected gene sets
######################

# Attach packages
library("Biobase")
library("circlize")
library("colorspace")
library("ComplexHeatmap")
library("DESeq2")
library("GSVA")
library("limma")
library("magrittr")
library("tidyverse")
library("viridis")

# Set color scheme and breaks
## For gene expression
col <- inferno(n = 100)
breaks <- seq(-2, 2, length.out = 100)

# Load data
dge_eset <- readRDS(
  file = here::here("output/gsva/dge_eset.rds")
)

## Subset to only non-responders between D4, 9, 17
dge_eset_subset <- dge_eset
dge_eset_subset <- dge_eset_subset[,
  dge_eset_subset$treatment_response == "NonResponder"
]
dge_eset_subset <- dge_eset_subset[,
  dge_eset_subset$day %in% c("D4", "D9", "D17")
]

## Calculate log expression
dge_vst <- dge_eset_subset %>%
  exprs() %>%
  vst()

# Load lists of genes for the gene sets of interest
# Hallmark MYC targets V2
# Hallmark interferon alpha
## Define gene sets
hallmark_sel <- c(
  "HALLMARK_MYC_TARGETS_V1",
  "HALLMARK_MYC_TARGETS_V2",
  "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION"
)

## Get list of ENTREZ IDs within gene set
genes_sel <- hallmark_sel %>%
  set_names(., .) %>%
  map(\(hallmark) {
    gsva_gmt %>%
      use_series(mh) %>%
      geneIds() %>%
      extract2(hallmark)
  })

# Subset data to GOIs
## Get Z-scores
dge_vst_sel <- genes_sel %>%
  map(\(entrezid) {
    # Filter
    dge_vst_sel <- dge_vst %>%
      # Filter for genes in provided gene list
      extract(rownames(.) %in% entrezid, ) %>%

      # get per-row Z-scores
      t() %>%
      scale() %>%
      t()

    # Recode row names
    rownames(dge_vst_sel) <- rownames(dge_vst_sel) %>%
      recode_values(
        from = featureData(dge_eset)$ENTREZID,
        to = featureData(dge_eset)$SYMBOL
      )

    # Return data
    return(dge_vst_sel)
  })

# Create extra factors for splitting columns
split_day <- dge_eset_subset$day

# Build annotation; include only necessary metadata
## Dataframe of annotation data
heatmap_anno_tibble <- dge_eset_subset@phenoData %$%
  tibble(
    # `Replicate` = .$biological_replicate,
    `Day` = .$day
  )

## Build ComplexHeatmap annotation object
heatmap_anno_object <- HeatmapAnnotation(
  df = heatmap_anno_tibble,
  which = "col",
  col = palette_heatmap,
  annotation_height = 0.6,
  annotation_legend_param = list(
    `Day` = list(
      nrow = 3,
      title = "Day",
      title_position = "lefttop-rot",
      legend_direction = "vertical",
      title_gp = gpar(fontsize = 8, fontface = "bold"),
      labels_gp = gpar(fontsize = 8, fontface = "bold")
    )
  )
)

# Try heatmap
heatmap_genes_sel <- dge_vst_sel %>%
  map(\(mat) {
    # Construct heatmap
    heatmap <- Heatmap(
      mat,
      name = "Expression Z-score",
      col = colorRamp2(breaks, col),
      border = FALSE,

      # parameters for the colour-bar that represents gradient of expression
      heatmap_legend_param = list(
        color_bar = "continuous",
        legend_direction = "vertical",
        legend_width = unit(8, "cm"),
        legend_height = unit(5.0, "cm"),
        legend_position = "left",
        title_position = "lefttop-rot",
        title_gp = gpar(fontsize = 8, fontface = "bold"),
        labels_gp = gpar(fontsize = 8, fontface = "bold")
      ),

      # row (gene) parameters
      row_km = 3,
      row_km_repeats = 5,
      cluster_rows = TRUE,
      cluster_row_slices = FALSE,
      show_row_dend = FALSE,
      row_title_side = "left",
      row_title_gp = gpar(fontsize = 10, fontface = "bold"),
      row_title_rot = 90,
      show_row_names = FALSE,
      row_names_gp = gpar(fontsize = 5),
      row_names_side = "left",
      row_names_max_width = max_text_width(
        rownames(dge_vst_sel[[1]]),
        gp = gpar(fontsize = 5)
      ),

      # column (sample) parameters
      column_split = split_day,
      # column_order = heatmap_column_order,
      column_title = NULL,
      cluster_column_slices = FALSE,
      cluster_columns = TRUE,
      show_column_dend = FALSE,
      show_column_names = FALSE,

      # specify top and bottom annotations
      top_annotation = heatmap_anno_object
    )

    # Return data
    return(heatmap)
  })

# Save heatmaps
iwalk(heatmap_genes_sel, \(heat, name) {
  png(
    file = here::here(
      "output",
      "dge",
      "plots_heatmap",
      str_c("dge_in_", name, "_split_day.png")
    ),
    width = 8,
    height = 8.5,
    units = "in",
    res = 600
  )

  draw(heat)

  dev.off()
})
