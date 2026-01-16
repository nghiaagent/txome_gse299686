here::i_am("r/06_gsva_02_heatmap.R")

# Attach packages
library("Biobase")
library("circlize")
library("ComplexHeatmap")
library("GSVA")
library("limma")
library("magrittr")
library("tidyverse")
library("viridis")

# Load data
## GSVA ES scores
gsva_quant <- readRDS(
  here::here("output/gsva/gsva_quant.rds")
)

## GSVA limma model fit
gsva_fit_contrasts <- readRDS(
  here::here("output/gsva/gsva_fit_contrasts.rds")
)

# Select relevant collections
gsva_quant_subset <- gsva_quant[c(
  "mh",
  "m2_cgp",
  "m2_cp",
  "m5_go"
)]

gsva_fit_contrasts_subset <- gsva_fit_contrasts[c(
  "mh",
  "m2_cgp",
  "m2_cp",
  "m5_go"
)]

# Get significant gene sets
gsva_genesets_significant <- gsva_fit_contrasts_subset %>%
  map(\(fit) {
    top <- fit %>%
      topTable(
        coef = NULL,
        number = Inf,
        sort.by = "none"
      ) %>%
      filter(adj.P.Val < 0.05) %>%
      rownames()

    # Return data
    return(top)
  })

## Subset GSVA objects to only significant gene sets
gsva_quant_heatmap <- map2(
  gsva_quant_subset,
  gsva_genesets_significant,
  \(quant, genesets) quant[genesets, ]
)

# Get Z-scores
gsva_exprs_heatmap <- gsva_quant_heatmap %>%
  map(\(x) {
    x %>%
      exprs() %>%
      t() %>%
      scale() %>%
      t()
  })

# Set color scheme and breaks
## For gene set activity expression (GSVA ES)
col <- inferno(n = 100)
breaks <- seq(-2, 2, length.out = 100)

# Create extra factors for splitting columns
split_day_response <- gsva_quant_subset[[1]]@phenoData$treatment_day

split_response <- gsva_quant_subset[[1]]@phenoData$treatment_response

split_day <- gsva_quant_subset[[1]]@phenoData$day

# Build annotation; include only necessary metadata
## Dataframe of annotation data
heatmap_anno_tibble <- gsva_quant_subset[[1]]@phenoData %$%
  tibble(
    # `Replicate` = .$biological_replicate,
    `Day` = .$day,
    `Response` = .$treatment_response %>%
      case_match(
        "Responder" ~ "Responder",
        "NonResponder" ~ "Non-Responder"
      )
  )

## Build ComplexHeatmap annotation object
heatmap_anno_object <- HeatmapAnnotation(
  df = heatmap_anno_tibble,
  which = "col",
  col = palette_heatmap,
  annotation_height = 0.6,
  annotation_legend_param = list(
    `Day` = list(
      nrow = 4,
      title = "Day",
      title_position = "topleft",
      legend_direction = "vertical",
      title_gp = gpar(fontsize = 8, fontface = "bold"),
      labels_gp = gpar(fontsize = 8, fontface = "bold")
    ),
    `Response` = list(
      nrow = 2,
      title = "Response",
      title_position = "topleft",
      legend_direction = "vertical",
      title_gp = gpar(fontsize = 8, fontface = "bold"),
      labels_gp = gpar(fontsize = 8, fontface = "bold")
    )
  )
)

# Format names of Hallmark gene sets
rownames(gsva_exprs_heatmap[[1]]) <- rownames(gsva_exprs_heatmap[[1]]) %>%
  str_remove("HALLMARK_") %>%
  str_replace_all("_", " ")

# Create a heatmap to extract column orders from
heatmap_column_order <- Heatmap(
  gsva_exprs_heatmap[[1]],
  name = "ES\nZ-\nscore",
  col = colorRamp2(breaks, col),
  border = FALSE,

  # parameters for the colour-bar that represents gradient of expression
  heatmap_legend_param = list(
    color_bar = "continuous",
    legend_direction = "vertical",
    legend_width = unit(8, "cm"),
    legend_height = unit(5.0, "cm"),
    legend_position = "left",
    title_position = "topcenter",
    title_gp = gpar(fontsize = 8, fontface = "bold"),
    labels_gp = gpar(fontsize = 8, fontface = "bold")
  ),

  # row (gene) parameters
  cluster_rows = TRUE,
  show_row_dend = FALSE,
  row_title_side = "left",
  row_title_gp = gpar(fontsize = 10, fontface = "bold"),
  row_title_rot = 90,
  show_row_names = TRUE,
  row_names_side = "left",

  # column (sample) parameters
  column_split = split_day_response,
  # column_order = order,
  column_title = NULL,
  cluster_column_slices = FALSE,
  cluster_columns = TRUE,
  show_column_dend = FALSE,
  show_column_names = FALSE,

  # specify top and bottom annotations
  top_annotation = heatmap_anno_object
) %>%
  plot() %>%
  column_order() %>%
  unlist()

# Try heatmap
## Only hallmark (with gene set name)
heatmap_hallmark <- Heatmap(
  gsva_exprs_heatmap[[1]],
  name = "ES\nZ-\nscore",
  col = colorRamp2(breaks, col),
  border = FALSE,

  # parameters for the colour-bar that represents gradient of expression
  heatmap_legend_param = list(
    color_bar = "continuous",
    legend_direction = "vertical",
    legend_width = unit(8, "cm"),
    legend_height = unit(5.0, "cm"),
    legend_position = "left",
    title_position = "topcenter",
    title_gp = gpar(fontsize = 8, fontface = "bold"),
    labels_gp = gpar(fontsize = 8, fontface = "bold")
  ),

  # row (gene) parameters
  cluster_rows = TRUE,
  show_row_dend = FALSE,
  row_title_side = "left",
  row_title_gp = gpar(fontsize = 10, fontface = "bold"),
  row_title_rot = 90,
  show_row_names = TRUE,
  row_names_gp = gpar(fontsize = 8),
  row_names_side = "left",
  row_names_max_width = max_text_width(
    rownames(gsva_exprs_heatmap[[1]]),
    gp = gpar(fontsize = 8)
  ),

  # column (sample) parameters
  column_split = split_day,
  column_order = heatmap_column_order,
  column_title = NULL,
  cluster_column_slices = FALSE,
  cluster_columns = FALSE,
  show_column_dend = FALSE,
  show_column_names = FALSE,

  # specify top and bottom annotations
  top_annotation = heatmap_anno_object
)

## Other collections (no gene set name)
heatmap_other_collections <- gsva_exprs_heatmap %>%
  map(\(exprs) {
    heatmap <- Heatmap(
      exprs,
      name = "ES\nZ-\nscore",
      col = colorRamp2(breaks, col),
      border = FALSE,

      # parameters for the colour-bar that represents gradient of expression
      heatmap_legend_param = list(
        color_bar = "continuous",
        legend_direction = "vertical",
        legend_width = unit(8, "cm"),
        legend_height = unit(5.0, "cm"),
        legend_position = "left",
        title_position = "topcenter",
        title_gp = gpar(fontsize = 8, fontface = "bold"),
        labels_gp = gpar(fontsize = 8, fontface = "bold")
      ),

      # row (gene) parameters
      cluster_rows = TRUE,
      show_row_dend = FALSE,
      row_title_side = "left",
      row_title_gp = gpar(fontsize = 10, fontface = "bold"),
      row_title_rot = 90,
      show_row_names = FALSE,
      row_names_side = "left",

      # column (sample) parameters
      column_split = split_day,
      column_order = heatmap_column_order,
      column_title = NULL,
      cluster_column_slices = FALSE,
      cluster_columns = FALSE,
      show_column_dend = FALSE,
      show_column_names = FALSE,

      # specify top and bottom annotations
      top_annotation = heatmap_anno_object
    )

    return(heatmap)
  })

# Export heatmaps
png(
  file = here::here(
    "output",
    "gsva",
    "plots_heatmap",
    str_c("gsva_mh_split_day_response_with_names.png")
  ),
  width = 12,
  height = 8,
  units = "in",
  res = 600
)

plot(heatmap_hallmark)

dev.off()

imap(
  heatmap_other_collections,
  \(heatmap, name) {
    png(
      file = here::here(
        "output",
        "gsva",
        "plots_heatmap",
        str_c("gsva_", name, "_split_day_response.png")
      ),
      width = 12,
      height = 6,
      units = "in",
      res = 600
    )

    plot(heatmap)

    dev.off()
  }
)
