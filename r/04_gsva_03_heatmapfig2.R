here::i_am("r/04_gsva_03_heatmapfig2.R")

######################
# Build GSVA heatmap only for desired comparison
######################

# Attach packages
library("Biobase")
library("circlize")
library("colorspace")
library("ComplexHeatmap")
library("GSVA")
library("limma")
library("magrittr")
library("nghiaagentrutils")
library("tidyverse")
library("viridis")

# Set color scheme and breaks
## For gene set activity expression (GSVA ES)
col <- mako(n = 100) %>% darken(0.15)
breaks <- seq(-2, 2, length.out = 100)

# Load data
## GSVA ES scores
gsva_quant <- readRDS(
  here::here("output/gsva/gsva_quant.rds")
)

## Subset GSVA objects to only MSigDB Hallmark gene sets
## Only non-responders between D4, 9, 17
gsva_quant_subset <- gsva_quant$mh
gsva_quant_subset <- gsva_quant_subset[
  ,
  gsva_quant_subset$treatment_response == "NonResponder"
]
gsva_quant_subset <- gsva_quant_subset[
  ,
  gsva_quant_subset$day %in% c("D4", "D9", "D17")
]

# Get Z-scores
gsva_h_exprs_heatmap <- gsva_quant_subset %>%
  exprs() %>%
  t() %>%
  scale() %>%
  t()

# Create extrac factors for splitting column
split_day <- gsva_quant_subset$day

# Build annotation; include only necessary metadata
## Dataframe of annotation data
heatmap_anno_tibble <- gsva_quant_subset@phenoData %$%
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
      title_position = "topleft",
      legend_direction = "vertical",
      title_gp = gpar(fontsize = 8, fontface = "bold"),
      labels_gp = gpar(fontsize = 8, fontface = "bold")
    )
  )
)

# Format names of Hallmark gene sets
rownames(gsva_h_exprs_heatmap) <- rownames(gsva_h_exprs_heatmap) %>%
  recode_msigdbh()

# Create a heatmap to extract column orders from
heatmap_column_order <- Heatmap(
  gsva_h_exprs_heatmap,
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
  column_split = split_day,
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
  gsva_h_exprs_heatmap,
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
  row_km = 3,
  row_km_repeats = 5,
  cluster_row_slices = FALSE,
  cluster_rows = TRUE,
  show_row_dend = FALSE,
  row_title_side = "left",
  row_title_gp = gpar(fontsize = 10, fontface = "bold"),
  row_title_rot = 90,
  show_row_names = TRUE,
  row_names_gp = gpar(fontsize = 8),
  row_names_side = "left",
  row_names_max_width = max_text_width(
    rownames(gsva_h_exprs_heatmap),
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

# Export heatmaps
png(
  file = here::here(
    "output",
    "gsva",
    "plots_heatmap",
    str_c("gsva_fig2_mh_split_day_response_with_names.png")
  ),
  width = 8,
  height = 7,
  units = "in",
  res = 600
)

plot(heatmap_hallmark)

dev.off()
