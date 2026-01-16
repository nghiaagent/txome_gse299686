here::i_am("r/00_define_colours.R")

# Attach necessary packages
library("colorspace")
library("grDevices")
library("tidyverse")

# Define colours for heatmaps
palette_heatmap <- list(
  `Day` = c(
    "D4" = palette.colors()[6],
    "D9" = palette.colors()[7],
    "D17" = palette.colors()[4],
    "D24" = palette.colors()[2]
  ),
  `Response` = c(
    "Responder" = hcl.colors(7, palette = "Purple-Brown")[2],
    "Non-Responder" = hcl.colors(7, palette = "Purple-Brown")[5]
  )
)

# Define colours for 3D volcano plots
palette_volcano3d <- c(
  palette.colors()[[9]],
  palette.colors()[[7]],
  palette.colors()[[5]],
  palette.colors()[[4]],
  palette.colors()[[3]],
  palette.colors()[[6]],
  palette.colors()[[8]]
)

message("Sourced ", here::here("r/00_define_colours.R"))