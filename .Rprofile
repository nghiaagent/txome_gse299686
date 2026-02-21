source("renv/activate.R")

# Source project setup files
library("tidyverse")

c(
  "r/00_load_genesets.R",
  "r/00_define_colours.R",
  "r/00_define_contrasts.R"
) %>%
  purrr::walk(\(x) {
    message(paste0("Sourcing ", here::here("scripts", x)))
    source(here::here(x), echo = FALSE, verbose = FALSE)
  })
