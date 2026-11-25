#!/usr/bin/env Rscript

cran_packages <- c(
  "readxl", "dplyr", "tidyr", "ggplot2", "readr", "matrixStats",
  "patchwork", "jsonlite", "scales", "statmod"
)
missing_cran <- cran_packages[
  !vapply(cran_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_cran)) {
  install.packages(missing_cran, repos = "https://cloud.r-project.org")
}

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager", repos = "https://cloud.r-project.org")
}
bioc_packages <- c("limma", "proDA")
missing_bioc <- bioc_packages[
  !vapply(bioc_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_bioc)) {
  BiocManager::install(missing_bioc, ask = FALSE, update = FALSE)
}

message("All FOT pipeline dependencies are available.")
