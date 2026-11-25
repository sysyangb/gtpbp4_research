options(stringsAsFactors = FALSE)

pipeline_root <- Sys.getenv("FOT_PIPELINE_ROOT", unset = "")
if (!nzchar(pipeline_root)) {
  stop("FOT_PIPELINE_ROOT is not set. Run scripts/run_all.R instead.")
}
pipeline_root <- normalizePath(pipeline_root, mustWork = TRUE)
project_root <- normalizePath(file.path(pipeline_root, ".."), mustWork = TRUE)

cfg <- list(
  pipeline_root = pipeline_root,
  project_root = project_root,
  input_xlsx = file.path(
    project_root,
    "cloud_report_312553",
    "2-Basic_analysis",
    "MS_identified_information.xlsx"
  ),
  input_sheet = "Protein_quant",
  report_pdf = file.path(project_root, "cloud_report_312553", "report.pdf"),
  annotation_xlsx = file.path(
    project_root,
    "cloud_report_312553",
    "3-Functional_annotation",
    "Annotation_combine.xlsx"
  ),
  old_result_csv = file.path(
    project_root,
    "reanalysis_exclude_polysome3_20260901",
    "results",
    "protein_statistics_all.csv"
  ),
  old_up_502_tsv = file.path(
    project_root,
    "polysome_up_GO_BP_plus_AMFs_total_human_background_20260902",
    "source_data",
    "Polysome_up_502_genes.tsv"
  ),
  sample_ids = c(
    "free_1", "free_2", "free_3",
    "polysome_1", "polysome_2", "polysome_3"
  ),
  main_samples = c("free_1", "free_2", "free_3", "polysome_1", "polysome_2"),
  excluded_sample = "polysome_3",
  groups = c(
    free_1 = "Free", free_2 = "Free", free_3 = "Free",
    polysome_1 = "Polysome", polysome_2 = "Polysome",
    polysome_3 = "Polysome"
  ),
  unique_peptides_min = 5L,
  min_observed_per_group = 2L,
  log2fc_threshold = 1,
  fdr_threshold = 0.05,
  raw_p_threshold = 0.05,
  ifot_scale = 1e5,
  seed = 20260928L,
  results_dir = file.path(pipeline_root, "results"),
  tables_dir = file.path(pipeline_root, "results", "tables"),
  figures_dir = file.path(pipeline_root, "results", "figures"),
  data_dir = file.path(pipeline_root, "results", "data"),
  logs_dir = file.path(pipeline_root, "results", "logs")
)

invisible(lapply(
  c(cfg$results_dir, cfg$tables_dir, cfg$figures_dir, cfg$data_dir, cfg$logs_dir),
  dir.create,
  recursive = TRUE,
  showWarnings = FALSE
))

required_packages <- c(
  "readxl", "dplyr", "tidyr", "ggplot2", "readr", "limma", "matrixStats",
  "patchwork", "jsonlite", "scales", "proDA"
)
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages)) {
  stop("Missing required R packages: ", paste(missing_packages, collapse = ", "))
}

theme_fot <- function(base_size = 11) {
  ggplot2::theme_classic(base_size = base_size) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", hjust = 0),
      axis.title = ggplot2::element_text(face = "bold"),
      legend.title = ggplot2::element_text(face = "bold"),
      strip.background = ggplot2::element_blank(),
      strip.text = ggplot2::element_text(face = "bold")
    )
}

write_csv_stable <- function(x, path) {
  readr::write_csv(x, path, na = "")
  invisible(path)
}

save_plot_both <- function(plot, stem, width, height, dpi = 300) {
  png_path <- file.path(cfg$figures_dir, paste0(stem, ".png"))
  pdf_path <- file.path(cfg$figures_dir, paste0(stem, ".pdf"))
  ggplot2::ggsave(png_path, plot, width = width, height = height, dpi = dpi, bg = "white")
  ggplot2::ggsave(pdf_path, plot, width = width, height = height, device = grDevices::cairo_pdf)
  invisible(c(png_path, pdf_path))
}
