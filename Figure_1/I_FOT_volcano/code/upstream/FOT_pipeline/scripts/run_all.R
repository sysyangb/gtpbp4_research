#!/usr/bin/env Rscript

locate_pipeline_root <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- args[grepl("^--file=", args)]
  if (!length(file_arg)) stop("Could not determine run_all.R path.")
  script_path <- normalizePath(sub("^--file=", "", file_arg[[1]]), mustWork = TRUE)
  dirname(dirname(script_path))
}

main <- function() {
  pipeline_root <- locate_pipeline_root()
  Sys.setenv(FOT_PIPELINE_ROOT = pipeline_root)
  source(file.path(pipeline_root, "scripts", "config.R"), local = .GlobalEnv)
  source(file.path(pipeline_root, "scripts", "analysis_helpers.R"), local = .GlobalEnv)

  log_path <- file.path(cfg$logs_dir, "run_all.log")
  log_con <- file(log_path, open = "wt", encoding = "UTF-8")
  sink(log_con, type = "output", split = TRUE)
  sink(log_con, type = "message")
  on.exit({
    while (sink.number(type = "message") > 2) sink(type = "message")
    while (sink.number(type = "output") > 0) sink(type = "output")
    close(log_con)
  }, add = TRUE)

  started <- Sys.time()
  message("FOT pipeline started at ", format(started, tz = "Asia/Shanghai", usetz = TRUE))
  set.seed(cfg$seed)

  parameter_table <- data.frame(
    parameter = c(
      "input_xlsx", "input_sheet", "main_samples", "excluded_sample",
      "FOT_numerator", "FOT_denominator", "iFOT_scale", "unique_peptides_min",
      "min_observed_per_group", "log2fc_threshold", "fdr_threshold",
      "raw_p_threshold", "primary_model", "primary_missing_value_policy",
      "dropout_model_sensitivity", "random_seed"
    ),
    value = c(
      cfg$input_xlsx,
      cfg$input_sheet,
      paste(cfg$main_samples, collapse = ","),
      cfg$excluded_sample,
      "vendor Normalized intensity (not iBAQ)",
      "sum of all positive non-decoy/non-contaminant normalized intensities per sample before Unique-peptide filtering",
      format(cfg$ifot_scale, scientific = FALSE),
      cfg$unique_peptides_min,
      cfg$min_observed_per_group,
      cfg$log2fc_threshold,
      cfg$fdr_threshold,
      cfg$raw_p_threshold,
      "limma eBayes(trend=TRUE, robust=TRUE) plus TREAT(lfc=1)",
      "No imputation; require at least two observed values per group; detection-only proteins reported separately",
      "proDA probabilistic intensity-dependent dropout model; UP>=5 and at least 2/5 observations",
      cfg$seed
    ),
    stringsAsFactors = FALSE
  )
  write_csv_stable(parameter_table, file.path(cfg$tables_dir, "analysis_parameters.csv"))

  source(file.path(pipeline_root, "scripts", "01_prepare_fot_qc.R"), local = .GlobalEnv)
  source(file.path(pipeline_root, "scripts", "02_differential_limma.R"), local = .GlobalEnv)
  source(file.path(pipeline_root, "scripts", "03_sensitivity_and_comparison.R"), local = .GlobalEnv)
  source(file.path(pipeline_root, "scripts", "04_proDA_missingness_model.R"), local = .GlobalEnv)
  source(file.path(pipeline_root, "scripts", "05_validate_outputs.R"), local = .GlobalEnv)

  package_names <- c(
    "R", "limma", "readxl", "dplyr", "tidyr", "readr", "ggplot2",
    "matrixStats", "patchwork", "jsonlite", "scales", "statmod", "proDA"
  )
  package_versions <- data.frame(
    package = package_names,
    version = c(
      paste(R.version$major, R.version$minor, sep = "."),
      vapply(package_names[-1], function(x) as.character(utils::packageVersion(x)), character(1))
    ),
    stringsAsFactors = FALSE
  )
  write_csv_stable(package_versions, file.path(cfg$tables_dir, "software_versions.csv"))
  writeLines(capture.output(utils::sessionInfo()), file.path(cfg$logs_dir, "sessionInfo.txt"))

  finished <- Sys.time()
  run_status <- list(
    status = "success",
    started = format(started, tz = "Asia/Shanghai", usetz = TRUE),
    finished = format(finished, tz = "Asia/Shanghai", usetz = TRUE),
    elapsed_seconds = as.numeric(difftime(finished, started, units = "secs")),
    pipeline_root = pipeline_root
  )
  jsonlite::write_json(
    run_status,
    file.path(cfg$logs_dir, "run_status.json"),
    pretty = TRUE,
    auto_unbox = TRUE
  )
  message("FOT pipeline completed successfully at ", run_status$finished)
  invisible(run_status)
}

main()
