#!/usr/bin/env Rscript

# Reproducible entry point for the FigureS1-C singscore panel.
# The shared plotting implementation is parameterized so Figure1-C and
# FigureS1-C remain visually identical apart from the scoring method.

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (!length(script_arg)) stop("Run this file with Rscript.")
script_path <- normalizePath(sub("^--file=", "", script_arg[[1L]]), mustWork = TRUE)
script_dir <- dirname(script_path)
project_dir <- dirname(script_dir)
shared_script <- file.path(project_dir, "Figure1-C", "plot_Figure1_C_ssGSEA.R")
if (!file.exists(shared_script)) stop("Shared plotting script not found: ", shared_script)

Sys.setenv(
  FIGURE_SCORE_METHOD = "singscore",
  FIGURE_SCORE_OUTDIR = script_dir,
  FIGURE_SCORE_PREFIX = "FigureS1-C"
)
source(shared_script, chdir = FALSE)
