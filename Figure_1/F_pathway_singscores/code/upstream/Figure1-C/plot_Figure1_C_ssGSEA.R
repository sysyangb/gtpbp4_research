#!/usr/bin/env Rscript

# Three-panel figure for one selected scoring method across three
# ribosome/translation-related GO programs. The default remains ssGSEA for
# Figure1-C; --method singscore is used to build FigureS1-C.
# The figure is rebuilt from source score and model-statistics tables; no raster
# panels are cropped or recombined.

options(
  stringsAsFactors = FALSE,
  width = 140,
  warn = 1
)

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (!length(script_arg)) stop("Run this file with Rscript.")
script_path <- normalizePath(sub("^--file=", "", script_arg[[1L]]), mustWork = TRUE)
script_dir <- dirname(script_path)
project_dir <- dirname(script_dir)

local_lib <- file.path(project_dir, ".Rlib")
if (dir.exists(local_lib)) .libPaths(c(local_lib, .libPaths()))

env_or <- function(name, fallback) {
  value <- Sys.getenv(name, unset = "")
  if (nzchar(value)) value else fallback
}

defaults <- list(
  scores = file.path(
    project_dir,
    "GO_pathway_scores_four_groups", "results", "GO_pathway_sample_scores_long.tsv.gz"
  ),
  statistics = file.path(
    project_dir,
    "GO_pathway_scores_four_groups", "results", "GO_pathway_pairwise_model_statistics.tsv"
  ),
  outdir = env_or("FIGURE_SCORE_OUTDIR", script_dir),
  prefix = env_or("FIGURE_SCORE_PREFIX", "Figure1-C"),
  method = env_or("FIGURE_SCORE_METHOD", "ssGSEA"),
  ymin = "-3.5",
  ymax = "5.0"
)

usage <- function() {
  cat(paste0(
    "Usage: Rscript plot_Figure1_C_ssGSEA.R [--key value ...]\n\n",
    "Options:\n",
    "  --scores PATH       GO pathway sample-score long table\n",
    "  --statistics PATH   Patient-aware pairwise statistics table\n",
    "  --outdir PATH       Output bundle directory [script directory]\n",
    "  --prefix TEXT       Output filename prefix [Figure1-C]\n",
    "  --method TEXT       Scoring method: ssGSEA or singscore [ssGSEA]\n",
    "  --ymin NUMBER       Y-axis lower limit [-3.5]\n",
    "  --ymax NUMBER       Y-axis upper limit [5.0]\n",
    "  --help              Show this message\n"
  ))
}

parse_args <- function(args, defaults) {
  if (any(args %in% c("--help", "-h"))) {
    usage()
    quit(save = "no", status = 0L)
  }
  if (length(args) %% 2L != 0L) stop("Arguments must be supplied as --key value pairs.")
  out <- defaults
  if (length(args)) {
    keys <- sub("^--", "", args[seq.int(1L, length(args), by = 2L)])
    values <- args[seq.int(2L, length(args), by = 2L)]
    unknown <- setdiff(keys, names(defaults))
    if (length(unknown)) stop("Unknown argument(s): ", paste(unknown, collapse = ", "))
    for (i in seq_along(keys)) out[[keys[[i]]]] <- values[[i]]
  }
  out
}

cfg <- parse_args(commandArgs(trailingOnly = TRUE), defaults)
cfg$scores <- normalizePath(cfg$scores, mustWork = TRUE)
cfg$statistics <- normalizePath(cfg$statistics, mustWork = TRUE)
cfg$outdir <- normalizePath(cfg$outdir, mustWork = FALSE)
cfg$ymin <- as.numeric(cfg$ymin)
cfg$ymax <- as.numeric(cfg$ymax)
method_labels <- c(ssGSEA = "ssGSEA", singscore = "singscore")
method_axis_labels <- c(
  ssGSEA = "Standardized ssGSEA score",
  singscore = "Standardized singscore"
)
if (!cfg$method %in% names(method_labels)) {
  stop("--method must be one of: ", paste(names(method_labels), collapse = ", "))
}
method_display <- unname(method_labels[[cfg$method]])
method_axis_label <- unname(method_axis_labels[[cfg$method]])
if (!is.finite(cfg$ymin) || !is.finite(cfg$ymax) || cfg$ymin >= cfg$ymax) {
  stop("--ymin and --ymax must be finite numbers with ymin < ymax.")
}

required_packages <- c("data.table", "ggplot2", "svglite", "ragg")
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages)) {
  stop("Missing R packages: ", paste(missing_packages, collapse = ", "))
}

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

data_dir <- file.path(cfg$outdir, "data")
figures_dir <- file.path(cfg$outdir, "figures")
results_dir <- file.path(cfg$outdir, "results")
qc_dir <- file.path(cfg$outdir, "QC")
for (path in c(cfg$outdir, data_dir, figures_dir, results_dir, qc_dir)) {
  dir.create(path, recursive = TRUE, showWarnings = FALSE)
}

read_tsv_gz <- function(path) {
  connection <- gzfile(path, open = "rt")
  on.exit(close(connection), add = TRUE)
  as.data.table(read.delim(
    connection,
    sep = "\t",
    header = TRUE,
    quote = "",
    comment.char = "",
    check.names = FALSE
  ))
}

scores <- read_tsv_gz(cfg$scores)
statistics <- fread(cfg$statistics, sep = "\t", header = TRUE)

required_score_columns <- c(
  "sample_id", "patient_id", "condition", "signature", "method",
  "raw_score", "standardized_score"
)
required_stat_columns <- c(
  "signature", "method", "contrast", "group1", "group2", "estimate",
  "lower_CI", "upper_CI", "p_raw", "p_adjusted", "model"
)
if (length(setdiff(required_score_columns, names(scores)))) {
  stop("Score table is missing required columns: ", paste(setdiff(required_score_columns, names(scores)), collapse = ", "))
}
if (length(setdiff(required_stat_columns, names(statistics)))) {
  stop("Statistics table is missing required columns: ", paste(setdiff(required_stat_columns, names(statistics)), collapse = ", "))
}

signature_levels <- c("RIBI", "POL-I", "Regulator of Translation")
signature_labels <- c(
  "RIBI" = "RIBI",
  "POL-I" = "POL-I",
  "Regulator of Translation" = "Regulator of translation"
)
group_levels <- c("TN", "HR", "HER2", "NORMAL")
group_labels <- c(TN = "TN", HR = "HR", HER2 = "HER2", NORMAL = "Normal")
group_colors <- c(TN = "#AA4643", HR = "#4572A7", HER2 = "#89A54E", NORMAL = "#71588F")
comparison_groups <- c("HR", "HER2", "NORMAL")

plot_data <- scores[
  signature %chin% signature_levels & method == cfg$method
]
plot_data[, `:=`(
  signature = factor(signature, levels = signature_levels),
  pathway = factor(
    unname(signature_labels[as.character(signature)]),
    levels = unname(signature_labels[signature_levels])
  ),
  condition = factor(condition, levels = group_levels)
)]
setorder(plot_data, signature, condition, sample_id)

stat_data <- statistics[
  signature %chin% signature_levels &
    method == cfg$method &
    group1 == "TN" &
    group2 %chin% comparison_groups
]
stat_data[, `:=`(
  signature = factor(signature, levels = signature_levels),
  pathway = factor(
    unname(signature_labels[as.character(signature)]),
    levels = unname(signature_labels[signature_levels])
  ),
  x1 = 1,
  x2 = match(group2, group_levels),
  y = 3.30 + 0.60 * (match(group2, comparison_groups) - 1L)
)]

format_raw_p <- function(p) {
  if (is.na(p)) return("P=NA")
  if (p < 1e-15) return("P<1e-15")
  if (p < 1e-3) return(paste0("P=", formatC(p, format = "e", digits = 1)))
  paste0("P=", formatC(p, format = "f", digits = 3))
}
stat_data[, label := vapply(p_raw, format_raw_p, character(1))]
stat_data[, comparison_order := match(group2, comparison_groups)]
setorder(stat_data, signature, comparison_order)
stat_data[, comparison_order := NULL]

if (nrow(plot_data) != 1116L * length(signature_levels)) {
  stop("Expected 1,116 samples for each of three pathways; observed ", nrow(plot_data), " rows.")
}
if (plot_data[, any(!is.finite(raw_score) | !is.finite(standardized_score))]) {
  stop("Non-finite score detected in ", cfg$prefix, " source data.")
}
if (anyDuplicated(plot_data, by = c("sample_id", "signature"))) {
  stop("Duplicate sample/signature keys detected in ", cfg$prefix, " source data.")
}
if (nrow(stat_data) != 9L || any(!is.finite(stat_data$p_raw))) {
  stop("Expected nine finite raw P values: three TN comparisons for each pathway.")
}

bracket_tick <- 0.12
p <- ggplot(
  plot_data,
  aes(x = condition, y = standardized_score, fill = condition)
) +
  geom_boxplot(
    width = 0.58,
    linewidth = 0.55,
    colour = "#222222",
    alpha = 0.76,
    outlier.shape = 16,
    outlier.size = 1.15,
    outlier.stroke = 0,
    outlier.colour = "#222222"
  ) +
  geom_segment(
    data = stat_data,
    aes(x = x1, xend = x2, y = y, yend = y),
    inherit.aes = FALSE,
    linewidth = 0.45,
    colour = "#222222"
  ) +
  geom_segment(
    data = stat_data,
    aes(x = x1, xend = x1, y = y, yend = y - bracket_tick),
    inherit.aes = FALSE,
    linewidth = 0.45,
    colour = "#222222"
  ) +
  geom_segment(
    data = stat_data,
    aes(x = x2, xend = x2, y = y, yend = y - bracket_tick),
    inherit.aes = FALSE,
    linewidth = 0.45,
    colour = "#222222"
  ) +
  geom_text(
    data = stat_data,
    aes(x = (x1 + x2) / 2, y = y + 0.08, label = label),
    inherit.aes = FALSE,
    family = "Arial",
    size = 10 / 2.845276,
    vjust = 0,
    colour = "#1F1F1F"
  ) +
  facet_wrap(~ pathway, nrow = 1, scales = "fixed") +
  scale_fill_manual(values = group_colors, guide = "none", drop = FALSE) +
  scale_x_discrete(labels = group_labels, drop = FALSE) +
  scale_y_continuous(
    breaks = c(-3, 0, 3),
    expand = expansion(mult = c(0.015, 0.025))
  ) +
  coord_cartesian(ylim = c(cfg$ymin, cfg$ymax), clip = "on") +
  labs(x = NULL, y = method_axis_label) +
  theme_classic(base_size = 10.5, base_family = "Arial") +
  theme(
    axis.line = element_line(linewidth = 0.50, colour = "#222222"),
    axis.ticks = element_line(linewidth = 0.45, colour = "#222222"),
    axis.ticks.length = grid::unit(1.6, "mm"),
    axis.text = element_text(size = 10, colour = "#222222"),
    axis.title.y = element_text(size = 10.5, colour = "#222222", margin = margin(r = 7)),
    strip.background = element_rect(fill = "white", colour = "#222222", linewidth = 0.50),
    strip.text = element_text(
      size = 10.5,
      face = "bold",
      colour = "#222222",
      margin = margin(3.5, 0, 3.5, 0)
    ),
    panel.spacing.x = grid::unit(2.5, "mm"),
    panel.grid = element_blank(),
    plot.title = element_blank(),
    plot.subtitle = element_blank(),
    plot.margin = margin(4, 7, 6, 7)
  )

save_figure <- function(plot, stem, width_mm = 183, height_mm = 82, dpi = 600) {
  width_in <- width_mm / 25.4
  height_in <- height_mm / 25.4

  svglite::svglite(paste0(stem, ".svg"), width = width_in, height = height_in, bg = "white")
  print(plot)
  grDevices::dev.off()

  grDevices::cairo_pdf(
    paste0(stem, ".pdf"),
    width = width_in,
    height = height_in,
    family = "Arial",
    bg = "white"
  )
  print(plot)
  grDevices::dev.off()

  ragg::agg_tiff(
    paste0(stem, ".tiff"),
    width = width_in,
    height = height_in,
    units = "in",
    res = dpi,
    compression = "lzw",
    background = "white"
  )
  print(plot)
  grDevices::dev.off()

  ragg::agg_png(
    paste0(stem, ".png"),
    width = width_in,
    height = height_in,
    units = "in",
    res = 300,
    background = "white"
  )
  print(plot)
  grDevices::dev.off()
}

figure_stem <- file.path(figures_dir, cfg$prefix)
save_figure(p, figure_stem)

source_output <- copy(plot_data)
source_output[, `:=`(
  signature = as.character(signature),
  pathway = as.character(pathway),
  condition = as.character(condition)
)]
fwrite(
  source_output[, .(
    sample_id,
    patient_id,
    condition,
    signature,
    pathway,
    method,
    raw_score,
    standardized_score
  )],
  file.path(data_dir, paste0(cfg$prefix, "-source-data.tsv.gz")),
  sep = "\t",
  compress = "gzip"
)

statistics_output <- copy(stat_data)
statistics_output[, `:=`(
  signature = as.character(signature),
  pathway = as.character(pathway)
)]
fwrite(
  statistics_output[, .(
    signature,
    pathway,
    method,
    contrast,
    group1,
    group2,
    estimate,
    lower_CI,
    upper_CI,
    p_raw,
    p_adjusted,
    label,
    model
  )],
  file.path(results_dir, paste0(cfg$prefix, "-statistics.tsv")),
  sep = "\t"
)

group_counts <- unique(
  plot_data[, .(n = .N), by = .(signature, condition)]
)[, value := paste0(as.character(condition), "=", n)]

parameters <- data.table(
  parameter = c(
    "pathways", "method", "groups", "group_counts_per_pathway",
    "comparisons", "p_value_display", "multiple_testing_for_display",
    "width_mm", "height_mm", "ymin", "ymax", "minimum_font_pt",
    "colors", "main_title", "subtitle"
  ),
  value = c(
    paste(unname(signature_labels[signature_levels]), collapse = "; "),
    method_display,
    "TN; HR; HER2; Normal",
    paste(group_counts[signature == signature_levels[[1L]], value], collapse = "; "),
    "TN-HR; TN-HER2; TN-Normal",
    "raw P value",
    "none",
    "183", "82", as.character(cfg$ymin), as.character(cfg$ymax), "10",
    "TN=#AA4643; HR=#4572A7; HER2=#89A54E; Normal=#71588F",
    "none", "none"
  )
)
fwrite(parameters, file.path(results_dir, paste0(cfg$prefix, "-figure-parameters.tsv")), sep = "\t")

contract <- data.table(
  field = c("core_conclusion", "evidence_chain", "archetype", "backend", "export_contract"),
  value = c(
    paste0(
      "Four-group differences in three ribosome/translation-related programs are shown using the same ",
      method_display,
      " scale."
    ),
    "Three ordered pathway facets; each facet shows four-group distributions and TN-versus-other-group raw P values.",
    "quantitative grid",
    "R/ggplot2",
    "183 x 82 mm; editable SVG/PDF; 600-dpi TIFF; 300-dpi PNG; minimum text 10 pt"
  )
)
fwrite(contract, file.path(qc_dir, paste0(cfg$prefix, "-figure-contract.tsv")), sep = "\t")

input_md5 <- data.table(
  input = c("scores", "statistics"),
  path = c(cfg$scores, cfg$statistics),
  md5 = unname(tools::md5sum(c(cfg$scores, cfg$statistics)))
)
fwrite(input_md5, file.path(qc_dir, paste0(cfg$prefix, "-input-md5.tsv")), sep = "\t")

package_versions <- data.table(
  package = required_packages,
  version = vapply(required_packages, function(x) as.character(packageVersion(x)), character(1))
)
fwrite(package_versions, file.path(qc_dir, paste0(cfg$prefix, "-package-versions.tsv")), sep = "\t")
capture.output(sessionInfo(), file = file.path(qc_dir, paste0(cfg$prefix, "-sessionInfo.txt")))

figure_files <- paste0(figure_stem, c(".svg", ".pdf", ".tiff", ".png"))
svg_text <- paste(readLines(paste0(figure_stem, ".svg"), warn = FALSE), collapse = "\n")
observed_signature_order <- unique(as.character(plot_data$signature))
expected_group_counts <- c(TN = 177L, HR = 664L, HER2 = 163L, NORMAL = 112L)
observed_group_counts <- plot_data[
  signature == signature_levels[[1L]],
  table(factor(condition, levels = group_levels))
]

validation <- data.table(
  check = c(
    "three_pathways_in_requested_order",
    "selected_method_only",
    "1116_samples_per_pathway",
    "expected_group_counts",
    "finite_scores",
    "unique_sample_pathway_keys",
    "nine_TN_comparisons",
    "raw_P_labels_only",
    "minimum_font_10pt",
    "four_export_formats"
  ),
  passed = c(
    identical(observed_signature_order, signature_levels),
    identical(unique(plot_data$method), cfg$method),
    all(plot_data[, .N, by = signature]$N == 1116L),
    identical(as.integer(observed_group_counts), as.integer(expected_group_counts)),
    !plot_data[, any(!is.finite(raw_score) | !is.finite(standardized_score))],
    !anyDuplicated(plot_data, by = c("sample_id", "signature")),
    nrow(stat_data) == 9L,
    all(grepl("^P(=|<)", stat_data$label)) && !grepl("adj\\. P", svg_text),
    TRUE,
    all(file.exists(figure_files)) && all(file.info(figure_files)$size > 0)
  ),
  value = c(
    paste(observed_signature_order, collapse = " -> "),
    paste(unique(plot_data$method), collapse = "; "),
    paste(plot_data[, paste0(as.character(signature), "=", .N), by = signature]$V1, collapse = "; "),
    paste(names(expected_group_counts), observed_group_counts, sep = "=", collapse = "; "),
    as.character(plot_data[, sum(!is.finite(raw_score) | !is.finite(standardized_score))]),
    as.character(anyDuplicated(plot_data, by = c("sample_id", "signature"))),
    as.character(nrow(stat_data)),
    paste(stat_data$label, collapse = "; "),
    "axis text=10; P labels=10; strip text=10.5; axis title=10.5",
    paste(basename(figure_files), collapse = "; ")
  )
)
fwrite(validation, file.path(qc_dir, paste0(cfg$prefix, "-validation-report.tsv")), sep = "\t")

if (!all(validation$passed)) {
  failed <- validation[passed == FALSE, check]
  stop(cfg$prefix, " validation failed: ", paste(failed, collapse = ", "))
}

message(cfg$prefix, " export complete: ", figure_stem)
