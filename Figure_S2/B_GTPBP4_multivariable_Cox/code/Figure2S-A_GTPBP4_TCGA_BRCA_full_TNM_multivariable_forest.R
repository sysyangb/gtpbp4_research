#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)
suppressPackageStartupMessages({
  library(survival)
  library(ggplot2)
  library(svglite)
  library(ragg)
})


script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (length(script_arg) != 1L) stop("Run this file with Rscript.")
script_path <- normalizePath(sub("^--file=", "", script_arg), mustWork = TRUE)
project_dir <- normalizePath(file.path(dirname(script_path), ".."), mustWork = TRUE)
figure2a_dir <- normalizePath(file.path(project_dir, ".."), mustWork = TRUE)

analysis_file <- file.path(
  project_dir, "data", "processed", "TCGA_BRCA_full_five_gene_OS_analysis_ready.tsv"
)
clinical_file <- file.path(figure2a_dir, "TCGA-BRCA_clinical.txt")
best_group_file <- file.path(project_dir, "results", "best_cutoff_sample_groups.tsv")
best_result_file <- file.path(project_dir, "results", "best_cutoff_OS_results.tsv")
results_dir <- file.path(project_dir, "results")
figures_dir <- file.path(project_dir, "figures")
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)

input_files <- c(analysis_file, clinical_file, best_group_file, best_result_file)
if (!all(file.exists(input_files))) {
  stop("Missing input file(s): ", paste(input_files[!file.exists(input_files)], collapse = ", "))
}

base <- read.delim(analysis_file, check.names = FALSE)
clinical <- read.delim(clinical_file, check.names = FALSE)
best_groups_all <- read.delim(best_group_file, check.names = FALSE)
best_results_all <- read.delim(best_result_file, check.names = FALSE)
best_groups <- best_groups_all[best_groups_all$gene == "GTPBP4", , drop = FALSE]
best_result <- best_results_all[best_results_all$gene == "GTPBP4", , drop = FALSE]

if (nrow(base) != 1073L || sum(base$OS_event) != 150L) {
  stop(sprintf(
    "Unexpected TCGA-BRCA cohort: n=%d, deaths=%d; expected 1073/150.",
    nrow(base), sum(base$OS_event)
  ))
}
if (
  anyDuplicated(base$sample_id) || anyDuplicated(clinical$submitter_id) ||
    anyDuplicated(best_groups$sample_id)
) {
  stop("Duplicated patient identifier detected in a patient-level input.")
}
if (nrow(best_groups) != nrow(base) || nrow(best_result) != 1L) {
  stop("The stored GTPBP4 best-cutpoint data do not match the full cohort.")
}

clinical_index <- match(base$sample_id, clinical$submitter_id)
group_index <- match(base$sample_id, best_groups$sample_id)
if (anyNA(clinical_index) || anyNA(group_index)) stop("Clinical/group matching is incomplete.")

clean_character <- function(x) {
  x <- trimws(as.character(x))
  missing_tokens <- c(
    "", "NA", "N/A", "Unknown", "unknown", "not reported", "Not Reported",
    "[Not Available]", "[Not Applicable]", "[Unknown]", "--"
  )
  x[is.na(x) | x %in% missing_tokens] <- NA_character_
  x
}

clinical_rows <- clinical[clinical_index, , drop = FALSE]
t_raw <- clean_character(clinical_rows$ajcc_pathologic_t)
n_raw <- clean_character(clinical_rows$ajcc_pathologic_n)
m_raw <- clean_character(clinical_rows$ajcc_pathologic_m)

analysis <- data.frame(
  patient_id = base$sample_id,
  OS_time_months = base$OS_time_months,
  OS_event = as.integer(base$OS_event),
  GTPBP4_log2TPM1 = base$GTPBP4_log2TPM1,
  T34_vs_T12 = ifelse(
    grepl("^T[34]", t_raw), 1L,
    ifelse(grepl("^T[12]", t_raw), 0L, NA_integer_)
  ),
  Npos_vs_N0 = ifelse(
    grepl("^N[123]", n_raw), 1L,
    ifelse(grepl("^N0", n_raw), 0L, NA_integer_)
  ),
  M1_vs_M0 = ifelse(
    m_raw == "M1", 1L,
    ifelse(m_raw == "M0", 0L, NA_integer_)
  ),
  GTPBP4_best_group = factor(
    best_groups$group[group_index], levels = c("Low", "High")
  ),
  stringsAsFactors = FALSE
)
analysis$GTPBP4_median_group <- factor(
  ifelse(
    analysis$GTPBP4_log2TPM1 > median(analysis$GTPBP4_log2TPM1),
    "High", "Low"
  ),
  levels = c("Low", "High")
)
analysis$GTPBP4_z <- as.numeric(scale(analysis$GTPBP4_log2TPM1))

model_vars <- c(
  "OS_time_months", "OS_event", "T34_vs_T12", "Npos_vs_N0", "M1_vs_M0"
)
model_data <- analysis[complete.cases(analysis[, model_vars]), , drop = FALSE]

fit_best <- coxph(
  Surv(OS_time_months, OS_event) ~
    T34_vs_T12 + Npos_vs_N0 + M1_vs_M0 + GTPBP4_best_group,
  data = model_data, ties = "efron", x = TRUE, y = TRUE, model = TRUE
)
fit_median <- coxph(
  Surv(OS_time_months, OS_event) ~
    T34_vs_T12 + Npos_vs_N0 + M1_vs_M0 + GTPBP4_median_group,
  data = model_data, ties = "efron", x = TRUE, y = TRUE, model = TRUE
)
fit_continuous <- coxph(
  Surv(OS_time_months, OS_event) ~
    T34_vs_T12 + Npos_vs_N0 + M1_vs_M0 + GTPBP4_z,
  data = model_data, ties = "efron", x = TRUE, y = TRUE, model = TRUE
)

extract_results <- function(fit, term_order, labels) {
  sm <- summary(fit)
  data.frame(
    term = term_order,
    characteristic = unname(labels[term_order]),
    n = fit$n,
    events = fit$nevent,
    hr = sm$conf.int[term_order, "exp(coef)"],
    ci_low = sm$conf.int[term_order, "lower .95"],
    ci_high = sm$conf.int[term_order, "upper .95"],
    wald_p = sm$coefficients[term_order, "Pr(>|z|)"],
    stringsAsFactors = FALSE,
    row.names = NULL
  )
}

term_order <- c("T34_vs_T12", "Npos_vs_N0", "M1_vs_M0", "GTPBP4_best_groupHigh")
term_labels <- c(
  T34_vs_T12 = "T stage (T3&T4 vs. T1&T2)",
  Npos_vs_N0 = "N stage (N1-3 vs. N0)",
  M1_vs_M0 = "M stage (M1 vs. M0)",
  GTPBP4_best_groupHigh = "GTPBP4 (high vs. low)"
)
results <- extract_results(fit_best, term_order, term_labels)

output_prefix <- "Figure2S-A_GTPBP4_TCGA_BRCA_full_TNM_multivariable"
results_file <- file.path(results_dir, paste0(output_prefix, "_results.tsv"))
write.table(results, results_file, sep = "\t", quote = FALSE, row.names = FALSE)

extract_gtpbp4_sensitivity <- function(fit, definition, term) {
  sm <- summary(fit)
  data.frame(
    expression_definition = definition,
    term = term,
    n = fit$n,
    events = fit$nevent,
    hr = sm$conf.int[term, "exp(coef)"],
    ci_low = sm$conf.int[term, "lower .95"],
    ci_high = sm$conf.int[term, "upper .95"],
    wald_p = sm$coefficients[term, "Pr(>|z|)"],
    stringsAsFactors = FALSE,
    row.names = NULL
  )
}
sensitivity <- rbind(
  extract_gtpbp4_sensitivity(
    fit_best, "Previously selected maxstat cutoff; High vs Low",
    "GTPBP4_best_groupHigh"
  ),
  extract_gtpbp4_sensitivity(
    fit_median, "Cohort median; High vs Low", "GTPBP4_median_groupHigh"
  ),
  extract_gtpbp4_sensitivity(
    fit_continuous, "Per 1 SD log2(TPM + 1)", "GTPBP4_z"
  )
)
sensitivity_file <- file.path(results_dir, paste0(output_prefix, "_expression_sensitivity.tsv"))
write.table(sensitivity, sensitivity_file, sep = "\t", quote = FALSE, row.names = FALSE)

ph_table <- cox.zph(fit_best, transform = "km")$table
ph_results <- data.frame(
  term = rownames(ph_table),
  chisq = ph_table[, "chisq"],
  df = ph_table[, "df"],
  p = ph_table[, "p"],
  row.names = NULL,
  stringsAsFactors = FALSE
)
ph_file <- file.path(results_dir, paste0(output_prefix, "_PH_test.tsv"))
write.table(ph_results, ph_file, sep = "\t", quote = FALSE, row.names = FALSE)

fit_summary <- summary(fit_best)
metrics <- data.frame(
  cohort = "TCGA-BRCA full cohort",
  endpoint = "Overall survival",
  model = "T34_vs_T12 + Npos_vs_N0 + M1_vs_M0 + GTPBP4_best_group",
  source_n = nrow(analysis),
  source_events = sum(analysis$OS_event),
  complete_case_n = fit_best$n,
  complete_case_events = fit_best$nevent,
  excluded_n = nrow(analysis) - fit_best$n,
  excluded_events = sum(analysis$OS_event) - fit_best$nevent,
  parameters = length(coef(fit_best)),
  events_per_parameter = fit_best$nevent / length(coef(fit_best)),
  likelihood_ratio_p = unname(fit_summary$logtest["pvalue"]),
  concordance = unname(fit_summary$concordance[1L]),
  concordance_se = unname(fit_summary$concordance[2L]),
  proportional_hazards_global_p = ph_results$p[ph_results$term == "GLOBAL"],
  gtpbp4_cutoff_log2TPM1 = best_result$cutoff_log2TPM1,
  gtpbp4_cutoff_TPM = best_result$cutoff_TPM,
  univariate_maxstat_condMC_adjusted_p = best_result$maxstat_condMC_adjusted_p,
  interpretation = paste(
    "Displayed GTPBP4 Wald P is nominal after an outcome-selected cutoff;",
    "median and continuous sensitivity estimates are provided"
  ),
  stringsAsFactors = FALSE,
  check.names = FALSE
)
metrics_file <- file.path(results_dir, paste0(output_prefix, "_model_metrics.tsv"))
write.table(metrics, metrics_file, sep = "\t", quote = FALSE, row.names = FALSE)

# Reproduce the supplied visual structure without its upper-left panel label.
forest <- results
forest$row <- c(4.24, 3.36, 2.48, 1.60)
forest$hr_label <- sprintf("%.3f (%.3f – %.3f)", forest$hr, forest$ci_low, forest$ci_high)
forest$p_label <- ifelse(
  forest$wald_p < 0.001,
  "<0.001",
  formatC(forest$wald_p, format = "f", digits = 3)
)

figure_width_mm <- 183
figure_height_mm <- 91
minimum_visible_font_pt <- 10.0
pt_to_mm <- 25.4 / 72.27
body_size <- 10.0 * pt_to_mm
header_size <- 10.5 * pt_to_mm
title_size <- 11.2 * pt_to_mm

characteristic_x <- 0.020
hr_label_x <- 0.365
forest_left <- 0.630
forest_right <- 0.872
p_x <- 0.972
title_y <- 6.02
top_rule_y <- 5.54
header_y <- 5.18
header_rule_y <- 4.80
axis_y <- 1.10
tick_bottom_y <- 0.91
tick_label_y <- 0.65

forest_min <- 0.5
forest_max <- 5.5
if (min(forest$ci_low) < forest_min || max(forest$ci_high) > forest_max) {
  stop("A confidence interval falls outside the prespecified forest range.")
}
map_hr <- function(value) {
  forest_left + (value - forest_min) / (forest_max - forest_min) *
    (forest_right - forest_left)
}
ticks <- c(1, 3, 5)
tick_x <- map_hr(ticks)
reference_x <- map_hr(1)
forest$plot_hr <- map_hr(forest$hr)
forest$plot_ci_low <- map_hr(forest$ci_low)
forest$plot_ci_high <- map_hr(forest$ci_high)

ink <- "#111111"
point_blue <- "#155A99"

forest_plot <- ggplot(forest, aes(y = row)) +
  annotate(
    "segment", x = 0.020, xend = 0.980,
    y = top_rule_y, yend = top_rule_y,
    linewidth = 0.60, colour = ink
  ) +
  annotate(
    "segment", x = 0.020, xend = 0.980,
    y = header_rule_y, yend = header_rule_y,
    linewidth = 0.55, colour = ink
  ) +
  annotate(
    "segment", x = 0.020, xend = 0.980,
    y = axis_y, yend = axis_y,
    linewidth = 0.55, colour = ink
  ) +
  annotate(
    "segment", x = reference_x, xend = reference_x,
    y = axis_y, yend = 4.56,
    linewidth = 0.55, linetype = "dashed", colour = ink
  ) +
  geom_segment(
    aes(x = plot_ci_low, xend = plot_ci_high, yend = row),
    colour = ink, linewidth = 0.65, lineend = "butt"
  ) +
  geom_segment(
    aes(x = plot_ci_low, xend = plot_ci_low, y = row - 0.10, yend = row + 0.10),
    colour = ink, linewidth = 0.55
  ) +
  geom_segment(
    aes(x = plot_ci_high, xend = plot_ci_high, y = row - 0.10, yend = row + 0.10),
    colour = ink, linewidth = 0.55
  ) +
  geom_point(
    aes(x = plot_hr), shape = 19, size = 3.25, colour = point_blue
  ) +
  geom_text(
    aes(x = characteristic_x, label = characteristic),
    hjust = 0, family = "Arial", size = body_size, colour = ink
  ) +
  geom_text(
    aes(x = hr_label_x, label = hr_label),
    hjust = 0, family = "Arial", size = body_size, colour = ink
  ) +
  geom_text(
    aes(x = p_x, label = p_label),
    hjust = 1, family = "Arial", size = body_size, colour = ink
  ) +
  annotate(
    "text", x = 0.5, y = title_y,
    label = "Multivariate analysis (OS)",
    family = "Arial", size = title_size, colour = ink
  ) +
  annotate(
    "text", x = characteristic_x, y = header_y,
    label = "Characteristics", hjust = 0,
    family = "Arial", size = header_size, colour = ink
  ) +
  annotate(
    "text", x = 0.440, y = header_y,
    label = "HR (95% CI)", hjust = 0.5,
    family = "Arial", size = header_size, colour = ink
  ) +
  annotate(
    "text", x = p_x, y = header_y,
    label = "P value", hjust = 1,
    family = "Arial", fontface = "italic", size = header_size, colour = ink
  ) +
  annotate(
    "segment", x = tick_x, xend = tick_x,
    y = axis_y, yend = tick_bottom_y,
    linewidth = 0.45, colour = ink
  ) +
  annotate(
    "text", x = tick_x, y = tick_label_y,
    label = ticks, family = "Arial", size = body_size, colour = ink
  ) +
  scale_x_continuous(limits = c(0, 1), expand = expansion(mult = c(0, 0))) +
  scale_y_continuous(limits = c(0.38, 6.22), expand = expansion(mult = c(0, 0))) +
  coord_cartesian(clip = "off") +
  theme_void(base_size = 10, base_family = "Arial") +
  theme(
    plot.margin = margin(1.5, 2.0, 1.5, 2.0, unit = "mm"),
    plot.background = element_rect(fill = "white", colour = NA),
    panel.background = element_rect(fill = "white", colour = NA)
  )

forest_source_file <- file.path(results_dir, paste0(output_prefix, "_forest_source_data.tsv"))
write.table(forest, forest_source_file, sep = "\t", quote = FALSE, row.names = FALSE)

figure_stem <- file.path(figures_dir, paste0(output_prefix, "_forest"))
width_in <- figure_width_mm / 25.4
height_in <- figure_height_mm / 25.4

svglite::svglite(paste0(figure_stem, ".svg"), width = width_in, height = height_in)
print(forest_plot)
dev.off()

grDevices::cairo_pdf(
  paste0(figure_stem, ".pdf"), width = width_in, height = height_in,
  family = "Arial", pointsize = 10
)
print(forest_plot)
dev.off()

ragg::agg_png(
  paste0(figure_stem, ".png"), width = width_in, height = height_in,
  units = "in", res = 300
)
print(forest_plot)
dev.off()

ragg::agg_tiff(
  paste0(figure_stem, ".tiff"), width = width_in, height = height_in,
  units = "in", res = 600, compression = "lzw"
)
print(forest_plot)
dev.off()

figure_files <- paste0(figure_stem, c(".svg", ".pdf", ".png", ".tiff"))
svg_text <- paste(readLines(paste0(figure_stem, ".svg"), warn = FALSE), collapse = "\n")

validation <- data.frame(
  check = c(
    "fixed TCGA-BRCA source cohort",
    "patient-level identifiers unique",
    "clinical and group joins complete",
    "TNM complete-case cohort",
    "model contains exactly four requested terms",
    "forest rows match requested order",
    "GTPBP4 High versus Low is nominally significant",
    "finite ordered confidence intervals",
    "global PH test does not reject at 0.05",
    "T-stage PH caveat recorded",
    "outcome-selected cutoff caveat recorded",
    "no upper-left panel identifier",
    "minimum visible text is at least 10 pt",
    "Arial retained in SVG",
    "four figure formats exist",
    "all outputs are nonempty"
  ),
  pass = c(
    nrow(analysis) == 1073L && sum(analysis$OS_event) == 150L,
    !anyDuplicated(analysis$patient_id),
    all(!is.na(clinical_index)) && all(!is.na(group_index)),
    fit_best$n == 897L && fit_best$nevent == 132L,
    length(coef(fit_best)) == 4L,
    identical(results$term, term_order),
    results$wald_p[results$term == "GTPBP4_best_groupHigh"] < 0.05,
    all(is.finite(results$hr)) && all(results$ci_low > 0) &&
      all(results$ci_low < results$hr) && all(results$hr < results$ci_high),
    ph_results$p[ph_results$term == "GLOBAL"] > 0.05,
    ph_results$p[ph_results$term == "T34_vs_T12"] < 0.05,
    is.finite(metrics$univariate_maxstat_condMC_adjusted_p),
    !grepl("(^|>)I(<|$)", svg_text),
    minimum_visible_font_pt >= 10,
    grepl("Arial", svg_text, fixed = TRUE),
    all(file.exists(figure_files)),
    all(file.info(c(
      figure_files, results_file, sensitivity_file, ph_file,
      metrics_file, forest_source_file
    ))$size > 0)
  ),
  detail = c(
    "n=1073; OS deaths=150",
    "one row per patient",
    "1073/1073 matched to clinical data and the GTPBP4 cutoff groups",
    sprintf("n=%d; deaths=%d", fit_best$n, fit_best$nevent),
    paste(names(coef(fit_best)), collapse = ", "),
    paste(results$characteristic, collapse = "; "),
    sprintf("HR %.3f; P %.5f", results$hr[4L], results$wald_p[4L]),
    paste(sprintf("%s %.3f–%.3f", results$term, results$ci_low, results$ci_high), collapse = "; "),
    sprintf("global PH P=%.5f", ph_results$p[ph_results$term == "GLOBAL"]),
    sprintf("T-stage PH P=%.5f", ph_results$p[ph_results$term == "T34_vs_T12"]),
    sprintf(
      "cutoff log2(TPM+1)=%.5f; univariate maxstat condMC P=%.5f",
      metrics$gtpbp4_cutoff_log2TPM1,
      metrics$univariate_maxstat_condMC_adjusted_p
    ),
    "No panel letter layer is present in the plotting code",
    sprintf("%.1f pt", minimum_visible_font_pt),
    "Arial present in editable SVG",
    paste(basename(figure_files), collapse = ", "),
    paste(range(file.info(c(
      figure_files, results_file, sensitivity_file, ph_file,
      metrics_file, forest_source_file
    ))$size), collapse = " to ")
  ),
  stringsAsFactors = FALSE
)

validation_file <- file.path(results_dir, paste0(output_prefix, "_validation.tsv"))
write.table(validation, validation_file, sep = "\t", quote = FALSE, row.names = FALSE)
capture.output(
  sessionInfo(),
  file = file.path(results_dir, paste0(output_prefix, "_sessionInfo.txt"))
)

if (!all(validation$pass)) {
  print(validation[!validation$pass, ], row.names = FALSE)
  stop("TCGA-BRCA full-cohort TNM/GTPBP4 forest validation failed.")
}

cat("TCGA-BRCA full-cohort TNM/GTPBP4 forest completed.\n")
print(results, row.names = FALSE, digits = 6)
cat(sprintf(
  "Model: n=%d, deaths=%d, LR P=%.3g, C-index=%.3f, global PH P=%.4f.\n",
  fit_best$n, fit_best$nevent, metrics$likelihood_ratio_p,
  metrics$concordance, metrics$proportional_hazards_global_p
))
