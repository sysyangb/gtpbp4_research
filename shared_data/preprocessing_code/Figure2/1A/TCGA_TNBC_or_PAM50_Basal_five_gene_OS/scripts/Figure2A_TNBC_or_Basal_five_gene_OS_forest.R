#!/usr/bin/env Rscript


options(stringsAsFactors = FALSE)
script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (length(script_arg) != 1L) stop("Run this script with Rscript.")
script_path <- normalizePath(sub("^--file=", "", script_arg), mustWork = TRUE)
project_dir <- normalizePath(file.path(dirname(script_path), ".."), mustWork = TRUE)
figure2a_dir <- normalizePath(file.path(project_dir, ".."), mustWork = TRUE)
workspace_dir <- normalizePath(file.path(figure2a_dir, "..", ".."), mustWork = TRUE)

# The existing Figure2/1A renv has the exact automatic-cutpoint packages.
renv_library <- file.path(
  figure2a_dir, "renv", "library", "macos", "R-4.5", "aarch64-apple-darwin20.0.0"
)
if (dir.exists(renv_library)) .libPaths(c(renv_library, .libPaths()))
suppressPackageStartupMessages({
  library(survival)
  library(survminer)
  library(maxstat)
  library(ggplot2)
  library(svglite)
  library(ragg)
})

genes <- c("DKC1", "GTPBP4", "ILF2", "NOP2", "RPP38")
# The historical figure used 0.20. Set FIGURE2A_MINPROP=0.10 to use the
# survminer::surv_cutpoint() default without overwriting the historical files.
minprop <- as.numeric(Sys.getenv("FIGURE2A_MINPROP", unset = "0.20"))
if (!is.finite(minprop) || minprop <= 0 || minprop >= 0.5) {
  stop("FIGURE2A_MINPROP must be between 0 and 0.5 (exclusive).")
}
condMC_B <- 99999L
base_seed <- 20261024L
prefix <- if (isTRUE(all.equal(minprop, 0.20))) {
  "Figure2A_TNBC_or_PAM50_Basal"
} else {
  paste0("Figure2A_TNBC_or_PAM50_Basal_minprop", sprintf("%02d", round(100 * minprop)))
}

basal_file <- file.path(
  figure2a_dir, "TCGA_PAM50_Basal_five_gene_OS", "data", "processed",
  "Figure2A_PAM50_Basal_sample_level_data.tsv"
)
tn_file <- file.path(
  workspace_dir, "01TCGA 筛选", "TN_five_gene_TPM_survival", "data",
  "TN_five_gene_TPM_survival_ready_data.tsv"
)
receptor_file <- file.path(
  workspace_dir, "TCGA", "TCGA-BRCA_receptor_subtype_by_sample.tsv"
)
os_file <- file.path(
  figure2a_dir, "PMC11391272_GTPBP4_reproduction", "source", "derived",
  "TCGA_BRCA_OS_selected_primary_samples.tsv"
)
source_files <- c(basal_file, tn_file, receptor_file, os_file)
if (!all(file.exists(source_files))) {
  stop("Missing input(s): ", paste(source_files[!file.exists(source_files)], collapse = ", "))
}

data_dir <- file.path(project_dir, "data", "processed")
results_dir <- file.path(project_dir, "results")
figures_dir <- file.path(project_dir, "figures")
for (path in c(data_dir, results_dir, figures_dir)) {
  dir.create(path, recursive = TRUE, showWarnings = FALSE)
}

basal <- read.delim(basal_file, check.names = FALSE)
tn <- read.delim(tn_file, check.names = FALSE)
receptor <- read.delim(receptor_file, check.names = FALSE)
os <- read.delim(os_file, check.names = FALSE)
required_expression <- c("patient_id", "sample_id", paste0(genes, "_TPM"),
                         paste0(genes, "_log2TPM1"))
stopifnot(
  all(required_expression %in% names(basal)),
  all(required_expression %in% names(tn)),
  all(c("patient_id", "sample_id", "sample_type", "label_recommended") %in% names(receptor)),
  all(c("patient_id", "sample_id", "event", "os_days") %in% names(os)),
  !anyDuplicated(basal$patient_id), !anyDuplicated(tn$patient_id),
  !anyDuplicated(os$patient_id),
  nrow(basal) == 139L, nrow(tn) == 169L,
  identical(unique(basal$PAM50), "Basal")
)

# Reconfirm strict receptor TNBC membership directly from the source labels.
strict_tn <- subset(receptor, label_recommended == "TN" & sample_type == "PRIMARY_TUMOR")
stopifnot(
  nrow(strict_tn) == 169L,
  !anyDuplicated(strict_tn$patient_id),
  setequal(strict_tn$patient_id, tn$patient_id),
  all(tn$sample_id == strict_tn$sample_id[match(tn$patient_id, strict_tn$patient_id)])
)

# Check that the two processed expression sources use the identical primary
# tumor and TPM values for every overlapping patient before deduplication.
overlap_ids <- intersect(tn$patient_id, basal$patient_id)
tn_overlap <- tn[match(overlap_ids, tn$patient_id), , drop = FALSE]
basal_overlap <- basal[match(overlap_ids, basal$patient_id), , drop = FALSE]
stopifnot(
  length(overlap_ids) == 101L,
  identical(tn_overlap$sample_id, basal_overlap$sample_id)
)
for (gene in genes) {
  for (suffix in c("_TPM", "_log2TPM1")) {
    field <- paste0(gene, suffix)
    if (!isTRUE(all.equal(tn_overlap[[field]], basal_overlap[[field]], tolerance = 1e-10))) {
      stop("Expression differs between TNBC and Basal sources for ", gene, suffix)
    }
  }
}

# One row/patient. The OS endpoint is the same TCGA selected-primary-sample
# table used in the original Basal forest, rather than mixing OS definitions.
cohort <- rbind(
  basal[, required_expression],
  tn[!tn$patient_id %in% basal$patient_id, required_expression]
)
cohort <- cohort[order(cohort$patient_id), , drop = FALSE]
rownames(cohort) <- NULL
cohort$TNBC <- cohort$patient_id %in% tn$patient_id
cohort$PAM50_Basal <- cohort$patient_id %in% basal$patient_id
os_match <- match(cohort$patient_id, os$patient_id)
stopifnot(
  nrow(cohort) == 207L,
  !anyDuplicated(cohort$patient_id),
  !anyNA(os_match),
  all(cohort$sample_id == os$sample_id[os_match])
)
cohort$OS_time_days <- os$os_days[os_match]
cohort$OS_event <- as.integer(os$event[os_match])
stopifnot(
  all(is.finite(cohort$OS_time_days) & cohort$OS_time_days >= 0),
  all(cohort$OS_event %in% c(0L, 1L)),
  sum(cohort$OS_event) == 36L
)
basal_match <- match(basal$patient_id, cohort$patient_id)
stopifnot(
  all(basal$OS_event == cohort$OS_event[basal_match]),
  isTRUE(all.equal(basal$OS_time_days, cohort$OS_time_days[basal_match]))
)
for (gene in genes) {
  tpm <- cohort[[paste0(gene, "_TPM")]]
  transformed <- cohort[[paste0(gene, "_log2TPM1")]]
  stopifnot(all(is.finite(tpm) & tpm >= 0),
            isTRUE(all.equal(log2(tpm + 1), transformed, tolerance = 1e-10)))
}

cohort_file <- file.path(data_dir, paste0(prefix, "_patient_level_union.tsv"))
write.table(cohort, cohort_file, sep = "\t", quote = FALSE, row.names = FALSE)

cohort_summary <- data.frame(
  membership = c("TNBC only", "PAM50 Basal only", "Both", "Union"),
  n = c(sum(cohort$TNBC & !cohort$PAM50_Basal),
        sum(!cohort$TNBC & cohort$PAM50_Basal),
        sum(cohort$TNBC & cohort$PAM50_Basal), nrow(cohort)),
  deaths = c(sum(cohort$OS_event[cohort$TNBC & !cohort$PAM50_Basal]),
             sum(cohort$OS_event[!cohort$TNBC & cohort$PAM50_Basal]),
             sum(cohort$OS_event[cohort$TNBC & cohort$PAM50_Basal]),
             sum(cohort$OS_event)),
  stringsAsFactors = FALSE
)
cohort_summary_file <- file.path(results_dir, paste0(prefix, "_cohort_summary.tsv"))
write.table(cohort_summary, cohort_summary_file, sep = "\t", quote = FALSE, row.names = FALSE)

cut_data <- data.frame(time = cohort$OS_time_days, event = cohort$OS_event)
for (gene in genes) cut_data[[gene]] <- cohort[[paste0(gene, "_log2TPM1")]]
cutpoint_object <- survminer::surv_cutpoint(
  cut_data, time = "time", event = "event", variables = genes,
  minprop = minprop, progressbar = FALSE
)
categorized <- survminer::surv_categorize(
  cutpoint_object, variables = genes, labels = c("Low", "High")
)

group_rows <- vector("list", length(genes))
result_rows <- vector("list", length(genes))
for (i in seq_along(genes)) {
  gene <- genes[i]
  group <- factor(categorized[[gene]], levels = c("Low", "High"))
  dat <- data.frame(time = cut_data$time, event = cut_data$event,
                    expression = cut_data[[gene]], group = group)
  cutoff <- unname(cutpoint_object$cutpoint[gene, "cutpoint"])
  fit <- survival::coxph(survival::Surv(time, event) ~ group, data = dat,
                         ties = "efron", x = TRUE)
  cox <- summary(fit)
  lr <- survival::survdiff(survival::Surv(time, event) ~ group, data = dat)
  set.seed(base_seed + i)
  corrected <- maxstat::maxstat.test(
    survival::Surv(time, event) ~ expression, data = dat,
    smethod = "LogRank", pmethod = "condMC", B = condMC_B,
    minprop = minprop, maxprop = 1 - minprop
  )
  group_rows[[i]] <- data.frame(
    patient_id = cohort$patient_id, gene = gene,
    expression_log2TPM1 = dat$expression, cutoff_log2TPM1 = cutoff,
    group = as.character(group), OS_time_days = dat$time,
    OS_event = dat$event, stringsAsFactors = FALSE
  )
  result_rows[[i]] <- data.frame(
    gene = gene, cohort = "TCGA-BRCA strict TNBC OR PAM50 Basal",
    n = nrow(dat), deaths = sum(dat$event),
    min_group_fraction = minprop,
    cutoff_log2TPM1 = cutoff, cutoff_TPM = 2^cutoff - 1,
    low_n = sum(group == "Low"), low_deaths = sum(dat$event[group == "Low"]),
    high_n = sum(group == "High"), high_deaths = sum(dat$event[group == "High"]),
    high_vs_low_hr = unname(cox$conf.int["groupHigh", "exp(coef)"]),
    ci_low = unname(cox$conf.int["groupHigh", "lower .95"]),
    ci_high = unname(cox$conf.int["groupHigh", "upper .95"]),
    cox_wald_p = unname(cox$coefficients["groupHigh", "Pr(>|z|)"]),
    selected_cutpoint_logrank_p = pchisq(lr$chisq, df = 1, lower.tail = FALSE),
    maxstat_condMC_adjusted_p = unname(corrected$p.value),
    maxstat_condMC_B = condMC_B, maxstat_condMC_seed = base_seed + i,
    cutoff_reproduced_by_condMC = isTRUE(all.equal(cutoff, unname(corrected$estimate))),
    stringsAsFactors = FALSE
  )
}
results <- do.call(rbind, result_rows)
results$BH_FDR_on_maxstat_adjusted_p <- p.adjust(results$maxstat_condMC_adjusted_p,
                                                method = "BH")
groups <- do.call(rbind, group_rows)
results_file <- file.path(results_dir, paste0(prefix, "_five_gene_results.tsv"))
groups_file <- file.path(results_dir, paste0(prefix, "_five_gene_groups.tsv"))
write.table(results, results_file, sep = "\t", quote = FALSE, row.names = FALSE)
write.table(groups, groups_file, sep = "\t", quote = FALSE, row.names = FALSE)
saveRDS(cutpoint_object, file.path(results_dir, paste0(prefix, "_surv_cutpoint.rds")))
writeLines(capture.output(sessionInfo()),
           file.path(results_dir, paste0(prefix, "_sessionInfo.txt")))

# Faithful compact four-column composition: Gene | HR (95% CI) | forest | P.
# Arial body/axis text is 10 pt; all plot geometry is kept inside the canvas.
forest <- results
forest$row <- rev(seq_along(genes))
forest$hr_ci_label <- sprintf("%.2f (%.2f–%.2f)", forest$high_vs_low_hr,
                              forest$ci_low, forest$ci_high)
forest$p_label <- ifelse(
  forest$selected_cutpoint_logrank_p < 0.001,
  formatC(forest$selected_cutpoint_logrank_p, format = "e", digits = 2),
  formatC(forest$selected_cutpoint_logrank_p, format = "f", digits = 4)
)
pt_to_mm <- 25.4 / 72.27
body_size <- 10 * pt_to_mm
header_size <- 10.2 * pt_to_mm
width_mm <- 150
height_mm <- 70
gene_x <- 0.032
hr_text_x <- 0.205
forest_left <- 0.535
forest_right <- 0.842
p_x <- 0.982
top_rule_y <- 6.00
header_y <- 5.66
header_rule_y <- 5.35
bottom_rule_y <- 0.72
tick_bottom_y <- 0.57
tick_label_y <- 0.30

# Expand the log-HR scale if the union produces a wider interval than Basal.
candidate_bounds <- c(0.1, 0.2, 0.25, 0.5, 1, 2, 5, 10, 20, 50, 100)
hr_min <- max(candidate_bounds[candidate_bounds <= min(c(0.5, forest$ci_low))])
hr_max <- min(candidate_bounds[candidate_bounds >= max(c(10, forest$ci_high))])
if (!is.finite(hr_min) || !is.finite(hr_max) || hr_min >= hr_max) {
  stop("Cannot determine a finite forest scale.")
}
map_hr <- function(x) forest_left +
  (log(x) - log(hr_min)) / (log(hr_max) - log(hr_min)) *
  (forest_right - forest_left)
ticks <- if (hr_max / hr_min >= 100) {
  # The broad CI range needs sparser labels to keep 10-pt Arial ticks apart.
  c(0.25, 1, 5, 20, 50)
} else {
  c(0.25, 0.5, 1, 2, 5, 10, 20)
}
ticks <- ticks[ticks >= hr_min & ticks <= hr_max]
forest$plot_hr <- map_hr(forest$high_vs_low_hr)
forest$plot_ci_low <- map_hr(forest$ci_low)
forest$plot_ci_high <- map_hr(forest$ci_high)
ink <- "#242424"
blue <- "#2E79B5"

plot <- ggplot(forest, aes(y = row)) +
  annotate("segment", x = 0.018, xend = 0.992,
           y = top_rule_y, yend = top_rule_y, linewidth = 0.52, colour = ink) +
  annotate("segment", x = 0.018, xend = 0.992,
           y = header_rule_y, yend = header_rule_y, linewidth = 0.45, colour = ink) +
  annotate("segment", x = 0.018, xend = 0.992,
           y = bottom_rule_y, yend = bottom_rule_y, linewidth = 0.52, colour = ink) +
  annotate("segment", x = map_hr(1), xend = map_hr(1),
           y = bottom_rule_y, yend = header_rule_y,
           linewidth = 0.40, linetype = "dashed", colour = "#555555") +
  geom_segment(aes(x = plot_ci_low, xend = plot_ci_high, yend = row),
               linewidth = 0.72, lineend = "butt", colour = ink) +
  geom_point(aes(x = plot_hr), shape = 17, size = 2.45, colour = blue) +
  geom_text(aes(x = gene_x, label = gene), hjust = 0, family = "Arial",
            size = body_size, colour = ink) +
  geom_text(aes(x = hr_text_x, label = hr_ci_label), hjust = 0,
            family = "Arial", size = body_size, colour = ink) +
  geom_text(aes(x = p_x, label = p_label), hjust = 1, family = "Arial",
            size = body_size, colour = ink) +
  annotate("text", x = gene_x, y = header_y, label = "Gene",
           family = "Arial", size = header_size, hjust = 0, colour = ink) +
  annotate("text", x = 0.355, y = header_y, label = "HR (95% CI)",
           family = "Arial", size = header_size, hjust = 0.5, colour = ink) +
  annotate("text", x = p_x, y = header_y, label = "P value",
           family = "Arial", fontface = "italic", size = header_size,
           hjust = 1, colour = ink) +
  annotate("segment", x = map_hr(ticks), xend = map_hr(ticks),
           y = bottom_rule_y, yend = tick_bottom_y,
           linewidth = 0.38, colour = ink) +
  annotate("text", x = map_hr(ticks), y = tick_label_y,
           label = as.character(ticks), family = "Arial",
           size = body_size, colour = ink) +
  scale_x_continuous(limits = c(0, 1), expand = expansion(mult = c(0, 0))) +
  scale_y_continuous(limits = c(0.08, 6.12), expand = expansion(mult = c(0, 0))) +
  coord_cartesian(clip = "off") +
  theme_void(base_family = "Arial") +
  theme(plot.margin = margin(1.6, 2, 1.6, 2, unit = "mm"),
        plot.background = element_rect(fill = "white", colour = NA),
        panel.background = element_rect(fill = "white", colour = NA))

figure_stem <- file.path(figures_dir, paste0(prefix, "_five_gene_OS_forest"))
width_in <- width_mm / 25.4
height_in <- height_mm / 25.4
svglite::svglite(paste0(figure_stem, ".svg"), width = width_in, height = height_in)
print(plot)
dev.off()
grDevices::cairo_pdf(paste0(figure_stem, ".pdf"), width = width_in,
                     height = height_in, family = "Arial")
print(plot)
dev.off()
ragg::agg_tiff(paste0(figure_stem, ".tiff"), width = width_in,
               height = height_in, units = "in", res = 600, compression = "lzw")
print(plot)
dev.off()
ragg::agg_png(paste0(figure_stem, ".png"), width = width_in,
              height = height_in, units = "in", res = 300)
print(plot)
dev.off()

source_data_file <- file.path(results_dir, paste0(prefix, "_forest_source_data.tsv"))
write.table(
  forest[, c("gene", "n", "deaths", "low_n", "low_deaths", "high_n", "high_deaths",
             "high_vs_low_hr", "ci_low", "ci_high", "selected_cutpoint_logrank_p",
             "maxstat_condMC_adjusted_p")],
  source_data_file, sep = "\t", quote = FALSE, row.names = FALSE
)
figure_files <- paste0(figure_stem, c(".svg", ".pdf", ".tiff", ".png"))
svg_text <- paste(readLines(paste0(figure_stem, ".svg"), warn = FALSE), collapse = "\n")
qa <- data.frame(
  check = c("patient IDs unique", "membership union correct", "overlap counted once",
            "fixed OS source", "five genes and total N retained", "group minimum proportion",
            "each split preserves 207 patients", "cutpoint independently reproduced",
            "finite HR and CI", "valid nominal and adjusted P", "all figures exist",
            "Arial in SVG", "minimum visible font at least 10 pt"),
  pass = c(!anyDuplicated(cohort$patient_id),
           all(cohort$TNBC | cohort$PAM50_Basal),
           sum(cohort$TNBC & cohort$PAM50_Basal) == length(overlap_ids),
           sum(cohort$OS_event) == 36L,
           identical(results$gene, genes) && all(results$n == 207L),
           # maxstat indexes the admissible range with floor(n * minprop).
           all(results$low_n >= floor(minprop * 207L)) &&
             all(results$high_n >= floor(minprop * 207L)),
           all(results$low_n + results$high_n == 207L),
           all(results$cutoff_reproduced_by_condMC),
           all(is.finite(results$high_vs_low_hr)) &&
             all(is.finite(results$ci_low)) && all(is.finite(results$ci_high)) &&
             all(results$ci_low < results$ci_high),
           all(results$selected_cutpoint_logrank_p >= 0 &
                 results$selected_cutpoint_logrank_p <= 1) &&
             all(results$maxstat_condMC_adjusted_p >= 0 &
                   results$maxstat_condMC_adjusted_p <= 1),
           all(file.exists(figure_files)) && all(file.info(figure_files)$size > 1000),
           grepl("Arial", svg_text, fixed = TRUE),
           body_size / pt_to_mm >= 10),
  stringsAsFactors = FALSE
)
qa_file <- file.path(results_dir, paste0(prefix, "_validation.tsv"))
write.table(qa, qa_file, sep = "\t", quote = FALSE, row.names = FALSE)
if (!all(qa$pass)) {
  print(qa[!qa$pass, ], row.names = FALSE)
  stop("Union analysis validation failed.")
}

cat("Cohort composition and events:\n")
print(cohort_summary, row.names = FALSE)
cat("Five-gene selected-cutpoint OS results:\n")
print(results[, c("gene", "low_n", "high_n", "high_vs_low_hr", "ci_low",
                  "ci_high", "selected_cutpoint_logrank_p",
                  "maxstat_condMC_adjusted_p", "BH_FDR_on_maxstat_adjusted_p")],
      row.names = FALSE, digits = 5)
