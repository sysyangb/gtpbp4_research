#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)
suppressPackageStartupMessages({
  library(survival)
  library(maxstat)
  library(ggplot2)
  library(patchwork)
  library(survminer)
  library(svglite)
  library(ragg)
})

genes <- c("DKC1", "GTPBP4", "ILF2", "NOP2", "RPP38")
minprop <- 0.20
maxprop <- 0.80
condMC_B <- 99999L
base_seed <- 20260904L

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (length(script_arg) != 1L) stop("Run this file with Rscript.")
script_path <- normalizePath(sub("^--file=", "", script_arg), mustWork = TRUE)
project_dir <- normalizePath(file.path(dirname(script_path), ".."), mustWork = TRUE)
figure_dir <- normalizePath(file.path(project_dir, ".."), mustWork = TRUE)
data_dir <- file.path(project_dir, "data", "processed")
results_dir <- file.path(project_dir, "results")
figures_dir <- file.path(project_dir, "figures")
dir.create(data_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)

expression_file <- file.path(figure_dir, "exprTrain.txt")
survival_file <- file.path(figure_dir, "status.txt")
stopifnot(file.exists(expression_file), file.exists(survival_file))

expression_table <- read.delim(expression_file, check.names = FALSE)
survival_table <- read.delim(survival_file, check.names = FALSE)
stopifnot(
  identical(names(expression_table)[1], "geneSymbol"),
  identical(names(survival_table), c("status", "tme")),
  ncol(expression_table) - 1L == nrow(survival_table),
  !anyDuplicated(expression_table$geneSymbol),
  !anyDuplicated(names(expression_table)[-1L]),
  all(survival_table$status %in% c(0L, 1L)),
  all(is.finite(survival_table$tme)),
  all(survival_table$tme >= 0)
)
missing_genes <- setdiff(genes, expression_table$geneSymbol)
if (length(missing_genes)) stop("Missing genes: ", paste(missing_genes, collapse = ", "))

expression_raw <- as.matrix(expression_table[match(genes, expression_table$geneSymbol), -1L, drop = FALSE])
storage.mode(expression_raw) <- "double"
rownames(expression_raw) <- genes
colnames(expression_raw) <- names(expression_table)[-1L]
if (any(!is.finite(expression_raw)) || any(expression_raw < 0)) stop("Invalid TPM expression values")

# Primary analysis excludes OS time = 0 because these records have no observable follow-up.
positive_followup <- survival_table$tme > 0
sample_ids <- colnames(expression_raw)[positive_followup]
status <- as.integer(survival_table$status[positive_followup])
time_days <- as.double(survival_table$tme[positive_followup])
expression_raw <- expression_raw[, positive_followup, drop = FALSE]
expression_log2 <- log2(expression_raw + 1)

analysis_data <- data.frame(
  sample_id = sample_ids,
  OS_time_days = time_days,
  OS_time_months = time_days / 30.4375,
  OS_event = status,
  stringsAsFactors = FALSE,
  check.names = FALSE
)
for (gene in genes) {
  analysis_data[[paste0(gene, "_TPM")]] <- expression_raw[gene, ]
  analysis_data[[paste0(gene, "_log2TPM1")]] <- expression_log2[gene, ]
}
write.table(
  analysis_data,
  file.path(data_dir, "TCGA_BRCA_full_five_gene_OS_analysis_ready.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

extract_cox <- function(fit, term) {
  sm <- summary(fit)
  c(
    hr = unname(sm$conf.int[term, "exp(coef)"]),
    ci_low = unname(sm$conf.int[term, "lower .95"]),
    ci_high = unname(sm$conf.int[term, "upper .95"]),
    wald_p = unname(sm$coefficients[term, "Pr(>|z|)"]),
    likelihood_ratio_p = unname(sm$logtest["pvalue"]),
    score_p = unname(sm$sctest["pvalue"])
  )
}

median_objects <- list()
single_results <- do.call(rbind, lapply(genes, function(gene) {
  dat <- data.frame(
    time = time_days,
    event = status,
    expression = expression_log2[gene, ],
    stringsAsFactors = FALSE
  )
  dat$expression_z <- as.numeric(scale(dat$expression))
  dat$group <- factor(ifelse(dat$expression > median(dat$expression), "High", "Low"), levels = c("Low", "High"))

  continuous_fit <- coxph(Surv(time, event) ~ expression_z, data = dat, ties = "efron", x = TRUE)
  continuous <- extract_cox(continuous_fit, "expression_z")
  ph_p <- unname(cox.zph(continuous_fit, transform = "km")$table["expression_z", "p"])
  group_fit <- coxph(Surv(time, event) ~ group, data = dat, ties = "efron", x = TRUE)
  group <- extract_cox(group_fit, "groupHigh")
  logrank <- survdiff(Surv(time, event) ~ group, data = dat)

  median_objects[[gene]] <<- list(data = dat, fit = survfit(Surv(time, event) ~ group, data = dat))
  data.frame(
    gene = gene,
    n = nrow(dat),
    events = sum(dat$event),
    expression_scale = "log2(TPM+1)",
    continuous_hr_per_1sd = continuous["hr"],
    continuous_ci_low = continuous["ci_low"],
    continuous_ci_high = continuous["ci_high"],
    continuous_wald_p = continuous["wald_p"],
    continuous_likelihood_ratio_p = continuous["likelihood_ratio_p"],
    proportional_hazards_p = ph_p,
    median_cutoff_log2TPM1 = median(dat$expression),
    median_cutoff_TPM = 2^median(dat$expression) - 1,
    low_n = sum(dat$group == "Low"),
    low_events = sum(dat$event[dat$group == "Low"]),
    high_n = sum(dat$group == "High"),
    high_events = sum(dat$event[dat$group == "High"]),
    median_high_vs_low_hr = group["hr"],
    median_high_vs_low_ci_low = group["ci_low"],
    median_high_vs_low_ci_high = group["ci_high"],
    median_high_vs_low_cox_p = group["wald_p"],
    median_split_logrank_p = 1 - pchisq(logrank$chisq, df = 1),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
}))
single_results$continuous_fdr_bh <- p.adjust(single_results$continuous_wald_p, method = "BH")
single_results$median_logrank_fdr_bh <- p.adjust(single_results$median_split_logrank_p, method = "BH")
single_results$continuous_significant_fdr_lt_0_05 <- single_results$continuous_fdr_bh < 0.05
single_results$median_significant_fdr_lt_0_05 <- single_results$median_logrank_fdr_bh < 0.05
write.table(single_results, file.path(results_dir, "single_gene_OS_results.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

cat("Running maxstat best-cutoff analyses (99,999 conditional Monte Carlo replicates each)...\n")
best_objects <- list()
best_results <- do.call(rbind, lapply(seq_along(genes), function(i) {
  gene <- genes[i]
  cat(sprintf("  %d/%d %s\n", i, length(genes), gene))
  dat <- data.frame(time = time_days, event = status, expression = expression_log2[gene, ])
  n <- nrow(dat)
  min_low_n <- ceiling(minprop * n)
  max_low_n <- floor(maxprop * n)
  epsilon <- 1e-8 / n
  set.seed(base_seed + i)
  maxstat_fit <- maxstat.test(
    Surv(time, event) ~ expression,
    data = dat,
    smethod = "LogRank",
    pmethod = "condMC",
    B = condMC_B,
    minprop = min_low_n / n + epsilon,
    maxprop = max_low_n / n + epsilon
  )
  cutoff <- unname(maxstat_fit$estimate)
  dat$group <- factor(ifelse(dat$expression <= cutoff, "Low", "High"), levels = c("Low", "High"))
  group_fit <- coxph(Surv(time, event) ~ group, data = dat, ties = "efron", x = TRUE)
  group <- extract_cox(group_fit, "groupHigh")
  logrank <- survdiff(Surv(time, event) ~ group, data = dat)
  best_objects[[gene]] <<- list(data = dat, fit = survfit(Surv(time, event) ~ group, data = dat))

  data.frame(
    gene = gene,
    n = n,
    events = sum(dat$event),
    cutoff_log2TPM1 = cutoff,
    cutoff_TPM = 2^cutoff - 1,
    low_n = sum(dat$group == "Low"),
    low_events = sum(dat$event[dat$group == "Low"]),
    high_n = sum(dat$group == "High"),
    high_events = sum(dat$event[dat$group == "High"]),
    high_vs_low_hr = group["hr"],
    high_vs_low_ci_low = group["ci_low"],
    high_vs_low_ci_high = group["ci_high"],
    cox_wald_p_at_selected_cutoff = group["wald_p"],
    naive_best_logrank_p = 1 - pchisq(logrank$chisq, df = 1),
    maxstat_condMC_adjusted_p = unname(maxstat_fit$p.value),
    maxstat_standardized_statistic = unname(maxstat_fit$statistic),
    condMC_B = condMC_B,
    condMC_seed = base_seed + i,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
}))
best_results$maxstat_adjusted_fdr_bh <- p.adjust(best_results$maxstat_condMC_adjusted_p, method = "BH")
best_results$naive_fdr_bh <- p.adjust(best_results$naive_best_logrank_p, method = "BH")
best_results$maxstat_significant_fdr_lt_0_05 <- best_results$maxstat_adjusted_fdr_bh < 0.05
write.table(best_results, file.path(results_dir, "best_cutoff_OS_results.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

best_groups <- do.call(rbind, lapply(genes, function(gene) {
  data.frame(
    sample_id = sample_ids,
    gene = gene,
    OS_time_days = time_days,
    OS_event = status,
    expression_log2TPM1 = best_objects[[gene]]$data$expression,
    group = best_objects[[gene]]$data$group,
    stringsAsFactors = FALSE
  )
}))
write.table(best_groups, file.path(results_dir, "best_cutoff_sample_groups.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

# Sensitivity analysis retaining zero-time records to show their impact explicitly.
zero_time_sensitivity <- do.call(rbind, lapply(genes, function(gene) {
  x <- log2(as.numeric(expression_table[expression_table$geneSymbol == gene, -1L]) + 1)
  make_row <- function(keep, rule) {
    dat <- data.frame(time = survival_table$tme[keep], event = survival_table$status[keep], expression = x[keep])
    dat$z <- as.numeric(scale(dat$expression))
    dat$group <- factor(ifelse(dat$expression > median(dat$expression), "High", "Low"), levels = c("Low", "High"))
    cf <- coxph(Surv(time, event) ~ z, data = dat, ties = "efron")
    gf <- coxph(Surv(time, event) ~ group, data = dat, ties = "efron")
    lr <- survdiff(Surv(time, event) ~ group, data = dat)
    data.frame(
      gene = gene, rule = rule, n = nrow(dat), events = sum(dat$event),
      continuous_hr_per_1sd = exp(coef(cf)), continuous_wald_p = summary(cf)$coefficients["z", "Pr(>|z|)"],
      median_high_vs_low_hr = exp(coef(gf)), median_logrank_p = 1 - pchisq(lr$chisq, 1),
      stringsAsFactors = FALSE
    )
  }
  rbind(
    make_row(survival_table$tme > 0, "Primary: OS time > 0"),
    make_row(survival_table$tme >= 0, "Sensitivity: include OS time = 0")
  )
}))
write.table(zero_time_sensitivity, file.path(results_dir, "zero_time_sensitivity.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

palette <- c(Low = "#2878B5", High = "#D9534F")
theme_pub <- function(base_size = 7) {
  theme_classic(base_size = base_size, base_family = "Arial") +
    theme(
      axis.line = element_line(linewidth = 0.35, colour = "black"),
      axis.ticks = element_line(linewidth = 0.35, colour = "black"),
      axis.text = element_text(colour = "black"),
      legend.title = element_blank(),
      legend.position = "top",
      legend.key.width = grid::unit(8, "mm"),
      plot.title = element_text(face = "bold", size = base_size + 0.7),
      plot.subtitle = element_text(size = base_size - 0.3, colour = "#444444"),
      plot.tag = element_text(face = "bold", size = 8.5)
    )
}

save_bundle <- function(plot, stem, width_mm, height_mm) {
  width_in <- width_mm / 25.4
  height_in <- height_mm / 25.4
  svglite::svglite(paste0(stem, ".svg"), width = width_in, height = height_in)
  print(plot); dev.off()
  grDevices::cairo_pdf(paste0(stem, ".pdf"), width = width_in, height = height_in, family = "Arial")
  print(plot); dev.off()
  ragg::agg_tiff(paste0(stem, ".tiff"), width = width_in, height = height_in, units = "in", res = 600, compression = "lzw")
  print(plot); dev.off()
  ragg::agg_png(paste0(stem, ".png"), width = width_in, height = height_in, units = "in", res = 300)
  print(plot); dev.off()
}

km_panel <- function(obj, gene, annotation, subtitle) {
  curve <- survminer::surv_summary(obj$fit, data = obj$data)
  curve$group <- factor(sub("^group=", "", curve$strata), levels = c("Low", "High"))
  ggplot(curve, aes(x = time / 30.4375, y = surv, colour = group)) +
    geom_step(linewidth = 0.65, direction = "hv") +
    geom_point(data = subset(curve, n.censor > 0), shape = 3, size = 0.6, stroke = 0.45) +
    annotate("label", x = Inf, y = 0.05, label = annotation, hjust = 1.03, vjust = 0,
             size = 2.1, family = "Arial", linewidth = 0.2, label.padding = grid::unit(1.1, "mm")) +
    scale_colour_manual(values = palette) +
    scale_y_continuous(limits = c(0, 1.02), breaks = c(0, 0.25, 0.5, 0.75, 1), expand = expansion(mult = c(0, 0.01))) +
    scale_x_continuous(expand = expansion(mult = c(0, 0.03))) +
    labs(x = "Overall survival (months)", y = "Survival probability", title = gene, subtitle = subtitle) +
    theme_pub(6.8)
}

median_panels <- lapply(genes, function(gene) {
  row <- single_results[single_results$gene == gene, ]
  annotation <- sprintf(
    "HR %.2f (%.2f–%.2f)\nlog-rank P = %.3g\nFDR = %.3g",
    row$median_high_vs_low_hr, row$median_high_vs_low_ci_low, row$median_high_vs_low_ci_high,
    row$median_split_logrank_p, row$median_logrank_fdr_bh
  )
  km_panel(median_objects[[gene]], gene, annotation, sprintf("Median split: Low n=%d; High n=%d", row$low_n, row$high_n))
})
median_figure <- wrap_plots(median_panels, ncol = 3, guides = "collect") +
  plot_annotation(
    title = "TCGA-BRCA full-cohort OS by five prespecified genes",
    subtitle = "Primary analysis: gene-wise median split; n=1,073, deaths=150",
    caption = "Censoring is marked by +. FDR: Benjamini–Hochberg across five genes.",
    tag_levels = "a"
  ) & theme(legend.position = "top")
save_bundle(median_figure, file.path(figures_dir, "five_gene_median_KM"), 183, 132)

forest <- single_results
forest$gene <- factor(forest$gene, levels = rev(genes))
forest$significant <- factor(ifelse(forest$median_logrank_fdr_bh < 0.05, "FDR < 0.05", "FDR ≥ 0.05"), levels = c("FDR < 0.05", "FDR ≥ 0.05"))
forest_plot <- ggplot(forest, aes(x = median_high_vs_low_hr, y = gene, colour = significant)) +
  geom_vline(xintercept = 1, linetype = "dashed", linewidth = 0.35, colour = "#777777") +
  geom_segment(aes(x = median_high_vs_low_ci_low, xend = median_high_vs_low_ci_high, yend = gene), linewidth = 0.65) +
  geom_point(size = 2.2) +
  scale_x_log10(limits = c(0.55, 2.2), breaks = c(0.6, 0.8, 1, 1.5, 2)) +
  scale_colour_manual(values = c("FDR < 0.05" = "#D9534F", "FDR ≥ 0.05" = "#555555")) +
  labs(
    x = "High versus Low hazard ratio (95% CI)", y = NULL,
    title = "Only GTPBP4 remains significant after five-gene correction",
    subtitle = "TCGA-BRCA full cohort; median split; univariate Cox OS",
    caption = "n=1,073; deaths=150. Significance is based on BH-adjusted log-rank P."
  ) + theme_pub(7.4)
save_bundle(forest_plot, file.path(figures_dir, "five_gene_median_forest"), 120, 80)

best_panels <- lapply(genes, function(gene) {
  row <- best_results[best_results$gene == gene, ]
  annotation <- sprintf(
    "HR %.2f (%.2f–%.2f)\nnaive P = %.3g\nmaxstat P = %.3g\nFDR = %.3g",
    row$high_vs_low_hr, row$high_vs_low_ci_low, row$high_vs_low_ci_high,
    row$naive_best_logrank_p, row$maxstat_condMC_adjusted_p, row$maxstat_adjusted_fdr_bh
  )
  km_panel(best_objects[[gene]], gene, annotation, sprintf("Exploratory best cutoff: Low n=%d; High n=%d", row$low_n, row$high_n))
})
best_figure <- wrap_plots(best_panels, ncol = 3, guides = "collect") +
  plot_annotation(
    title = "Exploratory best-cutoff OS analysis",
    subtitle = "20%–80% group constraint; 99,999 conditional Monte Carlo replicates",
    caption = "Use maxstat-adjusted P for inference; naive log-rank P is optimistic. FDR is across five adjusted P values.",
    tag_levels = "a"
  ) & theme(legend.position = "top")
save_bundle(best_figure, file.path(figures_dir, "five_gene_best_cutoff_KM"), 183, 132)

median_curve_source <- do.call(rbind, lapply(genes, function(gene) {
  curve <- survminer::surv_summary(median_objects[[gene]]$fit, data = median_objects[[gene]]$data)
  curve$gene <- gene
  curve
}))
write.table(median_curve_source, file.path(results_dir, "median_KM_curve_source_data.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

best_curve_source <- do.call(rbind, lapply(genes, function(gene) {
  curve <- survminer::surv_summary(best_objects[[gene]]$fit, data = best_objects[[gene]]$data)
  curve$gene <- gene
  curve
}))
write.table(best_curve_source, file.path(results_dir, "best_cutoff_KM_curve_source_data.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

data_quality <- data.frame(
  check = c(
    "input_samples", "input_events", "zero_time_records_excluded", "zero_time_events_excluded",
    "primary_samples", "primary_events", "missing_expression_values", "duplicated_gene_symbols",
    "duplicated_sample_ids", "genes_found", "OS_event_definition"
  ),
  result = c(
    nrow(survival_table), sum(survival_table$status), sum(survival_table$tme == 0),
    sum(survival_table$status[survival_table$tme == 0]), length(time_days), sum(status),
    sum(is.na(expression_raw)), anyDuplicated(expression_table$geneSymbol), anyDuplicated(sample_ids),
    paste(genes, collapse = ", "), "1=death; 0=censored"
  )
)
write.table(data_quality, file.path(results_dir, "data_quality_checks.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

hashes <- tools::md5sum(c(expression_file, survival_file))
manifest <- c(
  paste0("analysis_time: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
  "cohort: TCGA-BRCA full primary-tumor cohort represented by exprTrain/status",
  "primary_eligibility: finite OS time > 0; binary status; complete target-gene expression",
  "endpoint: overall survival; status 1=death, 0=censored",
  "expression: log2(TPM+1)",
  "prespecified_primary_split: within-cohort median, High > median",
  "continuous_model: per 1 SD log2(TPM+1)",
  "best_cutoff: maxstat LogRank; 20%-80%; condMC B=99999",
  "multiple_testing: BH across five genes",
  paste0(basename(names(hashes)), "_md5: ", unname(hashes))
)
writeLines(manifest, file.path(results_dir, "analysis_manifest.txt"))
writeLines(capture.output(sessionInfo()), file.path(results_dir, "sessionInfo.txt"))

cat("Analysis complete.\n\nMedian split results:\n")
print(single_results[, c("gene", "median_high_vs_low_hr", "median_high_vs_low_ci_low", "median_high_vs_low_ci_high", "median_split_logrank_p", "median_logrank_fdr_bh")], row.names = FALSE)
cat("\nContinuous Cox results:\n")
print(single_results[, c("gene", "continuous_hr_per_1sd", "continuous_ci_low", "continuous_ci_high", "continuous_wald_p", "continuous_fdr_bh")], row.names = FALSE)
cat("\nBest cutoff results:\n")
print(best_results[, c("gene", "low_n", "high_n", "high_vs_low_hr", "naive_best_logrank_p", "maxstat_condMC_adjusted_p", "maxstat_adjusted_fdr_bh")], row.names = FALSE)
