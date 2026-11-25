#!/usr/bin/env Rscript

# TCGA-BRCA full-cohort OS: five gene-specific survminer cutpoints, using
# its default 10% minimum group fraction. The displayed log-rank P is nominal
# after cutpoint selection; the companion table includes maxstat condMC P.
# This creates new Figure2A files and never overwrites the historical results.

options(stringsAsFactors = FALSE)
script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (length(script_arg) != 1L) stop("Run with Rscript.")
script_path <- normalizePath(sub("^--file=", "", script_arg), mustWork = TRUE)
project_dir <- normalizePath(file.path(dirname(script_path), ".."), mustWork = TRUE)
figure2a_dir <- normalizePath(file.path(project_dir, ".."), mustWork = TRUE)
renv_library <- file.path(figure2a_dir, "renv", "library", "macos", "R-4.5",
                          "aarch64-apple-darwin20.0.0")
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
minprop <- 0.10
condMC_B <- 99999L
base_seed <- 20261024L
prefix <- "Figure2A_TCGA_BRCA_full_minprop10"
old_file <- file.path(project_dir, "data", "processed",
                      "TCGA_BRCA_full_five_gene_OS_analysis_ready.tsv")
union_file <- file.path(figure2a_dir, "TCGA_TNBC_or_PAM50_Basal_five_gene_OS",
                        "data", "processed",
                        "Figure2A_TNBC_or_PAM50_Basal_minprop10_patient_level_union.tsv")
if (!all(file.exists(c(old_file, union_file)))) stop("Source data missing.")
data_dir <- file.path(project_dir, "data", "processed")
results_dir <- file.path(project_dir, "results")
figures_dir <- file.path(project_dir, "figures")
for (path in c(data_dir, results_dir, figures_dir))
  dir.create(path, recursive = TRUE, showWarnings = FALSE)

cohort <- read.delim(old_file, check.names = FALSE)
union <- read.delim(union_file, check.names = FALSE)
names(cohort)[names(cohort) == "sample_id"] <- "patient_id"
expr_fields <- as.vector(rbind(paste0(genes, "_TPM"),
                               paste0(genes, "_log2TPM1")))
stopifnot(nrow(cohort) == 1073L, sum(cohort$OS_event) == 150L,
          !anyDuplicated(cohort$patient_id), !anyDuplicated(union$patient_id),
          all(cohort$OS_time_days > 0), all(cohort$OS_event %in% 0:1),
          all(expr_fields %in% names(cohort)), all(expr_fields %in% names(union)))
shared <- intersect(cohort$patient_id, union$patient_id)
full_i <- match(shared, cohort$patient_id)
union_i <- match(shared, union$patient_id)
stopifnot(length(shared) == 204L,
          identical(cohort$OS_time_days[full_i], union$OS_time_days[union_i]),
          identical(cohort$OS_event[full_i], union$OS_event[union_i]))

# The old full-cohort file used a different primary-tumor aliquot for exactly
# one shared patient. Align that patient's expression to the previously used
# TNBC/Basal selected-primary-sample source; retain all 1073 full-cohort cases.
discordant <- rep(FALSE, length(shared))
for (field in expr_fields) {
  a <- cohort[[field]][full_i]
  b <- union[[field]][union_i]
  discordant <- discordant | abs(a - b) > 1e-6
}
discordant_ids <- shared[discordant]
stopifnot(identical(discordant_ids, "TCGA-A7-A13E"))
audit <- data.frame(patient_id = discordant_ids,
                    selected_primary_sample_id = union$sample_id[union_i[discordant]],
                    old_full_OS_days = cohort$OS_time_days[full_i[discordant]],
                    OS_event = cohort$OS_event[full_i[discordant]])
for (field in expr_fields) {
  audit[[paste0("old_", field)]] <- cohort[[field]][full_i[discordant]]
  audit[[paste0("aligned_", field)]] <- union[[field]][union_i[discordant]]
  cohort[[field]][full_i[discordant]] <- union[[field]][union_i[discordant]]
}
for (field in expr_fields)
  stopifnot(all(abs(cohort[[field]][full_i] - union[[field]][union_i]) < 1e-6))
for (gene in genes) {
  tpm <- cohort[[paste0(gene, "_TPM")]]
  logged <- cohort[[paste0(gene, "_log2TPM1")]]
  stopifnot(all(is.finite(tpm) & tpm >= 0),
            isTRUE(all.equal(log2(tpm + 1), logged, tolerance = 1e-8)))
}
write.table(cohort, file.path(data_dir, paste0(prefix, "_analysis_ready.tsv")),
            sep = "\t", quote = FALSE, row.names = FALSE)
write.table(audit, file.path(results_dir, paste0(prefix, "_sample_alignment_audit.tsv")),
            sep = "\t", quote = FALSE, row.names = FALSE)

cut_data <- data.frame(time = cohort$OS_time_days, event = cohort$OS_event)
for (gene in genes) cut_data[[gene]] <- cohort[[paste0(gene, "_log2TPM1")]]
cutpoint_object <- survminer::surv_cutpoint(
  cut_data, time = "time", event = "event", variables = genes,
  minprop = minprop, progressbar = FALSE
)
categorized <- survminer::surv_categorize(
  cutpoint_object, variables = genes, labels = c("Low", "High")
)

result_rows <- vector("list", length(genes))
group_rows <- vector("list", length(genes))
for (i in seq_along(genes)) {
  gene <- genes[i]
  dat <- data.frame(time = cut_data$time, event = cut_data$event,
                    expression = cut_data[[gene]],
                    group = factor(categorized[[gene]], levels = c("Low", "High")))
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
    patient_id = cohort$patient_id, gene = gene, expression_log2TPM1 = dat$expression,
    cutoff_log2TPM1 = cutoff, group = as.character(dat$group),
    OS_time_days = dat$time, OS_event = dat$event
  )
  result_rows[[i]] <- data.frame(
    gene = gene, cohort = "TCGA-BRCA full primary-tumor cohort",
    n = nrow(dat), deaths = sum(dat$event), min_group_fraction = minprop,
    cutoff_log2TPM1 = cutoff, cutoff_TPM = 2^cutoff - 1,
    low_n = sum(dat$group == "Low"), low_deaths = sum(dat$event[dat$group == "Low"]),
    high_n = sum(dat$group == "High"), high_deaths = sum(dat$event[dat$group == "High"]),
    high_vs_low_hr = unname(cox$conf.int["groupHigh", "exp(coef)"]),
    ci_low = unname(cox$conf.int["groupHigh", "lower .95"]),
    ci_high = unname(cox$conf.int["groupHigh", "upper .95"]),
    cox_wald_p = unname(cox$coefficients["groupHigh", "Pr(>|z|)"]),
    selected_cutpoint_logrank_p = pchisq(lr$chisq, 1, lower.tail = FALSE),
    maxstat_condMC_adjusted_p = unname(corrected$p.value),
    maxstat_condMC_B = condMC_B, maxstat_condMC_seed = base_seed + i,
    cutoff_reproduced_by_condMC = isTRUE(all.equal(cutoff, unname(corrected$estimate)))
  )
}
results <- do.call(rbind, result_rows)
results$BH_FDR_on_maxstat_adjusted_p <- p.adjust(results$maxstat_condMC_adjusted_p,
                                                method = "BH")
groups <- do.call(rbind, group_rows)
write.table(results, file.path(results_dir, paste0(prefix, "_five_gene_results.tsv")),
            sep = "\t", quote = FALSE, row.names = FALSE)
write.table(groups, file.path(results_dir, paste0(prefix, "_five_gene_groups.tsv")),
            sep = "\t", quote = FALSE, row.names = FALSE)
saveRDS(cutpoint_object,
        file.path(results_dir, paste0(prefix, "_surv_cutpoint.rds")))
writeLines(capture.output(sessionInfo()),
           file.path(results_dir, paste0(prefix, "_sessionInfo.txt")))

# Four-column reference style: gene | HR (95% CI) | forest | nominal P.
# All visible lettering is >=10 pt Arial; no panel label, title, footer, or FDR.
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

candidate_bounds <- c(0.1, 0.2, 0.25, 0.5, 1, 2, 5, 10, 20, 50, 100)
hr_min <- max(candidate_bounds[candidate_bounds <= min(c(0.5, forest$ci_low))])
hr_max <- min(candidate_bounds[candidate_bounds >= max(c(2, forest$ci_high))])
if (!is.finite(hr_min) || !is.finite(hr_max) || hr_min >= hr_max)
  stop("Cannot determine finite forest scale.")
map_hr <- function(x) forest_left +
  (log(x) - log(hr_min)) / (log(hr_max) - log(hr_min)) *
  (forest_right - forest_left)
tick_candidates <- c(0.25, 0.5, 1, 1.5, 2, 3, 5, 10, 20, 50)
ticks <- tick_candidates[tick_candidates >= hr_min & tick_candidates <= hr_max]
# Keep >=10 pt tick labels separated by at least 6 mm.
while (length(ticks) > 2L && any(diff(map_hr(ticks)) * width_mm < 6))
  ticks <- ticks[-which.min(diff(map_hr(ticks))) - 1L]
forest$plot_hr <- map_hr(forest$high_vs_low_hr)
forest$plot_ci_low <- map_hr(forest$ci_low)
forest$plot_ci_high <- map_hr(forest$ci_high)
ink <- "#242424"
blue <- "#2E79B5"

plot <- ggplot(forest, aes(y = row)) +
  annotate("segment", x = 0.018, xend = 0.992, y = top_rule_y,
           yend = top_rule_y, linewidth = 0.52, colour = ink) +
  annotate("segment", x = 0.018, xend = 0.992, y = header_rule_y,
           yend = header_rule_y, linewidth = 0.45, colour = ink) +
  annotate("segment", x = 0.018, xend = 0.992, y = bottom_rule_y,
           yend = bottom_rule_y, linewidth = 0.52, colour = ink) +
  annotate("segment", x = map_hr(1), xend = map_hr(1), y = bottom_rule_y,
           yend = header_rule_y, linewidth = 0.40, linetype = "dashed",
           colour = "#555555") +
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

write.table(
  forest[, c("gene", "n", "deaths", "low_n", "low_deaths", "high_n",
             "high_deaths", "high_vs_low_hr", "ci_low", "ci_high",
             "selected_cutpoint_logrank_p", "maxstat_condMC_adjusted_p")],
  file.path(results_dir, paste0(prefix, "_forest_source_data.tsv")),
  sep = "\t", quote = FALSE, row.names = FALSE
)
figure_files <- paste0(figure_stem, c(".svg", ".pdf", ".tiff", ".png"))
svg_text <- paste(readLines(paste0(figure_stem, ".svg"), warn = FALSE),
                  collapse = "\n")
qa <- data.frame(
  check = c("unique patient IDs", "full-cohort N and deaths",
            "shared patient outcomes consistent", "selected primary sample aligned",
            "five genes and group totals", "10% group floor",
            "cutpoint independently reproduced", "finite HR and confidence intervals",
            "valid P values", "four figures exist", "Arial in SVG", "text at least 10 pt"),
  pass = c(!anyDuplicated(cohort$patient_id),
           nrow(cohort) == 1073L && sum(cohort$OS_event) == 150L,
           length(shared) == 204L,
           identical(discordant_ids, "TCGA-A7-A13E"),
           identical(results$gene, genes) && all(results$low_n + results$high_n == 1073L),
           all(results$low_n >= floor(minprop * 1073L)) &&
             all(results$high_n >= floor(minprop * 1073L)),
           all(results$cutoff_reproduced_by_condMC),
           all(is.finite(results$high_vs_low_hr)) &&
             all(is.finite(results$ci_low)) && all(is.finite(results$ci_high)) &&
             all(results$ci_low < results$ci_high),
           all(results$selected_cutpoint_logrank_p >= 0 &
                 results$selected_cutpoint_logrank_p <= 1) &&
             all(results$maxstat_condMC_adjusted_p >= 0 &
                   results$maxstat_condMC_adjusted_p <= 1),
           all(file.exists(figure_files)) && all(file.info(figure_files)$size > 1000),
           grepl("Arial", svg_text, fixed = TRUE), body_size / pt_to_mm >= 10)
)
write.table(qa, file.path(results_dir, paste0(prefix, "_validation.tsv")),
            sep = "\t", quote = FALSE, row.names = FALSE)
if (!all(qa$pass)) stop("Full-cohort analysis validation failed.")
cat("TCGA-BRCA full-cohort results (n=1073, deaths=150):\n")
print(results[, c("gene", "low_n", "high_n", "high_vs_low_hr", "ci_low",
                  "ci_high", "selected_cutpoint_logrank_p",
                  "maxstat_condMC_adjusted_p", "BH_FDR_on_maxstat_adjusted_p")],
      row.names = FALSE, digits = 5)
