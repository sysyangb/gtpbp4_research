#!/usr/bin/env Rscript

# Figure 2A: the five smallest TNBC-vs-normal DESeq2 FDR values among the
# established 36-gene AMF set, aligned with their TCGA-BRCA full-cohort OS
# estimates. The two evidence types use different patient sets by design.

options(stringsAsFactors = FALSE)
script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (length(script_arg) != 1L) stop("Run this file with Rscript.")
script_path <- normalizePath(sub("^--file=", "", script_arg), mustWork = TRUE)
output_dir <- normalizePath(file.path(dirname(script_path), ".."), mustWork = TRUE)
figure2_dir <- normalizePath(file.path(output_dir, ".."), mustWork = TRUE)
workspace_dir <- normalizePath(file.path(figure2_dir, "..", ".."), mustWork = TRUE)
renv_library <- file.path(figure2_dir, "renv", "library", "macos", "R-4.5",
                          "aarch64-apple-darwin20.0.0")
if (dir.exists(renv_library)) .libPaths(c(renv_library, .libPaths()))
suppressPackageStartupMessages({
  library(survival)
  library(survminer)
  library(maxstat)
  library(ggplot2)
  library(readxl)
  library(svglite)
  library(ragg)
})

prefix <- "Figure2A_TNBC_FDR_top5_BRCA_OS"
minprop <- 0.10
condMC_B <- 99999L
fig_dir <- file.path(output_dir, "figures")
result_dir <- file.path(output_dir, "results")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)

# Evidence contract: rank the 36 named genes using the original genome-wide
# adjusted DESeq2 P values. Do not read the AMF polysome FDR as a TNBC FDR.
amf_file <- file.path(workspace_dir, "outputs",
  "01a1773c-9bfe-70a2-bd0c-558728c6497d",
  "AMFs36_FOT_FDR_sorted_checked.xlsx")
de_file <- file.path(workspace_dir, "outputs",
  "TNBC_DESeq2_RBP_20261026_FC1p5", "results",
  "TN_vs_NORMAL_all_tested_genes.tsv.gz")
group_count_file <- file.path(workspace_dir, "outputs",
  "TNBC_DESeq2_RBP_20261026_FC1p5", "metadata",
  "sample_group_counts.tsv")
stopifnot(all(file.exists(c(amf_file, de_file, group_count_file))))
amf36 <- readxl::read_excel(amf_file, sheet = "AMFs_36_FDR_sorted")
gene_pool <- as.character(amf36$Gene)
stopifnot(length(gene_pool) == 36L, !anyNA(gene_pool),
          !anyDuplicated(gene_pool))
de <- read.delim(gzfile(de_file), check.names = FALSE)
de <- de[de$gene_symbol %in% gene_pool, , drop = FALSE]
stopifnot(nrow(de) == 36L, !anyDuplicated(de$gene_symbol),
          setequal(de$gene_symbol, gene_pool),
          all(is.finite(de$padj)), all(de$padj > 0),
          all(is.finite(de$log2FoldChange)),
          all(de$contrast == "TN_vs_NORMAL"))
de <- de[order(de$padj, match(de$gene_symbol, gene_pool)), , drop = FALSE]
top5 <- de[seq_len(5L), , drop = FALSE]
genes <- top5$gene_symbol
stopifnot(identical(genes,
                    c("ILF2", "DKC1", "GTPBP4", "GNL2", "NOP58")))
groups_de <- read.delim(group_count_file)
stopifnot(groups_de$N[match("TN", groups_de$condition)] == 177L,
          groups_de$N[match("NORMAL", groups_de$condition)] == 112L)

# Use the exact full-BRCA patient set and aliquot alignment from the prior
# Figure 2A analysis. Read all five gene rows from the original TPM matrix;
# the older 51-gene wide expression file had a shifted header index.
prior_dir <- file.path(figure2_dir, "TCGA_BRCA_full_five_gene_OS")
cohort_file <- file.path(prior_dir, "data", "processed",
  "Figure2A_TCGA_BRCA_full_minprop10_analysis_ready.tsv")
alignment_file <- file.path(prior_dir, "results",
  "Figure2A_TCGA_BRCA_full_minprop10_sample_alignment_audit.tsv")
map_file <- file.path(figure2_dir, "fifty_one_gene_survival",
                      "sample_to_tpm_aliquot.tsv")
id_file <- file.path(figure2_dir, "fifty_one_gene_survival",
                     "gene_list_51.tsv")
tpm_file <- file.path(figure2_dir, "TCGA-BRCA_TPM.txt")
stopifnot(all(file.exists(c(cohort_file, alignment_file, map_file,
                            id_file, tpm_file))))
cohort <- read.delim(cohort_file, check.names = FALSE)
alignment <- read.delim(alignment_file, check.names = FALSE)
sample_map <- read.delim(map_file, check.names = FALSE)
gene_ids <- read.delim(id_file, check.names = FALSE)
stopifnot(nrow(cohort) == 1073L, sum(cohort$OS_event) == 150L,
          !anyDuplicated(cohort$patient_id),
          all(cohort$OS_time_days > 0), all(cohort$OS_event %in% 0:1),
          !anyDuplicated(sample_map$patient_id),
          !anyDuplicated(gene_ids$gene_symbol),
          all(genes %in% gene_ids$gene_symbol), nrow(alignment) == 1L,
          identical(alignment$patient_id, "TCGA-A7-A13E"))
chosen_aliquot <- sample_map$tpm_aliquot[
  match(cohort$patient_id, sample_map$patient_id)]
chosen_aliquot[cohort$patient_id == alignment$patient_id] <-
  alignment$selected_primary_sample_id
if (anyNA(chosen_aliquot)) stop("A cohort patient's TPM aliquot is missing.")

strip_quotes <- function(x) sub('"$', "", sub('^"', "", x))
con <- file(tpm_file, open = "r")
header <- strip_quotes(strsplit(readLines(con, n = 1L), "\t", fixed = TRUE)[[1L]])
stopifnot(!anyDuplicated(header), all(chosen_aliquot %in% header))
aliquot_index <- match(chosen_aliquot, header)
id_to_gene <- setNames(genes,
  gene_ids$ensembl_id[match(genes, gene_ids$gene_symbol)])
stopifnot(!anyNA(names(id_to_gene)), !anyDuplicated(names(id_to_gene)))
expression <- matrix(NA_real_, nrow = length(genes), ncol = nrow(cohort),
                     dimnames = list(genes, cohort$patient_id))
source_ensembl <- setNames(rep(NA_character_, length(genes)), genes)
repeat {
  lines <- readLines(con, n = 500L, warn = FALSE)
  if (!length(lines)) break
  row_ids <- strip_quotes(sub("\t.*$", "", lines))
  row_keys <- sub("\\..*$", "", row_ids)
  hits <- which(row_keys %in% names(id_to_gene))
  for (j in hits) {
    fields <- strsplit(lines[j], "\t", fixed = TRUE)[[1L]]
    if (length(fields) != length(header) + 1L)
      stop("Unexpected TPM row width for ", row_ids[j])
    gene <- id_to_gene[[row_keys[j]]]
    if (!is.na(source_ensembl[gene]))
      stop("Duplicate TPM row for ", gene)
    expression[gene, ] <- as.double(fields[aliquot_index + 1L])
    source_ensembl[gene] <- row_ids[j]
  }
}
close(con)
stopifnot(all(!is.na(source_ensembl)), all(is.finite(expression)),
          all(expression >= 0))
for (gene in intersect(genes, c("ILF2", "DKC1", "GTPBP4"))) {
  stopifnot(all(abs(expression[gene, ] -
                      cohort[[paste0(gene, "_TPM")]]) < 1e-8))
}
stopifnot(identical(unname(sub("\\..*$", "", source_ensembl[genes])),
                    unname(top5$gene_id)))

cut_data <- data.frame(time = cohort$OS_time_days,
                       event = cohort$OS_event)
for (gene in genes)
  cut_data[[gene]] <- log2(expression[gene, ] + 1)
cut_object <- survminer::surv_cutpoint(
  cut_data, time = "time", event = "event", variables = genes,
  minprop = minprop, progressbar = FALSE)
categorized <- survminer::surv_categorize(
  cut_object, variables = genes, labels = c("Low", "High"))

prior_results <- read.delim(file.path(prior_dir, "results",
  "Figure2A_TCGA_BRCA_full_minprop10_five_gene_results.tsv"))
seed_by_gene <- c(ILF2 = 20261027L, DKC1 = 20261025L,
                  GTPBP4 = 20261026L, GNL2 = 20261150L,
                  NOP58 = 20261151L)
rows <- vector("list", length(genes))
for (i in seq_along(genes)) {
  gene <- genes[i]
  dat <- data.frame(time = cut_data$time, event = cut_data$event,
                    expression = cut_data[[gene]],
                    group = factor(categorized[[gene]],
                                   levels = c("Low", "High")))
  cutoff <- unname(cut_object$cutpoint[gene, "cutpoint"])
  fit <- survival::coxph(survival::Surv(time, event) ~ group,
                         data = dat, ties = "efron", x = TRUE)
  cox <- summary(fit)
  lr <- survival::survdiff(survival::Surv(time, event) ~ group,
                           data = dat)
  set.seed(seed_by_gene[[gene]])
  corrected <- maxstat::maxstat.test(
    survival::Surv(time, event) ~ expression, data = dat,
    smethod = "LogRank", pmethod = "condMC", B = condMC_B,
    minprop = minprop, maxprop = 1 - minprop)
  rows[[i]] <- data.frame(
    rank_by_TNBC_FDR = i, gene = gene, TNBC_vs_normal_log2FC =
      top5$log2FoldChange[i], TNBC_vs_normal_FDR = top5$padj[i],
    TNBC_minus_log10_FDR = -log10(top5$padj[i]),
    source_ensembl_id = source_ensembl[gene],
    OS_n = nrow(dat), OS_deaths = sum(dat$event),
    cutoff_log2TPM1 = cutoff, cutoff_TPM = 2^cutoff - 1,
    low_n = sum(dat$group == "Low"),
    low_deaths = sum(dat$event[dat$group == "Low"]),
    high_n = sum(dat$group == "High"),
    high_deaths = sum(dat$event[dat$group == "High"]),
    high_vs_low_HR = unname(cox$conf.int["groupHigh", "exp(coef)"]),
    CI_low = unname(cox$conf.int["groupHigh", "lower .95"]),
    CI_high = unname(cox$conf.int["groupHigh", "upper .95"]),
    Cox_Wald_P = unname(cox$coefficients["groupHigh", "Pr(>|z|)"]),
    selected_cutpoint_logrank_P = pchisq(lr$chisq, 1, lower.tail = FALSE),
    maxstat_condMC_adjusted_P = unname(corrected$p.value),
    maxstat_B = condMC_B, maxstat_seed = seed_by_gene[[gene]],
    cutoff_reproduced_by_condMC = isTRUE(all.equal(
      cutoff, unname(corrected$estimate))))
  if (gene %in% prior_results$gene) {
    old <- prior_results[prior_results$gene == gene, ]
    stopifnot(nrow(old) == 1L,
      isTRUE(all.equal(rows[[i]]$cutoff_log2TPM1,
                       old$cutoff_log2TPM1, tolerance = 1e-10)),
      isTRUE(all.equal(rows[[i]]$selected_cutpoint_logrank_P,
                       old$selected_cutpoint_logrank_p, tolerance = 1e-10)),
      isTRUE(all.equal(rows[[i]]$high_vs_low_HR,
                       old$high_vs_low_hr, tolerance = 1e-10)))
  }
  message(sprintf("%d/5 %s: HR=%.3f, nominal log-rank P=%.4g, maxstat P=%.4g",
                  i, gene, rows[[i]]$high_vs_low_HR,
                  rows[[i]]$selected_cutpoint_logrank_P,
                  rows[[i]]$maxstat_condMC_adjusted_P))
}
results <- do.call(rbind, rows)
results$BH_FDR_on_five_maxstat_P <- p.adjust(
  results$maxstat_condMC_adjusted_P, method = "BH")
stopifnot(identical(results$gene, genes),
          all(results$low_n + results$high_n == 1073L),
          all(results$low_n >= ceiling(minprop * 1073L)),
          all(results$high_n >= ceiling(minprop * 1073L)),
          all(results$cutoff_reproduced_by_condMC),
          all(is.finite(results$high_vs_low_HR)),
          all(results$CI_low < results$CI_high))
write.table(results, file.path(result_dir, paste0(prefix, "_source_data.tsv")),
            sep = "\t", quote = FALSE, row.names = FALSE)
saveRDS(cut_object, file.path(result_dir,
                             paste0(prefix, "_surv_cutpoint.rds")))

# One coordinated coordinate system keeps the FDR bars and forest rows aligned.
# Size is deliberately compact for a single Figure 2A panel (172 x 84 mm).
pt_mm <- 25.4 / 72.27
body <- 10.0 * pt_mm
header_size <- 10.2 * pt_mm
small <- 10.0 * pt_mm
width_mm <- 172
height_mm <- 84
results$y <- c(4.85, 3.92, 2.99, 2.06, 1.13)

bar_left <- 50
bar_right <- 207
bar_axis_max <- ceiling(max(results$TNBC_minus_log10_FDR) / 25) * 25
results$bar_end <- bar_left + results$TNBC_minus_log10_FDR /
  bar_axis_max * (bar_right - bar_left)
fc_min <- min(results$TNBC_vs_normal_log2FC)
fc_max <- max(results$TNBC_vs_normal_log2FC)
palette <- grDevices::colorRampPalette(c("#E7BFBD", "#AA4643"))(101)
colour_for_fc <- function(x) {
  k <- 1L + round((x - fc_min) / (fc_max - fc_min) * 100)
  palette[pmax(1L, pmin(101L, k))]
}
results$bar_colour <- colour_for_fc(results$TNBC_vs_normal_log2FC)

forest_left <- 317
forest_right <- 394
hr_min <- 0.5
hr_max <- max(c(5, 2^ceiling(log2(max(results$CI_high)))))
map_hr <- function(x) forest_left +
  (log(x) - log(hr_min)) / (log(hr_max) - log(hr_min)) *
  (forest_right - forest_left)
results$point_x <- map_hr(results$high_vs_low_HR)
results$ci_left_x <- map_hr(results$CI_low)
results$ci_right_x <- map_hr(results$CI_high)
results$hr_label <- sprintf("%.2f (%.2f–%.2f)",
  results$high_vs_low_HR, results$CI_low, results$CI_high)
results$p_label <- ifelse(results$selected_cutpoint_logrank_P < 0.0001,
  formatC(results$selected_cutpoint_logrank_P, format = "e", digits = 2),
  formatC(results$selected_cutpoint_logrank_P, format = "f", digits = 4))

ticks_bar <- seq(0, bar_axis_max, by = 50)
ticks_bar_x <- bar_left + ticks_bar / bar_axis_max * (bar_right - bar_left)
ticks_hr <- c(0.5, 1, 2, 5, 10, 20)
ticks_hr <- ticks_hr[ticks_hr <= hr_max]
gradient_data <- data.frame(
  ymin = seq(1.08, 3.08, length.out = 101)[1:100],
  ymax = seq(1.08, 3.08, length.out = 101)[2:101],
  fill = palette[1:100])

fig <- ggplot() +
  geom_rect(data = results,
    aes(xmin = bar_left, xmax = bar_end, ymin = y - 0.25,
        ymax = y + 0.25), fill = results$bar_colour,
    colour = NA) +
  geom_rect(data = gradient_data,
    aes(xmin = 178, xmax = 183, ymin = ymin, ymax = ymax),
    fill = gradient_data$fill, colour = NA) +
  geom_segment(aes(x = 0, xend = 430, y = 5.86, yend = 5.86),
               linewidth = 0.45) +
  geom_segment(aes(x = 0, xend = 430, y = 5.36, yend = 5.36),
               linewidth = 0.40) +
  geom_segment(aes(x = map_hr(1), xend = map_hr(1),
                   y = 0.72, yend = 5.30),
               linewidth = 0.35, linetype = "22", colour = "#777777") +
  geom_segment(data = results,
    aes(x = ci_left_x, xend = ci_right_x, y = y, yend = y),
    linewidth = 0.62, colour = "#4572A7") +
  geom_segment(data = results,
    aes(x = ci_left_x, xend = ci_left_x,
        y = y - 0.075, yend = y + 0.075),
    linewidth = 0.52, colour = "#4572A7") +
  geom_segment(data = results,
    aes(x = ci_right_x, xend = ci_right_x,
        y = y - 0.075, yend = y + 0.075),
    linewidth = 0.52, colour = "#4572A7") +
  geom_point(data = results, aes(x = point_x, y = y),
             shape = 17, size = 2.9, colour = "#4572A7") +
  geom_text(data = results, aes(x = 4, y = y, label = gene),
            hjust = 0, family = "Arial", size = body) +
  geom_text(data = results, aes(x = 224, y = y, label = hr_label),
            hjust = 0, family = "Arial", size = body) +
  geom_text(data = results, aes(x = 425, y = y, label = p_label),
            hjust = 1, family = "Arial", size = body) +
  geom_text(aes(x = 168, y = 2.08, label = "Log₂FC"), angle = 90,
            family = "Arial", size = small) +
  geom_text(aes(x = 187, y = 1.08, label = sprintf("%.2f", fc_min)),
            hjust = 0, family = "Arial", size = small) +
  geom_text(aes(x = 187, y = 3.08, label = sprintf("%.2f", fc_max)),
            hjust = 0, family = "Arial", size = small) +
  geom_text(aes(x = 4, y = 5.60, label = "Gene"), hjust = 0,
            family = "Arial", size = header_size) +
  geom_text(aes(x = (bar_left + bar_right) / 2, y = 5.60,
                label = "Differential expression"), hjust = 0.5,
            family = "Arial", size = header_size) +
  geom_text(aes(x = 224, y = 5.60,
                label = "Overall survival, HR (95% CI)"),
            hjust = 0, family = "Arial", size = header_size) +
  geom_text(aes(x = 425, y = 5.60, label = "P value"),
            hjust = 1, family = "Arial", size = header_size) +
  geom_segment(aes(x = bar_left, xend = bar_right,
                   y = 0.62, yend = 0.62), linewidth = 0.44) +
  geom_segment(aes(x = forest_left, xend = forest_right,
                   y = 0.62, yend = 0.62), linewidth = 0.44) +
  geom_segment(data = data.frame(x = ticks_bar_x),
    aes(x = x, xend = x, y = 0.62, yend = 0.53), linewidth = 0.4) +
  geom_text(data = data.frame(x = ticks_bar_x, label = ticks_bar),
    aes(x = x, y = 0.35, label = label), family = "Arial", size = small) +
  geom_segment(data = data.frame(x = map_hr(ticks_hr)),
    aes(x = x, xend = x, y = 0.62, yend = 0.53), linewidth = 0.4) +
  geom_text(data = data.frame(x = map_hr(ticks_hr), label = ticks_hr),
    aes(x = x, y = 0.35, label = label), family = "Arial", size = small) +
  geom_text(aes(x = (bar_left + bar_right) / 2, y = 0.03,
                label = "−Log₁₀(FDR)"), family = "Arial", size = body) +
  geom_text(aes(x = (forest_left + forest_right) / 2, y = 0.03,
                label = "Hazard ratio"), family = "Arial", size = body) +
  geom_segment(aes(x = 0, xend = 430, y = -0.23, yend = -0.23),
               linewidth = 0.45) +
  coord_cartesian(xlim = c(0, 430), ylim = c(-0.33, 6.0),
                  expand = FALSE, clip = "off") +
  theme_void(base_family = "Arial", base_size = 10) +
  theme(plot.margin = margin(2, 2.5, 2, 2.5, unit = "mm"),
        plot.background = element_rect(fill = "white", colour = NA))

stem <- file.path(fig_dir, prefix)
w <- width_mm / 25.4
h <- height_mm / 25.4
svglite::svglite(paste0(stem, ".svg"), width = w, height = h)
print(fig)
dev.off()
grDevices::cairo_pdf(paste0(stem, ".pdf"), width = w, height = h,
                     family = "Arial")
print(fig)
dev.off()
ragg::agg_png(paste0(stem, ".png"), width = w, height = h,
              units = "in", res = 400)
print(fig)
dev.off()
ragg::agg_tiff(paste0(stem, ".tiff"), width = w, height = h,
               units = "in", res = 600, compression = "lzw")
print(fig)
dev.off()

files <- paste0(stem, c(".svg", ".pdf", ".png", ".tiff"))
svg <- paste(readLines(paste0(stem, ".svg"), warn = FALSE), collapse = "\n")
qa <- data.frame(
  check = c("36 unique AMF genes matched to TNBC DE results",
            "top five in ascending original genome-wide FDR",
            "one-to-one TPM aliquot alignment, original expression rows",
            "1073 patients and 150 OS events",
            "10 percent minimum per survival group",
            "cutoffs independently reproduced by maxstat",
            "forest order equals FDR order",
            "finite HR and 95 percent CI",
            "four figure formats exported",
            "Arial SVG and 10 pt body"),
  pass = c(length(gene_pool) == 36L && nrow(de) == 36L,
           !is.unsorted(top5$padj),
           all(!is.na(chosen_aliquot)) && all(is.finite(expression)),
           nrow(cohort) == 1073L && sum(cohort$OS_event) == 150L,
           all(results$low_n >= 108L & results$high_n >= 108L),
           all(results$cutoff_reproduced_by_condMC),
           identical(results$gene, top5$gene_symbol),
           all(is.finite(results$high_vs_low_HR)) &&
             all(results$CI_low < results$CI_high),
           all(file.exists(files)) && all(file.info(files)$size > 1000),
           grepl("Arial", svg, fixed = TRUE) &&
             min(body, header_size, small) / pt_mm >= 10)
)
write.table(qa, file.path(result_dir, paste0(prefix, "_validation.tsv")),
            sep = "\t", quote = FALSE, row.names = FALSE)
if (!all(qa$pass)) stop("Figure validation failed; inspect validation.tsv.")
write.table(data.frame(parameter=c("expression_transform","minprop","condMC_resamples"),
 value=c("log2(TPM+1)",as.character(minprop),as.character(condMC_B))),
 file.path(result_dir,paste0(prefix,"_analysis_parameters.tsv")),sep="\t",quote=FALSE,row.names=FALSE)
writeLines(capture.output(sessionInfo()),
           file.path(result_dir, paste0(prefix, "_sessionInfo.txt")))
print(results[, c("gene", "TNBC_vs_normal_FDR", "TNBC_vs_normal_log2FC",
                  "low_n", "high_n", "high_vs_low_HR", "CI_low", "CI_high",
                  "selected_cutpoint_logrank_P", "maxstat_condMC_adjusted_P")],
      row.names = FALSE, digits = 5)
cat("Figure 2A complete: ", stem, "\n", sep = "")
