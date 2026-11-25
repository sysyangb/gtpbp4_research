#!/usr/bin/env Rscript

# Supplementary Figure S2B: TCGA-BRCA full-cohort, five-gene OS KM grid.
# Default 10%-minimum groups exactly match the full-cohort five-gene forest.
# Set FIGURE2SB_MINPROP=0.20 to reproduce the historical S2A grouping.
# Display nominal log-rank P at each outcome-selected cutpoint;
# adjusted P values are retained in the companion results table, not the plot.

options(stringsAsFactors = FALSE)
script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (length(script_arg) != 1L) stop("Run with Rscript.")
script_path <- normalizePath(sub("^--file=", "", script_arg), mustWork = TRUE)
project_dir <- normalizePath(file.path(dirname(script_path), ".."), mustWork = TRUE)
suppressPackageStartupMessages({
  library(survival)
  library(svglite)
  library(ragg)
})

genes <- c("DKC1", "GTPBP4", "ILF2", "NOP2", "RPP38")
minprop <- as.numeric(Sys.getenv("FIGURE2SB_MINPROP", unset = "0.10"))
if (!minprop %in% c(0.10, 0.20)) stop("Only 0.10 and 0.20 are supported.")
is_forest_cutoff <- isTRUE(all.equal(minprop, 0.10))
prefix <- paste0("Figure2S-B_TCGA_BRCA_full_five_gene_",
                 if (is_forest_cutoff) "minprop10_" else "", "OS_KM")
data_file <- file.path(project_dir, "data", "processed",
                       if (is_forest_cutoff) {
                         "Figure2A_TCGA_BRCA_full_minprop10_analysis_ready.tsv"
                       } else {
                         "TCGA_BRCA_full_five_gene_OS_analysis_ready.tsv"
                       })
groups_file <- file.path(project_dir, "results",
                         if (is_forest_cutoff) {
                           "Figure2A_TCGA_BRCA_full_minprop10_five_gene_groups.tsv"
                         } else {
                           "best_cutoff_sample_groups.tsv"
                         })
reference_file <- file.path(project_dir, "results",
                            if (is_forest_cutoff) {
                              "Figure2A_TCGA_BRCA_full_minprop10_five_gene_results.tsv"
                            } else {
                              "best_cutoff_OS_results.tsv"
                            })
figures_dir <- file.path(project_dir, "figures")
results_dir <- file.path(project_dir, "results")
if (!all(file.exists(c(data_file, groups_file, reference_file))))
  stop("A required full-cohort analysis source is missing.")
dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)

base <- read.delim(data_file, check.names = FALSE)
groups <- read.delim(groups_file, check.names = FALSE)
reference <- read.delim(reference_file, check.names = FALSE)
base_id <- if (is_forest_cutoff) base$patient_id else base$sample_id
stopifnot(nrow(base) == 1073L, sum(base$OS_event) == 150L,
          all(base$OS_time_days > 0), !anyDuplicated(base_id),
          setequal(genes, unique(groups$gene)), setequal(genes, reference$gene),
          nrow(groups) == nrow(base) * length(genes))

km_objects <- vector("list", length(genes))
names(km_objects) <- genes
result_rows <- vector("list", length(genes))
source_rows <- vector("list", length(genes))
for (i in seq_along(genes)) {
  gene <- genes[i]
  ref <- reference[reference$gene == gene, , drop = FALSE]
  g <- groups[groups$gene == gene, , drop = FALSE]
  group_id <- if (is_forest_cutoff) g$patient_id else g$sample_id
  idx <- match(base_id, group_id)
  stopifnot(nrow(g) == nrow(base), !anyDuplicated(group_id), !anyNA(idx))
  g <- g[idx, , drop = FALSE]
  d <- data.frame(
    patient_id = base_id,
    months = base$OS_time_months,
    days = base$OS_time_days,
    event = as.integer(base$OS_event),
    expression = base[[paste0(gene, "_log2TPM1")]],
    group = factor(g$group, levels = c("Low", "High"))
  )
  stopifnot(
    all(d$event == g$OS_event),
    isTRUE(all.equal(d$days, g$OS_time_days)),
    isTRUE(all.equal(d$expression, g$expression_log2TPM1)),
    all(as.character(d$group) == ifelse(d$expression <= ref$cutoff_log2TPM1,
                                       "Low", "High")),
    sum(d$group == "Low") == ref$low_n,
    sum(d$group == "High") == ref$high_n,
    min(ref$low_n, ref$high_n) >=
      if (is_forest_cutoff) floor(minprop * nrow(d))
      else ceiling(minprop * nrow(d))
  )
  fit <- survfit(Surv(months, event) ~ group, data = d)
  logrank <- survdiff(Surv(months, event) ~ group, data = d)
  p_value <- pchisq(logrank$chisq, df = 1, lower.tail = FALSE)
  reference_p <- if (is_forest_cutoff) ref$selected_cutpoint_logrank_p
                 else ref$naive_best_logrank_p
  stopifnot(abs(p_value - reference_p) < 1e-10)
  km_objects[[gene]] <- list(fit = fit, data = d, ref = ref, p = p_value)
  result_rows[[i]] <- data.frame(
    gene = gene, cohort = "TCGA-BRCA full cohort", n = nrow(d),
    deaths = sum(d$event), min_group_fraction = minprop,
    cutoff_log2TPM1 = ref$cutoff_log2TPM1,
    low_n = ref$low_n,
    low_deaths = if (is_forest_cutoff) ref$low_deaths else ref$low_events,
    high_n = ref$high_n,
    high_deaths = if (is_forest_cutoff) ref$high_deaths else ref$high_events,
    high_vs_low_hr = ref$high_vs_low_hr,
    ci_low = if (is_forest_cutoff) ref$ci_low else ref$high_vs_low_ci_low,
    ci_high = if (is_forest_cutoff) ref$ci_high else ref$high_vs_low_ci_high,
    displayed_nominal_logrank_p = p_value,
    maxstat_condMC_adjusted_p = ref$maxstat_condMC_adjusted_p,
    BH_FDR_on_maxstat_adjusted_p = if (is_forest_cutoff) {
      ref$BH_FDR_on_maxstat_adjusted_p
    } else {
      ref$maxstat_adjusted_fdr_bh
    }
  )
  sf <- summary(fit, censored = TRUE)
  source_rows[[i]] <- data.frame(
    gene = gene, group = sub("^group=", "", as.character(sf$strata)),
    time_months = sf$time, survival = sf$surv, n_risk = sf$n.risk,
    n_event = sf$n.event, n_censor = sf$n.censor
  )
}
results <- do.call(rbind, result_rows)
source_data <- do.call(rbind, source_rows)
write.table(results, file.path(results_dir, paste0(prefix, "_results.tsv")),
            sep = "\t", quote = FALSE, row.names = FALSE)
write.table(source_data, file.path(results_dir, paste0(prefix, "_curve_source_data.tsv")),
            sep = "\t", quote = FALSE, row.names = FALSE)
writeLines(capture.output(sessionInfo()),
           file.path(results_dir, paste0(prefix, "_sessionInfo.txt")))

high_color <- "#B34A43"
low_color <- "#3F70A5"
curve_lwd <- 2.2
pointsize <- 10
width_mm <- 183
height_mm <- 138

draw_km_panel <- function(gene, compact = TRUE) {
  item <- km_objects[[gene]]
  d <- item$data
  fit <- item$fit
  if (compact) {
    par(family = "Arial", mar = c(3.35, 3.75, 1.85, 0.65),
        las = 1, bty = "n", mgp = c(2.25, 0.55, 0), tcl = -0.25,
        xaxs = "i", yaxs = "i", lend = "butt", ljoin = "mitre")
  } else {
    par(family = "Arial", mar = c(3.7, 4.2, 1.8, 0.8),
        las = 1, bty = "n", mgp = c(2.40, 0.65, 0), tcl = -0.26,
        xaxs = "i", yaxs = "i", lend = "butt", ljoin = "mitre")
  }
  plot(fit, col = c(low_color, high_color), lty = c(1, 1),
       lwd = c(curve_lwd, curve_lwd), mark.time = TRUE,
       pch = 3, cex = 0.48, xlim = c(0, 305), ylim = c(0, 1.045),
       xaxt = "n", yaxt = "n", xlab = "", ylab = "", conf.int = FALSE)
  axis(1, at = c(0, 100, 200, 300), labels = c("0", "100", "200", "300"),
       pos = 0, lwd = 1.45, lwd.ticks = 1.45, cex.axis = 1.00)
  axis(2, at = c(0, 0.25, 0.50, 0.75, 1),
       labels = c("0", "25%", "50%", "75%", "100%"),
       pos = 0, lwd = 1.45, lwd.ticks = 1.45, cex.axis = 1.00)
  mtext("Months", side = 1, line = 2.15, cex = 1.00)
  mtext("Probability of Survival", side = 2, line = 2.55,
        las = 0, cex = 1.00)
  title(main = gene, line = 0.35, cex.main = 1.10, font.main = 1)
  group_labels <- if (compact) {
    c(sprintf("High (n=%d)", sum(d$group == "High")),
      sprintf("Low (n=%d)", sum(d$group == "Low")))
  } else {
    c(sprintf("High expression (n=%d)", sum(d$group == "High")),
      sprintf("Low expression (n=%d)", sum(d$group == "Low")))
  }
  legend(15, 0.315, xjust = 0, yjust = 1,
         legend = group_labels, title = sprintf("p=%.4f", item$p),
         title.adj = 0, col = c(high_color, low_color), lwd = curve_lwd,
         lty = 1, pch = 3, pt.cex = 0.78, bty = "n", cex = 1.00,
         y.intersp = 0.92, x.intersp = 0.42, seg.len = 1.15)
}

draw_composite <- function() {
  old_par <- par(no.readonly = TRUE)
  on.exit(par(old_par))
  layout(matrix(c(1, 1, 2, 2, 3, 3,
                  0, 4, 4, 5, 5, 0), nrow = 2, byrow = TRUE),
         widths = rep(1, 6), heights = c(1, 1))
  par(oma = c(0, 0, 1.35, 0))
  for (gene in genes) draw_km_panel(gene, compact = TRUE)
  mtext("B", side = 3, outer = TRUE, line = -0.55, adj = 0,
        family = "Arial", cex = 1.10, font = 2)
}

stem <- file.path(figures_dir, prefix)
width_in <- width_mm / 25.4
height_in <- height_mm / 25.4
svglite::svglite(paste0(stem, ".svg"), width = width_in, height = height_in,
                 pointsize = pointsize)
draw_composite()
dev.off()
grDevices::cairo_pdf(paste0(stem, ".pdf"), width = width_in,
                     height = height_in, family = "Arial", pointsize = pointsize)
draw_composite()
dev.off()
ragg::agg_tiff(paste0(stem, ".tiff"), width = width_in,
               height = height_in, units = "in", res = 600,
               pointsize = pointsize, compression = "lzw")
draw_composite()
dev.off()
ragg::agg_png(paste0(stem, ".png"), width = width_in,
              height = height_in, units = "in", res = 300,
              pointsize = pointsize)
draw_composite()
dev.off()

for (gene in genes) {
  individual_stem <- file.path(figures_dir,
                               paste0("Figure2S-B_", gene, "_TCGA_BRCA_full_",
                                      if (is_forest_cutoff) "minprop10_" else "",
                                      "OS_KM"))
  grDevices::cairo_pdf(paste0(individual_stem, ".pdf"),
                       width = 88 / 25.4, height = 88 / 25.4,
                       family = "Arial", pointsize = pointsize)
  draw_km_panel(gene, compact = FALSE)
  dev.off()
  ragg::agg_png(paste0(individual_stem, ".png"),
                width = 88 / 25.4, height = 88 / 25.4,
                units = "in", res = 300, pointsize = pointsize)
  draw_km_panel(gene, compact = FALSE)
  dev.off()
}

figure_files <- c(paste0(stem, c(".svg", ".pdf", ".tiff", ".png")),
                  as.vector(outer(
                    file.path(figures_dir,
                              paste0("Figure2S-B_", genes, "_TCGA_BRCA_full_",
                                     if (is_forest_cutoff) "minprop10_" else "",
                                     "OS_KM")),
                    c(".pdf", ".png"), paste0)))
svg_text <- paste(readLines(paste0(stem, ".svg"), warn = FALSE),
                  collapse = "\n")
qa <- data.frame(
  check = c("patient-level full TCGA-BRCA cohort", "five genes and 150 deaths",
            "same cutpoint groups as chosen forest", "nominal P values reproduced",
            "identical axis ranges and curve widths", "10 pt Arial requested",
            "composite SVG contains Arial", "all composite and individual exports exist"),
  pass = c(nrow(base) == 1073L && !anyDuplicated(base_id),
           identical(results$gene, genes) && all(results$deaths == 150L),
           all(results$min_group_fraction == minprop) &&
             isTRUE(all.equal(results$cutoff_log2TPM1,
                              reference$cutoff_log2TPM1[match(genes, reference$gene)])),
           isTRUE(all.equal(results$displayed_nominal_logrank_p,
                            if (is_forest_cutoff) {
                              reference$selected_cutpoint_logrank_p[match(genes, reference$gene)]
                            } else {
                              reference$naive_best_logrank_p[match(genes, reference$gene)]
                            },
                            tolerance = 1e-10)),
           identical(curve_lwd, 2.2) && max(base$OS_time_months) <= 300,
           pointsize >= 10, grepl("Arial", svg_text, fixed = TRUE),
           all(file.exists(figure_files)) && all(file.info(figure_files)$size > 1000))
)
write.table(qa, file.path(results_dir, paste0(prefix, "_validation.tsv")),
            sep = "\t", quote = FALSE, row.names = FALSE)
if (!all(qa$pass)) stop("S2B full-cohort KM validation failed.")
cat("Supplementary Figure S2B: TCGA-BRCA full cohort (n=1073, deaths=150).\n")
print(results[, c("gene", "low_n", "high_n", "displayed_nominal_logrank_p",
                  "maxstat_condMC_adjusted_p")], row.names = FALSE, digits = 5)
