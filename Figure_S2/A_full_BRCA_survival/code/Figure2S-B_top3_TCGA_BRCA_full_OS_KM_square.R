#!/usr/bin/env Rscript

# Three separate square Supplementary Figure S2B OS Kaplan-Meier panels.
# Default 10%-minimum groups exactly match the TCGA-BRCA full-cohort forest.
# Set FIGURE2SB_MINPROP=0.20 to reproduce the historical S2A/S2B groups.
# The panel displays nominal log-rank P; corrected P remains in source results.

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

genes <- c("DKC1", "GTPBP4", "ILF2")
minprop <- as.numeric(Sys.getenv("FIGURE2SB_MINPROP", unset = "0.10"))
if (!minprop %in% c(0.10, 0.20)) stop("Only 0.10 and 0.20 are supported.")
is_forest_cutoff <- isTRUE(all.equal(minprop, 0.10))
data_name <- if (is_forest_cutoff) {
  "Figure2A_TCGA_BRCA_full_minprop10_analysis_ready.tsv"
} else {
  "TCGA_BRCA_full_five_gene_OS_analysis_ready.tsv"
}
groups_name <- if (is_forest_cutoff) {
  "Figure2A_TCGA_BRCA_full_minprop10_five_gene_groups.tsv"
} else {
  "best_cutoff_sample_groups.tsv"
}
reference_name <- if (is_forest_cutoff) {
  "Figure2A_TCGA_BRCA_full_minprop10_five_gene_results.tsv"
} else {
  "best_cutoff_OS_results.tsv"
}
base <- read.delim(file.path(project_dir, "data", "processed", data_name),
                   check.names = FALSE)
groups <- read.delim(file.path(project_dir, "results", groups_name),
                     check.names = FALSE)
reference <- read.delim(file.path(project_dir, "results", reference_name),
                        check.names = FALSE)
base_id <- if (is_forest_cutoff) base$patient_id else base$sample_id
figures_dir <- file.path(project_dir, "figures")
results_dir <- file.path(project_dir, "results")
dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)
stopifnot(nrow(base) == 1073L, sum(base$OS_event) == 150L,
          !anyDuplicated(base_id), all(base$OS_time_days > 0))

high_color <- "#B34A43"
low_color <- "#3F70A5"
pointsize <- 10
side_mm <- 88
curve_lwd <- 2.2
summary_rows <- vector("list", length(genes))
square_ratios <- numeric(length(genes))

for (i in seq_along(genes)) {
  gene <- genes[i]
  ref <- reference[reference$gene == gene, , drop = FALSE]
  g <- groups[groups$gene == gene, , drop = FALSE]
  group_id <- if (is_forest_cutoff) g$patient_id else g$sample_id
  idx <- match(base_id, group_id)
  stopifnot(nrow(ref) == 1L, nrow(g) == nrow(base),
            !anyDuplicated(group_id), !anyNA(idx))
  g <- g[idx, , drop = FALSE]
  d <- data.frame(
    months = base$OS_time_months,
    days = base$OS_time_days,
    event = as.integer(base$OS_event),
    expression = base[[paste0(gene, "_log2TPM1")]],
    group = factor(g$group, levels = c("Low", "High"))
  )
  stopifnot(all(d$event == g$OS_event),
            isTRUE(all.equal(d$days, g$OS_time_days)),
            isTRUE(all.equal(d$expression, g$expression_log2TPM1)),
            all(as.character(d$group) == ifelse(
              d$expression <= ref$cutoff_log2TPM1, "Low", "High")),
            sum(d$group == "Low") == ref$low_n,
            sum(d$group == "High") == ref$high_n,
            min(ref$low_n, ref$high_n) >=
              if (is_forest_cutoff) floor(minprop * nrow(d))
              else ceiling(minprop * nrow(d)))

  fit <- survfit(Surv(months, event) ~ group, data = d)
  lr <- survdiff(Surv(months, event) ~ group, data = d)
  p_value <- pchisq(lr$chisq, 1, lower.tail = FALSE)
  reference_p <- if (is_forest_cutoff) ref$selected_cutpoint_logrank_p
                 else ref$naive_best_logrank_p
  stopifnot(abs(p_value - reference_p) < 1e-10)
  low_n <- sum(d$group == "Low")
  high_n <- sum(d$group == "High")

  draw_panel <- function() {
    par(family = "Arial", mar = c(3.7, 4.2, 1.8, 0.8),
        pty = "s", las = 1, bty = "n", mgp = c(2.40, 0.65, 0),
        tcl = -0.26, xaxs = "i", yaxs = "i",
        lend = "butt", ljoin = "mitre")
    plot(fit, col = c(low_color, high_color), lty = c(1, 1),
         lwd = c(curve_lwd, curve_lwd), mark.time = TRUE,
         pch = 3, cex = 0.48, xlim = c(0, 305), ylim = c(0, 1.045),
         xaxt = "n", yaxt = "n", xlab = "", ylab = "", conf.int = FALSE)
    axis(1, at = c(0, 100, 200, 300),
         labels = c("0", "100", "200", "300"), pos = 0,
         lwd = 1.45, lwd.ticks = 1.45, cex.axis = 1)
    axis(2, at = c(0, 0.25, 0.50, 0.75, 1),
         labels = c("0", "25%", "50%", "75%", "100%"), pos = 0,
         lwd = 1.45, lwd.ticks = 1.45, cex.axis = 1)
    mtext("Months", side = 1, line = 2.15, cex = 1)
    mtext("Probability of Survival", side = 2, line = 2.55,
          las = 0, cex = 1)
    title(main = gene, line = 0.35, cex.main = 1.10, font.main = 1)
    legend(15, 0.315, xjust = 0, yjust = 1,
           legend = c(sprintf("High expression(n=%d)", high_n),
                      sprintf("Low expression(n=%d)", low_n)),
           title = if (is_forest_cutoff) sprintf("p=%.4f", p_value)
                   else sprintf("p=%.3f", p_value), title.adj = 0,
           col = c(high_color, low_color), lwd = curve_lwd,
           lty = 1, pch = 3, pt.cex = 0.78, bty = "n", cex = 1,
           y.intersp = 0.92, x.intersp = 0.42, seg.len = 1.15)
    pin <- par("pin")
    unname(pin[1] / pin[2])
  }

  stem <- file.path(figures_dir,
                    paste0("Figure2S-B_", gene, "_TCGA_BRCA_full_",
                           if (is_forest_cutoff) "minprop10_" else "",
                           "OS_KM_square"))
  side_in <- side_mm / 25.4
  svglite::svglite(paste0(stem, ".svg"), width = side_in,
                   height = side_in, pointsize = pointsize)
  ratio <- draw_panel()
  dev.off()
  grDevices::cairo_pdf(paste0(stem, ".pdf"), width = side_in,
                       height = side_in, family = "Arial", pointsize = pointsize)
  draw_panel()
  dev.off()
  ragg::agg_tiff(paste0(stem, ".tiff"), width = side_in,
                 height = side_in, units = "in", res = 600,
                 compression = "lzw", pointsize = pointsize)
  draw_panel()
  dev.off()
  ragg::agg_png(paste0(stem, ".png"), width = side_in,
                height = side_in, units = "in", res = 300,
                pointsize = pointsize)
  draw_panel()
  dev.off()

  square_ratios[i] <- ratio
  summary_rows[[i]] <- data.frame(
    gene = gene, n = nrow(d), deaths = sum(d$event), min_group_fraction = minprop,
    cutoff_log2TPM1 = ref$cutoff_log2TPM1,
    low_n = low_n, low_deaths = if (is_forest_cutoff) ref$low_deaths else ref$low_events,
    high_n = high_n, high_deaths = if (is_forest_cutoff) ref$high_deaths else ref$high_events,
    nominal_logrank_p = p_value,
    maxstat_condMC_adjusted_p = ref$maxstat_condMC_adjusted_p,
    plot_region_width_to_height = ratio,
    output_stem = basename(stem)
  )
}

summary_table <- do.call(rbind, summary_rows)
prefix <- paste0("Figure2S-B_top3_TCGA_BRCA_full_",
                 if (is_forest_cutoff) "minprop10_" else "",
                 "OS_KM_square")
write.table(summary_table, file.path(results_dir, paste0(prefix, "_results.tsv")),
            sep = "\t", quote = FALSE, row.names = FALSE)
figure_files <- as.vector(outer(
  file.path(figures_dir, summary_table$output_stem),
  c(".svg", ".pdf", ".tiff", ".png"), paste0))
qa <- data.frame(
  check = c("only top three genes", "same full cohort and selected-cutoff groups",
            "all canvases square", "all plot regions square",
            "Arial at least 10 pt", "nominal P reproduced", "all exports nonempty"),
  pass = c(identical(summary_table$gene, genes),
           all(summary_table$n == 1073L & summary_table$deaths == 150L),
           side_mm == 88 && pointsize >= 10,
           all(abs(square_ratios - 1) < 0.015),
           all(vapply(summary_table$output_stem, function(x) {
             svg <- paste(readLines(file.path(figures_dir, paste0(x, ".svg")),
                                    warn = FALSE), collapse = "\n")
             grepl("Arial", svg, fixed = TRUE)
           }, logical(1))),
           all(abs(summary_table$nominal_logrank_p -
                     if (is_forest_cutoff) {
                       reference$selected_cutpoint_logrank_p[match(genes, reference$gene)]
                     } else {
                       reference$naive_best_logrank_p[match(genes, reference$gene)]
                     }) < 1e-10),
           all(file.exists(figure_files)) && all(file.info(figure_files)$size > 1000))
)
write.table(qa, file.path(results_dir, paste0(prefix, "_validation.tsv")),
            sep = "\t", quote = FALSE, row.names = FALSE)
if (!all(qa$pass)) stop("Square three-panel KM validation failed.")
cat("Three separate 88 x 88 mm full-cohort KM figures created.\n")
print(summary_table[, c("gene", "low_n", "high_n", "nominal_logrank_p",
                        "plot_region_width_to_height")], row.names = FALSE)
