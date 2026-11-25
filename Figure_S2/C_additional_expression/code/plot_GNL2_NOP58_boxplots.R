#!/usr/bin/env Rscript


suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(ragg)
})
data.table::setDTthreads(1L)
options(width = 160)

script_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
script_path <- gsub("~\\+~", " ", sub("^--file=", "", script_arg[[1]]))
out_dir <- normalizePath(dirname(script_path), mustWork = TRUE)
project_dir <- normalizePath(file.path(out_dir, "..", ".."), mustWork = TRUE)
tcga_dir <- file.path(project_dir, "TCGA")
genes <- c("GNL2", "NOP58")
groups <- c("TN", "Non-TN", "Normal")
colors <- c(TN = "#AA4643", `Non-TN` = "#89A54E", Normal = "#4572A7")
font_pt <- 14
size_in <- 4
dpi <- 600L

metadata_path <- file.path(tcga_dir, "DESeq2_four_group_marker_analysis",
                           "metadata", "selected_samples.tsv")
tpm_path <- file.path(tcga_dir, "TCGA-BRCA_TPM.txt")
annotation_path <- file.path(tcga_dir, "DESeq2_TN_RBP_analysis",
  "KEGG_enrichment_top500_TN_all_genes", "reference_snapshot",
  "Ensembl_BioMart_hsapiens_gene_annotation.tsv")
stopifnot(all(file.exists(c(metadata_path, tpm_path, annotation_path))))

message("Checking sample metadata and gene identifiers...")
meta <- fread(metadata_path, showProgress = FALSE)[condition %chin% c("TN", "HR", "HER2", "NORMAL")]
stopifnot(!anyDuplicated(meta$sample_id), !anyNA(meta$patient_id))
expected <- c(TN = 177L, HR = 664L, HER2 = 163L, NORMAL = 112L)
counts <- meta[, .N, by = condition]
stopifnot(setequal(counts$condition, names(expected)), all(counts$N == expected[counts$condition]))
stopifnot(all(meta[condition != "NORMAL", sample_type] == "PRIMARY_TUMOR"),
          all(meta[condition == "NORMAL", sample_type] == "SOLID_TISSUE_NORMAL"))
meta[, display_group := fcase(condition == "TN", "TN",
  condition %chin% c("HR", "HER2"), "Non-TN", condition == "NORMAL", "Normal")]
stopifnot(!anyDuplicated(meta[, .(patient_id, display_group)]))

annotation <- fread(annotation_path, colClasses = "character", showProgress = FALSE)
gene_map <- unique(annotation[hgnc_symbol %chin% genes, .(gene = hgnc_symbol, ensembl_gene_id)])
stopifnot(nrow(gene_map) == 2L, uniqueN(gene_map$gene) == 2L,
          uniqueN(gene_map$ensembl_gene_id) == 2L)

# Read the exact two gene rows in a single pass, keeping all original columns.
extract_tpm <- function(path, ids) {
  con <- file(path, open = "r")
  on.exit(close(con))
  header <- readLines(con, n = 1L, warn = FALSE)
  pattern <- paste0('^"?(', paste(ids, collapse = "|"), ')([.]|"?\\t)')
  selected <- character()
  repeat {
    chunk <- readLines(con, n = 5000L, warn = FALSE)
    if (!length(chunk)) break
    selected <- c(selected, grep(pattern, chunk, value = TRUE))
  }
  stopifnot(length(selected) == length(ids))
  # The input's header omits the row-name field.
  fread(text = paste(c(paste("ensembl_gene_id_version", header, sep = "\t"), selected),
                     collapse = "\n"), check.names = FALSE, showProgress = FALSE)
}

message("Extracting GNL2 and NOP58 from the TPM matrix...")
tpm <- extract_tpm(tpm_path, gene_map$ensembl_gene_id)
stopifnot(all(meta$sample_id %in% names(tpm)))
tpm <- tpm[, c("ensembl_gene_id_version", meta$sample_id), with = FALSE]
tpm[, ensembl_gene_id := sub("\\..*$", "", ensembl_gene_id_version)]
tpm <- merge(tpm, gene_map, by = "ensembl_gene_id", all = FALSE, sort = FALSE)
stopifnot(nrow(tpm) == length(genes))
dat <- melt(tpm, id.vars = c("gene", "ensembl_gene_id", "ensembl_gene_id_version"),
            measure.vars = meta$sample_id, variable.name = "sample_id",
            value.name = "TPM", variable.factor = FALSE)
dat <- merge(dat, meta, by = "sample_id", all.x = TRUE, sort = FALSE)
stopifnot(nrow(dat) == length(genes) * nrow(meta), !anyNA(dat$display_group),
          all(is.finite(dat$TPM)), all(dat$TPM >= 0),
          !anyDuplicated(dat[, .(gene, sample_id)]))
dat[, log2TPM1 := log2(TPM + 1)]
dat[, display_group := factor(display_group, levels = groups)]
dat[, gene := factor(gene, levels = genes)]
setorder(dat, gene, display_group, patient_id)

summary <- dat[, .(n = .N, patients = uniqueN(patient_id),
  median_TPM = median(TPM), median_log2TPM1 = median(log2TPM1),
  Q1_log2TPM1 = unname(quantile(log2TPM1, .25)),
  Q3_log2TPM1 = unname(quantile(log2TPM1, .75)),
  minimum_log2TPM1 = min(log2TPM1), maximum_log2TPM1 = max(log2TPM1)),
  by = .(gene, display_group)]
stopifnot(all(summary$n == c(TN = 177L, `Non-TN` = 827L, Normal = 112L)[as.character(summary$display_group)]))

stats <- rbindlist(lapply(genes, function(g) {
  rbindlist(lapply(c("Non-TN", "Normal"), function(ref) {
    a <- dat[gene == g & display_group == "TN", log2TPM1]
    b <- dat[gene == g & display_group == ref, log2TPM1]
    fit <- wilcox.test(a, b, alternative = "two.sided", exact = FALSE, correct = FALSE)
    data.table(gene = g, comparison = paste("TN vs", ref), TN_n = length(a),
      reference_n = length(b), TN_median = median(a), reference_median = median(b),
      median_difference = median(a) - median(b), W = unname(fit$statistic),
      raw_P = fit$p.value, x_start = 1, x_end = match(ref, groups))
  }))
}))
stats[, BH_FDR := p.adjust(raw_P, method = "BH")]
stats[, significance := fcase(BH_FDR < 1e-4, "****", BH_FDR < 1e-3, "***",
  BH_FDR < 1e-2, "**", BH_FDR < .05, "*", default = "ns")]
stopifnot(nrow(stats) == 4L, all(is.finite(stats$raw_P)), all(is.finite(stats$BH_FDR)))

make_plot <- function(g) {
  d <- dat[gene == g]
  s <- copy(stats[gene == g])
  low <- min(d$log2TPM1)
  high <- max(d$log2TPM1)
  span <- max(high - low, 1)
  step <- max(span * .08, .28)
  s[, bracket_y := high + step * (0.45 + (x_end - 2) * .95)]
  s[, tip := step * .20]
  s[, label_y := bracket_y + step * .07]
  lower <- low - span * .035
  upper <- max(s$label_y) + max(step * .8, .24)
  breaks <- pretty(c(low, high), n = 4L)
  breaks <- breaks[breaks >= lower & breaks <= high]
  p <- ggplot(d, aes(display_group, log2TPM1, fill = display_group)) +
    geom_boxplot(width = .56, linewidth = .65, colour = "black",
                 outlier.shape = 16, outlier.size = .65, outlier.alpha = .35,
                 outlier.colour = "#666666") +
    geom_segment(data = s, aes(x = x_start, xend = x_end, y = bracket_y, yend = bracket_y),
                 inherit.aes = FALSE, linewidth = .55) +
    geom_segment(data = s, aes(x = x_start, xend = x_start, y = bracket_y, yend = bracket_y - tip),
                 inherit.aes = FALSE, linewidth = .55) +
    geom_segment(data = s, aes(x = x_end, xend = x_end, y = bracket_y, yend = bracket_y - tip),
                 inherit.aes = FALSE, linewidth = .55) +
    geom_text(data = s, aes(x = (x_start + x_end) / 2, y = label_y, label = significance),
              inherit.aes = FALSE, family = "Arial", fontface = "bold",
              size = font_pt / ggplot2::.pt, vjust = 0) +
    scale_fill_manual(values = colors, drop = FALSE) +
    scale_x_discrete(drop = FALSE, expand = expansion(add = .48)) +
    scale_y_continuous(breaks = breaks, limits = c(lower, upper), expand = expansion(mult = 0)) +
    labs(title = g, x = NULL, y = "Log\u2082(TPM + 1)") +
    theme_classic(base_size = font_pt, base_family = "Arial") +
    theme(text = element_text(size = font_pt, colour = "black"),
      plot.title = element_text(size = font_pt, face = "bold", hjust = .5, margin = margin(b = 6)),
      axis.title = element_text(size = font_pt),
      axis.title.y = element_text(size = font_pt, margin = margin(r = 6)),
      axis.text = element_text(size = font_pt, colour = "black"),
      axis.text.x = element_text(size = font_pt, margin = margin(t = 4)),
      axis.text.y = element_text(size = font_pt, margin = margin(r = 3)),
      axis.line = element_line(linewidth = .65, colour = "black"),
      axis.ticks = element_line(linewidth = .55, colour = "black"),
      axis.ticks.length = grid::unit(1.5, "mm"),
      legend.position = "none", plot.margin = margin(8, 9, 8, 8))
  # Verify the plotted box statistics against direct quantiles of all input values.
  built <- ggplot_build(p)$data[[1]]
  direct <- d[, .(lower = unname(quantile(log2TPM1, .25)), middle = median(log2TPM1),
                  upper = unname(quantile(log2TPM1, .75))), by = display_group]
  setorder(direct, display_group)
  stopifnot(isTRUE(all.equal(built$lower, direct$lower)),
            isTRUE(all.equal(built$middle, direct$middle)),
            isTRUE(all.equal(built$upper, direct$upper)))
  p
}

message("Exporting the two square PDF/PNG figures...")
for (g in genes) {
  p <- make_plot(g)
  basename <- file.path(out_dir, paste0(g, "_TN_Non-TN_Normal"))
  cairo_pdf(paste0(basename, ".pdf"), width = size_in, height = size_in, family = "Arial")
  print(p)
  dev.off()
  agg_png(paste0(basename, ".png"), width = size_in, height = size_in,
          units = "in", res = dpi, background = "white")
  print(p)
  dev.off()
}

fwrite(dat, file.path(out_dir, "source_data.tsv"), sep = "\t")
fwrite(summary, file.path(out_dir, "group_summary.tsv"), sep = "\t")
fwrite(stats, file.path(out_dir, "Wilcoxon_statistics.tsv"), sep = "\t")
fwrite(gene_map, file.path(out_dir, "gene_mapping.tsv"), sep = "\t")
quality <- meta[, .N, by = .(display_group, condition, condition_strict_qc, quality_flag)]
fwrite(quality, file.path(out_dir, "cohort_quality_flags.tsv"), sep = "\t")
capture.output(sessionInfo(), file = file.path(out_dir, "R_sessionInfo.txt"))
outputs <- unlist(lapply(genes, function(g) file.path(out_dir,
  paste0(g, "_TN_Non-TN_Normal", c(".pdf", ".png")))))
stopifnot(all(file.exists(outputs)), all(file.info(outputs)$size > 0))
print(summary)
print(stats[, .(gene, comparison, raw_P, BH_FDR, significance)])
message("Completed: 2 PDF files and 2 PNG files; all text 14 pt; 4 inch square.")
