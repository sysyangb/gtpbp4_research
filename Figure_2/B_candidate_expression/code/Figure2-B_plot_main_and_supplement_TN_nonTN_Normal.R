#!/usr/bin/env Rscript

# Figure2-B main and supplementary boxplots.
# Main: ILF2, DKC1 and GTPBP4.
# Supplementary: RPP38 and NOP2.
# Within each gene panel: TN (left), non-TN = HR + HER2 (middle), and Normal
# (right). Expression is log2(TPM + 1). Two-sided Wilcoxon rank-sum tests
# compare TN vs non-TN and TN vs Normal. Stars reflect BH FDR across all ten
# planned comparisons.

args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
script_path <- if (length(file_arg)) sub("^--file=", "", file_arg[[1]]) else "."
script_path <- gsub("~\\+~", " ", script_path)
output_dir <- normalizePath(dirname(script_path), mustWork = TRUE)
project_dir <- normalizePath(file.path(output_dir, "..", ".."), mustWork = TRUE)
tcga_dir <- file.path(project_dir, "TCGA")

local_lib <- file.path(tcga_dir, ".Rlib")
if (dir.exists(local_lib)) .libPaths(c(local_lib, .libPaths()))

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(ragg)
  library(svglite)
})

Sys.setenv(
  OMP_NUM_THREADS = "1",
  OPENBLAS_NUM_THREADS = "1",
  MKL_NUM_THREADS = "1",
  VECLIB_MAXIMUM_THREADS = "1"
)
options(mc.cores = 1L, width = 180)
data.table::setDTthreads(1L)

main_genes <- c("ILF2", "DKC1", "GTPBP4")
supplementary_genes <- c("RPP38", "NOP2")
genes <- c(main_genes, supplementary_genes)
source_conditions <- c("TN", "HR", "HER2", "NORMAL")
display_groups <- c("TN", "non-TN", "Normal")
prefix <- "Figure2-B"

tpm_path <- file.path(tcga_dir, "TCGA-BRCA_TPM.txt")
sample_meta_path <- file.path(
  tcga_dir, "DESeq2_four_group_marker_analysis", "metadata", "selected_samples.tsv"
)
annotation_path <- file.path(
  tcga_dir, "DESeq2_TN_RBP_analysis", "KEGG_enrichment_top500_TN_all_genes",
  "reference_snapshot", "Ensembl_BioMart_hsapiens_gene_annotation.tsv"
)
required_files <- c(tpm_path, sample_meta_path, annotation_path)
if (any(!file.exists(required_files))) {
  stop("Missing required input(s): ", paste(required_files[!file.exists(required_files)], collapse = ", "))
}

message("Reading the four-group sample metadata...")
sample_meta <- fread(sample_meta_path, showProgress = FALSE)[condition %chin% source_conditions]
required_meta <- c("sample_id", "patient_id", "condition", "sample_type")
if (!all(required_meta %in% names(sample_meta))) stop("Sample metadata lacks required columns")
if (anyDuplicated(sample_meta$sample_id)) stop("Duplicated sample_id in selected metadata")

expected_source_counts <- c(TN = 177L, HR = 664L, HER2 = 163L, NORMAL = 112L)
source_counts <- sample_meta[, .N, by = condition]
if (!all(expected_source_counts[source_counts$condition] == source_counts$N)) {
  stop("Unexpected four-group sample counts")
}

message("Mapping genes and extracting TPM values...")
annotation <- fread(annotation_path, colClasses = "character", showProgress = FALSE)
gene_map <- unique(annotation[
  toupper(hgnc_symbol) %chin% genes,
  .(gene = toupper(hgnc_symbol), ensembl_gene_id)
])
if (nrow(gene_map) != length(genes) || uniqueN(gene_map$gene) != length(genes)) {
  stop("The five genes did not map one-to-one")
}

tpm <- suppressWarnings(fread(
  tpm_path,
  select = c("V1", sample_meta$sample_id),
  check.names = FALSE,
  showProgress = FALSE
))
setnames(tpm, "V1", "ensembl_gene_id_version")
tpm[, ensembl_gene_id := sub("\\..*$", "", ensembl_gene_id_version)]
tpm <- merge(tpm, gene_map, by = "ensembl_gene_id", all = FALSE, sort = FALSE)
if (nrow(tpm) != length(genes) || uniqueN(tpm$gene) != length(genes)) {
  stop("Expected one TPM row per gene")
}

plot_data <- melt(
  tpm,
  id.vars = c("ensembl_gene_id", "ensembl_gene_id_version", "gene"),
  measure.vars = sample_meta$sample_id,
  variable.name = "sample_id",
  value.name = "TPM",
  variable.factor = FALSE
)
plot_data <- merge(
  plot_data,
  sample_meta[, .(sample_id, patient_id, condition, sample_type)],
  by = "sample_id",
  all.x = TRUE,
  sort = FALSE
)
if (anyNA(plot_data$condition)) stop("Expression-to-metadata merge lost samples")
if (any(!is.finite(plot_data$TPM)) || any(plot_data$TPM < 0)) stop("Invalid TPM detected")

plot_data[, log2TPM1 := log2(TPM + 1)]
plot_data[, display_group := fcase(
  condition == "TN", "TN",
  condition %chin% c("HR", "HER2"), "non-TN",
  condition == "NORMAL", "Normal",
  default = NA_character_
)]
if (anyNA(plot_data$display_group)) stop("Unexpected condition-to-display-group mapping")
plot_data[, display_group := factor(display_group, levels = display_groups)]
plot_data[, figure_set := fifelse(gene %chin% main_genes, "Main", "Supplementary")]
plot_data[, gene := factor(gene, levels = genes)]
setorder(plot_data, figure_set, gene, display_group, patient_id, sample_id)

message("Performing TN versus non-TN and TN versus Normal tests...")
comparison_definitions <- list(
  list(label = "TN_vs_non-TN", reference = "non-TN", x_end = 2L, level = 1L),
  list(label = "TN_vs_Normal", reference = "Normal", x_end = 3L, level = 2L)
)

statistics <- rbindlist(lapply(genes, function(gene_name) {
  rbindlist(lapply(comparison_definitions, function(definition) {
    tn_values <- plot_data[gene == gene_name & display_group == "TN", log2TPM1]
    reference_values <- plot_data[
      gene == gene_name & display_group == definition$reference,
      log2TPM1
    ]
    fit <- wilcox.test(
      tn_values,
      reference_values,
      alternative = "two.sided",
      exact = FALSE,
      correct = FALSE
    )
    data.table(
      gene = gene_name,
      figure_set = ifelse(gene_name %chin% main_genes, "Main", "Supplementary"),
      comparison = definition$label,
      reference_group = definition$reference,
      TN_n = length(tn_values),
      reference_n = length(reference_values),
      TN_median_TPM = median(2^tn_values - 1),
      reference_median_TPM = median(2^reference_values - 1),
      TN_median_log2TPM1 = median(tn_values),
      reference_median_log2TPM1 = median(reference_values),
      median_difference_log2TPM1 = median(tn_values) - median(reference_values),
      direction = ifelse(median(tn_values) > median(reference_values), "TN higher", "TN lower"),
      Wilcoxon_W = unname(fit$statistic),
      raw_P = fit$p.value,
      x_start = 1L,
      x_end = definition$x_end,
      annotation_level = definition$level
    )
  }))
}))
statistics[, BH_FDR := p.adjust(raw_P, method = "BH")]
statistics[, significance := fcase(
  BH_FDR < 1e-4, "****",
  BH_FDR < 1e-3, "***",
  BH_FDR < 1e-2, "**",
  BH_FDR < 0.05, "*",
  default = "ns"
)]
statistics[, gene := factor(gene, levels = genes)]

panel_ranges <- plot_data[, .(
  panel_min = min(log2TPM1),
  panel_max = max(log2TPM1),
  panel_range = max(log2TPM1) - min(log2TPM1)
), by = gene]
statistics <- merge(statistics, panel_ranges, by = "gene", all.x = TRUE, sort = FALSE)
statistics[, annotation_step := pmax(panel_range * 0.11, 0.14)]
statistics[, bracket_y := panel_max + annotation_step * (0.48 + 0.88 * (annotation_level - 1L))]
statistics[, tip_height := annotation_step * 0.19]
statistics[, label_y := bracket_y + annotation_step * 0.07]
setorder(statistics, gene, annotation_level)

group_summary <- plot_data[, .(
  n = .N,
  n_unique_patients = uniqueN(patient_id),
  median_TPM = median(TPM),
  Q1_TPM = quantile(TPM, 0.25),
  Q3_TPM = quantile(TPM, 0.75),
  median_log2TPM1 = median(log2TPM1),
  Q1_log2TPM1 = quantile(log2TPM1, 0.25),
  Q3_log2TPM1 = quantile(log2TPM1, 0.75)
), by = .(figure_set, gene, display_group)]
setorder(group_summary, figure_set, gene, display_group)

fill_colors <- c(TN = "#AA4643", `non-TN` = "#89A54E", Normal = "#4572A7")

build_plot <- function(genes_for_plot) {
  data_subset <- plot_data[gene %chin% genes_for_plot]
  data_subset[, gene := factor(as.character(gene), levels = genes_for_plot)]
  stats_subset <- statistics[gene %chin% genes_for_plot]
  stats_subset[, gene := factor(as.character(gene), levels = genes_for_plot)]

  ggplot(data_subset, aes(x = display_group, y = log2TPM1, fill = display_group)) +
    geom_boxplot(
      width = 0.56,
      linewidth = 0.50,
      colour = "black",
      outlier.shape = 16,
      outlier.size = 0.58,
      outlier.alpha = 0.50
    ) +
    geom_segment(
      data = stats_subset,
      aes(x = x_start, xend = x_end, y = bracket_y, yend = bracket_y),
      inherit.aes = FALSE,
      linewidth = 0.42,
      colour = "black"
    ) +
    geom_segment(
      data = stats_subset,
      aes(x = x_start, xend = x_start, y = bracket_y, yend = bracket_y - tip_height),
      inherit.aes = FALSE,
      linewidth = 0.42,
      colour = "black"
    ) +
    geom_segment(
      data = stats_subset,
      aes(x = x_end, xend = x_end, y = bracket_y, yend = bracket_y - tip_height),
      inherit.aes = FALSE,
      linewidth = 0.42,
      colour = "black"
    ) +
    geom_text(
      data = stats_subset,
      aes(x = (x_start + x_end) / 2, y = label_y, label = significance),
      inherit.aes = FALSE,
      family = "Arial",
      fontface = "bold",
      size = 10 / ggplot2::.pt,
      vjust = 0
    ) +
    facet_wrap(~gene, nrow = 1, scales = "free_y") +
    scale_fill_manual(values = fill_colors, drop = FALSE) +
    scale_x_discrete(drop = FALSE) +
    scale_y_continuous(
      name = expression(log[2] * "(TPM + 1)"),
      expand = expansion(mult = c(0.03, 0.055))
    ) +
    labs(x = NULL) +
    coord_cartesian(clip = "off") +
    theme_classic(base_size = 10, base_family = "Arial") +
    theme(
      axis.line = element_line(linewidth = 0.45, colour = "black"),
      axis.ticks = element_line(linewidth = 0.42, colour = "black"),
      axis.title.y = element_text(size = 11, margin = margin(r = 7)),
      axis.text.x = element_text(size = 10, colour = "black", margin = margin(t = 4)),
      axis.text.y = element_text(size = 9.5, colour = "black"),
      strip.background = element_rect(fill = "white", colour = "black", linewidth = 0.50),
      strip.text = element_text(size = 10.5, face = "bold", margin = margin(t = 4, b = 4)),
      panel.spacing.x = grid::unit(3.0, "mm"),
      legend.position = "none",
      plot.margin = margin(t = 5, r = 5, b = 5, l = 5)
    )
}

save_plot_bundle <- function(plot_object, basename, width_mm, height_mm) {
  width_in <- width_mm / 25.4
  height_in <- height_mm / 25.4

  grDevices::cairo_pdf(
    file.path(output_dir, paste0(basename, ".pdf")),
    width = width_in,
    height = height_in,
    family = "Arial"
  )
  print(plot_object)
  dev.off()

  svglite::svglite(
    file.path(output_dir, paste0(basename, ".svg")),
    width = width_in,
    height = height_in,
    system_fonts = list(sans = "Arial")
  )
  print(plot_object)
  dev.off()

  ragg::agg_png(
    file.path(output_dir, paste0(basename, ".png")),
    width = width_in,
    height = height_in,
    units = "in",
    res = 600,
    background = "white"
  )
  print(plot_object)
  dev.off()

  ragg::agg_tiff(
    file.path(output_dir, paste0(basename, ".tiff")),
    width = width_in,
    height = height_in,
    units = "in",
    res = 600,
    compression = "lzw",
    background = "white"
  )
  print(plot_object)
  dev.off()
}

message("Building and exporting the three-gene main Figure2-B...")
main_plot <- build_plot(main_genes)
save_plot_bundle(main_plot, prefix, width_mm = 183, height_mm = 84)

message("Building and exporting the two-gene supplementary figure...")
supplementary_basename <- paste0(prefix, "_supplementary_RPP38_NOP2")
supplementary_plot <- build_plot(supplementary_genes)
save_plot_bundle(supplementary_plot, supplementary_basename, width_mm = 122, height_mm = 84)

statistics_export <- statistics[, !c(
  "panel_min", "panel_max", "panel_range", "annotation_step",
  "bracket_y", "tip_height", "label_y", "x_start", "x_end", "annotation_level"
)]
fwrite(
  statistics_export,
  file.path(output_dir, paste0(prefix, "_Wilcoxon_TN_vs_non-TN_and_Normal.tsv")),
  sep = "\t"
)
fwrite(
  statistics_export[figure_set == "Main"],
  file.path(output_dir, paste0(prefix, "_main_statistics.tsv")),
  sep = "\t"
)
fwrite(
  statistics_export[figure_set == "Supplementary"],
  file.path(output_dir, paste0(prefix, "_supplementary_statistics.tsv")),
  sep = "\t"
)
fwrite(
  plot_data,
  file.path(output_dir, paste0(prefix, "_source_data_all_five_genes.tsv.gz")),
  sep = "\t"
)
fwrite(
  plot_data[gene %chin% main_genes],
  file.path(output_dir, paste0(prefix, "_main_source_data.tsv.gz")),
  sep = "\t"
)
fwrite(
  plot_data[gene %chin% supplementary_genes],
  file.path(output_dir, paste0(prefix, "_supplementary_source_data.tsv.gz")),
  sep = "\t"
)
fwrite(
  group_summary,
  file.path(output_dir, paste0(prefix, "_group_summary.tsv")),
  sep = "\t"
)

specifications <- data.table(
  item = c(
    "figure_claim", "archetype", "main_genes", "supplementary_genes",
    "expression_transform", "within_panel_order", "non-TN_definition",
    "planned_comparisons", "test", "multiple_testing", "star_thresholds",
    "boxplot_definition", "group_colors", "sample_counts", "font",
    "main_dimensions_mm", "supplementary_dimensions_mm", "resolution"
  ),
  value = c(
    "Compare gene-matched expression in TN tumors against non-TN tumors and Normal breast tissue",
    "quantitative grid with a three-gene main panel and two-gene supplementary panel",
    paste(main_genes, collapse = ", "),
    paste(supplementary_genes, collapse = ", "),
    "log2(TPM+1)",
    "TN, non-TN, Normal",
    "HR and HER2 primary tumors pooled at the sample level",
    "TN vs non-TN and TN vs Normal within every gene",
    "two-sided Wilcoxon rank-sum test; exact=FALSE; continuity correction=FALSE",
    "Benjamini-Hochberg across all ten planned comparisons",
    "ns: FDR>=0.05; *: <0.05; **: <0.01; ***: <0.001; ****: <0.0001",
    "median, IQR, whiskers to 1.5xIQR; points are boxplot outliers",
    "TN #AA4643; non-TN #89A54E; Normal #4572A7",
    "TN n=177, non-TN n=827 (HR 664 + HER2 163), Normal n=112 per gene",
    "Arial; >=9.5 pt tick labels and >=10 pt principal text",
    "183 x 84",
    "122 x 84",
    "PNG/TIFF 600 dpi; PDF/SVG vector"
  )
)
fwrite(specifications, file.path(output_dir, paste0(prefix, "_figure_specifications.tsv")), sep = "\t")

expected_display_counts <- c(TN = 177L, `non-TN` = 827L, Normal = 112L)
display_count_check <- plot_data[, .N, by = .(gene = as.character(gene), display_group = as.character(display_group))]
display_count_check[, expected_N := expected_display_counts[display_group]]
main_outputs <- file.path(output_dir, paste0(prefix, c(".pdf", ".svg", ".png", ".tiff")))
supp_outputs <- file.path(output_dir, paste0(supplementary_basename, c(".pdf", ".svg", ".png", ".tiff")))
validation <- data.table(
  check = c(
    "five genes mapped one-to-one", "expected sample counts in every gene/group",
    "group order is TN, non-TN, Normal", "ten planned tests completed",
    "each gene has both required comparisons", "all P/FDR values finite",
    "main figure contains three requested genes", "supplement contains two requested genes",
    "all main and supplementary outputs created"
  ),
  observed = c(
    as.character(nrow(gene_map) == 5L && uniqueN(gene_map$gene) == 5L),
    as.character(all(display_count_check$N == display_count_check$expected_N)),
    as.character(identical(levels(plot_data$display_group), display_groups)),
    as.character(nrow(statistics)),
    as.character(all(statistics[, uniqueN(comparison), by = gene]$V1 == 2L)),
    as.character(all(is.finite(statistics$raw_P)) && all(is.finite(statistics$BH_FDR))),
    as.character(setequal(as.character(unique(plot_data[figure_set == "Main", gene])), main_genes)),
    as.character(setequal(as.character(unique(plot_data[figure_set == "Supplementary", gene])), supplementary_genes)),
    as.character(all(file.exists(c(main_outputs, supp_outputs))))
  ),
  expected = c("TRUE", "TRUE", "TRUE", "10", "TRUE", "TRUE", "TRUE", "TRUE", "TRUE")
)
validation[, pass := observed == expected]
fwrite(validation, file.path(output_dir, paste0(prefix, "_validation_checks.tsv")), sep = "\t")

if (!all(validation$pass)) stop("One or more Figure2-B validation checks failed")
message("Figure2-B main and supplementary figures are complete.")
print(statistics_export[, .(
  figure_set, gene, comparison, TN_n, reference_n,
  TN_median_log2TPM1, reference_median_log2TPM1,
  median_difference_log2TPM1, direction, raw_P, BH_FDR, significance
)])
