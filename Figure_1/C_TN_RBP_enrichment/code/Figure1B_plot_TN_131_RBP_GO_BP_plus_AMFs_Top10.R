#!/usr/bin/env Rscript


args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
script_path <- if (length(file_arg)) sub("^--file=", "", file_arg[[1]]) else "."
script_path <- gsub("~\\+~", " ", script_path)
bundle_dir <- normalizePath(dirname(script_path), mustWork = TRUE)
project_dir <- normalizePath(file.path(bundle_dir, ".."), mustWork = TRUE)

local_lib_candidates <- c(
  file.path(bundle_dir, ".Rlib"),
  file.path(project_dir, "..", "TCGA", ".Rlib")
)
local_lib_candidates <- local_lib_candidates[dir.exists(local_lib_candidates)]
if (length(local_lib_candidates)) .libPaths(c(local_lib_candidates, .libPaths()))

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(scales)
  library(stringr)
  library(grid)
  library(ragg)
  library(svglite)
})

Sys.setenv(
  OMP_NUM_THREADS = "1",
  OPENBLAS_NUM_THREADS = "1",
  MKL_NUM_THREADS = "1",
  VECLIB_MAXIMUM_THREADS = "1"
)
options(mc.cores = 1L, width = 170)
data.table::setDTthreads(1L)
set.seed(20260829)

input_dir <- file.path(project_dir, "TN_131_RBP_GO_BP_plus_AMFs_enrichment")
input_results_dir <- file.path(input_dir, "results")
all_pathways_path <- file.path(input_results_dir, "TN_131_RBP_GO_BP_plus_AMFs_ORA_all_pathways.tsv")
amf_hits_path <- file.path(input_results_dir, "TN_131_RBP_AMFs_hit_genes.tsv")
required <- c(all_pathways_path, amf_hits_path)
if (any(!file.exists(required))) {
  stop("Missing input(s): ", paste(required[!file.exists(required)], collapse = ", "))
}

message("Loading the jointly corrected GO BP + AMFs enrichment results...")
result <- fread(all_pathways_path, showProgress = FALSE)
if (sum(result$term_id == "CUSTOM_AMFS") != 1L) stop("AMFs term missing or duplicated")
if (!identical(result$rank_raw_P, seq_len(nrow(result)))) stop("Input results are not in raw-P rank order")

plot_dt <- copy(result[candidate_count >= 1L & is.finite(pvalue)])
setorder(plot_dt, pvalue, p_adjust_BH, -fold_enrichment, -candidate_count, term_id)
plot_dt <- head(plot_dt, 10L)
plot_dt[, Figure1B_Top10_rank := seq_len(.N)]
plot_dt[, minus_log10_pvalue := -log10(pmax(pvalue, .Machine$double.xmin))]
plot_dt[, label := str_wrap(term_name, width = 43)]
plot_dt[, label := make.unique(label, sep = " ")]
plot_dt[, pathway_factor := factor(label, levels = rev(label))]

min_count <- min(plot_dt$candidate_count)
max_count <- max(plot_dt$candidate_count)
max_point_size_mm <- 7.2
size_breaks <- pretty(c(min_count, max_count), n = 3)
size_breaks <- size_breaks[size_breaks >= min_count & size_breaks <= max_count]
if (!length(size_breaks)) size_breaks <- sort(unique(c(min_count, max_count)))

p <- ggplot(
  plot_dt,
  aes(
    x = minus_log10_pvalue,
    y = pathway_factor,
    size = candidate_count,
    fill = fold_enrichment
  )
) +
  geom_point(shape = 21, colour = "#843B38", stroke = 0.75, alpha = 1) +
  scale_size_area(
    max_size = max_point_size_mm,
    breaks = size_breaks,
    limits = c(0, max_count),
    name = "Gene count"
  ) +
  scale_fill_gradient(
    low = "#F8E8E5",
    high = "#B34D49",
    name = "Enrichment score",
    breaks = breaks_pretty(n = 4)
  ) +
  scale_x_continuous(
    name = expression(-log[10](italic(P)~value)),
    breaks = breaks_pretty(n = 5),
    expand = expansion(mult = c(0.10, 0.08))
  ) +
  scale_y_discrete(name = NULL, expand = expansion(add = c(0.42, 0.42))) +
  guides(
    size = guide_legend(
      order = 1,
      title.position = "top",
      override.aes = list(fill = "#B34D49", colour = "#843B38", alpha = 1)
    ),
    fill = guide_colourbar(
      order = 2,
      title.position = "top",
      barheight = grid::unit(28, "mm"),
      barwidth = grid::unit(6, "mm"),
      ticks = TRUE,
      frame.colour = NA
    )
  ) +
  theme_classic(base_size = 10.5, base_family = "Arial") +
  theme(
    axis.line = element_line(linewidth = 0.7, colour = "#111111"),
    axis.ticks = element_line(linewidth = 0.65, colour = "#111111"),
    axis.ticks.length = grid::unit(2.0, "mm"),
    axis.text.x = element_text(size = 10.2, colour = "#222222"),
    axis.text.y = element_text(size = 10.2, colour = "#222222", lineheight = 0.88),
    axis.title.x = element_text(size = 10.8, margin = margin(t = 7)),
    legend.position = "right",
    legend.box = "vertical",
    legend.box.just = "left",
    legend.title = element_text(size = 10.5, face = "bold"),
    legend.text = element_text(size = 10.0),
    legend.key.height = grid::unit(6.5, "mm"),
    legend.spacing.y = grid::unit(3.5, "mm"),
    legend.margin = margin(l = 3),
    plot.margin = margin(7, 8, 9, 7, unit = "mm"),
    panel.grid = element_blank()
  )

collect_font_sizes <- function(grob) {
  out <- numeric()
  if (inherits(grob, "text") && !is.null(grob$gp$fontsize)) {
    out <- c(out, as.numeric(grob$gp$fontsize))
  }
  if (!is.null(grob$children)) {
    out <- c(out, unlist(lapply(grob$children, collect_font_sizes), use.names = FALSE))
  }
  if (!is.null(grob$grobs)) {
    out <- c(out, unlist(lapply(grob$grobs, collect_font_sizes), use.names = FALSE))
  }
  out
}
grDevices::pdf(NULL)
plot_grob <- ggplotGrob(p)
grDevices::dev.off()
font_sizes <- collect_font_sizes(plot_grob)
if (!length(font_sizes) || min(font_sizes) < 10) {
  stop("Detected plotted text below 10 pt; minimum was ", min(font_sizes))
}

stem <- file.path(bundle_dir, "Figure1B_TN_131_RBP_GO_BP_plus_AMFs_Top10_dotplot")
width_mm <- 185
height_mm <- 110
ggsave(paste0(stem, ".png"), p, device = ragg::agg_png, width = width_mm, height = height_mm, units = "mm", dpi = 600, bg = "white")
ggsave(paste0(stem, ".tiff"), p, device = ragg::agg_tiff, width = width_mm, height = height_mm, units = "mm", dpi = 600, compression = "lzw", bg = "white")
ggsave(paste0(stem, ".pdf"), p, device = cairo_pdf, width = width_mm, height = height_mm, units = "mm", family = "Arial", bg = "white")
ggsave(paste0(stem, ".svg"), p, device = svglite::svglite, width = width_mm, height = height_mm, units = "mm", bg = "white")

message("Writing Figure1B-prefixed source data and result files...")
fwrite(
  plot_dt[, setdiff(names(plot_dt), "pathway_factor"), with = FALSE],
  file.path(bundle_dir, "Figure1B_GO_BP_plus_AMFs_Top10_source_data.tsv"),
  sep = "\t",
  na = "NA"
)
fwrite(result, file.path(bundle_dir, "Figure1B_GO_BP_plus_AMFs_all_pathways.tsv"), sep = "\t", na = "NA")
amf_hits <- fread(amf_hits_path, showProgress = FALSE)
fwrite(amf_hits, file.path(bundle_dir, "Figure1B_AMFs_hit_genes.tsv"), sep = "\t", na = "NA")
size_scale_spec <- data.table(
  property = c(
    "mapping",
    "minimum_observed_gene_count",
    "maximum_observed_gene_count",
    "legend_breaks",
    "maximum_point_size_mm",
    "special_AMFs_point_style"
  ),
  value = c(
    "point area proportional to gene count",
    min_count,
    max_count,
    paste(size_breaks, collapse = "/"),
    max_point_size_mm,
    "none"
  )
)
fwrite(size_scale_spec, file.path(bundle_dir, "Figure1B_point_size_scale.tsv"), sep = "\t")

amf_result <- result[term_id == "CUSTOM_AMFS"]
summary_dt <- data.table(
  metric = c(
    "Input TN RBP genes",
    "Mapped TN RBP genes",
    "Combined mapped background genes",
    "Eligible GO BP plus AMFs terms",
    "Displayed pathways",
    "Observed Top10 gene-count range",
    "Maximum point size mm",
    "AMFs natural rank",
    "AMFs hits",
    "AMFs fold enrichment",
    "AMFs raw P",
    "AMFs BH FDR",
    "Minimum plotted text size pt",
    "Figure width mm",
    "Figure height mm"
  ),
  value = c(
    amf_result$query_input_count,
    amf_result$query_mapped_count,
    amf_result$mapped_background_count,
    nrow(result),
    nrow(plot_dt),
    paste0(min_count, "-", max_count),
    max_point_size_mm,
    amf_result$rank_raw_P,
    amf_result$candidate_count,
    signif(amf_result$fold_enrichment, 8),
    format(amf_result$pvalue, scientific = TRUE, digits = 8),
    format(amf_result$p_adjust_BH, scientific = TRUE, digits = 8),
    min(font_sizes),
    width_mm,
    height_mm
  )
)
fwrite(summary_dt, file.path(bundle_dir, "Figure1B_analysis_summary.tsv"), sep = "\t")

qa <- data.table(
  check = c(
    "top10_contains_exactly_10_terms",
    "top10_is_natural_raw_P_order",
    "AMFs_is_not_forced_and_remains_rank_1",
    "AMFs_has_no_special_point_style",
    "point_area_scale_uses_observed_counts_and_smaller_maximum",
    "BH_FDR_is_retained_from_joint_analysis",
    "all_drawn_text_is_at_least_10_pt",
    "no_extra_plot_title",
    "all_exports_exist",
    "all_bundle_files_start_with_Figure1B"
  ),
  passed = c(
    nrow(plot_dt) == 10L,
    identical(plot_dt$term_id, head(result[candidate_count >= 1L]$term_id, 10L)),
    amf_result$rank_raw_P == 1L && plot_dt$term_id[1L] == "CUSTOM_AMFS",
    length(p$layers) == 1L,
    max_point_size_mm == 7.2 && all(size_breaks >= min_count & size_breaks <= max_count),
    all(plot_dt$p_adjust_BH == result$p_adjust_BH[match(plot_dt$term_id, result$term_id)]),
    min(font_sizes) >= 10,
    is.null(p$labels$title),
    all(file.exists(paste0(stem, c(".png", ".tiff", ".pdf", ".svg")))),
    all(startsWith(basename(list.files(bundle_dir, full.names = TRUE)), "Figure1B"))
  )
)
fwrite(qa, file.path(bundle_dir, "Figure1B_validation_checks.tsv"), sep = "\t")
if (!all(qa$passed)) stop("Figure1B validation failed")


message("Completed Figure1B: ", stem)
print(plot_dt[, .(Figure1B_Top10_rank, term_id, term_name, candidate_count, fold_enrichment, pvalue, p_adjust_BH)])
print(summary_dt)
print(qa)
