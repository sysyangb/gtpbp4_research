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
  library(pheatmap)
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
options(mc.cores = 1L, width = 160)
data.table::setDTthreads(1L)
set.seed(20260831)

input_dir <- file.path(project_dir, "TCGA_BRCA_four_group_RBP_heatmap_ESR1_ERBB2_final")
input_source_dir <- file.path(input_dir, "source_data")
marker_path <- file.path(input_source_dir, "selected_marker_row_metadata_clustered.tsv")
sample_path <- file.path(input_source_dir, "sample_column_metadata.tsv")
zscore_path <- file.path(input_source_dir, "selected_marker_row_zscore_matrix_clustered.tsv.gz")
required <- c(marker_path, sample_path, zscore_path)
if (any(!file.exists(required))) {
  stop("Missing input(s): ", paste(required[!file.exists(required)], collapse = ", "))
}

marker_groups <- c("TN", "HR", "HER2")
sample_groups <- c("TN", "HR", "HER2", "NORMAL")
group_colors <- c(TN = "#89A54E", HR = "#71588F", HER2 = "#4198AF", NORMAL = "#DB843D")

message("Loading the validated four-group heatmap source data...")
marker_meta_all <- fread(marker_path, showProgress = FALSE)
sample_meta <- fread(sample_path, showProgress = FALSE)
setorder(sample_meta, heatmap_column_order)

z_connection <- gzfile(zscore_path, open = "rt")
z_dt <- as.data.table(read.delim(z_connection, check.names = FALSE, stringsAsFactors = FALSE))
close(z_connection)
gene_ids <- z_dt[[1L]]
heatmap_z_all <- as.matrix(z_dt[, -1L])
storage.mode(heatmap_z_all) <- "double"
rownames(heatmap_z_all) <- gene_ids

# Requested revision: remove NORMAL-specific marker genes, while retaining every
# NORMAL sample column as the fourth comparison group.
marker_meta <- marker_meta_all[marker_group %chin% marker_groups]
heatmap_z <- heatmap_z_all[marker_meta$Gene_ID, sample_meta$sample_id, drop = FALSE]

if (!identical(unique(marker_meta$marker_group), marker_groups)) stop("Marker group order changed")
if (!identical(unique(sample_meta$condition), sample_groups)) stop("Sample group order changed")
if (any(marker_meta$marker_group == "NORMAL")) stop("NORMAL-specific marker genes were not removed")
if (any(!is.finite(heatmap_z))) stop("Non-finite heatmap values found")
if (!setequal(marker_meta$Gene_ID, rownames(heatmap_z))) stop("Marker metadata and matrix genes differ")
if (!identical(sample_meta$sample_id, colnames(heatmap_z))) stop("Sample metadata and matrix columns differ")

message("Reconstructing complete-linkage row clustering within fixed TN, HR and HER2 blocks...")
group_dendrograms <- lapply(marker_groups, function(group) {
  ids <- marker_meta[marker_group == group, Gene_ID]
  as.dendrogram(hclust(dist(heatmap_z[ids, , drop = FALSE]), method = "complete"))
})
combine_dendrograms <- function(left, right) {
  merge(left, right, height = max(attr(left, "height"), attr(right, "height")) * 1.05 + 1e-8)
}
row_cluster <- as.hclust(Reduce(combine_dendrograms, group_dendrograms))
clustered_gene_ids <- row_cluster$labels[row_cluster$order]
clustered_groups <- marker_meta$marker_group[match(clustered_gene_ids, marker_meta$Gene_ID)]
if (!identical(unique(clustered_groups), marker_groups)) stop("Row clustering changed fixed marker-group order")

row_group_counts <- vapply(marker_groups, function(g) sum(marker_meta$marker_group == g), integer(1))
column_group_counts <- vapply(sample_groups, function(g) sum(sample_meta$condition == g), integer(1))
row_gaps <- head(cumsum(row_group_counts), -1L)
column_gaps <- head(cumsum(column_group_counts), -1L)
cluster_gaps <- as.integer(pheatmap:::find_gaps(row_cluster, length(marker_groups)))
if (!identical(cluster_gaps, as.integer(row_gaps))) stop("Cluster cuts do not reproduce marker-group boundaries")

annotation_col <- data.frame(
  `Sample type` = factor(sample_meta$condition, levels = sample_groups),
  row.names = sample_meta$sample_id,
  check.names = FALSE
)
annotation_colors <- list(`Sample type` = group_colors)
heatmap_palette <- colorRampPalette(c("#2166AC", "#F7F7F7", "#B2182B"))(101)
heatmap_breaks <- seq(-2.5, 2.5, length.out = length(heatmap_palette) + 1L)

# Keep narrow 2-bigpoint separators on both axes.
two_point_find_coordinates <- local({
  function(n, gaps, m = 1:n) {
    if (length(gaps) == 0) {
      return(list(coord = grid::unit(m / n, "npc"), size = grid::unit(1 / n, "npc")))
    }
    if (max(gaps) > n) stop("Gaps do not match matrix size")
    gap_width <- 2
    size <- (1 / n) * (grid::unit(1, "npc") - length(gaps) * grid::unit(gap_width, "bigpts"))
    gaps2 <- apply(sapply(gaps, function(gap, x) x > gap, m), 1, sum)
    coord <- m * size + gaps2 * grid::unit(gap_width, "bigpts")
    list(coord = coord, size = size)
  }
})
assignInNamespace("find_coordinates", two_point_find_coordinates, ns = "pheatmap")

message("Building Figure1A pheatmap at 190 x 90 mm with 10-pt text...")
ph <- pheatmap(
  heatmap_z,
  color = heatmap_palette,
  breaks = heatmap_breaks,
  cluster_rows = row_cluster,
  cutree_rows = length(marker_groups),
  cluster_cols = FALSE,
  gaps_col = column_gaps,
  annotation_col = annotation_col,
  annotation_colors = annotation_colors,
  annotation_names_col = FALSE,
  show_rownames = FALSE,
  show_colnames = FALSE,
  border_color = NA,
  legend = TRUE,
  legend_breaks = c(-2.5, 0, 2.5),
  legend_labels = c("−2.5", "0", "2.5"),
  annotation_legend = TRUE,
  drop_levels = FALSE,
  treeheight_row = 34,
  treeheight_col = 0,
  fontsize = 10,
  fontsize_row = 10,
  fontsize_col = 10,
  fontfamily = "Arial",
  main = NA,
  silent = TRUE
)

collect_font_sizes <- function(grob) {
  out <- numeric()
  if (inherits(grob, "text") && !is.null(grob$gp$fontsize)) out <- c(out, as.numeric(grob$gp$fontsize))
  if (!is.null(grob$children)) out <- c(out, unlist(lapply(grob$children, collect_font_sizes), use.names = FALSE))
  if (!is.null(grob$grobs)) out <- c(out, unlist(lapply(grob$grobs, collect_font_sizes), use.names = FALSE))
  out
}
font_sizes <- collect_font_sizes(ph$gtable)
if (!length(font_sizes) || min(font_sizes) < 10) stop("Detected plotted text below 10 pt")

draw_heatmap <- function() {
  grid.newpage()
  grid.draw(ph$gtable)
}

stem <- file.path(bundle_dir, "Figure1A_TCGA_BRCA_TN_HR_HER2_markers_with_NORMAL_reference_heatmap")
width_mm <- 190
height_mm <- 90
width_in <- width_mm / 25.4
height_in <- height_mm / 25.4

ragg::agg_png(paste0(stem, ".png"), width = width_in, height = height_in, units = "in", res = 600, background = "white")
draw_heatmap()
dev.off()
ragg::agg_tiff(paste0(stem, ".tiff"), width = width_in, height = height_in, units = "in", res = 600, compression = "lzw", background = "white")
draw_heatmap()
dev.off()
grDevices::cairo_pdf(paste0(stem, ".pdf"), width = width_in, height = height_in, family = "Arial")
draw_heatmap()
dev.off()
svglite::svglite(paste0(stem, ".svg"), width = width_in, height = height_in, bg = "white")
draw_heatmap()
dev.off()

# Short aliases for manuscript assembly.
file.copy(paste0(stem, ".png"), file.path(bundle_dir, "Figure1A.png"), overwrite = TRUE)
file.copy(paste0(stem, ".pdf"), file.path(bundle_dir, "Figure1A.pdf"), overwrite = TRUE)

message("Writing Figure1A-prefixed source data and QA files...")
clustered_meta <- marker_meta[match(clustered_gene_ids, Gene_ID)]
clustered_meta[, Figure1A_heatmap_row_order := seq_len(.N)]
setcolorder(clustered_meta, c("Figure1A_heatmap_row_order", setdiff(names(clustered_meta), "Figure1A_heatmap_row_order")))
fwrite(clustered_meta, file.path(bundle_dir, "Figure1A_marker_row_metadata.tsv"), sep = "\t", na = "NA")
fwrite(sample_meta, file.path(bundle_dir, "Figure1A_sample_column_metadata.tsv"), sep = "\t", na = "NA")
z_out <- as.data.table(heatmap_z[clustered_gene_ids, , drop = FALSE], keep.rownames = "Gene_ID")
fwrite(z_out, file.path(bundle_dir, "Figure1A_heatmap_row_zscore_matrix.tsv.gz"), sep = "\t", compress = "gzip")

representative_summary <- rbindlist(lapply(c("ESR1", "ERBB2"), function(symbol) {
  gene_id <- marker_meta[Gene_symbol == symbol, Gene_ID]
  data.table(condition = sample_meta$condition, row_z = as.numeric(heatmap_z[gene_id, sample_meta$sample_id]))[
    , .(mean_row_z = mean(row_z), median_row_z = median(row_z)), by = condition
  ][, Gene_symbol := symbol]
}))
setcolorder(representative_summary, c("Gene_symbol", "condition", "mean_row_z", "median_row_z"))
fwrite(representative_summary, file.path(bundle_dir, "Figure1A_ESR1_ERBB2_group_row_z_summary.tsv"), sep = "\t")

marker_summary <- rbindlist(lapply(marker_groups, function(group) {
  data.table(
    marker_group = group,
    RBP_gene_count = marker_meta[marker_group == group & Is_RBP == "Yes", .N],
    representative_non_RBP_count = marker_meta[marker_group == group & Is_RBP != "Yes", .N],
    total_heatmap_gene_count = marker_meta[marker_group == group, .N]
  )
}))
sample_summary <- data.table(
  sample_group = sample_groups,
  sample_count = vapply(sample_groups, function(group) sample_meta[condition == group, .N], integer(1))
)
fwrite(marker_summary, file.path(bundle_dir, "Figure1A_marker_group_counts.tsv"), sep = "\t")
fwrite(sample_summary, file.path(bundle_dir, "Figure1A_sample_group_counts.tsv"), sep = "\t")

spec <- data.table(
  property = c("width_mm", "height_mm", "aspect_ratio", "minimum_text_pt", "row_gap_bigpoints", "column_gap_bigpoints", "png_dpi", "tiff_dpi"),
  value = c(width_mm, height_mm, width_mm / height_mm, min(font_sizes), 2, 2, 600, 600)
)
fwrite(spec, file.path(bundle_dir, "Figure1A_figure_specifications.tsv"), sep = "\t")

expected_marker_counts <- c(TN = 131L, HR = 7L, HER2 = 21L)
expected_sample_counts <- c(TN = 177L, HR = 664L, HER2 = 163L, NORMAL = 112L)
export_paths <- paste0(stem, c(".png", ".tiff", ".pdf", ".svg"))
qa <- data.table(
  check = c(
    "NORMAL_specific_marker_genes_are_excluded",
    "marker_gene_count_is_159",
    "marker_group_counts_are_131_7_21",
    "all_1116_samples_are_retained",
    "sample_group_counts_are_177_664_163_112",
    "marker_group_order_is_TN_HR_HER2",
    "sample_group_order_is_TN_HR_HER2_NORMAL",
    "within_group_row_clustering_preserves_blocks",
    "two_row_and_three_column_separators_are_present",
    "sample_type_annotation_legend_is_present",
    "requested_group_colors_are_applied",
    "all_drawn_text_is_at_least_10_pt",
    "A4_full_width_flat_dimensions_are_used",
    "samples_and_genes_are_not_labeled",
    "all_exports_exist",
    "all_new_bundle_files_start_with_Figure1A"
  ),
  passed = c(
    !any(marker_meta$marker_group == "NORMAL"),
    nrow(marker_meta) == 159L,
    identical(unname(row_group_counts), unname(expected_marker_counts)),
    nrow(sample_meta) == 1116L,
    identical(unname(column_group_counts), unname(expected_sample_counts)),
    identical(unique(marker_meta$marker_group), marker_groups),
    identical(unique(sample_meta$condition), sample_groups),
    identical(unique(clustered_groups), marker_groups),
    length(row_gaps) == 2L && length(column_gaps) == 3L,
    sum(ph$gtable$layout$name == "annotation_legend") == 1L,
    identical(unname(group_colors), c("#89A54E", "#71588F", "#4198AF", "#DB843D")),
    min(font_sizes) >= 10,
    width_mm == 190 && height_mm == 90,
    TRUE,
    all(file.exists(export_paths)),
    all(startsWith(basename(list.files(bundle_dir, full.names = TRUE, all.files = FALSE)), "Figure1A") |
          basename(list.files(bundle_dir, full.names = TRUE, all.files = FALSE)) == "version0")
  )
)
fwrite(qa, file.path(bundle_dir, "Figure1A_validation_checks.tsv"), sep = "\t")
if (!all(qa$passed)) stop("Figure1A validation failed: ", paste(qa[passed == FALSE, check], collapse = ", "))


message("Completed revised Figure1A: ", stem)
print(marker_summary)
print(sample_summary)
print(spec)
print(qa)
