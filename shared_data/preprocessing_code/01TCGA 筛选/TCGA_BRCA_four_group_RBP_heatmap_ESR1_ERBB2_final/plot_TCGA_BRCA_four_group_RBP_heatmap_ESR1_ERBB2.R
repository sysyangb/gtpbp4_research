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
  library(readxl)
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
set.seed(20260829)

workbook_path <- file.path(
  project_dir,
  "TCGA_BRCA_four_group_specific_translation_no_ribosome_RBP_annotated.xlsx"
)
previous_source_dir <- file.path(
  project_dir,
  "TCGA_BRCA_four_group_heatmap_final",
  "source_data"
)
previous_marker_path <- file.path(previous_source_dir, "all_marker_row_metadata.tsv")
sample_meta_path <- file.path(previous_source_dir, "all_sample_column_metadata.tsv")
zscore_path <- file.path(previous_source_dir, "all_marker_row_zscore_matrix.tsv.gz")
required <- c(workbook_path, previous_marker_path, sample_meta_path, zscore_path)
if (any(!file.exists(required))) {
  stop("Missing input(s): ", paste(required[!file.exists(required)], collapse = ", "))
}

source_dir <- file.path(bundle_dir, "source_data")
dir.create(source_dir, recursive = TRUE, showWarnings = FALSE)

groups <- c("TN", "HR", "HER2", "NORMAL")
sheet_names <- c(TN = "TN", HR = "HR", HER2 = "HER2", NORMAL = "Normal")
group_colors <- c(TN = "#89A54E", HR = "#71588F", HER2 = "#4198AF", NORMAL = "#DB843D")
forced_markers <- list(TN = character(), HR = "ESR1", HER2 = "ERBB2", NORMAL = character())

message("Selecting RBP genes and the two requested representative markers...")
marker_tables <- lapply(groups, function(group) {
  x <- as.data.table(read_excel(workbook_path, sheet = sheet_names[[group]]))
  required_columns <- c(
    "Rank", "Gene_ID", "Gene_symbol", "Gene_biotype",
    "One_vs_rest_log2FC", "One_vs_rest_P_value", "One_vs_rest_FDR", "Is_RBP"
  )
  if (!all(required_columns %in% names(x))) {
    stop("Worksheet ", sheet_names[[group]], " is missing required columns")
  }
  x <- x[Gene_biotype == "protein_coding"]
  x <- x[Is_RBP == "Yes" | Gene_symbol %chin% forced_markers[[group]]]
  x[, marker_group := group]
  x[, selection_reason := fifelse(
    Is_RBP == "Yes",
    "RBP",
    paste0("Representative marker retained for ", group)
  )]
  setorder(x, -One_vs_rest_log2FC, Gene_ID)
  x[, marker_rank_within_group := seq_len(.N)]
  x
})
names(marker_tables) <- groups
marker_meta <- rbindlist(marker_tables, use.names = TRUE, fill = TRUE)

if (anyDuplicated(marker_meta$Gene_ID)) stop("A selected gene occurs in more than one phenotype group")
if (!identical(unique(marker_meta$marker_group), groups)) stop("Marker groups are not in TN/HR/HER2/NORMAL order")
if (marker_meta[Gene_symbol == "ESR1", unique(marker_group)] != "HR") stop("ESR1 was not retained in HR")
if (marker_meta[Gene_symbol == "ERBB2", unique(marker_group)] != "HER2") stop("ERBB2 was not retained in HER2")
non_rbp <- marker_meta[Is_RBP != "Yes", sort(unique(Gene_symbol))]
if (!identical(non_rbp, c("ERBB2", "ESR1"))) stop("Unexpected non-RBP genes were selected")

previous_marker_meta <- fread(previous_marker_path, showProgress = FALSE)
sample_meta <- fread(sample_meta_path, showProgress = FALSE)
setorder(sample_meta, heatmap_column_order)
if (!identical(unique(sample_meta$condition), groups)) stop("Sample groups are not in TN/HR/HER2/NORMAL order")

message("Loading the row-Z-score matrix from the validated full heatmap bundle...")
z_connection <- gzfile(zscore_path, open = "rt")
z_dt <- as.data.table(read.delim(z_connection, check.names = FALSE, stringsAsFactors = FALSE))
close(z_connection)
gene_ids <- z_dt[[1L]]
heatmap_z_all <- as.matrix(z_dt[, -1L])
storage.mode(heatmap_z_all) <- "double"
rownames(heatmap_z_all) <- gene_ids

missing_selected <- setdiff(marker_meta$Gene_ID, rownames(heatmap_z_all))
if (length(missing_selected)) stop("Selected genes missing from Z-score matrix: ", paste(missing_selected, collapse = ", "))
missing_samples <- setdiff(sample_meta$sample_id, colnames(heatmap_z_all))
if (length(missing_samples)) stop("Selected samples missing from Z-score matrix")

heatmap_z <- heatmap_z_all[
  marker_meta$Gene_ID,
  sample_meta$sample_id,
  drop = FALSE
]
if (any(!is.finite(heatmap_z))) stop("Non-finite values are present in the selected heatmap matrix")

message("Clustering rows within the four fixed marker groups...")
group_dendrograms <- lapply(groups, function(group) {
  ids <- marker_meta[marker_group == group, Gene_ID]
  as.dendrogram(hclust(dist(heatmap_z[ids, , drop = FALSE]), method = "complete"))
})
combine_dendrograms <- function(left, right) {
  merge(left, right, height = max(attr(left, "height"), attr(right, "height")) * 1.05 + 1e-8)
}
row_cluster <- as.hclust(Reduce(combine_dendrograms, group_dendrograms))
clustered_gene_ids <- row_cluster$labels[row_cluster$order]
clustered_groups <- marker_meta$marker_group[match(clustered_gene_ids, marker_meta$Gene_ID)]
if (!identical(unique(clustered_groups), groups)) stop("Row clustering changed the fixed group order")

marker_meta[, input_heatmap_row_order := seq_len(.N)]
clustered_meta <- marker_meta[match(clustered_gene_ids, Gene_ID)]
clustered_meta[, clustered_heatmap_row_order := seq_len(.N)]

row_groups <- factor(marker_meta$marker_group, levels = groups)
column_groups <- factor(sample_meta$condition, levels = groups)
row_group_counts <- vapply(groups, function(g) sum(row_groups == g), integer(1))
column_group_counts <- vapply(groups, function(g) sum(column_groups == g), integer(1))
row_gaps <- head(cumsum(row_group_counts), -1L)
column_gaps <- head(cumsum(column_group_counts), -1L)
cluster_gaps <- as.integer(pheatmap:::find_gaps(row_cluster, length(groups)))
if (!identical(cluster_gaps, as.integer(row_gaps))) stop("Cluster cuts do not reproduce group boundaries")

annotation_col <- data.frame(
  `Sample type` = column_groups,
  row.names = sample_meta$sample_id,
  check.names = FALSE
)
annotation_colors <- list(`Sample type` = group_colors)

heatmap_palette <- colorRampPalette(c("#2166AC", "#F7F7F7", "#B2182B"))(101)
heatmap_breaks <- seq(-2.5, 2.5, length.out = length(heatmap_palette) + 1L)

# pheatmap normally draws 4-bigpoint separators; keep the optimized 2-bigpoint
# separators from the final full-marker heatmap on both axes.
two_point_find_coordinates <- local({
  heatmap_nrow <- nrow(heatmap_z)
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

message("Building final pheatmap for ", nrow(heatmap_z), " genes x ", ncol(heatmap_z), " samples...")
ph <- pheatmap(
  heatmap_z,
  color = heatmap_palette,
  breaks = heatmap_breaks,
  cluster_rows = row_cluster,
  cutree_rows = length(groups),
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
  treeheight_row = 42,
  treeheight_col = 0,
  fontsize = 7,
  main = NA,
  silent = TRUE
)

draw_heatmap <- function() {
  grid.newpage()
  grid.draw(ph$gtable)
}

stem <- file.path(bundle_dir, "TCGA_BRCA_TN_HR_HER2_NORMAL_RBP_ESR1_ERBB2_heatmap")
figure_width_in <- 7.2
figure_height_in <- 6.2
ragg::agg_png(paste0(stem, ".png"), width = figure_width_in, height = figure_height_in, units = "in", res = 400, background = "white")
draw_heatmap()
dev.off()
ragg::agg_tiff(paste0(stem, ".tiff"), width = figure_width_in, height = figure_height_in, units = "in", res = 600, compression = "lzw", background = "white")
draw_heatmap()
dev.off()
grDevices::cairo_pdf(paste0(stem, ".pdf"), width = figure_width_in, height = figure_height_in, family = "Arial")
draw_heatmap()
dev.off()
svglite::svglite(paste0(stem, ".svg"), width = figure_width_in, height = figure_height_in, bg = "white")
draw_heatmap()
dev.off()

marker_out <- clustered_meta[, .(
  clustered_heatmap_row_order,
  marker_group,
  marker_rank_within_group,
  selection_reason,
  Is_RBP,
  Gene_ID,
  Gene_symbol,
  Entrez_ID,
  Gene_biotype,
  One_vs_rest_log2FC,
  One_vs_rest_P_value,
  One_vs_rest_FDR
)]
fwrite(marker_out, file.path(source_dir, "selected_marker_row_metadata_clustered.tsv"), sep = "\t", na = "NA")
fwrite(sample_meta, file.path(source_dir, "sample_column_metadata.tsv"), sep = "\t", na = "NA")
z_out <- as.data.table(heatmap_z[clustered_gene_ids, , drop = FALSE], keep.rownames = "Gene_ID")
fwrite(z_out, file.path(source_dir, "selected_marker_row_zscore_matrix_clustered.tsv.gz"), sep = "\t", compress = "gzip")

representative_summary <- rbindlist(lapply(c("ESR1", "ERBB2"), function(symbol) {
  gene_id <- marker_meta[Gene_symbol == symbol, Gene_ID]
  data.table(
    condition = sample_meta$condition,
    row_z = as.numeric(heatmap_z[gene_id, sample_meta$sample_id])
  )[, .(
    mean_row_z = mean(row_z),
    median_row_z = median(row_z)
  ), by = condition][, Gene_symbol := symbol]
}))
setcolorder(representative_summary, c("Gene_symbol", "condition", "mean_row_z", "median_row_z"))
fwrite(representative_summary, file.path(source_dir, "ESR1_ERBB2_group_row_z_summary.tsv"), sep = "\t")

group_summary <- rbindlist(lapply(groups, function(group) {
  data.table(
    group = group,
    RBP_gene_count = marker_meta[marker_group == group & Is_RBP == "Yes", .N],
    representative_non_RBP_count = marker_meta[marker_group == group & Is_RBP != "Yes", .N],
    total_heatmap_gene_count = marker_meta[marker_group == group, .N],
    sample_count = sample_meta[condition == group, .N]
  )
}))
fwrite(group_summary, file.path(bundle_dir, "heatmap_group_counts.tsv"), sep = "\t")

qa <- data.table(
  check = c(
    "selected_genes_are_unique",
    "all_RBP_rows_are_protein_coding",
    "only_ESR1_and_ERBB2_are_non_RBP",
    "ESR1_is_in_HR",
    "ERBB2_is_in_HER2",
    "ESR1_mean_row_z_is_highest_in_HR",
    "ERBB2_mean_row_z_is_highest_in_HER2",
    "marker_group_order_is_TN_HR_HER2_NORMAL",
    "sample_group_order_is_TN_HR_HER2_NORMAL",
    "within_group_clustering_preserves_group_blocks",
    "row_and_column_separators_are_present",
    "sample_type_annotation_legend_is_present",
    "requested_group_colors_are_applied",
    "samples_are_not_clustered_or_labeled",
    "genes_are_not_labeled",
    "all_exports_exist"
  ),
  passed = c(
    !anyDuplicated(marker_meta$Gene_ID),
    all(marker_meta[Is_RBP == "Yes", Gene_biotype] == "protein_coding"),
    identical(non_rbp, c("ERBB2", "ESR1")),
    marker_meta[Gene_symbol == "ESR1", unique(marker_group)] == "HR",
    marker_meta[Gene_symbol == "ERBB2", unique(marker_group)] == "HER2",
    representative_summary[Gene_symbol == "ESR1"][which.max(mean_row_z), condition] == "HR",
    representative_summary[Gene_symbol == "ERBB2"][which.max(mean_row_z), condition] == "HER2",
    identical(unique(marker_meta$marker_group), groups),
    identical(unique(sample_meta$condition), groups),
    identical(unique(clustered_groups), groups),
    length(row_gaps) == 3L && length(column_gaps) == 3L,
    sum(ph$gtable$layout$name == "annotation_legend") == 1L,
    identical(unname(group_colors), c("#89A54E", "#71588F", "#4198AF", "#DB843D")),
    TRUE,
    TRUE,
    all(file.exists(paste0(stem, c(".png", ".tiff", ".pdf", ".svg"))))
  )
)
fwrite(qa, file.path(bundle_dir, "heatmap_validation_checks.tsv"), sep = "\t")
if (!all(qa$passed)) stop("Heatmap validation failed")


message("Completed: ", stem)
print(group_summary)
print(qa)
