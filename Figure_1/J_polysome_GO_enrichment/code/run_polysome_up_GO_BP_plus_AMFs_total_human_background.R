#!/usr/bin/env Rscript


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
options(mc.cores = 1L, width = 180, scipen = 999, warn = 1)
data.table::setDTthreads(1L)
set.seed(20260902)

args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
script_path <- if (length(file_arg)) sub("^--file=", "", file_arg[[1]]) else "."
script_path <- normalizePath(script_path, mustWork = TRUE)
scripts_dir <- dirname(script_path)
analysis_dir <- normalizePath(file.path(scripts_dir, ".."), mustWork = TRUE)
project_dir <- normalizePath(file.path(analysis_dir, ".."), mustWork = TRUE)
dongyiming_dir <- normalizePath(file.path(project_dir, ".."), mustWork = TRUE)

results_dir <- file.path(analysis_dir, "results")
source_dir <- file.path(analysis_dir, "source_data")
figures_dir <- file.path(analysis_dir, "figures")
qc_dir <- file.path(analysis_dir, "QC")
for (path in c(results_dir, source_dir, figures_dir, qc_dir)) {
  dir.create(path, recursive = TRUE, showWarnings = FALSE)
}

protein_statistics_path <- file.path(
  project_dir,
  "reanalysis_exclude_polysome3_20260901",
  "results",
  "protein_statistics_all.csv"
)
amf_path <- file.path(
  dongyiming_dir,
  "01TCGA 筛选",
  "TN_131_RBP_GO_BP_plus_AMFs_enrichment",
  "source_data",
  "AMFs_215_original_rows_with_clean_symbols.tsv"
)
go_reference_dir <- file.path(
  dongyiming_dir,
  "TCGA",
  "DESeq2_TN_RBP_analysis",
  "GO_BP_enrichment_two_gene_sets",
  "reference_snapshot"
)
go_map_path <- file.path(go_reference_dir, "GO_BP_propagated_symbol_to_term.tsv.gz")
go_names_path <- file.path(go_reference_dir, "GO_BP_term_names.tsv")
required_inputs <- c(protein_statistics_path, amf_path, go_map_path, go_names_path)
if (any(!file.exists(required_inputs))) {
  stop("Missing input(s): ", paste(required_inputs[!file.exists(required_inputs)], collapse = ", "))
}

clean_symbol <- function(values) {
  cleaned <- toupper(trimws(as.character(values)))
  cleaned[is.na(cleaned) | cleaned == "" | cleaned == "NA"] <- NA_character_
  cleaned
}

message("Preparing the 502 polysome-up genes; no experiment-specific background is used...")
protein_statistics <- fread(protein_statistics_path, showProgress = FALSE)
required_columns <- c(
  "Protein accession", "Gene name", "Unique peptides",
  "p_value_excluding_polysome_3", "q_value_excluding_polysome_3",
  "log2FC_excluding_polysome_3"
)
if (!all(required_columns %in% names(protein_statistics))) {
  stop("Protein statistics file is missing required columns")
}
protein_statistics[, gene_symbol := clean_symbol(`Gene name`)]
query_rows <- protein_statistics[
  !is.na(p_value_excluding_polysome_3) &
    p_value_excluding_polysome_3 < 0.05 &
    !is.na(log2FC_excluding_polysome_3) &
    log2FC_excluding_polysome_3 > 1 &
    `Unique peptides` >= 5 &
    !is.na(gene_symbol)
]
query_symbols <- sort(unique(query_rows$gene_symbol))
if (length(query_symbols) != 502L) {
  stop("Expected 502 polysome-up genes, recovered ", length(query_symbols))
}

fwrite(
  unique(query_rows[, .(
    gene_symbol,
    protein_accession = `Protein accession`,
    unique_peptides = `Unique peptides`,
    log2FC = log2FC_excluding_polysome_3,
    p_value = p_value_excluding_polysome_3,
    BH_q_value = q_value_excluding_polysome_3
  )]),
  file.path(source_dir, "Polysome_up_502_genes.tsv"),
  sep = "\t",
  na = "NA"
)

message("Building the complete human GO Biological Process membership library...")
go_map <- fread(
  cmd = paste("gzip -dc", shQuote(go_map_path)),
  colClasses = "character",
  showProgress = FALSE
)
go_names <- fread(go_names_path, colClasses = "character", showProgress = FALSE)
go_names <- go_names[
  namespace == "biological_process",
  .(term_id = go_id, term_name)
]
go_membership <- unique(go_map[, .(
  term_id = go_id,
  gene = clean_symbol(gene_symbol)
)][!is.na(gene)])
go_membership <- merge(go_membership, go_names, by = "term_id", all = FALSE)
go_membership <- unique(go_membership[, .(
  term_id,
  term_name,
  source_type = "GO Biological Process",
  gene
)])

message("Appending all 215 AMFs genes as a custom BP-like term...")
amf_source <- fread(amf_path, colClasses = "character", showProgress = FALSE)
if (!"gene_symbol_clean" %in% names(amf_source)) {
  stop("AMFs source is missing gene_symbol_clean")
}
amf_source[, gene_symbol_clean := clean_symbol(gene_symbol_clean)]
amf_symbols <- sort(unique(na.omit(amf_source$gene_symbol_clean)))
if (length(amf_symbols) != 215L) {
  stop("Expected 215 unique AMFs symbols, recovered ", length(amf_symbols))
}
amf_membership <- data.table(
  term_id = "CUSTOM_AMFS",
  term_name = "Assembly and maturation factors (AMFs)",
  source_type = "Custom AMFs",
  gene = amf_symbols
)
membership <- unique(rbindlist(list(go_membership, amf_membership), use.names = TRUE))

# The effective total-human background is the union of all human genes that can
# be selected by at least one tested GO BP or AMFs term. No proteomics
# detectability, P-value, fold-change, or peptide-count filter is applied here.
go_background <- sort(unique(go_membership$gene))
background_symbols <- sort(unique(membership$gene))
query_mapped <- intersect(query_symbols, background_symbols)
query_unmapped <- setdiff(query_symbols, background_symbols)
N <- length(background_symbols)
n <- length(query_mapped)
if (N == 0L || n == 0L) stop("No mapped background or query genes")
if (N < 10000L) stop("The total-human annotated background is unexpectedly small: ", N)

background_annotation <- data.table(gene_symbol = background_symbols)
background_annotation[, in_GO_BP := gene_symbol %chin% go_background]
background_annotation[, in_AMFs := gene_symbol %chin% amf_symbols]
background_annotation[, annotation_source := fifelse(
  in_GO_BP & in_AMFs, "GO BP + AMFs",
  fifelse(in_GO_BP, "GO BP only", "AMFs only")
)]
fwrite(
  background_annotation,
  file.path(source_dir, "Total_human_GO_BP_plus_AMFs_annotated_background.tsv"),
  sep = "\t"
)
fwrite(
  data.table(gene_symbol = query_mapped),
  file.path(source_dir, "Polysome_up_mapped_to_total_human_library.tsv"),
  sep = "\t"
)
fwrite(
  data.table(gene_symbol = query_unmapped),
  file.path(source_dir, "Polysome_up_unmapped_to_total_human_library.tsv"),
  sep = "\t"
)
fwrite(amf_source, file.path(source_dir, "AMFs_215_original_rows.tsv"), sep = "\t", na = "NA")

message("Running one-sided hypergeometric ORA against the total-human annotated background...")
term_background <- membership[, .(
  background_count = uniqueN(gene)
), by = .(term_id, term_name, source_type)]
term_background <- term_background[background_count >= 10L & background_count <= 500L]

hits <- membership[
  term_id %chin% term_background$term_id & gene %chin% query_mapped,
  .(
    candidate_count = uniqueN(gene),
    candidate_symbols = paste(sort(unique(gene)), collapse = "/")
  ),
  by = .(term_id, term_name, source_type)
]
result <- merge(
  term_background,
  hits,
  by = c("term_id", "term_name", "source_type"),
  all.x = TRUE
)
result[is.na(candidate_count), candidate_count := 0L]
result[is.na(candidate_symbols), candidate_symbols := ""]
result[, pvalue := phyper(
  candidate_count - 1L,
  background_count,
  N - background_count,
  n,
  lower.tail = FALSE
)]
result[, p_adjust_BH := p.adjust(pvalue, method = "BH")]
result[, expected_count := n * background_count / N]
result[, fold_enrichment := fifelse(
  expected_count > 0,
  candidate_count / expected_count,
  NA_real_
)]
result[, gene_ratio := candidate_count / n]
result[, background_ratio := background_count / N]
result[, GeneRatio := paste0(candidate_count, "/", n)]
result[, BgRatio := paste0(background_count, "/", N)]
result[, `:=`(
  query_input_count = length(query_symbols),
  query_mapped_count = n,
  total_human_GO_BP_genes = length(go_background),
  total_human_combined_background = N,
  background_definition = "All human genes represented in the GO BP snapshot or AMFs list",
  pathway_url = fifelse(
    source_type == "GO Biological Process",
    paste0("https://amigo.geneontology.org/amigo/term/", term_id),
    NA_character_
  )
)]
setorder(result, pvalue, p_adjust_BH, -fold_enrichment, -candidate_count, term_id)
result[, rank_raw_P := seq_len(.N)]
setcolorder(result, c(
  "rank_raw_P", "term_id", "term_name", "source_type", "candidate_count",
  "background_count", "GeneRatio", "BgRatio", "gene_ratio", "background_ratio",
  "expected_count", "fold_enrichment", "pvalue", "p_adjust_BH",
  "candidate_symbols", "query_input_count", "query_mapped_count",
  "total_human_GO_BP_genes", "total_human_combined_background",
  "background_definition", "pathway_url"
))

all_pathways_path <- file.path(results_dir, "Polysome_up_GO_BP_plus_AMFs_total_human_ORA_all_pathways.tsv")
fwrite(result, all_pathways_path, sep = "\t", na = "NA")

amf_result <- result[term_id == "CUSTOM_AMFS"]
if (nrow(amf_result) != 1L) stop("The custom AMFs term was not recovered exactly once")
fwrite(
  amf_result,
  file.path(results_dir, "Polysome_up_AMFs_total_human_background_enrichment.tsv"),
  sep = "\t",
  na = "NA"
)
amf_hit_symbols <- sort(intersect(query_symbols, amf_symbols))
amf_hit_detail <- merge(
  data.table(gene_symbol_clean = amf_hit_symbols),
  amf_source,
  by = "gene_symbol_clean",
  all.x = TRUE,
  sort = TRUE
)
fwrite(
  amf_hit_detail,
  file.path(results_dir, "Polysome_up_AMFs_hit_genes.tsv"),
  sep = "\t",
  na = "NA"
)

plot_dt <- copy(result[candidate_count >= 1L & is.finite(pvalue)])
setorder(plot_dt, pvalue, p_adjust_BH, -fold_enrichment, -candidate_count, term_id)
plot_dt <- head(plot_dt, 10L)
plot_dt[, top10_rank_raw_P := seq_len(.N)]
plot_dt[, minus_log10_pvalue := -log10(pmax(pvalue, .Machine$double.xmin))]
plot_dt[, label := str_wrap(term_name, width = 43)]
plot_dt[, label := make.unique(label, sep = " ")]
plot_dt[, pathway_factor := factor(label, levels = rev(label))]
fwrite(
  plot_dt[, setdiff(names(plot_dt), "pathway_factor"), with = FALSE],
  file.path(results_dir, "Polysome_up_GO_BP_plus_AMFs_total_human_Top10_source_data.tsv"),
  sep = "\t",
  na = "NA"
)

min_count <- min(plot_dt$candidate_count)
max_count <- max(plot_dt$candidate_count)
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
    max_size = 7.2,
    breaks = size_breaks,
    limits = c(0, max_count),
    name = "Gene count"
  ) +
  scale_fill_gradient(
    low = "#F8E8E5",
    high = "#B34D49",
    name = "Fold enrichment",
    breaks = breaks_pretty(n = 4)
  ) +
  scale_x_continuous(
    name = expression(-log[10](italic(P)~value)),
    breaks = breaks_pretty(n = 5),
    expand = expansion(mult = c(0.08, 0.08))
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
  theme_classic(base_size = 10.5, base_family = "Helvetica") +
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
  output <- numeric()
  if (inherits(grob, "text") && !is.null(grob$gp$fontsize)) {
    output <- c(output, as.numeric(grob$gp$fontsize))
  }
  if (!is.null(grob$children)) {
    output <- c(output, unlist(lapply(grob$children, collect_font_sizes), use.names = FALSE))
  }
  if (!is.null(grob$grobs)) {
    output <- c(output, unlist(lapply(grob$grobs, collect_font_sizes), use.names = FALSE))
  }
  output
}
grDevices::pdf(NULL)
plot_grob <- ggplotGrob(p)
grDevices::dev.off()
font_sizes <- collect_font_sizes(plot_grob)
if (!length(font_sizes) || min(font_sizes) < 10) {
  stop("Detected plotted text below 10 pt; minimum was ", min(font_sizes))
}

figure_stem <- file.path(figures_dir, "Polysome_up_GO_BP_plus_AMFs_total_human_Top10_dotplot")
width_mm <- 185
height_mm <- 118
ggsave(
  paste0(figure_stem, ".png"), p, device = ragg::agg_png,
  width = width_mm, height = height_mm, units = "mm", dpi = 600, bg = "white"
)
ggsave(
  paste0(figure_stem, ".tiff"), p, device = ragg::agg_tiff,
  width = width_mm, height = height_mm, units = "mm", dpi = 600,
  compression = "lzw", bg = "white"
)
ggsave(
  paste0(figure_stem, ".pdf"), p, device = cairo_pdf,
  width = width_mm, height = height_mm, units = "mm", family = "Helvetica", bg = "white"
)
ggsave(
  paste0(figure_stem, ".svg"), p, device = svglite::svglite,
  width = width_mm, height = height_mm, units = "mm", bg = "white"
)

summary_dt <- data.table(
  metric = c(
    "Polysome-up input genes",
    "Polysome-up genes mapped to combined library",
    "Polysome-up genes not mapped",
    "Complete human GO BP annotated genes",
    "AMFs source genes",
    "AMFs genes already represented in GO BP",
    "Combined total-human GO BP plus AMFs background genes",
    "Eligible standard GO BP terms",
    "Eligible combined terms including AMFs",
    "AMFs polysome-up hits",
    "AMFs raw-P rank",
    "AMFs raw P",
    "AMFs BH FDR",
    "AMFs fold enrichment",
    "Top 10 terms with BH FDR < 0.05",
    "Minimum plotted text size pt",
    "Figure width mm",
    "Figure height mm"
  ),
  value = c(
    length(query_symbols),
    n,
    length(query_unmapped),
    length(go_background),
    length(amf_symbols),
    length(intersect(amf_symbols, go_background)),
    N,
    sum(result$source_type == "GO Biological Process"),
    nrow(result),
    amf_result$candidate_count,
    amf_result$rank_raw_P,
    format(amf_result$pvalue, scientific = TRUE, digits = 8),
    format(amf_result$p_adjust_BH, scientific = TRUE, digits = 8),
    signif(amf_result$fold_enrichment, 8),
    sum(plot_dt$p_adjust_BH < 0.05),
    min(font_sizes),
    width_mm,
    height_mm
  )
)
fwrite(summary_dt, file.path(results_dir, "analysis_summary.tsv"), sep = "\t")


package_names <- c("data.table", "ggplot2", "scales", "stringr", "ragg", "svglite")
package_versions <- data.table(
  package = package_names,
  version = vapply(package_names, function(pkg) as.character(packageVersion(pkg)), character(1))
)
fwrite(package_versions, file.path(qc_dir, "R_package_versions.tsv"), sep = "\t")
capture.output(sessionInfo(), file = file.path(qc_dir, "sessionInfo.txt"))

input_md5 <- data.table(
  input = c(
    "protein_statistics_all.csv",
    "AMFs_215_original_rows_with_clean_symbols.tsv",
    "GO_BP_propagated_symbol_to_term.tsv.gz",
    "GO_BP_term_names.tsv"
  ),
  md5 = unname(tools::md5sum(required_inputs))
)
fwrite(input_md5, file.path(qc_dir, "input_md5.tsv"), sep = "\t")

export_paths <- paste0(figure_stem, c(".png", ".tiff", ".pdf", ".svg"))
required_result_paths <- c(
  all_pathways_path,
  file.path(results_dir, "Polysome_up_GO_BP_plus_AMFs_total_human_Top10_source_data.tsv"),
  file.path(results_dir, "Polysome_up_AMFs_total_human_background_enrichment.tsv"),
  file.path(results_dir, "Polysome_up_AMFs_hit_genes.tsv"),
  file.path(results_dir, "analysis_summary.tsv")
)
qa <- data.table(
  check = c(
    "protein_statistics_contains_7347_rows",
    "query_contains_502_unique_polysome_up_genes",
    "query_filter_is_log2FC_gt_1_raw_P_lt_0.05_peptides_ge_5",
    "experiment_specific_2674_gene_background_is_not_used",
    "full_GO_BP_background_contains_16384_genes",
    "AMFs_contains_215_unique_source_genes",
    "combined_background_is_exact_GO_BP_AMFs_union",
    "combined_total_human_background_contains_16396_genes",
    "query_maps_484_genes_to_combined_library",
    "combined_term_size_is_10_to_500",
    "AMFs_is_present_once_in_combined_results",
    "BH_FDR_is_jointly_recomputed_for_all_terms",
    "top10_contains_exactly_10_terms",
    "top10_is_natural_raw_P_order",
    "AMFs_is_not_forced_into_top10",
    "all_top10_terms_have_query_hits",
    "all_drawn_text_is_at_least_10_pt",
    "figure_has_no_main_title",
    "all_four_figure_exports_exist_and_are_nonempty",
    "all_required_result_tables_exist"
  ),
  passed = c(
    nrow(protein_statistics) == 7347L,
    length(query_symbols) == 502L,
    all(query_rows$log2FC_excluding_polysome_3 > 1 &
          query_rows$p_value_excluding_polysome_3 < 0.05 &
          query_rows$`Unique peptides` >= 5),
    N != 2674L && N > 10000L,
    length(go_background) == 16384L,
    length(amf_symbols) == 215L,
    identical(background_symbols, sort(union(go_background, amf_symbols))),
    N == 16396L,
    n == 484L,
    all(result$background_count >= 10L & result$background_count <= 500L),
    sum(result$term_id == "CUSTOM_AMFS") == 1L,
    isTRUE(all.equal(result$p_adjust_BH, p.adjust(result$pvalue, method = "BH"))),
    nrow(plot_dt) == 10L,
    identical(plot_dt$term_id, head(result[candidate_count >= 1L]$term_id, 10L)),
    (amf_result$rank_raw_P <= 10L) == ("CUSTOM_AMFS" %chin% plot_dt$term_id),
    all(plot_dt$candidate_count >= 1L),
    min(font_sizes) >= 10,
    is.null(p$labels$title),
    all(file.exists(export_paths) & file.info(export_paths)$size > 0),
    all(file.exists(required_result_paths))
  )
)
fwrite(qa, file.path(qc_dir, "analysis_and_figure_validation.tsv"), sep = "\t")
if (!all(qa$passed)) {
  stop("Validation failed: ", paste(qa[passed == FALSE, check], collapse = ", "))
}


message("Completed: ", figure_stem)
print(summary_dt)
print(plot_dt[, .(
  top10_rank_raw_P,
  term_id,
  term_name,
  source_type,
  candidate_count,
  background_count,
  fold_enrichment,
  pvalue,
  p_adjust_BH
)])
print(qa)
