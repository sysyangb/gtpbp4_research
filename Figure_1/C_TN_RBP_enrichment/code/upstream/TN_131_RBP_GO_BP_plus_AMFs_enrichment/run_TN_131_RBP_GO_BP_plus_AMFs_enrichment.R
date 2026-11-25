#!/usr/bin/env Rscript


args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
script_path <- if (length(file_arg)) sub("^--file=", "", file_arg[[1]]) else "."
script_path <- gsub("~\\+~", " ", script_path)
out_dir <- normalizePath(dirname(script_path), mustWork = TRUE)
project_dir <- normalizePath(file.path(out_dir, ".."), mustWork = TRUE)
tcga_root <- normalizePath(file.path(project_dir, "..", "TCGA"), mustWork = TRUE)

suppressPackageStartupMessages({
  library(data.table)
  library(readxl)
  library(ggplot2)
  library(scales)
  library(stringr)
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

workbook_path <- file.path(
  project_dir,
  "TCGA_BRCA_four_group_specific_translation_no_ribosome_RBP_annotated.xlsx"
)
analysis_root <- file.path(tcga_root, "DESeq2_TN_RBP_analysis")
tested_path <- file.path(analysis_root, "results", "DESeq2_TN_vs_NON_TN_all_genes.tsv")
annotation_path <- file.path(
  analysis_root,
  "KEGG_enrichment_top500_TN_all_genes",
  "reference_snapshot",
  "Ensembl_BioMart_hsapiens_gene_annotation.tsv"
)
go_dir <- file.path(analysis_root, "GO_BP_enrichment_two_gene_sets", "reference_snapshot")
go_map_path <- file.path(go_dir, "GO_BP_propagated_symbol_to_term.tsv.gz")
go_names_path <- file.path(go_dir, "GO_BP_term_names.tsv")
amf_path <- file.path(
  analysis_root,
  "KEGG_plus_AMFs_enrichment",
  "reference_background",
  "AMFs_215_gene_list_from_workbook.tsv"
)
required <- c(
  workbook_path, tested_path, annotation_path, go_map_path, go_names_path, amf_path
)
if (any(!file.exists(required))) {
  stop("Missing input(s): ", paste(required[!file.exists(required)], collapse = ", "))
}

results_dir <- file.path(out_dir, "results")
source_dir <- file.path(out_dir, "source_data")
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(source_dir, recursive = TRUE, showWarnings = FALSE)

clean_symbol <- function(x) {
  x <- toupper(trimws(as.character(x)))
  x[is.na(x) | x == "" | x == "NA"] <- NA_character_
  x
}

message("Reading the 131 TN RBP genes...")
tn <- as.data.table(read_excel(workbook_path, sheet = "TN"))
required_tn <- c("Gene_ID", "Gene_symbol", "Gene_biotype", "Is_RBP")
if (!all(required_tn %in% names(tn))) stop("TN worksheet is missing required columns")
tn_rbp <- tn[Is_RBP == "Yes" & Gene_biotype == "protein_coding"]
tn_rbp[, gene_symbol_clean := clean_symbol(Gene_symbol)]
query_symbols <- sort(unique(na.omit(tn_rbp$gene_symbol_clean)))
if (length(query_symbols) != 131L) {
  stop("Expected 131 unique TN RBP symbols, recovered ", length(query_symbols))
}
fwrite(tn_rbp, file.path(source_dir, "TN_131_RBP_rows_used.tsv"), sep = "\t", na = "NA")

message("Preparing the full 34,291-gene tested background...")
tested <- fread(tested_path, select = "gene_id", colClasses = "character", showProgress = FALSE)
if (nrow(tested) != 34291L || uniqueN(tested$gene_id) != 34291L) {
  stop("The expected 34,291-gene DESeq2-tested background was not recovered")
}
annotation <- fread(
  annotation_path,
  colClasses = list(
    character = c("ensembl_gene_id", "hgnc_symbol", "entrezgene_id", "gene_biotype")
  ),
  showProgress = FALSE
)
annotation[, gene_symbol := clean_symbol(hgnc_symbol)]
annotation_tested <- unique(annotation[ensembl_gene_id %chin% tested$gene_id])
background_symbols <- sort(unique(na.omit(annotation_tested$gene_symbol)))

message("Preparing standard GO BP memberships...")
go_map <- fread(
  cmd = paste("gzip -dc", shQuote(go_map_path)),
  colClasses = "character",
  showProgress = FALSE
)
go_names <- fread(go_names_path, colClasses = "character", showProgress = FALSE)
go_names <- go_names[namespace == "biological_process", .(term_id = go_id, term_name)]
go_membership <- unique(go_map[, .(
  term_id = go_id,
  gene = clean_symbol(gene_symbol)
)][!is.na(gene)])
go_membership <- merge(go_membership, go_names, by = "term_id", all.x = TRUE)
go_membership[is.na(term_name) | term_name == "", term_name := term_id]
go_membership <- unique(go_membership[
  gene %chin% background_symbols,
  .(
    term_id,
    term_name,
    source_type = "GO Biological Process",
    gene
  )
])

message("Appending the custom AMFs term to the GO BP membership library...")
amf_source <- fread(amf_path, colClasses = "character", showProgress = FALSE)
amf_source[, gene_symbol_clean := clean_symbol(workbook_gene_symbol)]
amf_symbols <- sort(unique(na.omit(amf_source$gene_symbol_clean)))
if (length(amf_symbols) != 215L) {
  stop("Expected 215 unique AMFs symbols, recovered ", length(amf_symbols))
}
amf_in_tested <- intersect(amf_symbols, background_symbols)
amf_membership <- data.table(
  term_id = "CUSTOM_AMFS",
  term_name = "Assembly and maturation factors (AMFs)",
  source_type = "Custom AMFs",
  gene = amf_in_tested
)
membership <- unique(rbindlist(list(go_membership, amf_membership), use.names = TRUE))

# One combined mapped universe is used for every standard GO BP and AMFs term.
# This prevents the custom AMFs term from being evaluated against a different denominator.
background_mapped <- sort(unique(membership$gene))
query_in_tested <- intersect(query_symbols, background_symbols)
query_mapped <- intersect(query_in_tested, background_mapped)
N <- length(background_mapped)
n <- length(query_mapped)
if (!N || !n) stop("No mapped background or query genes")

fwrite(
  amf_source,
  file.path(source_dir, "AMFs_215_original_rows_with_clean_symbols.tsv"),
  sep = "\t",
  na = "NA"
)
fwrite(
  data.table(gene_symbol = background_mapped),
  file.path(source_dir, "combined_GO_BP_plus_AMFs_mapped_background_symbols.tsv"),
  sep = "\t"
)
fwrite(
  data.table(gene_symbol = query_mapped),
  file.path(source_dir, "TN_RBP_mapped_to_combined_library.tsv"),
  sep = "\t"
)

message("Running one-sided hypergeometric ORA on the combined library...")
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
  query_in_tested_background = length(query_in_tested),
  query_mapped_count = n,
  mapped_background_count = N,
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
  "candidate_symbols", "query_input_count", "query_in_tested_background",
  "query_mapped_count", "mapped_background_count", "pathway_url"
))
fwrite(
  result,
  file.path(results_dir, "TN_131_RBP_GO_BP_plus_AMFs_ORA_all_pathways.tsv"),
  sep = "\t",
  na = "NA"
)

amf_result <- result[term_id == "CUSTOM_AMFS"]
if (nrow(amf_result) != 1L) stop("The custom AMFs term was not recovered exactly once")
fwrite(
  amf_result,
  file.path(results_dir, "TN_131_RBP_AMFs_enrichment_result.tsv"),
  sep = "\t",
  na = "NA"
)
amf_hit_symbols <- sort(intersect(query_symbols, amf_in_tested))
amf_hit_detail <- merge(
  data.table(gene_symbol_clean = amf_hit_symbols),
  amf_source,
  by = "gene_symbol_clean",
  all.x = TRUE,
  sort = TRUE
)
fwrite(
  amf_hit_detail,
  file.path(results_dir, "TN_131_RBP_AMFs_hit_genes.tsv"),
  sep = "\t",
  na = "NA"
)

plot_dt <- copy(result[candidate_count >= 1L & is.finite(pvalue)])
setorder(plot_dt, pvalue, p_adjust_BH, -fold_enrichment, -candidate_count, term_id)
plot_dt <- head(plot_dt, 20L)
plot_dt[, top20_rank_raw_P := seq_len(.N)]
plot_dt[, minus_log10_pvalue := -log10(pmax(pvalue, .Machine$double.xmin))]
plot_dt[, label := str_wrap(term_name, width = 48)]
plot_dt[, label := make.unique(label, sep = " ")]
plot_dt[, pathway_factor := factor(label, levels = rev(label))]
fwrite(
  plot_dt[, setdiff(names(plot_dt), "pathway_factor"), with = FALSE],
  file.path(results_dir, "TN_131_RBP_GO_BP_plus_AMFs_Top20_rawP_source.tsv"),
  sep = "\t",
  na = "NA"
)

max_count <- max(plot_dt$candidate_count)
size_breaks <- pretty(c(0, max_count), n = 3)
size_breaks <- size_breaks[size_breaks > 0 & size_breaks <= max_count]
if (!length(size_breaks)) size_breaks <- max_count

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
  geom_point(
    data = plot_dt[term_id == "CUSTOM_AMFS"],
    shape = 21,
    colour = "#156082",
    stroke = 1.25,
    alpha = 1,
    show.legend = FALSE
  ) +
  scale_size_area(
    max_size = 10.2,
    breaks = size_breaks,
    limits = c(0, max_count),
    name = "Gene count"
  ) +
  scale_fill_gradient(
    low = "#F8E8E5",
    high = "#B34D49",
    name = "Enrichment score",
    breaks = breaks_pretty(n = 5)
  ) +
  scale_x_continuous(
    name = expression(-log[10](italic(P)~value)),
    breaks = breaks_pretty(n = 5),
    expand = expansion(mult = c(0.03, 0.10))
  ) +
  scale_y_discrete(name = NULL, expand = expansion(add = c(0.45, 0.45))) +
  guides(
    size = guide_legend(
      order = 1,
      title.position = "top",
      override.aes = list(fill = "#B34D49", colour = "#843B38", alpha = 1)
    ),
    fill = guide_colourbar(
      order = 2,
      title.position = "top",
      barheight = grid::unit(34, "mm"),
      barwidth = grid::unit(6, "mm"),
      ticks = TRUE,
      frame.colour = NA
    )
  ) +
  theme_classic(base_size = 10.8, base_family = "Arial") +
  theme(
    axis.line = element_line(linewidth = 0.7, colour = "#111111"),
    axis.ticks = element_line(linewidth = 0.65, colour = "#111111"),
    axis.ticks.length = grid::unit(2.2, "mm"),
    axis.text.x = element_text(size = 11, colour = "#222222"),
    axis.text.y = element_text(size = 10.5, colour = "#222222", lineheight = 0.86),
    axis.title.x = element_text(size = 11.8, margin = margin(t = 9)),
    legend.position = "right",
    legend.box = "vertical",
    legend.box.just = "left",
    legend.title = element_text(size = 11.2, face = "bold"),
    legend.text = element_text(size = 10.6),
    legend.key.height = grid::unit(7, "mm"),
    legend.spacing.y = grid::unit(4, "mm"),
    legend.margin = margin(l = 4),
    plot.margin = margin(9, 11, 13, 9, unit = "mm"),
    panel.grid = element_blank()
  )

stem <- file.path(out_dir, "TN_131_RBP_GO_BP_plus_AMFs_Top20_rawP_dotplot")
width_mm <- 245
height_mm <- 260
ggsave(
  paste0(stem, ".png"), p, device = ragg::agg_png,
  width = width_mm, height = height_mm, units = "mm", dpi = 400, bg = "white"
)
ggsave(
  paste0(stem, ".tiff"), p, device = ragg::agg_tiff,
  width = width_mm, height = height_mm, units = "mm", dpi = 600,
  compression = "lzw", bg = "white"
)
ggsave(
  paste0(stem, ".pdf"), p, device = cairo_pdf,
  width = width_mm, height = height_mm, units = "mm", family = "Arial", bg = "white"
)
ggsave(
  paste0(stem, ".svg"), p, device = svglite::svglite,
  width = width_mm, height = height_mm, units = "mm", bg = "white"
)

qa <- data.table(
  check = c(
    "query_contains_131_unique_TN_RBP_symbols",
    "background_contains_34291_DESeq2_tested_genes",
    "AMFs_contains_215_unique_source_symbols",
    "combined_term_size_is_10_to_500",
    "AMFs_is_in_same_combined_result_table",
    "BH_FDR_was_computed_jointly_for_all_terms",
    "top20_ranked_by_raw_P",
    "plot_has_20_terms",
    "AMFs_natural_rank_not_forced",
    "all_exports_exist"
  ),
  passed = c(
    length(query_symbols) == 131L,
    nrow(tested) == 34291L,
    length(amf_symbols) == 215L,
    all(result$background_count >= 10L & result$background_count <= 500L),
    sum(result$term_id == "CUSTOM_AMFS") == 1L,
    isTRUE(all.equal(result$p_adjust_BH, p.adjust(result$pvalue, method = "BH"))),
    identical(plot_dt$term_id, head(result[candidate_count >= 1L]$term_id, 20L)),
    nrow(plot_dt) == 20L,
    (amf_result$rank_raw_P <= 20L) == ("CUSTOM_AMFS" %chin% plot_dt$term_id),
    all(file.exists(paste0(stem, c(".png", ".tiff", ".pdf", ".svg"))))
  )
)
fwrite(qa, file.path(out_dir, "analysis_and_figure_validation.tsv"), sep = "\t")
if (!all(qa$passed)) stop("Validation failed")

summary_dt <- data.table(
  metric = c(
    "TN RBP input symbols",
    "TN RBP symbols in tested background",
    "TN RBP symbols mapped to combined library",
    "Combined mapped background symbols",
    "Standard GO BP mapped background symbols",
    "AMFs source symbols",
    "AMFs symbols in tested background",
    "Eligible standard GO BP terms",
    "Eligible combined terms including AMFs",
    "AMFs TN-RBP hits",
    "AMFs raw-P rank",
    "AMFs raw P",
    "AMFs BH FDR",
    "AMFs fold enrichment",
    "Top 20 terms with BH FDR < 0.05"
  ),
  value = c(
    length(query_symbols), length(query_in_tested), length(query_mapped), N,
    uniqueN(go_membership$gene), length(amf_symbols), length(amf_in_tested),
    sum(result$source_type == "GO Biological Process"), nrow(result),
    amf_result$candidate_count, amf_result$rank_raw_P,
    format(amf_result$pvalue, scientific = TRUE, digits = 8),
    format(amf_result$p_adjust_BH, scientific = TRUE, digits = 8),
    format(amf_result$fold_enrichment, scientific = FALSE, digits = 8),
    sum(plot_dt$p_adjust_BH < 0.05)
  )
)
fwrite(summary_dt, file.path(out_dir, "analysis_summary.tsv"), sep = "\t")


message("Completed: ", stem)
print(summary_dt)
print(plot_dt[, .(
  top20_rank_raw_P, term_id, term_name, source_type, candidate_count,
  fold_enrichment, pvalue, p_adjust_BH
)])
