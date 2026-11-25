#!/usr/bin/env Rscript

# TCGA-BRCA sample-level scoring for four prespecified GO Biological Process
# gene sets. The cohort, preprocessing, scoring methods, patient-aware models,
# and Figure 1D styling follow the existing AMF workflow.

options(
  stringsAsFactors = FALSE,
  width = 140,
  warn = 1,
  contrasts = c("contr.treatment", "contr.poly")
)
Sys.setenv(
  OMP_NUM_THREADS = "1",
  OPENBLAS_NUM_THREADS = "1",
  MKL_NUM_THREADS = "1",
  VECLIB_MAXIMUM_THREADS = "1",
  NUMEXPR_NUM_THREADS = "1"
)

root <- normalizePath(getwd())
local_lib <- file.path(root, ".Rlib")
if (dir.exists(local_lib)) .libPaths(c(local_lib, .libPaths()))

defaults <- list(
  tpm = file.path("..", "TCGA", "TCGA-BRCA_TPM.txt"),
  selected_samples = file.path("AMF_scores_four_groups", "QC", "selected_samples.tsv"),
  go_mapping = file.path(
    "..", "TCGA", "DESeq2_TN_RBP_analysis", "GO_BP_enrichment_two_gene_sets",
    "reference_snapshot", "GO_BP_propagated_symbol_to_term.tsv.gz"
  ),
  go_direct_mapping = file.path(
    "..", "TCGA", "DESeq2_TN_RBP_analysis", "GO_BP_enrichment_two_gene_sets",
    "reference_snapshot", "GO_BP_direct_symbol_to_term.tsv.gz"
  ),
  go_terms = file.path(
    "..", "TCGA", "DESeq2_TN_RBP_analysis", "GO_BP_enrichment_two_gene_sets",
    "reference_snapshot", "GO_BP_term_names.tsv"
  ),
  annotation = file.path(
    "..", "TCGA", "DESeq2_TN_RBP_analysis", "KEGG_enrichment_top500_TN_all_genes",
    "reference_snapshot", "Ensembl_BioMart_hsapiens_gene_annotation.tsv"
  ),
  outdir = "GO_pathway_scores_four_groups",
  seed = "20260831",
  ymin = "-3.5",
  ymax = "5.0"
)

usage <- function() {
  cat(paste0(
    "Usage: Rscript run_GO_pathway_scores.R [--key value ...]\n\n",
    "Options:\n",
    "  --tpm PATH               TCGA-BRCA TPM matrix\n",
    "  --selected_samples PATH  Selected 1,116-sample cohort from the AMF workflow\n",
    "  --go_mapping PATH        Propagated GO BP symbol-to-term mapping\n",
    "  --go_direct_mapping PATH Direct GO BP symbol-to-term mapping\n",
    "  --go_terms PATH          GO term-name table\n",
    "  --annotation PATH        Ensembl-symbol-biotype annotation\n",
    "  --outdir PATH            Output directory [GO_pathway_scores_four_groups]\n",
    "  --seed INTEGER           Random seed [20260831]\n",
    "  --ymin NUMBER            Figure lower display limit [-3.5]\n",
    "  --ymax NUMBER            Figure upper display limit [5.0]\n",
    "  --help                   Show this message\n"
  ))
}

parse_args <- function(args, defaults) {
  if (any(args %in% c("--help", "-h"))) {
    usage()
    quit(save = "no", status = 0L)
  }
  if (length(args) %% 2L != 0L) stop("Arguments must be supplied as --key value pairs.")
  out <- defaults
  if (length(args)) {
    keys <- sub("^--", "", args[seq.int(1L, length(args), by = 2L)])
    values <- args[seq.int(2L, length(args), by = 2L)]
    unknown <- setdiff(keys, names(defaults))
    if (length(unknown)) stop("Unknown argument(s): ", paste(unknown, collapse = ", "))
    for (i in seq_along(keys)) out[[keys[[i]]]] <- values[[i]]
  }
  out
}

cfg <- parse_args(commandArgs(trailingOnly = TRUE), defaults)
cfg$seed <- as.integer(cfg$seed)
cfg$ymin <- as.numeric(cfg$ymin)
cfg$ymax <- as.numeric(cfg$ymax)
if (!is.finite(cfg$seed)) stop("--seed must be a finite integer.")
if (!is.finite(cfg$ymin) || !is.finite(cfg$ymax) || cfg$ymin >= cfg$ymax) {
  stop("--ymin and --ymax must be finite numbers with ymin < ymax.")
}

path_keys <- c("tpm", "selected_samples", "go_mapping", "go_direct_mapping", "go_terms", "annotation")
for (key in path_keys) cfg[[key]] <- normalizePath(cfg[[key]], mustWork = FALSE)
cfg$outdir <- normalizePath(cfg$outdir, mustWork = FALSE)
missing_inputs <- path_keys[!vapply(path_keys, function(x) file.exists(cfg[[x]]), logical(1))]
if (length(missing_inputs)) {
  stop(
    "Missing input(s): ",
    paste(sprintf("%s=%s", missing_inputs, vapply(missing_inputs, function(x) cfg[[x]], character(1))), collapse = ", ")
  )
}

required_packages <- c(
  "data.table", "matrixStats", "ggplot2", "svglite", "ragg",
  "GSVA", "singscore", "BiocParallel", "GSEABase",
  "lme4", "lmerTest", "emmeans", "sandwich", "lmtest"
)
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages)) {
  stop("Missing R packages: ", paste(missing_packages, collapse = ", "), ". Run setup_R_packages.R first.")
}

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})
data.table::setDTthreads(1L)
set.seed(cfg$seed)

results_dir <- file.path(cfg$outdir, "results")
figures_dir <- file.path(cfg$outdir, "figures")
qc_dir <- file.path(cfg$outdir, "QC")
logs_dir <- file.path(cfg$outdir, "logs")
for (d in c(cfg$outdir, results_dir, figures_dir, qc_dir, logs_dir)) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

log_file <- file.path(logs_dir, "analysis.log")
log_con <- file(log_file, open = "wt")
sink(log_con, type = "output", split = TRUE)
sink(log_con, type = "message", append = TRUE)
on.exit({
  try(sink(type = "message"), silent = TRUE)
  try(sink(type = "output"), silent = TRUE)
  try(close(log_con), silent = TRUE)
}, add = TRUE)

read_table <- function(path) {
  if (grepl("\\.gz$", path, ignore.case = TRUE)) {
    data.table::fread(cmd = paste("gzip -dc", shQuote(path)))
  } else {
    data.table::fread(path)
  }
}

strip_version <- function(x) sub("\\..*$", "", gsub('^"|"$', "", x))

safe_z <- function(x) {
  s <- stats::sd(x)
  if (!is.finite(s) || s == 0) stop("Cannot standardize a score with zero or non-finite variance.")
  (x - mean(x)) / s
}

sha256_file <- function(path) {
  ans <- system2("shasum", c("-a", "256", shQuote(path)), stdout = TRUE, stderr = TRUE)
  if (!length(ans)) return(NA_character_)
  strsplit(ans[[1L]], "[[:space:]]+")[[1L]][[1L]]
}

groups <- c("TN", "HR", "HER2", "NORMAL")
group_labels <- c(TN = "TN", HR = "HR", HER2 = "HER2", NORMAL = "Normal")
group_colors <- c(TN = "#AA4643", HR = "#4572A7", HER2 = "#89A54E", NORMAL = "#71588F")
methods <- c("singscore", "ssGSEA", "mean_z")
method_labels <- c(singscore = "singscore", ssGSEA = "ssGSEA", mean_z = "Mean z-score")

pathways <- data.table(
  signature = c("RIBI", "POL-I", "Translation", "Regulator of Translation"),
  go_id = c("GO:0042254", "GO:0006360", "GO:0006412", "GO:0006417"),
  membership_source = c("propagated", "propagated", "direct", "propagated"),
  expected_term_name = c(
    "ribosome biogenesis",
    "transcription by RNA polymerase I",
    "translation",
    "regulation of translation"
  ),
  figure_prefix = c(
    "Figure1-D-RIBI-score",
    "Figure1-D-POL-I-score",
    "Figure1-D-Translation-score",
    "Figure1-D-Regulator-of-Translation-score"
  ),
  y_label = c(
    "Standardized RIBI score",
    "Standardized POL-I score",
    "Standardized translation score",
    "Standardized regulation-of-\ntranslation score"
  )
)

message("Analysis started: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"))
message("R version: ", as.character(getRversion()))

input_checksums <- data.table(
  input = path_keys,
  path = vapply(path_keys, function(x) cfg[[x]], character(1)),
  sha256 = vapply(path_keys, function(x) sha256_file(cfg[[x]]), character(1))
)
fwrite(input_checksums, file.path(qc_dir, "input_sha256.tsv"), sep = "\t")

run_parameters <- data.table(
  parameter = c(names(cfg), "background_rule", "gene_mapping_rule", "statistics", "figure_comparisons"),
  value = c(
    vapply(cfg, as.character, character(1)),
    "TPM > 0 in >=20% selected samples and nonzero variance",
    paste0(
      "GO membership: Translation (GO:0006412) uses direct GOA annotations; RIBI, POL-I, and Regulator of Translation use ",
      "is_a/part_of-propagated annotations. Unique protein-coding HGNC symbols map to unversioned Ensembl IDs; ambiguous symbols excluded"
    ),
    "score ~ condition + (1|patient_id), Tukey pairwise; cluster-robust fallback for singular/nonconvergent fits",
    "TN-HR, TN-HER2, TN-NORMAL"
  )
)
fwrite(run_parameters, file.path(qc_dir, "run_parameters.tsv"), sep = "\t")

# -----------------------------------------------------------------------------
# Fixed cohort
# -----------------------------------------------------------------------------

selected_meta <- fread(cfg$selected_samples, na.strings = c("", "NA"))
required_sample_columns <- c("sample_id", "patient_id", "condition")
missing_sample_columns <- setdiff(required_sample_columns, names(selected_meta))
if (length(missing_sample_columns)) {
  stop("Selected-sample table is missing: ", paste(missing_sample_columns, collapse = ", "))
}
if (anyDuplicated(selected_meta$sample_id)) stop("Selected sample IDs are not unique.")
if (!all(selected_meta$condition %in% groups)) stop("Unexpected condition labels in selected-sample table.")
selected_meta[, condition := factor(condition, levels = groups)]
selected_ids <- selected_meta$sample_id

observed_counts <- selected_meta[, .N, by = condition]
setorder(observed_counts, condition)
expected_counts <- c(TN = 177L, HR = 664L, HER2 = 163L, NORMAL = 112L)
observed_vector <- setNames(observed_counts$N, as.character(observed_counts$condition))[groups]
if (!identical(as.integer(observed_vector), as.integer(expected_counts))) {
  stop("Selected cohort differs from the validated AMF cohort: ", paste(groups, observed_vector, sep = "=", collapse = ", "))
}
fwrite(selected_meta, file.path(qc_dir, "selected_samples.tsv"), sep = "\t")
fwrite(observed_counts, file.path(qc_dir, "selected_sample_counts.tsv"), sep = "\t")
message("Fixed cohort: ", paste(groups, observed_vector, sep = "=", collapse = ", "))

# -----------------------------------------------------------------------------
# TPM matrix and fixed background
# -----------------------------------------------------------------------------

tpm_header <- gsub(
  '^"|"$', "",
  strsplit(readLines(cfg$tpm, n = 1L, warn = FALSE), "\t", fixed = TRUE)[[1L]]
)
if (!all(selected_ids %in% tpm_header)) stop("Selected sample IDs are missing from the TPM matrix.")

message("Reading TPM matrix for ", length(selected_ids), " selected samples...")
tpm_dt <- suppressWarnings(fread(
  cfg$tpm,
  select = c("V1", selected_ids),
  check.names = FALSE,
  showProgress = TRUE,
  nThread = 1L
))
setnames(tpm_dt, "V1", "gene_id_version")
gene_id <- strip_version(tpm_dt$gene_id_version)
tpm_versioned <- as.matrix(tpm_dt[, ..selected_ids])
storage.mode(tpm_versioned) <- "double"
rm(tpm_dt)
invisible(gc())
if (anyNA(tpm_versioned) || any(!is.finite(tpm_versioned)) || any(tpm_versioned < 0)) {
  stop("TPM matrix contains missing, non-finite, or negative values.")
}

message("Collapsing duplicate unversioned Ensembl IDs by summed TPM...")
tpm_matrix <- rowsum(tpm_versioned, group = gene_id, reorder = FALSE)
storage.mode(tpm_matrix) <- "double"
rm(tpm_versioned)
invisible(gc())
if (!identical(colnames(tpm_matrix), selected_ids)) stop("TPM sample order changed unexpectedly.")

gene_detection_fraction <- 1 - matrixStats::rowCounts(tpm_matrix, value = 0) / ncol(tpm_matrix)
gene_variance <- matrixStats::rowVars(tpm_matrix)
names(gene_detection_fraction) <- rownames(tpm_matrix)
names(gene_variance) <- rownames(tpm_matrix)
background_keep <- gene_detection_fraction >= 0.20 & is.finite(gene_variance) & gene_variance > 0
background_genes <- rownames(tpm_matrix)[background_keep]
message("Background genes retained: ", length(background_genes), "/", nrow(tpm_matrix))

# -----------------------------------------------------------------------------
# GO membership and gene-mapping QC
# -----------------------------------------------------------------------------

go_terms <- fread(cfg$go_terms, colClasses = "character")
if (!all(c("go_id", "term_name", "namespace") %in% names(go_terms))) {
  stop("GO term table must contain go_id, term_name, and namespace.")
}
term_check <- merge(pathways, go_terms, by = "go_id", all.x = TRUE, sort = FALSE)
term_check[, pathway_order := match(signature, pathways$signature)]
setorder(term_check, pathway_order)
term_check[, pathway_order := NULL]
if (term_check[, anyNA(term_name) | anyNA(namespace)] || term_check[, any(namespace != "biological_process")]) {
  stop("One or more prespecified GO BP terms are missing from the local snapshot.")
}
if (term_check[, any(tolower(term_name) != tolower(expected_term_name))]) {
  stop("GO term-name mismatch in the local snapshot.")
}
fwrite(term_check, file.path(qc_dir, "GO_pathway_definitions.tsv"), sep = "\t")

go_map_propagated <- read_table(cfg$go_mapping)
if (!all(c("gene_symbol", "go_id") %in% names(go_map_propagated))) {
  stop("Propagated GO mapping must contain gene_symbol and go_id.")
}
go_map_propagated <- unique(go_map_propagated[, .(
  gene_symbol = toupper(trimws(gene_symbol)),
  go_id
)])

go_map_direct <- read_table(cfg$go_direct_mapping)
if (!all(c("gene_symbol", "direct_go_id") %in% names(go_map_direct))) {
  stop("Direct GO mapping must contain gene_symbol and direct_go_id.")
}
go_map_direct <- unique(go_map_direct[, .(
  gene_symbol = toupper(trimws(gene_symbol)),
  go_id = direct_go_id
)])

go_members <- rbindlist(lapply(seq_len(nrow(pathways)), function(i) {
  pathway_row <- pathways[i]
  source_map <- if (pathway_row$membership_source == "direct") go_map_direct else go_map_propagated
  unique(source_map[go_id == pathway_row$go_id, .(
    go_id,
    gene_symbol,
    signature = pathway_row$signature,
    expected_term_name = pathway_row$expected_term_name,
    membership_source = pathway_row$membership_source
  )])
}))
if (go_members[, anyNA(signature)] || go_members[, any(!nzchar(gene_symbol))]) {
  stop("Invalid target GO membership rows were found.")
}

annotation <- fread(cfg$annotation, colClasses = "character", na.strings = c("", "NA"))
required_ann <- c("ensembl_gene_id", "hgnc_symbol", "gene_biotype")
missing_ann <- setdiff(required_ann, names(annotation))
if (length(missing_ann)) stop("Annotation table is missing: ", paste(missing_ann, collapse = ", "))
annotation[, `:=`(
  ensembl_gene_id = strip_version(ensembl_gene_id),
  hgnc_symbol = toupper(trimws(hgnc_symbol))
)]
ann_candidates <- unique(annotation[
  gene_biotype == "protein_coding" & nzchar(hgnc_symbol) & ensembl_gene_id %in% rownames(tpm_matrix),
  .(gene_symbol = hgnc_symbol, ensembl_gene_id)
])

symbol_resolution <- ann_candidates[, {
  ids <- sort(unique(ensembl_gene_id))
  list(
    candidate_count = length(ids),
    candidate_ensembl_ids = paste(ids, collapse = ";"),
    resolved_gene_id = if (length(ids) == 1L) ids[[1L]] else NA_character_
  )
}, by = gene_symbol]

mapping_qc <- merge(go_members, symbol_resolution, by = "gene_symbol", all.x = TRUE, sort = FALSE)
mapping_qc[, mapping_status := fcase(
  is.na(candidate_count), "unresolved_symbol",
  candidate_count == 1L, "unique_symbol",
  default = "ambiguous_symbol"
)]
mapping_qc[, `:=`(
  detection_fraction = NA_real_,
  mean_TPM = NA_real_,
  median_TPM = NA_real_,
  row_variance = NA_real_,
  inclusion_status = "excluded_unresolved_or_ambiguous"
)]
mapped_index <- which(!is.na(mapping_qc$resolved_gene_id))
if (length(mapped_index)) {
  ids <- mapping_qc$resolved_gene_id[mapped_index]
  mapping_qc[mapped_index, detection_fraction := gene_detection_fraction[ids]]
  mapping_qc[mapped_index, mean_TPM := matrixStats::rowMeans2(tpm_matrix[ids, , drop = FALSE])]
  mapping_qc[mapped_index, median_TPM := matrixStats::rowMedians(tpm_matrix[ids, , drop = FALSE])]
  mapping_qc[mapped_index, row_variance := gene_variance[ids]]
  mapping_qc[mapped_index, inclusion_status := fcase(
    detection_fraction == 0, "excluded_zero_expression",
    detection_fraction < 0.20, "excluded_low_prevalence",
    !is.finite(row_variance) | row_variance <= 0, "excluded_zero_variance",
    default = "included"
  )]
}
mapping_qc[, pathway_order := match(signature, pathways$signature)]
setorder(mapping_qc, pathway_order, gene_symbol)
mapping_qc[, pathway_order := NULL]
fwrite(mapping_qc, file.path(qc_dir, "GO_pathway_gene_mapping_QC.tsv"), sep = "\t")

gene_sets <- setNames(vector("list", nrow(pathways)), pathways$signature)
for (signature_name in pathways$signature) {
  genes <- unique(mapping_qc[signature == signature_name & inclusion_status == "included", resolved_gene_id])
  gene_sets[[signature_name]] <- genes
}

gene_set_qc <- mapping_qc[, .(
  GO_members = uniqueN(gene_symbol),
  uniquely_resolved = uniqueN(gene_symbol[mapping_status == "unique_symbol"]),
  included_genes = uniqueN(resolved_gene_id[inclusion_status == "included"]),
  unresolved_symbols = sum(mapping_status == "unresolved_symbol"),
  ambiguous_symbols = sum(mapping_status == "ambiguous_symbol"),
  excluded_low_prevalence = sum(inclusion_status == "excluded_low_prevalence"),
  excluded_zero_expression = sum(inclusion_status == "excluded_zero_expression"),
  excluded_zero_variance = sum(inclusion_status == "excluded_zero_variance")
), by = .(signature, go_id, term_name = expected_term_name, membership_source)]
gene_set_qc[, mapping_coverage := included_genes / GO_members]
gene_set_qc[, pathway_order := match(signature, pathways$signature)]
setorder(gene_set_qc, pathway_order)
gene_set_qc[, pathway_order := NULL]
fwrite(gene_set_qc, file.path(qc_dir, "GO_pathway_gene_set_summary.tsv"), sep = "\t")
if (any(lengths(gene_sets) < 10L)) {
  stop("At least one effective GO gene set contains fewer than 10 genes.")
}
message(
  "Effective GO gene sets: ",
  paste(names(gene_sets), lengths(gene_sets), sep = "=", collapse = ", ")
)

tpm_background <- tpm_matrix[background_genes, , drop = FALSE]
log_tpm <- log2(tpm_background + 1)
rm(tpm_matrix, tpm_background, go_map_propagated, go_map_direct, annotation, ann_candidates)
invisible(gc())

# -----------------------------------------------------------------------------
# Three scoring methods
# -----------------------------------------------------------------------------

message("Ranking genes and calculating singscore for four GO pathways...")
ranked <- singscore::rankGenes(log_tpm)

extract_singscore <- function(scored) {
  score_col <- intersect(c("TotalScore", "UpScore"), names(scored))
  if (!length(score_col)) score_col <- grep("Score$", names(scored), value = TRUE)[1L]
  if (!length(score_col) || is.na(score_col)) stop("Could not identify the singscore score column.")
  ans <- as.numeric(scored[[score_col[[1L]]]])
  names(ans) <- rownames(scored)
  ans
}

score_singscore <- function(genes, signature_name) {
  up_gene_set <- GSEABase::GeneSet(genes, setName = signature_name)
  scored <- singscore::simpleScore(
    ranked,
    upSet = up_gene_set,
    centerScore = TRUE,
    knownDirection = TRUE
  )
  extract_singscore(scored)[selected_ids]
}

singscore_values <- lapply(names(gene_sets), function(signature_name) {
  score_singscore(gene_sets[[signature_name]], signature_name)
})
names(singscore_values) <- names(gene_sets)

message("Calculating normalized ssGSEA for four GO pathways...")
ss_param <- GSVA::ssgseaParam(
  exprData = log_tpm,
  geneSets = gene_sets,
  minSize = 10L,
  maxSize = Inf,
  alpha = 0.25,
  normalize = TRUE
)
ssgsea_matrix <- GSVA::gsva(
  ss_param,
  verbose = TRUE,
  BPPARAM = BiocParallel::SerialParam(progressbar = FALSE)
)

message("Calculating mean gene z-scores for four GO pathways...")
mean_z_score <- function(genes) {
  zmat <- t(scale(t(log_tpm[genes, , drop = FALSE])))
  ans <- colMeans(zmat)
  ans[selected_ids]
}
mean_z_values <- lapply(gene_sets, mean_z_score)

score_tables <- list()
counter <- 1L
for (signature_name in pathways$signature) {
  method_vectors <- list(
    singscore = singscore_values[[signature_name]],
    ssGSEA = setNames(as.numeric(ssgsea_matrix[signature_name, selected_ids]), selected_ids),
    mean_z = mean_z_values[[signature_name]]
  )
  for (method_name in methods) {
    values <- method_vectors[[method_name]][selected_ids]
    score_tables[[counter]] <- data.table(
      sample_id = selected_ids,
      patient_id = selected_meta$patient_id,
      condition = as.character(selected_meta$condition),
      signature = signature_name,
      method = method_name,
      raw_score = as.numeric(values),
      standardized_score = safe_z(as.numeric(values))
    )
    counter <- counter + 1L
  }
}
scores <- rbindlist(score_tables)
if (scores[, any(!is.finite(raw_score) | !is.finite(standardized_score))]) {
  stop("Non-finite values were found in the score table.")
}
if (anyDuplicated(scores, by = c("sample_id", "signature", "method"))) {
  stop("Duplicate sample/signature/method rows were found in the score table.")
}
expected_score_rows <- length(selected_ids) * nrow(pathways) * length(methods)
if (nrow(scores) != expected_score_rows) stop("Unexpected number of score rows.")
fwrite(scores, file.path(results_dir, "GO_pathway_sample_scores_long.tsv.gz"), sep = "\t")

group_summaries <- scores[, .(
  n = .N,
  mean = mean(raw_score),
  sd = sd(raw_score),
  median = median(raw_score),
  q1 = quantile(raw_score, 0.25),
  q3 = quantile(raw_score, 0.75)
), by = .(signature, method, group = condition)]
group_summaries[, group := factor(group, levels = groups)]
group_summaries[, `:=`(
  pathway_order = match(signature, pathways$signature),
  method_order = match(method, methods)
)]
setorder(group_summaries, pathway_order, method_order, group)
group_summaries[, c("pathway_order", "method_order") := NULL]
fwrite(group_summaries, file.path(results_dir, "GO_pathway_group_score_summaries.tsv"), sep = "\t")

# -----------------------------------------------------------------------------
# Patient-aware models and Tukey-adjusted pairwise comparisons
# -----------------------------------------------------------------------------

fit_score_model <- function(score_dt, signature_name, method_name) {
  d <- copy(score_dt)
  d[, condition_model := factor(condition, levels = groups)]
  score_sd <- sd(d$raw_score)
  if (!is.finite(score_sd) || score_sd == 0) stop("Zero score variance for ", signature_name, "/", method_name)

  mixed_fit <- suppressWarnings(lmerTest::lmer(
    raw_score ~ condition_model + (1 | patient_id),
    data = d,
    REML = TRUE,
    control = lme4::lmerControl(
      optimizer = "bobyqa",
      optCtrl = list(maxfun = 200000),
      check.conv.singular = "ignore"
    )
  ))
  convergence_messages <- mixed_fit@optinfo$conv$lme4$messages
  singular <- lme4::isSingular(mixed_fit, tol = 1e-6)
  use_fallback <- singular || length(convergence_messages)

  if (!use_fallback) {
    model_type <- "linear_mixed_model_patient_random_intercept"
    anova_tab <- suppressMessages(anova(mixed_fit, ddf = "Satterthwaite"))
    global_p <- anova_tab["condition_model", "Pr(>F)"]
    emm <- emmeans::emmeans(mixed_fit, ~ condition_model, lmer.df = "satterthwaite")
  } else {
    model_type <- "cluster_robust_linear_model_fallback"
    linear_fit <- lm(raw_score ~ condition_model, data = d)
    robust_vcov <- sandwich::vcovCL(linear_fit, cluster = d$patient_id, type = "HC2")
    wald <- suppressMessages(lmtest::waldtest(linear_fit, . ~ 1, vcov = robust_vcov, test = "Chisq"))
    global_col <- grep("Pr", names(wald), value = TRUE)[1L]
    global_p <- wald[[global_col]][nrow(wald)]
    emm <- emmeans::emmeans(linear_fit, ~ condition_model, vcov. = robust_vcov, df = Inf)
  }

  contrast_raw <- summary(
    emmeans::contrast(emm, method = "pairwise", adjust = "none"),
    infer = c(TRUE, TRUE), adjust = "none"
  )
  contrast_adj <- summary(
    emmeans::contrast(emm, method = "pairwise", adjust = "tukey"),
    infer = c(TRUE, TRUE), adjust = "tukey"
  )
  raw_dt <- as.data.table(contrast_raw)
  adj_dt <- as.data.table(contrast_adj)
  lower_col <- intersect(c("lower.CL", "asymp.LCL"), names(adj_dt))[[1L]]
  upper_col <- intersect(c("upper.CL", "asymp.UCL"), names(adj_dt))[[1L]]
  split_groups <- tstrsplit(adj_dt$contrast, " - ", fixed = TRUE)

  pairwise <- data.table(
    signature = signature_name,
    method = method_name,
    contrast = adj_dt$contrast,
    group1 = split_groups[[1L]],
    group2 = split_groups[[2L]],
    estimate = adj_dt$estimate,
    SE = adj_dt$SE,
    df = adj_dt$df,
    lower_CI = adj_dt[[lower_col]],
    upper_CI = adj_dt[[upper_col]],
    p_raw = raw_dt$p.value[match(adj_dt$contrast, raw_dt$contrast)],
    p_adjusted = adj_dt$p.value,
    standardized_estimate = adj_dt$estimate / score_sd,
    standardized_lower_CI = adj_dt[[lower_col]] / score_sd,
    standardized_upper_CI = adj_dt[[upper_col]] / score_sd,
    global_p = as.numeric(global_p),
    model = model_type,
    singular_mixed_model = singular,
    convergence_message = if (length(convergence_messages)) paste(convergence_messages, collapse = " | ") else ""
  )

  global <- data.table(
    signature = signature_name,
    method = method_name,
    n_samples = nrow(d),
    n_patients = uniqueN(d$patient_id),
    global_p = as.numeric(global_p),
    model = model_type,
    singular_mixed_model = singular,
    convergence_message = if (length(convergence_messages)) paste(convergence_messages, collapse = " | ") else ""
  )
  list(pairwise = pairwise, global = global)
}

message("Fitting 12 patient-aware four-group models...")
model_pairwise <- list()
model_global <- list()
counter <- 1L
for (signature_name in pathways$signature) {
  for (method_name in methods) {
    ans <- fit_score_model(
      scores[signature == signature_name & method == method_name],
      signature_name,
      method_name
    )
    model_pairwise[[counter]] <- ans$pairwise
    model_global[[counter]] <- ans$global
    counter <- counter + 1L
  }
}
pairwise_stats <- rbindlist(model_pairwise)
global_stats <- rbindlist(model_global)
fwrite(pairwise_stats, file.path(results_dir, "GO_pathway_pairwise_model_statistics.tsv"), sep = "\t")
fwrite(global_stats, file.path(results_dir, "GO_pathway_global_model_statistics.tsv"), sep = "\t")

# Cross-method concordance is retained as a compact method-sensitivity check.
correlations <- rbindlist(lapply(pathways$signature, function(signature_name) {
  wide <- dcast(
    scores[signature == signature_name],
    sample_id ~ method,
    value.var = "raw_score"
  )
  cor_matrix <- cor(as.matrix(wide[, ..methods]), method = "spearman")
  out <- as.data.table(as.table(cor_matrix))
  setnames(out, c("method1", "method2", "spearman_rho"))
  out[, signature := signature_name]
  out[]
}))
setcolorder(correlations, c("signature", "method1", "method2", "spearman_rho"))
fwrite(correlations, file.path(results_dir, "GO_pathway_method_spearman_correlations.tsv"), sep = "\t")

# -----------------------------------------------------------------------------
# Figure 1D-style exports
# -----------------------------------------------------------------------------

format_raw_p <- function(p) {
  if (is.na(p)) return("P=NA")
  if (p < 1e-15) return("P<1e-15")
  if (p < 1e-3) return(paste0("P=", formatC(p, format = "e", digits = 1)))
  paste0("P=", formatC(p, format = "f", digits = 3))
}

save_figure <- function(plot, stem, width_mm = 183, height_mm = 82, dpi = 600) {
  width_in <- width_mm / 25.4
  height_in <- height_mm / 25.4

  svglite::svglite(paste0(stem, ".svg"), width = width_in, height = height_in, bg = "white")
  print(plot)
  grDevices::dev.off()

  grDevices::cairo_pdf(
    paste0(stem, ".pdf"), width = width_in, height = height_in,
    family = "Arial", bg = "white"
  )
  print(plot)
  grDevices::dev.off()

  ragg::agg_tiff(
    paste0(stem, ".tiff"), width = width_in, height = height_in,
    units = "in", res = dpi, compression = "lzw", background = "white"
  )
  print(plot)
  grDevices::dev.off()

  ragg::agg_png(
    paste0(stem, ".png"), width = width_in, height = height_in,
    units = "in", res = 300, background = "white"
  )
  print(plot)
  grDevices::dev.off()
}

comparison_groups <- c("HR", "HER2", "NORMAL")
bracket_tick <- 0.12

make_figure <- function(signature_name) {
  pathway_row <- pathways[signature == signature_name]
  plot_dt <- copy(scores[signature == signature_name])
  plot_dt[, `:=`(
    condition_factor = factor(condition, levels = groups),
    method_factor = factor(method, levels = methods, labels = unname(method_labels[methods]))
  )]

  stat_dt <- copy(pairwise_stats[
    signature == signature_name & group1 == "TN" & group2 %chin% comparison_groups
  ])
  if (nrow(stat_dt) != length(methods) * length(comparison_groups)) {
    stop("Expected three TN comparisons for each method in ", signature_name)
  }
  stat_dt[, `:=`(
    method_factor = factor(method, levels = methods, labels = unname(method_labels[methods])),
    x1 = 1,
    x2 = match(group2, groups),
    y = 3.30 + 0.60 * (match(group2, comparison_groups) - 1),
    label = vapply(p_raw, format_raw_p, character(1))
  )]

  p <- ggplot(
    plot_dt,
    aes(x = condition_factor, y = standardized_score, fill = condition_factor)
  ) +
    geom_boxplot(
      width = 0.58,
      linewidth = 0.55,
      colour = "#222222",
      alpha = 0.76,
      outlier.shape = 16,
      outlier.size = 1.15,
      outlier.stroke = 0,
      outlier.colour = "#222222"
    ) +
    geom_segment(
      data = stat_dt,
      aes(x = x1, xend = x2, y = y, yend = y),
      inherit.aes = FALSE,
      linewidth = 0.45,
      colour = "#222222"
    ) +
    geom_segment(
      data = stat_dt,
      aes(x = x1, xend = x1, y = y, yend = y - bracket_tick),
      inherit.aes = FALSE,
      linewidth = 0.45,
      colour = "#222222"
    ) +
    geom_segment(
      data = stat_dt,
      aes(x = x2, xend = x2, y = y, yend = y - bracket_tick),
      inherit.aes = FALSE,
      linewidth = 0.45,
      colour = "#222222"
    ) +
    geom_text(
      data = stat_dt,
      aes(x = (x1 + x2) / 2, y = y + 0.08, label = label),
      inherit.aes = FALSE,
      family = "Arial",
      size = 10 / 2.845276,
      vjust = 0,
      colour = "#1F1F1F"
    ) +
    facet_wrap(~ method_factor, nrow = 1) +
    scale_fill_manual(values = group_colors, guide = "none", drop = FALSE) +
    scale_x_discrete(labels = group_labels, drop = FALSE) +
    scale_y_continuous(breaks = c(-3, 0, 3), expand = expansion(mult = c(0.015, 0.025))) +
    coord_cartesian(ylim = c(cfg$ymin, cfg$ymax), clip = "on") +
    labs(x = NULL, y = pathway_row$y_label) +
    theme_classic(base_size = 10.5, base_family = "Arial") +
    theme(
      axis.line = element_line(linewidth = 0.50, colour = "#222222"),
      axis.ticks = element_line(linewidth = 0.45, colour = "#222222"),
      axis.ticks.length = grid::unit(1.6, "mm"),
      axis.text = element_text(size = 10, colour = "#222222"),
      axis.title.y = element_text(size = 10.5, colour = "#222222", margin = margin(r = 7)),
      strip.background = element_rect(fill = "white", colour = "#222222", linewidth = 0.50),
      strip.text = element_text(
        size = 10.5, face = "bold", colour = "#222222",
        margin = margin(3.5, 0, 3.5, 0)
      ),
      panel.spacing.x = grid::unit(2.5, "mm"),
      panel.grid = element_blank(),
      plot.title = element_blank(),
      plot.subtitle = element_blank(),
      plot.margin = margin(4, 7, 6, 7)
    )

  stem <- file.path(figures_dir, pathway_row$figure_prefix)
  save_figure(p, stem)

  fwrite(
    plot_dt[, .(
      sample_id,
      patient_id,
      condition = as.character(condition_factor),
      signature,
      method = as.character(method_factor),
      raw_score,
      standardized_score
    )],
    file.path(results_dir, paste0(pathway_row$figure_prefix, "-source-data.tsv.gz")),
    sep = "\t"
  )
  fwrite(
    stat_dt[, .(
      signature,
      method,
      contrast,
      group1,
      group2,
      estimate,
      lower_CI,
      upper_CI,
      p_raw,
      p_adjusted,
      label,
      model
    )],
    file.path(results_dir, paste0(pathway_row$figure_prefix, "-statistics.tsv")),
    sep = "\t"
  )
  invisible(stem)
}

message("Exporting four Figure 1D-style pathway-score figures...")
figure_stems <- vapply(pathways$signature, make_figure, character(1))

figure_parameters <- merge(
  pathways[, .(signature, go_id, expected_term_name, membership_source, figure_prefix, y_label)],
  gene_set_qc[, .(signature, GO_members, included_genes, mapping_coverage)],
  by = "signature",
  sort = FALSE
)
figure_parameters[, `:=`(
  width_mm = 183,
  height_mm = 82,
  y_display_min = cfg$ymin,
  y_display_max = cfg$ymax,
  minimum_font_pt = 10,
  comparisons = "TN-HR;TN-HER2;TN-NORMAL",
  p_value_display = "raw (unadjusted)",
  multiple_testing = "none for displayed P values",
  TN_color = group_colors[["TN"]],
  HR_color = group_colors[["HR"]],
  HER2_color = group_colors[["HER2"]],
  Normal_color = group_colors[["NORMAL"]]
)]
figure_parameters[, pathway_order := match(signature, pathways$signature)]
setorder(figure_parameters, pathway_order)
figure_parameters[, pathway_order := NULL]
fwrite(figure_parameters, file.path(results_dir, "GO_pathway_figure_parameters.tsv"), sep = "\t")

figure_files <- list.files(figures_dir, pattern = "^Figure1-D-.*-score\\.(svg|pdf|tiff|png)$", full.names = TRUE)
figure_manifest <- data.table(
  file = basename(figure_files),
  bytes = file.info(figure_files)$size,
  format = toupper(sub(".*\\.", "", figure_files))
)
setorder(figure_manifest, file)
fwrite(figure_manifest, file.path(qc_dir, "figure_export_manifest.tsv"), sep = "\t")

package_manifest <- data.table(
  package = required_packages,
  version = vapply(required_packages, function(x) as.character(packageVersion(x)), character(1)),
  library_path = vapply(required_packages, function(x) dirname(find.package(x)), character(1))
)
fwrite(package_manifest, file.path(qc_dir, "package_versions.tsv"), sep = "\t")
capture.output(sessionInfo(), file = file.path(qc_dir, "sessionInfo.txt"))

validation <- data.table(
  check = c(
    "fixed_sample_count_1116",
    "fixed_group_counts",
    "four_GO_definitions_resolved",
    "GO_membership_source_policy",
    "all_gene_sets_at_least_10_genes",
    "all_12_scores_per_sample_finite",
    "score_key_unique",
    "twelve_patient_aware_models",
    "nine_TN_comparisons_per_pathway_figure_total_36",
    "four_figures_four_formats"
  ),
  passed = c(
    length(selected_ids) == 1116L,
    identical(as.integer(observed_vector), as.integer(expected_counts)),
    nrow(term_check) == 4L && !anyNA(term_check$term_name),
    identical(pathways$membership_source, c("propagated", "propagated", "direct", "propagated")),
    all(lengths(gene_sets) >= 10L),
    nrow(scores) == expected_score_rows && !scores[, any(!is.finite(raw_score) | !is.finite(standardized_score))],
    !anyDuplicated(scores, by = c("sample_id", "signature", "method")),
    nrow(global_stats) == 12L,
    pairwise_stats[group1 == "TN" & group2 %chin% comparison_groups, .N] == 36L,
    nrow(figure_manifest) == 16L
  ),
  value = c(
    length(selected_ids),
    paste(groups, observed_vector, sep = "=", collapse = ";"),
    paste(pathways$go_id, term_check$term_name, sep = "=", collapse = ";"),
    paste(pathways$signature, pathways$membership_source, sep = "=", collapse = ";"),
    paste(names(gene_sets), lengths(gene_sets), sep = "=", collapse = ";"),
    nrow(scores),
    anyDuplicated(scores, by = c("sample_id", "signature", "method")),
    nrow(global_stats),
    pairwise_stats[group1 == "TN" & group2 %chin% comparison_groups, .N],
    nrow(figure_manifest)
  )
)
fwrite(validation, file.path(qc_dir, "validation_report.tsv"), sep = "\t")
if (!all(validation$passed)) stop("One or more validation checks failed. Inspect QC/validation_report.tsv.")


message("Analysis completed: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"))
message("Validation: PASS (", nrow(validation), "/", nrow(validation), ")")
message("Figure prefixes: ", paste(basename(figure_stems), collapse = ", "))
