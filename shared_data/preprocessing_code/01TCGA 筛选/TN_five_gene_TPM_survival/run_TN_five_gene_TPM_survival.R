#!/usr/bin/env Rscript

# Patient-level survival analysis of DKC1, GTPBP4, ILF2, NOP2 and RPP38
# in the receptor-defined TCGA-BRCA TN cohort.
#
# Primary analysis:
#   OS ~ continuous log2(TPM + 1), scaled per one cohort SD.
# Secondary analyses:
#   DSS, PFS and DFS; median-split Kaplan-Meier/log-rank; and an exploratory
#   Cox model adjusted for age and collapsed pathologic stage.

args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
script_path <- if (length(file_arg)) sub("^--file=", "", file_arg[[1]]) else "."
script_path <- gsub("~\\+~", " ", script_path)
bundle_dir <- normalizePath(dirname(script_path), mustWork = TRUE)
project_dir <- normalizePath(file.path(bundle_dir, ".."), mustWork = TRUE)
workspace_dir <- normalizePath(file.path(project_dir, ".."), mustWork = TRUE)
tcga_dir <- file.path(workspace_dir, "TCGA")

local_lib <- file.path(tcga_dir, ".Rlib")
if (dir.exists(local_lib)) .libPaths(c(local_lib, .libPaths()))

suppressPackageStartupMessages({
  library(data.table)
  library(survival)
})

Sys.setenv(
  OMP_NUM_THREADS = "1",
  OPENBLAS_NUM_THREADS = "1",
  MKL_NUM_THREADS = "1",
  VECLIB_MAXIMUM_THREADS = "1"
)
options(mc.cores = 1L, width = 180)
data.table::setDTthreads(1L)
set.seed(20260903)

genes <- c("DKC1", "GTPBP4", "ILF2", "NOP2", "RPP38")
primary_endpoint <- "OS"
endpoint_order <- c("OS", "DSS", "PFS", "DFS")

data_dir <- file.path(bundle_dir, "data")
results_dir <- file.path(bundle_dir, "results")
qc_dir <- file.path(bundle_dir, "QC")
for (d in c(data_dir, results_dir, qc_dir)) dir.create(d, recursive = TRUE, showWarnings = FALSE)

tpm_path <- file.path(tcga_dir, "TCGA-BRCA_TPM.txt")
sample_meta_path <- file.path(tcga_dir, "DESeq2_four_group_marker_analysis", "metadata", "selected_samples.tsv")
strict_receptor_path <- file.path(tcga_dir, "TCGA-BRCA_receptor_subtype_by_sample.tsv")
clinical_path <- file.path(
  tcga_dir,
  "external_receptor_sources",
  "cbioportal_brca_tcga_pan_can_atlas_2018_data_clinical_patient.txt"
)
annotation_path <- file.path(
  tcga_dir,
  "DESeq2_TN_RBP_analysis",
  "KEGG_enrichment_top500_TN_all_genes",
  "reference_snapshot",
  "Ensembl_BioMart_hsapiens_gene_annotation.tsv"
)

required <- c(tpm_path, sample_meta_path, strict_receptor_path, clinical_path, annotation_path)
if (any(!file.exists(required))) {
  stop("Missing required input(s): ", paste(required[!file.exists(required)], collapse = ", "))
}

strip_ensembl_version <- function(x) sub("\\..*$", "", as.character(x))

status_to_event <- function(x) {
  x <- trimws(as.character(x))
  out <- rep(NA_integer_, length(x))
  out[grepl("^0:", x)] <- 0L
  out[grepl("^1:", x)] <- 1L
  out
}

collapse_stage <- function(x) {
  x <- toupper(trimws(as.character(x)))
  fcase(
    grepl("^STAGE IV", x), "IV",
    grepl("^STAGE III", x), "III",
    grepl("^STAGE II", x), "II",
    grepl("^STAGE I($|A|B|C)", x), "I",
    default = NA_character_
  )
}

message("Reading the patient-deduplicated four-group sample metadata...")
sample_meta <- fread(sample_meta_path, showProgress = FALSE)
required_meta <- c("sample_id", "patient_id", "condition", "sample_type", "library_size")
if (!all(required_meta %in% names(sample_meta))) stop("Sample metadata is missing required columns")

strict_receptor <- fread(strict_receptor_path, showProgress = FALSE)
required_receptor <- c(
  "sample_id", "patient_id", "sample_type", "ER_status_selected",
  "PR_status_selected", "HER2_status_selected", "label_recommended"
)
if (!all(required_receptor %in% names(strict_receptor))) stop("Strict receptor metadata is missing required columns")
strict_receptor <- strict_receptor[
  sample_type == "PRIMARY_TUMOR" & label_recommended == "TN"
]
if (
  nrow(strict_receptor) != 169L || uniqueN(strict_receptor$patient_id) != 169L ||
    uniqueN(strict_receptor$sample_id) != 169L ||
    !all(strict_receptor$ER_status_selected == "Negative") ||
    !all(strict_receptor$PR_status_selected == "Negative") ||
    !all(strict_receptor$HER2_status_selected == "Negative")
) {
  stop("Expected 169 unique strict ER-/PR-/HER2- primary-tumor patients/samples")
}

tn_meta <- copy(sample_meta[sample_id %chin% strict_receptor$sample_id & sample_type == "PRIMARY_TUMOR"])
setorder(tn_meta, patient_id, -library_size, sample_id)
tn_meta <- tn_meta[, .SD[1L], by = patient_id]
setcolorder(tn_meta, c("patient_id", "sample_id", setdiff(names(tn_meta), c("patient_id", "sample_id"))))
if (
  nrow(tn_meta) != 169L || uniqueN(tn_meta$patient_id) != 169L ||
    uniqueN(tn_meta$sample_id) != 169L || !setequal(tn_meta$sample_id, strict_receptor$sample_id)
) {
  stop("Expected 169 unique strict TN patients/samples, recovered ", nrow(tn_meta))
}

message("Mapping the five gene symbols to Ensembl identifiers...")
annotation <- fread(annotation_path, colClasses = "character", showProgress = FALSE)
gene_annotation <- unique(annotation[
  toupper(hgnc_symbol) %chin% genes,
  .(
    gene_symbol = toupper(hgnc_symbol),
    ensembl_gene_id,
    entrezgene_id,
    gene_biotype
  )
])
gene_annotation[, requested_gene_order := match(gene_symbol, genes)]
setorder(gene_annotation, requested_gene_order)
gene_annotation[, requested_gene_order := NULL]
if (nrow(gene_annotation) != length(genes) || !setequal(gene_annotation$gene_symbol, genes)) {
  stop("The five requested genes did not map one-to-one in the local annotation snapshot")
}
fwrite(gene_annotation, file.path(data_dir, "TN_five_gene_annotation.tsv"), sep = "\t", na = "NA")

message("Extracting TPM values for the 169 strict TN tumor samples...")
tpm_subset <- suppressWarnings(fread(
  tpm_path,
  select = c("V1", tn_meta$sample_id),
  check.names = FALSE,
  showProgress = FALSE
))
setnames(tpm_subset, "V1", "ensembl_gene_id_version")
tpm_subset[, ensembl_gene_id := strip_ensembl_version(ensembl_gene_id_version)]
tpm_subset <- tpm_subset[ensembl_gene_id %chin% gene_annotation$ensembl_gene_id]
if (nrow(tpm_subset) != length(genes) || anyDuplicated(tpm_subset$ensembl_gene_id)) {
  stop("Expected one TPM row for each requested Ensembl gene")
}

tpm_long <- melt(
  tpm_subset,
  id.vars = c("ensembl_gene_id_version", "ensembl_gene_id"),
  measure.vars = tn_meta$sample_id,
  variable.name = "sample_id",
  value.name = "TPM",
  variable.factor = FALSE
)
tpm_long <- merge(tpm_long, gene_annotation, by = "ensembl_gene_id", all.x = TRUE, sort = FALSE)
if (any(!is.finite(tpm_long$TPM)) || any(tpm_long$TPM < 0)) stop("Invalid TPM value detected")
if (nrow(tpm_long) != nrow(tn_meta) * length(genes)) stop("Unexpected TPM long-table dimensions")

tpm_wide <- dcast(tpm_long, sample_id ~ gene_symbol, value.var = "TPM")
setcolorder(tpm_wide, c("sample_id", genes))
setnames(tpm_wide, genes, paste0(genes, "_TPM"))
for (gene in genes) {
  raw_col <- paste0(gene, "_TPM")
  log_col <- paste0(gene, "_log2TPM1")
  z_col <- paste0(gene, "_log2TPM1_z")
  group_col <- paste0(gene, "_median_group")
  tpm_wide[, (log_col) := log2(get(raw_col) + 1)]
  tpm_wide[, (z_col) := as.numeric(scale(get(log_col)))]
  gene_median <- median(tpm_wide[[log_col]], na.rm = TRUE)
  tpm_wide[, (group_col) := factor(
    ifelse(get(log_col) >= gene_median, "High", "Low"),
    levels = c("Low", "High")
  )]
}
tpm_wide <- merge(tn_meta, tpm_wide, by = "sample_id", all.x = TRUE, sort = FALSE)
setorder(tpm_wide, patient_id)
if (anyNA(tpm_wide[, paste0(genes, "_TPM"), with = FALSE])) stop("Missing TPM after TN sample merge")

message("Reading standardized TCGA PanCancer Atlas survival endpoints...")
clinical <- fread(clinical_path, skip = 4L, header = TRUE, fill = TRUE, showProgress = FALSE)
if (anyDuplicated(clinical$PATIENT_ID)) stop("Clinical source contains duplicated PATIENT_ID values")

clinical_columns <- c(
  "PATIENT_ID", "AGE", "SEX", "AJCC_PATHOLOGIC_TUMOR_STAGE",
  "PATH_T_STAGE", "PATH_N_STAGE", "PATH_M_STAGE",
  "OS_STATUS", "OS_MONTHS", "DSS_STATUS", "DSS_MONTHS",
  "PFS_STATUS", "PFS_MONTHS", "DFS_STATUS", "DFS_MONTHS"
)
if (!all(clinical_columns %in% names(clinical))) stop("Clinical source is missing survival/covariate columns")
clinical <- clinical[, ..clinical_columns]
setnames(clinical, "PATIENT_ID", "patient_id")

survival_data <- merge(tpm_wide, clinical, by = "patient_id", all.x = TRUE, sort = FALSE)
setorder(survival_data, patient_id)
survival_data[, clinical_match := !is.na(OS_MONTHS)]
survival_data[, pathologic_stage_group := collapse_stage(AJCC_PATHOLOGIC_TUMOR_STAGE)]
survival_data[, pathologic_stage_binary := factor(
  fcase(
    pathologic_stage_group %chin% c("I", "II"), "I-II",
    pathologic_stage_group %chin% c("III", "IV"), "III-IV",
    default = NA_character_
  ),
  levels = c("I-II", "III-IV")
)]
survival_data[, age_per_10y := AGE / 10]

for (endpoint in endpoint_order) {
  status_col <- paste0(endpoint, "_STATUS")
  months_col <- paste0(endpoint, "_MONTHS")
  event_col <- paste0(endpoint, "_event")
  time_col <- paste0(endpoint, "_time_months")
  survival_data[, (event_col) := status_to_event(get(status_col))]
  survival_data[, (time_col) := suppressWarnings(as.numeric(get(months_col)))]
  survival_data[get(time_col) < 0, (time_col) := NA_real_]
}

fwrite(tn_meta, file.path(data_dir, "TN_169_selected_primary_tumor_samples.tsv"), sep = "\t", na = "NA")
fwrite(tpm_long, file.path(data_dir, "TN_five_gene_TPM_long.tsv"), sep = "\t", na = "NA")
fwrite(tpm_wide, file.path(data_dir, "TN_five_gene_TPM_patient_wide.tsv"), sep = "\t", na = "NA")
fwrite(survival_data, file.path(data_dir, "TN_five_gene_TPM_survival_ready_data.tsv"), sep = "\t", na = "NA")
fwrite(survival_data, file.path(data_dir, "TN_five_gene_TPM_survival_ready_data.csv"), na = "NA")

expression_summary <- rbindlist(lapply(genes, function(gene) {
  x <- survival_data[[paste0(gene, "_TPM")]]
  y <- survival_data[[paste0(gene, "_log2TPM1")]]
  data.table(
    gene = gene,
    n = sum(is.finite(x)),
    TPM_min = min(x, na.rm = TRUE),
    TPM_Q1 = quantile(x, 0.25, na.rm = TRUE),
    TPM_median = median(x, na.rm = TRUE),
    TPM_mean = mean(x, na.rm = TRUE),
    TPM_Q3 = quantile(x, 0.75, na.rm = TRUE),
    TPM_max = max(x, na.rm = TRUE),
    log2TPM1_mean = mean(y, na.rm = TRUE),
    log2TPM1_SD = sd(y, na.rm = TRUE)
  )
}))
fwrite(expression_summary, file.path(data_dir, "TN_five_gene_TPM_expression_summary.tsv"), sep = "\t")

endpoint_completeness <- rbindlist(lapply(endpoint_order, function(endpoint) {
  time_col <- paste0(endpoint, "_time_months")
  event_col <- paste0(endpoint, "_event")
  ok <- is.finite(survival_data[[time_col]]) & !is.na(survival_data[[event_col]])
  data.table(
    endpoint = endpoint,
    eligible_patients = sum(ok),
    events = sum(survival_data[[event_col]][ok]),
    censored = sum(survival_data[[event_col]][ok] == 0L),
    missing_or_invalid = nrow(survival_data) - sum(ok),
    median_observed_time_months = median(survival_data[[time_col]][ok], na.rm = TRUE)
  )
}))
fwrite(endpoint_completeness, file.path(data_dir, "TN_survival_endpoint_completeness.tsv"), sep = "\t")

fit_continuous_cox <- function(dt, endpoint, gene) {
  time_col <- paste0(endpoint, "_time_months")
  event_col <- paste0(endpoint, "_event")
  expr_col <- paste0(gene, "_log2TPM1_z")
  model_dt <- dt[, .(
    time = get(time_col),
    event = get(event_col),
    expression_z = get(expr_col),
    age_per_10y,
    stage_binary = pathologic_stage_binary
  )]
  model_dt <- model_dt[is.finite(time) & time >= 0 & !is.na(event) & is.finite(expression_z)]
  if (nrow(model_dt) < 20L || sum(model_dt$event) < 5L) stop("Insufficient endpoint data for ", endpoint, " / ", gene)

  fit <- coxph(Surv(time, event) ~ expression_z, data = model_dt, ties = "efron", x = TRUE)
  coef_row <- summary(fit)$coefficients["expression_z", ]
  ci_row <- summary(fit)$conf.int["expression_z", ]
  ph_p <- tryCatch(
    as.numeric(cox.zph(fit, transform = "km")$table["expression_z", "p"]),
    error = function(e) NA_real_
  )

  adjusted_dt <- model_dt[complete.cases(age_per_10y, stage_binary)]
  adjusted <- list(n = nrow(adjusted_dt), events = sum(adjusted_dt$event), hr = NA_real_, lower = NA_real_, upper = NA_real_, p = NA_real_, ph_p = NA_real_)
  if (nrow(adjusted_dt) >= 20L && sum(adjusted_dt$event) >= 10L) {
    adjusted_fit <- coxph(
      Surv(time, event) ~ expression_z + age_per_10y + stage_binary,
      data = adjusted_dt,
      ties = "efron",
      x = TRUE
    )
    adjusted_coef <- summary(adjusted_fit)$coefficients["expression_z", ]
    adjusted_ci <- summary(adjusted_fit)$conf.int["expression_z", ]
    adjusted <- list(
      n = nrow(adjusted_dt),
      events = sum(adjusted_dt$event),
      hr = unname(adjusted_ci["exp(coef)"]),
      lower = unname(adjusted_ci["lower .95"]),
      upper = unname(adjusted_ci["upper .95"]),
      p = unname(adjusted_coef["Pr(>|z|)"]),
      ph_p = tryCatch(
        as.numeric(cox.zph(adjusted_fit, transform = "km")$table["expression_z", "p"]),
        error = function(e) NA_real_
      )
    )
  }

  data.table(
    endpoint = endpoint,
    gene = gene,
    expression_scale = "per 1 SD of log2(TPM+1)",
    n = nrow(model_dt),
    events = sum(model_dt$event),
    HR = unname(ci_row["exp(coef)"]),
    CI95_lower = unname(ci_row["lower .95"]),
    CI95_upper = unname(ci_row["upper .95"]),
    Wald_P = unname(coef_row["Pr(>|z|)"]),
    PH_assumption_P = ph_p,
    adjusted_covariates = "age per 10 years + pathologic stage I-II vs III-IV",
    adjusted_n = adjusted$n,
    adjusted_events = adjusted$events,
    adjusted_HR = adjusted$hr,
    adjusted_CI95_lower = adjusted$lower,
    adjusted_CI95_upper = adjusted$upper,
    adjusted_Wald_P = adjusted$p,
    adjusted_PH_assumption_P = adjusted$ph_p
  )
}

fit_median_split <- function(dt, endpoint, gene) {
  time_col <- paste0(endpoint, "_time_months")
  event_col <- paste0(endpoint, "_event")
  expr_col <- paste0(gene, "_log2TPM1")
  group_col <- paste0(gene, "_median_group")
  model_dt <- dt[, .(
    time = get(time_col),
    event = get(event_col),
    expression = get(expr_col),
    expression_group = factor(as.character(get(group_col)), levels = c("Low", "High"))
  )]
  model_dt <- model_dt[is.finite(time) & time >= 0 & !is.na(event) & is.finite(expression) & !is.na(expression_group)]
  if (uniqueN(model_dt$expression_group) != 2L) stop("Median split did not produce two groups for ", gene)

  fit <- coxph(Surv(time, event) ~ expression_group, data = model_dt, ties = "efron", x = TRUE)
  coef_name <- "expression_groupHigh"
  coef_row <- summary(fit)$coefficients[coef_name, ]
  ci_row <- summary(fit)$conf.int[coef_name, ]
  ph_p <- tryCatch(
    as.numeric(cox.zph(fit, transform = "km")$table["expression_group", "p"]),
    error = function(e) NA_real_
  )
  lr <- survdiff(Surv(time, event) ~ expression_group, data = model_dt, rho = 0)
  logrank_p <- pchisq(lr$chisq, df = length(lr$n) - 1L, lower.tail = FALSE)
  sf <- survfit(Surv(time, event) ~ expression_group, data = model_dt)
  sf_tab <- as.data.table(as.data.frame(summary(sf)$table), keep.rownames = "stratum")
  low_med <- sf_tab[grepl("Low$", stratum), get("median")]
  high_med <- sf_tab[grepl("High$", stratum), get("median")]

  data.table(
    endpoint = endpoint,
    gene = gene,
    cutoff_rule = "High >= TN-cohort median log2(TPM+1); Low < median",
    cutoff_TPM = median(dt[[paste0(gene, "_TPM")]], na.rm = TRUE),
    cutoff_log2TPM1 = median(dt[[expr_col]], na.rm = TRUE),
    n = nrow(model_dt),
    events = sum(model_dt$event),
    low_n = sum(model_dt$expression_group == "Low"),
    low_events = sum(model_dt$event[model_dt$expression_group == "Low"]),
    high_n = sum(model_dt$expression_group == "High"),
    high_events = sum(model_dt$event[model_dt$expression_group == "High"]),
    HR_high_vs_low = unname(ci_row["exp(coef)"]),
    CI95_lower = unname(ci_row["lower .95"]),
    CI95_upper = unname(ci_row["upper .95"]),
    Cox_Wald_P = unname(coef_row["Pr(>|z|)"]),
    logrank_P = logrank_p,
    PH_assumption_P = ph_p,
    low_median_survival_months = if (length(low_med)) low_med else NA_real_,
    high_median_survival_months = if (length(high_med)) high_med else NA_real_
  )
}

message("Running gene-wise continuous-expression Cox models...")
continuous_results <- rbindlist(lapply(endpoint_order, function(endpoint) {
  rbindlist(lapply(genes, function(gene) fit_continuous_cox(survival_data, endpoint, gene)))
}))
continuous_results[, Wald_FDR_BH_5genes := p.adjust(Wald_P, method = "BH"), by = endpoint]
continuous_results[, adjusted_Wald_FDR_BH_5genes := p.adjust(adjusted_Wald_P, method = "BH"), by = endpoint]
continuous_results[, significant_raw_P_lt_0.05 := Wald_P < 0.05]
continuous_results[, significant_BH_FDR_lt_0.05 := Wald_FDR_BH_5genes < 0.05]
continuous_results[, adjusted_significant_raw_P_lt_0.05 := adjusted_Wald_P < 0.05]
continuous_results[, adjusted_significant_BH_FDR_lt_0.05 := adjusted_Wald_FDR_BH_5genes < 0.05]
continuous_results[, PH_assumption_violation_P_lt_0.05 := PH_assumption_P < 0.05]
continuous_results[, adjusted_PH_assumption_violation_P_lt_0.05 := adjusted_PH_assumption_P < 0.05]
continuous_results[, `:=`(
  endpoint_order_index = match(endpoint, endpoint_order),
  gene_order_index = match(gene, genes)
)]
setorder(continuous_results, endpoint_order_index, gene_order_index)
continuous_results[, c("endpoint_order_index", "gene_order_index") := NULL]

message("Running fixed median-split Kaplan-Meier/log-rank analyses...")
median_results <- rbindlist(lapply(endpoint_order, function(endpoint) {
  rbindlist(lapply(genes, function(gene) fit_median_split(survival_data, endpoint, gene)))
}))
median_results[, logrank_FDR_BH_5genes := p.adjust(logrank_P, method = "BH"), by = endpoint]
median_results[, Cox_Wald_FDR_BH_5genes := p.adjust(Cox_Wald_P, method = "BH"), by = endpoint]
median_results[, significant_logrank_raw_P_lt_0.05 := logrank_P < 0.05]
median_results[, significant_logrank_BH_FDR_lt_0.05 := logrank_FDR_BH_5genes < 0.05]
median_results[, `:=`(
  endpoint_order_index = match(endpoint, endpoint_order),
  gene_order_index = match(gene, genes)
)]
setorder(median_results, endpoint_order_index, gene_order_index)
median_results[, c("endpoint_order_index", "gene_order_index") := NULL]

fwrite(continuous_results, file.path(results_dir, "TN_five_gene_all_endpoints_continuous_Cox.tsv"), sep = "\t", na = "NA")
fwrite(median_results, file.path(results_dir, "TN_five_gene_all_endpoints_median_split.tsv"), sep = "\t", na = "NA")
fwrite(continuous_results[endpoint == primary_endpoint], file.path(results_dir, "TN_five_gene_OS_continuous_Cox.tsv"), sep = "\t", na = "NA")
fwrite(median_results[endpoint == primary_endpoint], file.path(results_dir, "TN_five_gene_OS_median_split_KM_logrank.tsv"), sep = "\t", na = "NA")

os_summary <- merge(
  continuous_results[endpoint == primary_endpoint, .(
    gene, OS_n = n, OS_events = events,
    continuous_HR_per_SD = HR,
    continuous_CI95_lower = CI95_lower,
    continuous_CI95_upper = CI95_upper,
    continuous_Wald_P = Wald_P,
    continuous_BH_FDR = Wald_FDR_BH_5genes,
    continuous_PH_assumption_P = PH_assumption_P,
    adjusted_HR_per_SD = adjusted_HR,
    adjusted_CI95_lower = adjusted_CI95_lower,
    adjusted_CI95_upper = adjusted_CI95_upper,
    adjusted_Wald_P = adjusted_Wald_P,
    adjusted_BH_FDR = adjusted_Wald_FDR_BH_5genes,
    adjusted_PH_assumption_P = adjusted_PH_assumption_P
  )],
  median_results[endpoint == primary_endpoint, .(
    gene, median_cutoff_TPM = cutoff_TPM,
    low_n, low_events, high_n, high_events,
    median_split_HR_high_vs_low = HR_high_vs_low,
    median_split_CI95_lower = CI95_lower,
    median_split_CI95_upper = CI95_upper,
    median_split_logrank_P = logrank_P,
    median_split_logrank_BH_FDR = logrank_FDR_BH_5genes
  )],
  by = "gene",
  all = TRUE,
  sort = FALSE
)
os_summary[, requested_gene_order := match(gene, genes)]
setorder(os_summary, requested_gene_order)
os_summary[, requested_gene_order := NULL]
os_summary[, primary_OS_significant_after_BH := continuous_BH_FDR < 0.05]
os_summary[, interpretation := fifelse(
  primary_OS_significant_after_BH,
  "Significant continuous-expression OS association after BH correction",
  "No significant continuous-expression OS association after BH correction"
)]
fwrite(os_summary, file.path(results_dir, "TN_five_gene_OS_significance_summary.tsv"), sep = "\t", na = "NA")
fwrite(os_summary, file.path(results_dir, "TN_five_gene_OS_significance_summary.csv"), na = "NA")

analysis_parameters <- data.table(
  parameter = c(
    "analysis_population", "TN_definition_source", "sample_grain", "expression_source",
    "expression_transformation", "primary_endpoint", "secondary_endpoints",
    "primary_statistical_model", "secondary_grouping_rule", "secondary_group_test",
    "multiple_testing", "adjusted_sensitivity_model", "clinical_endpoint_source"
  ),
  value = c(
    "Receptor-defined TCGA-BRCA TN primary tumors",
    strict_receptor_path,
    "One primary-tumor RNA aliquot per patient; highest library size resolves duplicates",
    tpm_path,
    "log2(TPM+1), then z-standardized across all 169 strict TN patients",
    "Overall survival (OS)",
    "DSS, PFS and DFS",
    "Univariable Cox proportional-hazards model; HR per 1 SD log2(TPM+1)",
    "Fixed TN-cohort median for each gene; no optimized cutpoint search",
    "Kaplan-Meier estimate with two-sided log-rank test; binary Cox HR also reported",
    "Benjamini-Hochberg across the five genes separately within each endpoint/model",
    "Continuous-expression Cox adjusted for age and pathologic stage I-II vs III-IV",
    clinical_path
  )
)
fwrite(analysis_parameters, file.path(qc_dir, "TN_five_gene_survival_analysis_parameters.tsv"), sep = "\t")

validation <- data.table(
  check = c(
    "five_genes_map_one_to_one",
    "exactly_169_unique_strict_TN_patients",
    "all_TN_samples_exist_in_TPM_matrix",
    "exactly_five_TPM_rows_extracted",
    "no_missing_or_negative_TPM_values",
    "patient_level_data_has_unique_patient_and_sample_ids",
    "clinical_source_has_unique_patient_ids",
    "168_of_169_strict_TN_patients_match_survival_source",
    "OS_has_at_least_20_events",
    "event_variables_are_binary_or_missing",
    "continuous_result_has_5_genes_x_4_endpoints",
    "median_split_result_has_5_genes_x_4_endpoints",
    "OS_significance_summary_has_exactly_five_genes",
    "OS_is_primary_endpoint",
    "median_cutpoints_are_fixed_before_outcome_analysis",
    "BH_FDR_is_computed_within_each_endpoint_across_five_genes",
    "all_univariable_Cox_HRs_and_P_values_are_finite"
  ),
  passed = c(
    nrow(gene_annotation) == 5L && !anyDuplicated(gene_annotation$gene_symbol),
    nrow(tn_meta) == 169L && uniqueN(tn_meta$patient_id) == 169L &&
      setequal(tn_meta$sample_id, strict_receptor$sample_id),
    all(tn_meta$sample_id %chin% names(tpm_subset)),
    nrow(tpm_subset) == 5L,
    !anyNA(survival_data[, paste0(genes, "_TPM"), with = FALSE]) && all(unlist(survival_data[, paste0(genes, "_TPM"), with = FALSE], use.names = FALSE) >= 0),
    !anyDuplicated(survival_data$patient_id) && !anyDuplicated(survival_data$sample_id),
    !anyDuplicated(clinical$patient_id),
    sum(survival_data$clinical_match) == 168L,
    endpoint_completeness[endpoint == "OS", events] >= 20L,
    all(unlist(survival_data[, paste0(endpoint_order, "_event"), with = FALSE], use.names = FALSE) %in% c(0L, 1L, NA_integer_)),
    nrow(continuous_results) == 20L && setequal(unique(continuous_results$gene), genes),
    nrow(median_results) == 20L && setequal(unique(median_results$gene), genes),
    nrow(os_summary) == 5L && identical(os_summary$gene, genes),
    primary_endpoint == "OS",
    all(grepl("TN-cohort median", median_results$cutoff_rule, fixed = TRUE)),
    all(continuous_results[, isTRUE(all.equal(Wald_FDR_BH_5genes, p.adjust(Wald_P, method = "BH"))), by = endpoint]$V1) &&
      all(median_results[, isTRUE(all.equal(logrank_FDR_BH_5genes, p.adjust(logrank_P, method = "BH"))), by = endpoint]$V1),
    all(is.finite(continuous_results$HR)) && all(is.finite(continuous_results$Wald_P))
  )
)
fwrite(validation, file.path(qc_dir, "TN_five_gene_survival_validation_checks.tsv"), sep = "\t")
if (!all(validation$passed)) {
  stop("Validation failed: ", paste(validation[passed == FALSE, check], collapse = ", "))
}

capture.output(sessionInfo(), file = file.path(qc_dir, "R_sessionInfo.txt"))

os_cont <- continuous_results[endpoint == "OS"]
os_med <- median_results[endpoint == "OS"]

message("Completed: ", bundle_dir)
print(endpoint_completeness)
print(os_summary)
print(validation)
