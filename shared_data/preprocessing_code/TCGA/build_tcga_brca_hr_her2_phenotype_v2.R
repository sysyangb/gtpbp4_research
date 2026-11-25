#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(data.table))

root <- normalizePath(dirname(sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)[1])), mustWork = TRUE)

patient_input <- file.path(root, "TCGA-BRCA_receptor_subtype_by_patient.tsv")
sample_input <- file.path(root, "TCGA-BRCA_receptor_subtype_by_sample.tsv")
counts_input <- file.path(root, "TCGA-BRCA_counts.txt")
tpm_input <- file.path(root, "TCGA-BRCA_TPM.txt")

patient_output <- file.path(root, "TCGA-BRCA_HR_HER2_phenotype_by_patient_v2.tsv")
sample_output <- file.path(root, "TCGA-BRCA_HR_HER2_phenotype_by_sample_v2.tsv")
condition_output <- file.path(root, "TCGA-BRCA_HR_HER2_NORMAL.condition_v2.tsv")
qc_output <- file.path(root, "TCGA-BRCA_HR_HER2_phenotype_v2_QC.txt")

required_inputs <- c(patient_input, sample_input, counts_input, tpm_input)
if (!all(file.exists(required_inputs))) {
  stop("Missing input files: ", paste(basename(required_inputs[!file.exists(required_inputs)]), collapse = ", "))
}

normalize_binary <- function(x) {
  value <- tolower(trimws(fcoalesce(as.character(x), "")))
  fifelse(value == "positive", "Positive", fifelse(value == "negative", "Negative", "Unknown"))
}

derive_hr <- function(er, pr) {
  fifelse(
    er == "Positive" | pr == "Positive", "Positive",
    fifelse(er == "Negative" & pr == "Negative", "Negative", "Unknown")
  )
}

derive_phenotype <- function(hr, her2) {
  fifelse(
    her2 == "Positive", "HER2",
    fifelse(
      her2 == "Negative" & hr == "Positive", "HR",
      fifelse(her2 == "Negative" & hr == "Negative", "TN", "UNKNOWN")
    )
  )
}

phenotype_display <- function(x) {
  map <- c(
    HR = "HR (HER2-; ER+ and/or PR+)",
    HER2 = "HER2+",
    TN = "TN (HR-/HER2-)",
    NORMAL = "Normal",
    UNKNOWN = "Unknown"
  )
  unname(map[x])
}

read_expression_header <- function(path) {
  values <- strsplit(readLines(path, n = 1L, warn = FALSE), "\t", fixed = TRUE)[[1]]
  gsub("^\"|\"$", "", values)
}

append_quality_flags <- function(original, subtype_source_conflict, her2_assay_conflict,
                                 selected_phenotype, recommended_phenotype) {
  vapply(seq_along(original), function(i) {
    old <- strsplit(fcoalesce(original[i], ""), "\\|", perl = TRUE)[[1]]
    old <- old[nzchar(old)]
    add <- character()
    if (isTRUE(subtype_source_conflict[i])) add <- c(add, "SUBTYPE_SOURCE_CONFLICT")
    if (isTRUE(her2_assay_conflict[i])) add <- c(add, "HER2_ASSAY_CONFLICT_SUBTYPE_UNRESOLVED")
    if (selected_phenotype[i] == "UNKNOWN") add <- c(add, "INCOMPLETE_FOR_HR_HER2_SUBTYPE")
    if (recommended_phenotype[i] == "UNKNOWN" && selected_phenotype[i] != "UNKNOWN" &&
        !isTRUE(subtype_source_conflict[i]) && !isTRUE(her2_assay_conflict[i])) {
      add <- c(add, "ORIGINAL_TN_NON_TN_CLASS_CONFLICT")
    }
    if (length(add)) old <- setdiff(old, "OK")
    flags <- unique(c(old, add))
    if (!length(flags)) flags <- "OK"
    paste(flags, collapse = "|")
  }, character(1))
}

patient <- fread(patient_input, na.strings = NULL)
sample <- fread(sample_input, colClasses = c(sample_type_code = "character"), na.strings = NULL)

required_patient_columns <- c(
  "patient_id", "ER_status_selected", "PR_status_selected", "HER2_status_selected",
  "label_recommended", "label_best_available", "label_source", "confidence", "quality_flag",
  "gdc_er_status_by_ihc", "gdc_pr_status_by_ihc", "gdc_her2_status_by_ihc",
  "gdc_her2_fish_status", "nature2012_er_status", "nature2012_pr_status",
  "nature2012_her2_final_status"
)
missing_columns <- setdiff(required_patient_columns, names(patient))
if (length(missing_columns)) stop("Missing patient columns: ", paste(missing_columns, collapse = ", "))

# Selected receptor state: Nature 2012 final status first, with GDC fill, as in the
# audited receptor-label input. HR positive means ER or PR positive.
patient[, `:=`(
  HR_status_selected = derive_hr(ER_status_selected, PR_status_selected),
  HER2_status_selected_v2 = HER2_status_selected
)]
patient[, phenotype_best_available := derive_phenotype(HR_status_selected, HER2_status_selected_v2)]

# Source-specific states are retained to detect conflicts that are invisible at
# the coarse TN/non-TN level but change HR/HER2 subtype assignment.
patient[, `:=`(
  gdc_ER_status_v2 = normalize_binary(gdc_er_status_by_ihc),
  gdc_PR_status_v2 = normalize_binary(gdc_pr_status_by_ihc),
  gdc_HER2_IHC_status_v2 = normalize_binary(gdc_her2_status_by_ihc),
  gdc_HER2_FISH_status_v2 = normalize_binary(gdc_her2_fish_status),
  nature2012_ER_status_v2 = normalize_binary(nature2012_er_status),
  nature2012_PR_status_v2 = normalize_binary(nature2012_pr_status),
  nature2012_HER2_status_v2 = normalize_binary(nature2012_her2_final_status)
)]

patient[, gdc_HER2_assay_conflict :=
          gdc_HER2_IHC_status_v2 %in% c("Positive", "Negative") &
          gdc_HER2_FISH_status_v2 %in% c("Positive", "Negative") &
          gdc_HER2_IHC_status_v2 != gdc_HER2_FISH_status_v2]
patient[, gdc_HER2_status_v2 := fifelse(
  gdc_HER2_IHC_status_v2 == "Positive" | gdc_HER2_FISH_status_v2 == "Positive", "Positive",
  fifelse(gdc_HER2_IHC_status_v2 == "Negative" | gdc_HER2_FISH_status_v2 == "Negative", "Negative", "Unknown")
)]
patient[, `:=`(
  gdc_HR_status_v2 = derive_hr(gdc_ER_status_v2, gdc_PR_status_v2),
  nature2012_HR_status_v2 = derive_hr(nature2012_ER_status_v2, nature2012_PR_status_v2)
)]
patient[, `:=`(
  gdc_phenotype_v2 = derive_phenotype(gdc_HR_status_v2, gdc_HER2_status_v2),
  nature2012_phenotype_v2 = derive_phenotype(nature2012_HR_status_v2, nature2012_HER2_status_v2)
)]
patient[, subtype_source_conflict :=
          gdc_phenotype_v2 != "UNKNOWN" & nature2012_phenotype_v2 != "UNKNOWN" &
          gdc_phenotype_v2 != nature2012_phenotype_v2]

# Strict recommended phenotype: do not force a fine subtype when the two sources
# disagree at subtype level, when GDC IHC/FISH conflict, or when a required
# receptor state is missing. The best-available column preserves the Nature-first
# call for sensitivity analyses.
patient[, phenotype_recommended := fifelse(
  label_recommended == "TN" & phenotype_best_available == "TN", "TN",
  fifelse(
    label_recommended == "NON_TN" & phenotype_best_available %in% c("HR", "HER2") &
      !subtype_source_conflict & !gdc_HER2_assay_conflict,
    phenotype_best_available,
    "UNKNOWN"
  )
)]
patient[, `:=`(
  phenotype_recommended_display = phenotype_display(phenotype_recommended),
  phenotype_best_available_display = phenotype_display(phenotype_best_available)
)]
patient[, quality_flag_v2 := append_quality_flags(
  quality_flag, subtype_source_conflict, gdc_HER2_assay_conflict,
  phenotype_best_available, phenotype_recommended
)]
patient[, `:=`(
  phenotype_rule_based = phenotype_best_available,
  phenotype_rule_based_display = phenotype_best_available_display,
  phenotype_strict_qc = phenotype_recommended,
  phenotype_strict_qc_display = phenotype_recommended_display
)]

patient_fields <- c(
  "patient_id", "phenotype_rule_based", "phenotype_rule_based_display",
  "phenotype_strict_qc", "phenotype_strict_qc_display",
  "HR_status_selected", "HER2_status_selected_v2",
  "ER_status_selected", "PR_status_selected", "HER2_status_selected",
  "subtype_source_conflict", "gdc_HER2_assay_conflict",
  "gdc_HR_status_v2", "gdc_HER2_status_v2", "gdc_phenotype_v2",
  "nature2012_HR_status_v2", "nature2012_HER2_status_v2", "nature2012_phenotype_v2",
  "label_source", "confidence", "quality_flag_v2", "quality_flag"
)
patient_out <- patient[, ..patient_fields]

# Sample-level phenotype follows the tumor patient call, while sample type 11 is
# always NORMAL. This file remains in the exact counts/TPM column order.
patient_index <- match(sample$patient_id, patient$patient_id)
if (anyNA(patient_index)) stop("Sample table contains patients absent from patient table")
sample_out <- cbind(
  sample[, .(sample_id, patient_id, sample_type_code, sample_type)],
  patient_out[patient_index, -"patient_id"]
)
normal <- sample_out$sample_type_code == "11"
sample_out[normal, `:=`(
  phenotype_rule_based = "NORMAL",
  phenotype_rule_based_display = "Normal",
  phenotype_strict_qc = "NORMAL",
  phenotype_strict_qc_display = "Normal",
  HR_status_selected = "Not_Applicable",
  HER2_status_selected_v2 = "Not_Applicable",
  quality_flag_v2 = "NORMAL_NOT_RECEPTOR_CLASSIFIED"
)]

condition_out <- sample_out[, .(
  sample_id, patient_id, sample_type_code, sample_type,
  condition = phenotype_rule_based,
  condition_display = phenotype_rule_based_display,
  condition_strict_qc = phenotype_strict_qc,
  condition_strict_qc_display = phenotype_strict_qc_display,
  HR_status = HR_status_selected,
  HER2_status = HER2_status_selected_v2,
  quality_flag = quality_flag_v2
)]

counts_samples <- read_expression_header(counts_input)
tpm_samples <- read_expression_header(tpm_input)
allowed <- c("HR", "HER2", "TN", "NORMAL", "UNKNOWN")

validation <- c(
  counts_tpm_headers_identical = identical(counts_samples, tpm_samples),
  sample_order_matches_expression = identical(condition_out$sample_id, counts_samples),
  sample_ids_unique = !anyDuplicated(condition_out$sample_id),
  patient_ids_unique = !anyDuplicated(patient_out$patient_id),
  phenotype_values_valid = all(condition_out$condition %in% allowed),
  normal_only_sample_type_11 = all(condition_out$condition[normal] == "NORMAL") &&
    all(condition_out$condition[!normal] != "NORMAL"),
  tn_receptor_rule_valid = all(patient[phenotype_recommended == "TN",
    ER_status_selected == "Negative" & PR_status_selected == "Negative" & HER2_status_selected == "Negative"]),
  hr_rule_valid = all(patient[phenotype_recommended == "HR",
    HR_status_selected == "Positive" & HER2_status_selected_v2 == "Negative"]),
  her2_rule_valid = all(patient[phenotype_recommended == "HER2",
    HER2_status_selected_v2 == "Positive"]),
  strict_groups_exclude_subtype_conflicts = all(patient[subtype_source_conflict == TRUE, phenotype_recommended == "UNKNOWN"]),
  strict_groups_exclude_her2_assay_conflicts = all(patient[gdc_HER2_assay_conflict == TRUE, phenotype_recommended == "UNKNOWN"])
)
if (!all(validation)) stop("Validation failed: ", paste(names(validation)[!validation], collapse = ", "))

fwrite(patient_out, patient_output, sep = "\t", quote = FALSE, na = "")
fwrite(sample_out, sample_output, sep = "\t", quote = FALSE, na = "")
fwrite(condition_out, condition_output, sep = "\t", quote = FALSE, na = "")

format_counts <- function(dt, column) {
  tab <- dt[, .N, by = column][order(-N)]
  paste(sprintf("  %s: %d", tab[[column]], tab$N), collapse = "\n")
}

qc_lines <- c(
  "TCGA-BRCA HR/HER2 phenotype v2 QC",
  "",
  "Intended grain:",
  sprintf("  Patient phenotype rows: %d", nrow(patient_out)),
  sprintf("  Expression sample phenotype rows: %d", nrow(sample_out)),
  "",
  "Rule-based patient phenotype counts (main condition):",
  format_counts(patient_out, "phenotype_rule_based"),
  "",
  "Strict-QC patient phenotype counts (sensitivity condition):",
  format_counts(patient_out, "phenotype_strict_qc"),
  "",
  "Rule-based expression-sample condition counts (main condition):",
  format_counts(condition_out, "condition"),
  "",
  "Strict-QC expression-sample condition counts (sensitivity condition):",
  format_counts(condition_out, "condition_strict_qc"),
  "",
  sprintf("Subtype-changing GDC vs Nature 2012 conflicts: %d patients", sum(patient$subtype_source_conflict)),
  sprintf("GDC HER2 IHC/FISH conflicts: %d patients", sum(patient$gdc_HER2_assay_conflict)),
  sprintf("Original NON_TN patients incomplete for exact HR/HER2 subtype: %d",
          sum(patient$label_recommended == "NON_TN" & patient$phenotype_best_available == "UNKNOWN")),
  "",
  "Validation:",
  sprintf("  %s: %s", names(validation), ifelse(validation, "YES", "NO")),
  "",
  "Rules:",
  "  HR positive = ER positive OR PR positive.",
  "  HR negative = ER negative AND PR negative.",
  "  TN = HER2 negative AND ER negative AND PR negative.",
  "  HR = HER2 negative AND (ER positive OR PR positive).",
  "  HER2 = HER2 positive, regardless of ER/PR status.",
  "  Classification priority is HER2, then HR/TN among HER2-negative tumors.",
  "  Recommended subtypes require all receptor states needed by the rule to be conclusive.",
  "  Recommended fine subtypes exclude subtype-changing source conflicts and GDC HER2 IHC/FISH conflicts.",
  "  Best-available subtype retains the audited Nature-2012-first, GDC-fill receptor call for sensitivity analysis.",
  "  Sample type code 11 is NORMAL; NORMAL is not assigned at patient tumor-phenotype level.",
  "",
  "Input MD5:",
  sprintf("  %s: %s", basename(c(patient_input, sample_input, counts_input, tpm_input)),
          unname(tools::md5sum(c(patient_input, sample_input, counts_input, tpm_input))))
)
writeLines(qc_lines, qc_output, useBytes = TRUE)

cat("Wrote:", basename(patient_output), "-", nrow(patient_out), "rows\n")
cat("Wrote:", basename(sample_output), "-", nrow(sample_out), "rows\n")
cat("Wrote:", basename(condition_output), "-", nrow(condition_out), "rows\n")
cat("Rule-based patient phenotypes (main):\n")
print(patient_out[, .N, by = phenotype_rule_based][order(-N)])
cat("Rule-based expression-sample conditions (main):\n")
print(condition_out[, .N, by = condition][order(-N)])
cat("Strict-QC expression-sample conditions (sensitivity):\n")
print(condition_out[, .N, by = condition_strict_qc][order(-N)])
