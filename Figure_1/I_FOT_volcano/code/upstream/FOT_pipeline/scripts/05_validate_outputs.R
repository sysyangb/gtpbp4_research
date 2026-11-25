message("[05] Validating generated tables, calls, and figure exports")

prepared <- readRDS(file.path(cfg$data_dir, "prepared_FOT_data.rds"))
primary <- readr::read_csv(
  file.path(cfg$tables_dir, "FOT_limma_all_tested.csv"),
  show_col_types = FALSE
)
strict_up <- readr::read_csv(
  file.path(cfg$tables_dir, "FOT_TREAT_BH05_absFCgt1_UP.csv"),
  show_col_types = FALSE
)
strict_down <- readr::read_csv(
  file.path(cfg$tables_dir, "FOT_TREAT_BH05_absFCgt1_DOWN.csv"),
  show_col_types = FALSE
)
consensus_up <- readr::read_csv(
  file.path(cfg$tables_dir, "recommended_high_confidence_FOT_UP_consensus.csv"),
  show_col_types = FALSE
)
consensus_down <- readr::read_csv(
  file.path(cfg$tables_dir, "recommended_high_confidence_FOT_DOWN_consensus.csv"),
  show_col_types = FALSE
)
detection_up <- readr::read_csv(
  file.path(cfg$tables_dir, "detection_only_Polysome_2of2_Free_0of3.csv"),
  show_col_types = FALSE
)
proda <- readr::read_csv(
  file.path(cfg$tables_dir, "proDA_all_UP5_atleast2of5.csv"),
  show_col_types = FALSE
)

png_files <- list.files(cfg$figures_dir, pattern = "[.]png$", full.names = TRUE)
pdf_files <- list.files(cfg$figures_dir, pattern = "[.]pdf$", full.names = TRUE)
png_stems <- tools::file_path_sans_ext(basename(png_files))
pdf_stems <- tools::file_path_sans_ext(basename(pdf_files))

checks <- list(
  source_accessions_unique = !anyDuplicated(prepared$metadata$Protein_accession),
  fot_columns_sum_to_one = max(abs(colSums(prepared$fot, na.rm = TRUE) - 1)) < 1e-10,
  ifot_columns_sum_to_scale = max(abs(colSums(prepared$ifot, na.rm = TRUE) - cfg$ifot_scale)) < 1e-5,
  primary_accessions_unique = !anyDuplicated(primary$Protein_accession),
  primary_unique_peptides_filter = all(primary$Unique_peptides >= cfg$unique_peptides_min),
  primary_free_observation_filter = all(primary$n_observed_Free >= cfg$min_observed_per_group),
  primary_polysome_observation_filter = all(primary$n_observed_Polysome >= cfg$min_observed_per_group),
  primary_coefficients_finite = all(is.finite(primary$log2FC_Polysome_vs_Free)),
  primary_pvalues_valid = all(primary$P_Value >= 0 & primary$P_Value <= 1),
  primary_treat_pvalues_valid = all(primary$TREAT_P_Value >= 0 & primary$TREAT_P_Value <= 1),
  strict_up_call_consistent = all(strict_up$TREAT_call == "Up" & strict_up$log2FC_Polysome_vs_Free > cfg$log2fc_threshold),
  strict_down_call_consistent = all(strict_down$TREAT_call == "Down" & strict_down$log2FC_Polysome_vs_Free < -cfg$log2fc_threshold),
  strict_sets_partition_primary_calls =
    nrow(strict_up) == sum(primary$TREAT_call == "Up") &&
    nrow(strict_down) == sum(primary$TREAT_call == "Down"),
  consensus_up_is_subset_of_strict_up = all(consensus_up$Protein_accession %in% strict_up$Protein_accession),
  consensus_down_is_subset_of_strict_down = all(consensus_down$Protein_accession %in% strict_down$Protein_accession),
  consensus_up_call_consistent = all(consensus_up$consensus_call == "Up"),
  consensus_down_call_consistent = all(consensus_down$consensus_call == "Down"),
  detection_only_pattern_consistent = all(
    detection_up$Unique_peptides >= cfg$unique_peptides_min &
      detection_up$n_observed_Free == 0 &
      detection_up$n_observed_Polysome == 2
  ),
  proda_pvalues_valid = all(proda$proDA_P_Value >= 0 & proda$proDA_P_Value <= 1),
  proda_fdr_valid = all(proda$proDA_BH_FDR >= 0 & proda$proDA_BH_FDR <= 1),
  figure_png_pdf_stems_match = setequal(png_stems, pdf_stems),
  figure_files_nonempty = all(file.info(c(png_files, pdf_files))$size > 1000),
  required_docs_exist = all(file.exists(c(
    file.path(cfg$pipeline_root, "README.md"),
    file.path(cfg$pipeline_root, "docs", "workflow_exploration.md"),
    file.path(cfg$pipeline_root, "docs", "results_summary.md"),
    file.path(cfg$pipeline_root, "docs", "literature_search_record.md")
  )))
)

failed <- names(checks)[!unlist(checks, use.names = FALSE)]
qa <- list(
  status = if (length(failed)) "failed" else "passed",
  checked_at = format(Sys.time(), tz = "Asia/Shanghai", usetz = TRUE),
  checks = checks,
  failed_checks = failed,
  counts = list(
    source_proteins = nrow(prepared$metadata),
    primary_tested = nrow(primary),
    strict_up = nrow(strict_up),
    strict_down = nrow(strict_down),
    consensus_up = nrow(consensus_up),
    consensus_down = nrow(consensus_down),
    polysome_detection_only = nrow(detection_up),
    proDA_tested = nrow(proda),
    figure_pairs = length(png_files)
  )
)
jsonlite::write_json(
  qa,
  file.path(cfg$logs_dir, "final_QA.json"),
  pretty = TRUE,
  auto_unbox = TRUE,
  digits = NA
)
if (length(failed)) stop("Final QA failed: ", paste(failed, collapse = ", "))

message("[05] All ", length(checks), " final QA checks passed")
