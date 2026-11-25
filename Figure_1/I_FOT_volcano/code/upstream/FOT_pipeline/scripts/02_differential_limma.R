message("[02] Running primary FOT differential analysis with limma")

prepared_path <- file.path(cfg$data_dir, "prepared_FOT_data.rds")
stopifnot(file.exists(prepared_path))
dat <- readRDS(prepared_path)

metadata <- dat$metadata
log2_ifot <- dat$log2_ifot
ifot <- dat$ifot

free_samples <- cfg$main_samples[cfg$groups[cfg$main_samples] == "Free"]
poly_samples <- cfg$main_samples[cfg$groups[cfg$main_samples] == "Polysome"]

metadata$n_observed_Free <- count_observed(log2_ifot, free_samples)
metadata$n_observed_Polysome <- count_observed(log2_ifot, poly_samples)
metadata$n_observed_main5 <- count_observed(log2_ifot, cfg$main_samples)
metadata$complete_main5 <- metadata$n_observed_main5 == length(cfg$main_samples)

eligible <-
  !metadata$explicit_decoy_or_contaminant &
  metadata$Unique_peptides >= cfg$unique_peptides_min &
  metadata$n_observed_Free >= cfg$min_observed_per_group &
  metadata$n_observed_Polysome >= cfg$min_observed_per_group

analysis_meta <- metadata[eligible, , drop = FALSE]
analysis_matrix <- log2_ifot[eligible, cfg$main_samples, drop = FALSE]
if (!nrow(analysis_matrix)) stop("No proteins passed the primary analysis filter.")

fits <- fit_limma_two_group(
  log2_matrix = analysis_matrix,
  sample_ids = cfg$main_samples,
  group_labels = cfg$groups,
  lfc_threshold = cfg$log2fc_threshold
)
res <- extract_limma_results(fits, analysis_meta)
res$mean_log2_iFOT_Free <- rowMeans(analysis_matrix[, free_samples, drop = FALSE], na.rm = TRUE)
res$mean_log2_iFOT_Polysome <- rowMeans(analysis_matrix[, poly_samples, drop = FALSE], na.rm = TRUE)
res$mean_iFOT_Free <- rowMeans(ifot[eligible, free_samples, drop = FALSE], na.rm = TRUE)
res$mean_iFOT_Polysome <- rowMeans(ifot[eligible, poly_samples, drop = FALSE], na.rm = TRUE)
res <- add_result_calls(res, cfg)
res <- res[order(res$TREAT_P_Value, res$P_Value, -abs(res$log2FC_Polysome_vs_Free)), ]

write_csv_stable(res, file.path(cfg$tables_dir, "FOT_limma_all_tested.csv"))

write_result_subsets <- function(call_col, prefix) {
  up <- res[res[[call_col]] == "Up", , drop = FALSE]
  down <- res[res[[call_col]] == "Down", , drop = FALSE]
  write_csv_stable(up, file.path(cfg$tables_dir, paste0(prefix, "_UP.csv")))
  write_csv_stable(down, file.path(cfg$tables_dir, paste0(prefix, "_DOWN.csv")))
  c(Up = nrow(up), Down = nrow(down), Not_significant = nrow(res) - nrow(up) - nrow(down))
}

count_treat <- write_result_subsets("TREAT_call", "FOT_TREAT_BH05_absFCgt1")
count_bhfc <- write_result_subsets("BH_plus_FC_call", "FOT_eBayes_BH05_absFCgt1")
count_legacy <- write_result_subsets("legacy_rawP_plus_FC_call", "FOT_legacy_rawP05_absFCgt1")

call_counts <- dplyr::bind_rows(
  data.frame(method = "TREAT_lfc1_BH05", direction = names(count_treat), n = as.integer(count_treat)),
  data.frame(method = "eBayes_BH05_and_absFCgt1", direction = names(count_bhfc), n = as.integer(count_bhfc)),
  data.frame(method = "legacy_rawP05_and_absFCgt1", direction = names(count_legacy), n = as.integer(count_legacy))
)
write_csv_stable(call_counts, file.path(cfg$tables_dir, "FOT_differential_call_counts.csv"))

# Complete-case sensitivity: all five retained samples must be observed.
complete_eligible <- eligible & metadata$complete_main5
cc_meta <- metadata[complete_eligible, , drop = FALSE]
cc_matrix <- log2_ifot[complete_eligible, cfg$main_samples, drop = FALSE]
cc_fits <- fit_limma_two_group(
  log2_matrix = cc_matrix,
  sample_ids = cfg$main_samples,
  group_labels = cfg$groups,
  lfc_threshold = cfg$log2fc_threshold
)
cc_res <- extract_limma_results(cc_fits, cc_meta)
cc_res <- add_result_calls(cc_res, cfg)
cc_res <- cc_res[order(cc_res$TREAT_P_Value, cc_res$P_Value), ]
write_csv_stable(cc_res, file.path(cfg$tables_dir, "FOT_limma_complete_case_sensitivity.csv"))

# Proteins that cannot enter a continuous abundance model are reported separately.
presence_base <-
  !metadata$explicit_decoy_or_contaminant &
  metadata$Unique_peptides >= cfg$unique_peptides_min
poly_detection_only <- metadata[
  presence_base & metadata$n_observed_Free == 0 & metadata$n_observed_Polysome == length(poly_samples),
  , drop = FALSE
]
free_detection_only_strict <- metadata[
  presence_base & metadata$n_observed_Free == length(free_samples) & metadata$n_observed_Polysome == 0,
  , drop = FALSE
]
free_detection_only_lenient <- metadata[
  presence_base & metadata$n_observed_Free >= cfg$min_observed_per_group & metadata$n_observed_Polysome == 0,
  , drop = FALSE
]
write_csv_stable(poly_detection_only, file.path(cfg$tables_dir, "detection_only_Polysome_2of2_Free_0of3.csv"))
write_csv_stable(free_detection_only_strict, file.path(cfg$tables_dir, "detection_only_Free_3of3_Polysome_0of2.csv"))
write_csv_stable(free_detection_only_lenient, file.path(cfg$tables_dir, "detection_only_Free_atleast2of3_Polysome_0of2.csv"))

filter_summary <- data.frame(
  metric = c(
    "source_proteins",
    "unique_peptides_ge_5",
    "free_ge_2_and_polysome_2",
    "primary_tested_both_filters",
    "complete_case_primary_filter",
    "polysome_detection_only_2of2_vs_0of3",
    "free_detection_only_3of3_vs_0of2",
    "free_detection_only_atleast2of3_vs_0of2"
  ),
  n = c(
    nrow(metadata),
    sum(metadata$Unique_peptides >= cfg$unique_peptides_min, na.rm = TRUE),
    sum(metadata$n_observed_Free >= 2 & metadata$n_observed_Polysome >= 2),
    sum(eligible),
    sum(complete_eligible),
    nrow(poly_detection_only),
    nrow(free_detection_only_strict),
    nrow(free_detection_only_lenient)
  )
)
write_csv_stable(filter_summary, file.path(cfg$tables_dir, "analysis_filter_summary.csv"))

# Volcano and MA figures ----------------------------------------------------
plot_res <- res
plot_res$minus_log10_TREAT_P <- -log10(pmax(plot_res$TREAT_P_Value, .Machine$double.xmin))
plot_res$TREAT_call <- factor(plot_res$TREAT_call, levels = c("Not_significant", "Down", "Up"))
volcano_colors <- c(Not_significant = "#BDBDBD", Down = "#4C78A8", Up = "#E45756")
p_volcano <- ggplot2::ggplot(
  plot_res,
  ggplot2::aes(log2FC_Polysome_vs_Free, minus_log10_TREAT_P, color = TREAT_call)
) +
  ggplot2::geom_point(size = 1.25, alpha = 0.78) +
  ggplot2::geom_vline(xintercept = c(-cfg$log2fc_threshold, cfg$log2fc_threshold), linetype = 2, color = "grey45") +
  ggplot2::scale_color_manual(values = volcano_colors, drop = FALSE) +
  ggplot2::labs(
    x = "log2 fold change (Polysome / Free)",
    y = "-log10 TREAT P value",
    color = "TREAT + BH call",
    title = "FOT differential abundance",
    caption = "Colors denote TREAT BH-FDR < 0.05; vertical lines mark log2FC = +/-1"
  ) +
  theme_fot()
save_plot_both(p_volcano, "FOT_volcano_TREAT", 7.2, 5.6)

p_ma <- ggplot2::ggplot(
  plot_res,
  ggplot2::aes(AveExpr, log2FC_Polysome_vs_Free, color = TREAT_call)
) +
  ggplot2::geom_point(size = 1.2, alpha = 0.72) +
  ggplot2::geom_hline(yintercept = c(-cfg$log2fc_threshold, 0, cfg$log2fc_threshold), linetype = c(2, 1, 2), color = "grey45") +
  ggplot2::scale_color_manual(values = volcano_colors, drop = FALSE) +
  ggplot2::labs(
    x = "Average log2 iFOT",
    y = "log2 fold change (Polysome / Free)",
    color = "TREAT + BH call",
    title = "FOT mean-difference plot"
  ) +
  theme_fot()
save_plot_both(p_ma, "FOT_MA_TREAT", 7.2, 5.6)

saveRDS(
  list(
    eligible = eligible,
    complete_eligible = complete_eligible,
    primary_results = res,
    complete_case_results = cc_res,
    primary_fits = fits,
    complete_case_fits = cc_fits,
    call_counts = call_counts,
    filter_summary = filter_summary
  ),
  file.path(cfg$data_dir, "FOT_differential_objects.rds"),
  version = 3
)

validation <- list(
  tested_proteins = nrow(res),
  complete_case_tested_proteins = nrow(cc_res),
  all_tested_have_unique_peptides_ge_5 = all(res$Unique_peptides >= cfg$unique_peptides_min),
  all_tested_have_free_ge_2 = all(res$n_observed_Free >= cfg$min_observed_per_group),
  all_tested_have_polysome_ge_2 = all(res$n_observed_Polysome >= cfg$min_observed_per_group),
  finite_primary_coefficients = all(is.finite(res$log2FC_Polysome_vs_Free)),
  p_values_in_range = all(res$P_Value >= 0 & res$P_Value <= 1),
  treat_p_values_in_range = all(res$TREAT_P_Value >= 0 & res$TREAT_P_Value <= 1)
)
jsonlite::write_json(
  validation,
  file.path(cfg$logs_dir, "differential_validation.json"),
  pretty = TRUE,
  auto_unbox = TRUE
)

message(
  "[02] Tested ", nrow(res), " proteins; TREAT calls: ",
  count_treat[["Up"]], " up, ", count_treat[["Down"]], " down"
)
