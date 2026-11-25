message("[04] Fitting proDA probabilistic dropout model (no single-value imputation)")

dat <- readRDS(file.path(cfg$data_dir, "prepared_FOT_data.rds"))
diff_obj <- readRDS(file.path(cfg$data_dir, "FOT_differential_objects.rds"))
meta <- dat$metadata
y_all <- dat$log2_ifot[, cfg$main_samples, drop = FALSE]

# proDA explicitly models intensity-dependent dropout. Unlike the limma primary
# subset, it can retain proteins with uneven group-wise observation patterns.
# We still require identification reliability and at least two observations in
# the five retained samples so the model has empirical abundance information.
observed_total <- rowSums(!is.na(y_all))
proda_keep <-
  !meta$explicit_decoy_or_contaminant &
  meta$Unique_peptides >= cfg$unique_peptides_min &
  observed_total >= 2

y <- y_all[proda_keep, , drop = FALSE]
proda_meta <- meta[proda_keep, , drop = FALSE]
proda_meta$n_observed_main5 <- observed_total[proda_keep]
proda_meta$n_observed_Free <- rowSums(!is.na(y[, c("free_1", "free_2", "free_3"), drop = FALSE]))
proda_meta$n_observed_Polysome <- rowSums(!is.na(y[, c("polysome_1", "polysome_2"), drop = FALSE]))

col_data <- data.frame(
  sample_id = cfg$main_samples,
  group = factor(unname(cfg$groups[cfg$main_samples]), levels = c("Free", "Polysome")),
  row.names = cfg$main_samples,
  stringsAsFactors = FALSE
)

proda_fit <- proDA::proDA(
  y,
  design = ~ group,
  col_data = col_data,
  reference_level = "Free",
  data_is_log_transformed = TRUE,
  moderate_location = TRUE,
  moderate_variance = TRUE,
  verbose = FALSE
)
proda_test <- proDA::test_diff(
  proda_fit,
  contrast = "groupPolysome",
  alternative = "two.sided",
  pval_adjust_method = "BH",
  sort_by = NULL
)

idx <- match(proda_test$name, rownames(y))
if (anyNA(idx)) stop("Could not match all proDA result rows back to protein metadata.")
proda_res <- cbind(
  proda_meta[idx, , drop = FALSE],
  data.frame(
    proDA_log2FC_Polysome_vs_Free = proda_test$diff,
    proDA_P_Value = proda_test$pval,
    proDA_BH_FDR = proda_test$adj_pval,
    proDA_t = proda_test$t_statistic,
    proDA_SE = proda_test$se,
    proDA_df = proda_test$df,
    proDA_average_abundance = proda_test$avg_abundance,
    proDA_effective_n = proda_test$n_approx,
    proDA_observed_n = proda_test$n_obs,
    stringsAsFactors = FALSE
  )
)
proda_res$proDA_call <- dplyr::case_when(
  proda_res$proDA_BH_FDR < cfg$fdr_threshold &
    proda_res$proDA_log2FC_Polysome_vs_Free > cfg$log2fc_threshold ~ "Up",
  proda_res$proDA_BH_FDR < cfg$fdr_threshold &
    proda_res$proDA_log2FC_Polysome_vs_Free < -cfg$log2fc_threshold ~ "Down",
  TRUE ~ "Not_significant"
)
proda_res <- proda_res[order(proda_res$proDA_P_Value, -abs(proda_res$proDA_log2FC_Polysome_vs_Free)), ]

write_csv_stable(proda_res, file.path(cfg$tables_dir, "proDA_all_UP5_atleast2of5.csv"))
write_csv_stable(
  proda_res[proda_res$proDA_call == "Up", ],
  file.path(cfg$tables_dir, "proDA_BH05_absFCgt1_UP.csv")
)
write_csv_stable(
  proda_res[proda_res$proDA_call == "Down", ],
  file.path(cfg$tables_dir, "proDA_BH05_absFCgt1_DOWN.csv")
)

primary <- diff_obj$primary_results
common <- dplyr::inner_join(
  primary[, c(
    "Protein_accession", "Gene_name", "Protein_description", "Unique_peptides",
    "n_observed_Free", "n_observed_Polysome", "log2FC_Polysome_vs_Free",
    "TREAT_BH_FDR", "TREAT_call", "BH_FDR", "BH_plus_FC_call"
  )],
  proda_res[, c(
    "Protein_accession", "proDA_log2FC_Polysome_vs_Free",
    "proDA_P_Value", "proDA_BH_FDR", "proDA_call"
  )],
  by = "Protein_accession"
)
common$direction_agree <-
  sign(common$log2FC_Polysome_vs_Free) == sign(common$proDA_log2FC_Polysome_vs_Free)
common$consensus_call <- dplyr::case_when(
  common$TREAT_call == "Up" & common$proDA_call == "Up" & common$direction_agree ~ "Up",
  common$TREAT_call == "Down" & common$proDA_call == "Down" & common$direction_agree ~ "Down",
  TRUE ~ "Not_consensus"
)
write_csv_stable(common, file.path(cfg$tables_dir, "limma_proDA_common_protein_comparison.csv"))
write_csv_stable(
  common[common$consensus_call == "Up", ],
  file.path(cfg$tables_dir, "recommended_high_confidence_FOT_UP_consensus.csv")
)
write_csv_stable(
  common[common$consensus_call == "Down", ],
  file.path(cfg$tables_dir, "recommended_high_confidence_FOT_DOWN_consensus.csv")
)

proda_summary <- data.frame(
  metric = c(
    "proDA_input_proteins",
    "proDA_up_BH05_absFCgt1",
    "proDA_down_BH05_absFCgt1",
    "limma_proDA_common_proteins",
    "common_log2FC_Pearson",
    "common_log2FC_Spearman",
    "common_direction_agreement",
    "consensus_up_limmaTREAT_and_proDA",
    "consensus_down_limmaTREAT_and_proDA"
  ),
  value = c(
    nrow(proda_res),
    sum(proda_res$proDA_call == "Up"),
    sum(proda_res$proDA_call == "Down"),
    nrow(common),
    stats::cor(common$log2FC_Polysome_vs_Free, common$proDA_log2FC_Polysome_vs_Free, method = "pearson"),
    stats::cor(common$log2FC_Polysome_vs_Free, common$proDA_log2FC_Polysome_vs_Free, method = "spearman"),
    mean(common$direction_agree),
    sum(common$consensus_call == "Up"),
    sum(common$consensus_call == "Down")
  ),
  stringsAsFactors = FALSE
)
write_csv_stable(proda_summary, file.path(cfg$tables_dir, "proDA_summary.csv"))

# Missingness diagnostics: observed abundance drops as the number of missing
# measurements increases, supporting an intensity-dependent dropout component.
observed_mean_log2_intensity <- rowMeans(log2(dat$intensity[, cfg$main_samples, drop = FALSE]), na.rm = TRUE)
missing_diag <- data.frame(
  Protein_accession = meta$Protein_accession,
  Unique_peptides = meta$Unique_peptides,
  n_missing_main5 = rowSums(is.na(dat$intensity[, cfg$main_samples, drop = FALSE])),
  mean_observed_log2_intensity = observed_mean_log2_intensity,
  stringsAsFactors = FALSE
)
missing_diag$mean_observed_log2_intensity[!is.finite(missing_diag$mean_observed_log2_intensity)] <- NA_real_
missing_summary <- missing_diag |>
  dplyr::group_by(n_missing_main5) |>
  dplyr::summarise(
    proteins = dplyr::n(),
    median_observed_log2_intensity = stats::median(mean_observed_log2_intensity, na.rm = TRUE),
    mean_observed_log2_intensity = mean(mean_observed_log2_intensity, na.rm = TRUE),
    .groups = "drop"
  )
write_csv_stable(missing_diag, file.path(cfg$tables_dir, "protein_missingness_diagnostics.csv"))
write_csv_stable(missing_summary, file.path(cfg$tables_dir, "missingness_by_number_missing.csv"))

p_missing <- ggplot2::ggplot(
  missing_diag[!is.na(missing_diag$mean_observed_log2_intensity), ],
  ggplot2::aes(factor(n_missing_main5), mean_observed_log2_intensity)
) +
  ggplot2::geom_violin(fill = "#9ECAE1", color = "#3182BD", scale = "width", trim = TRUE) +
  ggplot2::geom_boxplot(width = 0.18, outlier.shape = NA, fill = "white") +
  ggplot2::labs(
    x = "Number of missing values among retained five samples",
    y = "Mean observed log2 normalized intensity",
    title = "Missingness is associated with lower abundance"
  ) +
  theme_fot()
save_plot_both(p_missing, "QC_missingness_vs_abundance", 7.2, 5.4)

p_proda <- ggplot2::ggplot(
  common,
  ggplot2::aes(log2FC_Polysome_vs_Free, proDA_log2FC_Polysome_vs_Free, color = consensus_call)
) +
  ggplot2::geom_point(alpha = 0.58, size = 1.15) +
  ggplot2::geom_abline(slope = 1, intercept = 0, linetype = 2, color = "grey40") +
  ggplot2::scale_color_manual(values = c(Up = "#E45756", Down = "#4C78A8", Not_consensus = "#BDBDBD")) +
  ggplot2::labs(
    x = "limma log2FC on FOT",
    y = "proDA log2FC on FOT",
    color = "Consensus",
    title = "Missingness-model sensitivity"
  ) +
  theme_fot()
save_plot_both(p_proda, "comparison_limma_vs_proDA_log2FC", 6.6, 5.7)

saveRDS(
  list(
    fit = proda_fit,
    all_results = proda_res,
    common_comparison = common,
    summary = proda_summary,
    missingness_summary = missing_summary
  ),
  file.path(cfg$data_dir, "proDA_missingness_model.rds"),
  version = 3
)

message(
  "[04] proDA completed: ", sum(proda_res$proDA_call == "Up"), " up, ",
  sum(proda_res$proDA_call == "Down"), " down; consensus with limma TREAT: ",
  sum(common$consensus_call == "Up"), " up, ", sum(common$consensus_call == "Down"), " down"
)
