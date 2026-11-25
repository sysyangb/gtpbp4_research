message("[03] Running normalization and model sensitivity analyses")

dat <- readRDS(file.path(cfg$data_dir, "prepared_FOT_data.rds"))
diff_obj <- readRDS(file.path(cfg$data_dir, "FOT_differential_objects.rds"))
eligible <- diff_obj$eligible
metadata <- dat$metadata[eligible, , drop = FALSE]
metadata$n_observed_Free <- count_observed(dat$intensity[eligible, , drop = FALSE], c("free_1", "free_2", "free_3"))
metadata$n_observed_Polysome <- count_observed(dat$intensity[eligible, , drop = FALSE], c("polysome_1", "polysome_2"))

all_log2_intensity <- log2(dat$intensity[, cfg$main_samples, drop = FALSE])
log2_intensity <- all_log2_intensity[eligible, , drop = FALSE]
# Estimate median-scaling factors from all quantified proteins before applying
# the Unique-peptide and testability filters, mirroring the FOT order of operations.
sample_medians <- apply(all_log2_intensity, 2, stats::median, na.rm = TRUE)
grand_median <- stats::median(sample_medians)
log2_intensity_median_centered <- sweep(
  log2_intensity,
  2,
  sample_medians - grand_median,
  FUN = "-"
)

fit_and_extract <- function(matrix, label) {
  fits <- fit_limma_two_group(
    log2_matrix = matrix,
    sample_ids = cfg$main_samples,
    group_labels = cfg$groups,
    lfc_threshold = cfg$log2fc_threshold
  )
  out <- extract_limma_results(fits, metadata)
  out <- add_result_calls(out, cfg)
  out$analysis_method <- label
  out
}

raw_res <- fit_and_extract(log2_intensity, "vendor_normalized_intensity_no_additional_scaling")
med_res <- fit_and_extract(log2_intensity_median_centered, "vendor_normalized_intensity_median_centered")
write_csv_stable(
  raw_res[order(raw_res$TREAT_P_Value, raw_res$P_Value), ],
  file.path(cfg$tables_dir, "sensitivity_limma_log2_normalized_intensity.csv")
)
write_csv_stable(
  med_res[order(med_res$TREAT_P_Value, med_res$P_Value), ],
  file.path(cfg$tables_dir, "sensitivity_limma_log2_intensity_median_centered.csv")
)

# Optional pairing sensitivity using only replicate labels 1 and 2. This is not
# the primary analysis because biological pairing was not documented in metadata.
paired_res <- fit_paired_two_pair_sensitivity(
  dat$log2_ifot[eligible, , drop = FALSE],
  metadata
)
paired_res$BH_plus_FC_call <- dplyr::case_when(
  paired_res$BH_FDR < cfg$fdr_threshold & paired_res$log2FC_Polysome_vs_Free > cfg$log2fc_threshold ~ "Up",
  paired_res$BH_FDR < cfg$fdr_threshold & paired_res$log2FC_Polysome_vs_Free < -cfg$log2fc_threshold ~ "Down",
  TRUE ~ "Not_significant"
)
write_csv_stable(
  paired_res[order(paired_res$P_Value), ],
  file.path(cfg$tables_dir, "sensitivity_FOT_paired_two_pairs_unconfirmed_design.csv")
)

primary <- diff_obj$primary_results
comparison <- dplyr::full_join(
  primary[, c(
    "Protein_accession", "Gene_name", "Unique_peptides",
    "log2FC_Polysome_vs_Free", "P_Value", "BH_FDR", "TREAT_BH_FDR",
    "TREAT_call", "BH_plus_FC_call", "legacy_rawP_plus_FC_call"
  )] |>
    dplyr::rename_with(~ paste0(.x, "_FOT"), -c(Protein_accession, Gene_name, Unique_peptides)),
  raw_res[, c(
    "Protein_accession", "Gene_name", "Unique_peptides",
    "log2FC_Polysome_vs_Free", "P_Value", "BH_FDR", "TREAT_BH_FDR",
    "TREAT_call", "BH_plus_FC_call", "legacy_rawP_plus_FC_call"
  )] |>
    dplyr::rename_with(~ paste0(.x, "_rawIntensity"), -c(Protein_accession, Gene_name, Unique_peptides)),
  by = c("Protein_accession", "Gene_name", "Unique_peptides")
)
comparison <- dplyr::full_join(
  comparison,
  med_res[, c(
    "Protein_accession", "Gene_name", "Unique_peptides",
    "log2FC_Polysome_vs_Free", "P_Value", "BH_FDR", "TREAT_BH_FDR",
    "TREAT_call", "BH_plus_FC_call", "legacy_rawP_plus_FC_call"
  )] |>
    dplyr::rename_with(~ paste0(.x, "_medianCentered"), -c(Protein_accession, Gene_name, Unique_peptides)),
  by = c("Protein_accession", "Gene_name", "Unique_peptides")
)
write_csv_stable(comparison, file.path(cfg$tables_dir, "normalization_method_protein_comparison.csv"))

safe_cor <- function(x, y, method = "pearson") {
  stats::cor(x, y, use = "pairwise.complete.obs", method = method)
}
method_summary <- data.frame(
  comparison = c("FOT_vs_rawIntensity", "FOT_vs_medianCentered", "rawIntensity_vs_medianCentered"),
  pearson_log2FC = c(
    safe_cor(comparison$log2FC_Polysome_vs_Free_FOT, comparison$log2FC_Polysome_vs_Free_rawIntensity),
    safe_cor(comparison$log2FC_Polysome_vs_Free_FOT, comparison$log2FC_Polysome_vs_Free_medianCentered),
    safe_cor(comparison$log2FC_Polysome_vs_Free_rawIntensity, comparison$log2FC_Polysome_vs_Free_medianCentered)
  ),
  spearman_log2FC = c(
    safe_cor(comparison$log2FC_Polysome_vs_Free_FOT, comparison$log2FC_Polysome_vs_Free_rawIntensity, "spearman"),
    safe_cor(comparison$log2FC_Polysome_vs_Free_FOT, comparison$log2FC_Polysome_vs_Free_medianCentered, "spearman"),
    safe_cor(comparison$log2FC_Polysome_vs_Free_rawIntensity, comparison$log2FC_Polysome_vs_Free_medianCentered, "spearman")
  ),
  median_log2FC_shift_first_minus_second = c(
    stats::median(comparison$log2FC_Polysome_vs_Free_FOT - comparison$log2FC_Polysome_vs_Free_rawIntensity, na.rm = TRUE),
    stats::median(comparison$log2FC_Polysome_vs_Free_FOT - comparison$log2FC_Polysome_vs_Free_medianCentered, na.rm = TRUE),
    stats::median(comparison$log2FC_Polysome_vs_Free_rawIntensity - comparison$log2FC_Polysome_vs_Free_medianCentered, na.rm = TRUE)
  ),
  direction_agreement = c(
    mean(sign(comparison$log2FC_Polysome_vs_Free_FOT) == sign(comparison$log2FC_Polysome_vs_Free_rawIntensity), na.rm = TRUE),
    mean(sign(comparison$log2FC_Polysome_vs_Free_FOT) == sign(comparison$log2FC_Polysome_vs_Free_medianCentered), na.rm = TRUE),
    mean(sign(comparison$log2FC_Polysome_vs_Free_rawIntensity) == sign(comparison$log2FC_Polysome_vs_Free_medianCentered), na.rm = TRUE)
  )
)
write_csv_stable(method_summary, file.path(cfg$tables_dir, "normalization_method_summary.csv"))

method_call_counts <- dplyr::bind_rows(
  data.frame(method = "FOT", call = primary$TREAT_call),
  data.frame(method = "rawIntensity", call = raw_res$TREAT_call),
  data.frame(method = "medianCentered", call = med_res$TREAT_call)
) |>
  dplyr::count(method, call, name = "n") |>
  tidyr::complete(
    method,
    call = c("Up", "Down", "Not_significant"),
    fill = list(n = 0L)
  )
write_csv_stable(method_call_counts, file.path(cfg$tables_dir, "normalization_method_TREAT_call_counts.csv"))

# Compare against the historical 502-protein list without redefining it.
if (file.exists(cfg$old_up_502_tsv)) {
  old502 <- readr::read_tsv(cfg$old_up_502_tsv, show_col_types = FALSE)
  old_accessions <- unique(stats::na.omit(as.character(old502$protein_accession)))
  new_sets <- list(
    FOT_TREAT_BH05 = unique(primary$Protein_accession[primary$TREAT_call == "Up"]),
    FOT_eBayes_BH05_FC1 = unique(primary$Protein_accession[primary$BH_plus_FC_call == "Up"]),
    FOT_legacy_rawP05_FC1 = unique(primary$Protein_accession[primary$legacy_rawP_plus_FC_call == "Up"]),
    rawIntensity_TREAT_BH05 = unique(raw_res$Protein_accession[raw_res$TREAT_call == "Up"])
  )
  overlap_summary <- dplyr::bind_rows(lapply(names(new_sets), function(nm) {
    new <- new_sets[[nm]]
    data.frame(
      comparison_set = nm,
      historical_502_n = length(old_accessions),
      new_set_n = length(new),
      overlap_n = length(intersect(old_accessions, new)),
      historical_only_n = length(setdiff(old_accessions, new)),
      new_only_n = length(setdiff(new, old_accessions))
    )
  }))
  write_csv_stable(overlap_summary, file.path(cfg$tables_dir, "historical_502_overlap_summary.csv"))

  union_accessions <- Reduce(union, c(list(old_accessions), new_sets))
  overlap_membership <- data.frame(Protein_accession = union_accessions)
  overlap_membership$in_historical_502 <- overlap_membership$Protein_accession %in% old_accessions
  for (nm in names(new_sets)) {
    overlap_membership[[paste0("in_", nm)]] <- overlap_membership$Protein_accession %in% new_sets[[nm]]
  }
  gene_lookup <- dat$metadata[, c("Protein_accession", "Gene_name")]
  overlap_membership <- dplyr::left_join(overlap_membership, gene_lookup, by = "Protein_accession")
  write_csv_stable(overlap_membership, file.path(cfg$tables_dir, "historical_502_overlap_membership.csv"))
}

# Composition diagnostics ---------------------------------------------------
top_n <- c(1L, 10L, 20L, 50L, 100L, 250L, 500L)
cumulative_share <- dplyr::bind_rows(lapply(cfg$sample_ids, function(s) {
  x <- sort(dat$intensity[, s], decreasing = TRUE, na.last = NA)
  data.frame(
    sample_id = s,
    group = unname(cfg$groups[s]),
    top_n = top_n,
    cumulative_intensity_share = vapply(top_n, function(n) sum(utils::head(x, n)) / sum(x), numeric(1))
  )
}))
write_csv_stable(cumulative_share, file.path(cfg$tables_dir, "sample_topN_cumulative_intensity_share.csv"))

shift_table <- data.frame(
  sample_id = cfg$main_samples,
  group = unname(cfg$groups[cfg$main_samples]),
  total_normalized_intensity = colSums(dat$intensity[, cfg$main_samples, drop = FALSE], na.rm = TRUE),
  log2_total_normalized_intensity = log2(colSums(dat$intensity[, cfg$main_samples, drop = FALSE], na.rm = TRUE)),
  median_log2_intensity = apply(log2(dat$intensity[, cfg$main_samples, drop = FALSE]), 2, stats::median, na.rm = TRUE),
  stringsAsFactors = FALSE
)
write_csv_stable(shift_table, file.path(cfg$tables_dir, "FOT_sample_scaling_shifts.csv"))

composition_plot_data <- cumulative_share[cumulative_share$sample_id %in% cfg$main_samples, ]
composition_plot_data$sample_id <- factor(composition_plot_data$sample_id, levels = cfg$main_samples)
composition_colors <- c(
  free_1 = "#4C78A8", free_2 = "#72A0CF", free_3 = "#2F5597",
  polysome_1 = "#E45756", polysome_2 = "#F28E8C"
)
p_composition <- ggplot2::ggplot(
  composition_plot_data,
  ggplot2::aes(top_n, cumulative_intensity_share, color = sample_id, linetype = group)
) +
  ggplot2::geom_line(linewidth = 0.8) +
  ggplot2::geom_point(size = 2) +
  ggplot2::scale_color_manual(values = composition_colors, drop = FALSE) +
  ggplot2::scale_x_log10(breaks = top_n) +
  ggplot2::scale_y_continuous(labels = scales::label_percent(accuracy = 1), limits = c(0, 1)) +
  ggplot2::labs(
    x = "Top N proteins by normalized intensity (log scale)",
    y = "Cumulative share of sample signal",
    color = "Sample",
    linetype = "Group",
    title = "Sample composition concentration"
  ) +
  theme_fot()
save_plot_both(p_composition, "QC_topN_composition_share", 7.6, 5.5)

scatter_df <- comparison[stats::complete.cases(
  comparison$log2FC_Polysome_vs_Free_FOT,
  comparison$log2FC_Polysome_vs_Free_rawIntensity
), ]
p_scatter <- ggplot2::ggplot(
  scatter_df,
  ggplot2::aes(log2FC_Polysome_vs_Free_rawIntensity, log2FC_Polysome_vs_Free_FOT)
) +
  ggplot2::geom_point(alpha = 0.55, size = 1.15, color = "#3B6FB6") +
  ggplot2::geom_abline(slope = 1, intercept = 0, linetype = 2, color = "grey40") +
  ggplot2::geom_smooth(method = "lm", se = FALSE, color = "#D1495B", linewidth = 0.8) +
  ggplot2::labs(
    x = "log2FC from vendor normalized intensity",
    y = "log2FC from intensity-derived FOT",
    title = "Normalization sensitivity of fold changes"
  ) +
  theme_fot()
save_plot_both(p_scatter, "comparison_FOT_vs_intensity_log2FC", 6.5, 5.7)

saveRDS(
  list(
    raw_intensity_results = raw_res,
    median_centered_results = med_res,
    paired_two_pair_results = paired_res,
    comparison = comparison,
    method_summary = method_summary,
    method_call_counts = method_call_counts,
    cumulative_share = cumulative_share,
    sample_scaling_shifts = shift_table
  ),
  file.path(cfg$data_dir, "sensitivity_analysis_objects.rds"),
  version = 3
)

message("[03] Sensitivity analyses and historical-set comparison completed")
