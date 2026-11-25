message("[01] Reading authoritative protein table and preparing FOT/QC data")

stopifnot(file.exists(cfg$input_xlsx))
raw <- as.data.frame(readxl::read_xlsx(cfg$input_xlsx, sheet = cfg$input_sheet))

required_meta <- c("Protein accession", "Protein description", "Gene name", "Unique peptides")
intensity_cols <- paste0("Intensity ", cfg$sample_ids)
relative_cols <- cfg$sample_ids
required_cols <- c(required_meta, intensity_cols, relative_cols)
missing_cols <- setdiff(required_cols, names(raw))
if (length(missing_cols)) {
  stop("Required columns absent from source workbook: ", paste(missing_cols, collapse = ", "))
}

metadata <- data.frame(
  row_id = seq_len(nrow(raw)),
  Protein_accession = as.character(raw[["Protein accession"]]),
  Gene_name = as.character(raw[["Gene name"]]),
  Protein_description = as.character(raw[["Protein description"]]),
  Unique_peptides = as.numeric(raw[["Unique peptides"]]),
  stringsAsFactors = FALSE
)
metadata$Gene_name[metadata$Gene_name %in% c("", "--", "NA")] <- NA_character_
metadata$explicit_decoy_or_contaminant <-
  grepl("^(CON__|REV__|DECOY)", metadata$Protein_accession, ignore.case = TRUE) |
  grepl("contaminant", metadata$Protein_description, ignore.case = TRUE)

protein_keys <- make.unique(ifelse(
  is.na(metadata$Protein_accession) | !nzchar(metadata$Protein_accession),
  paste0("row_", metadata$row_id),
  metadata$Protein_accession
))

intensity <- clean_positive_matrix(raw[, intensity_cols, drop = FALSE])
colnames(intensity) <- cfg$sample_ids
rownames(intensity) <- protein_keys

reported_relative <- clean_positive_matrix(raw[, relative_cols, drop = FALSE])
colnames(reported_relative) <- cfg$sample_ids
rownames(reported_relative) <- protein_keys

# The vendor report defines the unprefixed columns as R_ij = I_ij / mean_i(I_ij).
row_mean_intensity <- rowMeans(intensity, na.rm = TRUE)
row_mean_intensity[!is.finite(row_mean_intensity) | row_mean_intensity <= 0] <- NA_real_
reconstructed_relative <- sweep(intensity, 1, row_mean_intensity, "/")
relative_abs_error <- abs(reconstructed_relative - reported_relative)
relative_max_abs_error <- max(relative_abs_error, na.rm = TRUE)

# FOT denominator uses all positive normalized intensities before applying the
# Unique-peptide analysis filter. Only explicit decoy/contaminant accessions are excluded.
denominator_keep <- !metadata$explicit_decoy_or_contaminant
denominators <- colSums(intensity[denominator_keep, , drop = FALSE], na.rm = TRUE)
if (any(!is.finite(denominators) | denominators <= 0)) {
  stop("At least one FOT denominator is non-positive or non-finite.")
}
fot <- sweep(intensity, 2, denominators, "/")
fot[!denominator_keep, ] <- NA_real_
ifot <- fot * cfg$ifot_scale
log2_ifot <- log2(ifot)

top_share <- function(x, n) {
  x <- sort(x[is.finite(x) & x > 0], decreasing = TRUE)
  if (!length(x)) return(NA_real_)
  sum(utils::head(x, n)) / sum(x)
}

sample_qc <- data.frame(
  sample_id = cfg$sample_ids,
  group = unname(cfg$groups[cfg$sample_ids]),
  included_in_main_analysis = cfg$sample_ids %in% cfg$main_samples,
  detected_proteins = colSums(!is.na(intensity)),
  missing_proteins = colSums(is.na(intensity)),
  missing_fraction = colMeans(is.na(intensity)),
  total_normalized_intensity = colSums(intensity, na.rm = TRUE),
  median_positive_intensity = apply(intensity, 2, stats::median, na.rm = TRUE),
  FOT_denominator = denominators,
  FOT_sum = colSums(fot, na.rm = TRUE),
  iFOT_sum = colSums(ifot, na.rm = TRUE),
  top1_intensity_share = apply(intensity, 2, top_share, n = 1),
  top10_intensity_share = apply(intensity, 2, top_share, n = 10),
  top100_intensity_share = apply(intensity, 2, top_share, n = 100),
  stringsAsFactors = FALSE
)

cor_mat <- stats::cor(log2(intensity), use = "pairwise.complete.obs", method = "pearson")
detected <- !is.na(intensity)
jaccard_mat <- outer(seq_along(cfg$sample_ids), seq_along(cfg$sample_ids), Vectorize(
  function(i, j) sum(detected[, i] & detected[, j]) / sum(detected[, i] | detected[, j])
))
dimnames(jaccard_mat) <- list(cfg$sample_ids, cfg$sample_ids)

sample_metadata <- data.frame(
  sample_id = cfg$sample_ids,
  group = unname(cfg$groups[cfg$sample_ids]),
  replicate = sub(".*_", "", cfg$sample_ids),
  included_in_main_analysis = cfg$sample_ids %in% cfg$main_samples,
  exclusion_reason = ifelse(
    cfg$sample_ids == cfg$excluded_sample,
    "Pre-specified QC exclusion: only 970 proteins detected and low correlation with polysome_1/2",
    ""
  ),
  stringsAsFactors = FALSE
)

write_csv_stable(sample_qc, file.path(cfg$tables_dir, "sample_QC_metrics.csv"))
write_csv_stable(sample_metadata, file.path(cfg$tables_dir, "sample_metadata.csv"))
write_csv_stable(
  cbind(sample_id = rownames(cor_mat), as.data.frame(cor_mat, check.names = FALSE)),
  file.path(cfg$tables_dir, "sample_log2_intensity_Pearson_correlation.csv")
)
write_csv_stable(
  cbind(sample_id = rownames(jaccard_mat), as.data.frame(jaccard_mat, check.names = FALSE)),
  file.path(cfg$tables_dir, "sample_detection_Jaccard.csv")
)

wide_values <- metadata
for (s in cfg$sample_ids) {
  wide_values[[paste0("NormalizedIntensity_", s)]] <- intensity[, s]
  wide_values[[paste0("FOT_", s)]] <- fot[, s]
  wide_values[[paste0("iFOT_", s)]] <- ifot[, s]
  wide_values[[paste0("log2_iFOT_", s)]] <- log2_ifot[, s]
}
write_csv_stable(wide_values, file.path(cfg$data_dir, "FOT_values_all_proteins.csv"))

source_info <- file.info(c(cfg$input_xlsx, cfg$report_pdf))
source_manifest <- data.frame(
  role = c("authoritative_quantification_input", "vendor_method_report"),
  path = rownames(source_info),
  md5 = unname(tools::md5sum(rownames(source_info))),
  size_bytes = source_info$size,
  modified_time = format(source_info$mtime, tz = "Asia/Shanghai", usetz = TRUE),
  stringsAsFactors = FALSE
)
write_csv_stable(source_manifest, file.path(cfg$tables_dir, "source_manifest.csv"))

prepared <- list(
  metadata = metadata,
  protein_keys = protein_keys,
  intensity = intensity,
  reported_relative = reported_relative,
  fot = fot,
  ifot = ifot,
  log2_ifot = log2_ifot,
  denominator_keep = denominator_keep,
  sample_qc = sample_qc,
  sample_metadata = sample_metadata,
  correlation = cor_mat,
  jaccard = jaccard_mat,
  source_manifest = source_manifest,
  relative_max_abs_error = relative_max_abs_error,
  config_snapshot = cfg
)
saveRDS(prepared, file.path(cfg$data_dir, "prepared_FOT_data.rds"), version = 3)

# QC visualizations ---------------------------------------------------------
sample_levels <- cfg$sample_ids
qc_plot_data <- sample_qc
qc_plot_data$sample_id <- factor(qc_plot_data$sample_id, levels = sample_levels)

p_detect <- ggplot2::ggplot(qc_plot_data, ggplot2::aes(sample_id, detected_proteins, fill = group)) +
  ggplot2::geom_col(width = 0.72) +
  ggplot2::geom_text(ggplot2::aes(label = scales::comma(detected_proteins)), vjust = -0.3, size = 3) +
  ggplot2::scale_fill_manual(values = c(Free = "#4C78A8", Polysome = "#E45756")) +
  ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0, 0.10))) +
  ggplot2::labs(x = NULL, y = "Detected proteins", title = "Protein detection depth") +
  theme_fot() +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 35, hjust = 1), legend.position = "none")

p_total <- ggplot2::ggplot(qc_plot_data, ggplot2::aes(sample_id, total_normalized_intensity, fill = group)) +
  ggplot2::geom_col(width = 0.72) +
  ggplot2::scale_fill_manual(values = c(Free = "#4C78A8", Polysome = "#E45756")) +
  ggplot2::scale_y_continuous(labels = scales::label_scientific()) +
  ggplot2::labs(x = NULL, y = "Total normalized intensity", title = "Pre-FOT sample totals") +
  theme_fot() +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 35, hjust = 1), legend.position = "none")
save_plot_both(p_detect + p_total, "QC_detection_and_total_intensity", 10, 4.8)

intensity_long <- as.data.frame(intensity)
intensity_long$row_id <- seq_len(nrow(intensity_long))
intensity_long <- tidyr::pivot_longer(
  intensity_long,
  cols = -row_id,
  names_to = "sample_id",
  values_to = "intensity"
)
intensity_long <- intensity_long[is.finite(intensity_long$intensity) & intensity_long$intensity > 0, ]
intensity_long$sample_id <- factor(intensity_long$sample_id, levels = sample_levels)
intensity_long$group <- unname(cfg$groups[as.character(intensity_long$sample_id)])
p_density <- ggplot2::ggplot(
  intensity_long,
  ggplot2::aes(log2(intensity), color = sample_id, linetype = group)
) +
  ggplot2::geom_density(linewidth = 0.75, adjust = 1.1) +
  ggplot2::labs(x = "log2 normalized intensity", y = "Density", color = "Sample", linetype = "Group") +
  theme_fot()
save_plot_both(p_density, "QC_log2_intensity_density", 8, 5)

cor_long <- as.data.frame(as.table(cor_mat), stringsAsFactors = FALSE)
names(cor_long) <- c("sample_x", "sample_y", "correlation")
cor_long$sample_x <- factor(cor_long$sample_x, levels = rev(sample_levels))
cor_long$sample_y <- factor(cor_long$sample_y, levels = sample_levels)
p_cor <- ggplot2::ggplot(cor_long, ggplot2::aes(sample_y, sample_x, fill = correlation)) +
  ggplot2::geom_tile(color = "white", linewidth = 0.5) +
  ggplot2::geom_text(ggplot2::aes(label = sprintf("%.2f", correlation)), size = 3.2) +
  ggplot2::scale_fill_gradient2(low = "#3B4CC0", mid = "white", high = "#B40426", midpoint = 0.8, limits = c(0.5, 1)) +
  ggplot2::coord_equal() +
  ggplot2::labs(x = NULL, y = NULL, fill = "Pearson r", title = "Pairwise log2-intensity correlation") +
  theme_fot() +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 35, hjust = 1))
save_plot_both(p_cor, "QC_sample_correlation", 7.2, 6.2)

complete_all <- stats::complete.cases(log2_ifot)
pca_input <- log2_ifot[complete_all, , drop = FALSE]
if (nrow(pca_input) >= 3) {
  vars <- matrixStats::rowVars(pca_input)
  top_idx <- order(vars, decreasing = TRUE)[seq_len(min(1000L, length(vars)))]
  pca <- stats::prcomp(t(pca_input[top_idx, , drop = FALSE]), center = TRUE, scale. = FALSE)
  pca_var <- 100 * pca$sdev^2 / sum(pca$sdev^2)
  pca_df <- data.frame(
    sample_id = rownames(pca$x),
    PC1 = pca$x[, 1],
    PC2 = pca$x[, 2],
    group = unname(cfg$groups[rownames(pca$x)]),
    excluded = rownames(pca$x) == cfg$excluded_sample,
    stringsAsFactors = FALSE
  )
  pca_label_offsets <- data.frame(
    sample_id = cfg$sample_ids,
    dx = c(-4.5, -4.5, -4.5, -4.5, -5.5, 0),
    dy = c(1.2, -1.4, 2.2, 2.6, -2.6, 3.2),
    stringsAsFactors = FALSE
  )
  pca_df <- dplyr::left_join(pca_df, pca_label_offsets, by = "sample_id")
  pca_df$label_x <- pca_df$PC1 + pca_df$dx
  pca_df$label_y <- pca_df$PC2 + pca_df$dy
  write_csv_stable(pca_df, file.path(cfg$tables_dir, "PCA_sample_coordinates_all6.csv"))
  p_pca <- ggplot2::ggplot(pca_df, ggplot2::aes(PC1, PC2, color = group, shape = excluded)) +
    ggplot2::geom_point(size = 3.5) +
    ggplot2::geom_text(
      ggplot2::aes(x = label_x, y = label_y, label = sample_id),
      show.legend = FALSE,
      size = 3.3
    ) +
    ggplot2::scale_color_manual(values = c(Free = "#4C78A8", Polysome = "#E45756")) +
    ggplot2::scale_shape_manual(values = c(`FALSE` = 16, `TRUE` = 8), labels = c(`FALSE` = "Included", `TRUE` = "Excluded")) +
    ggplot2::labs(
      x = sprintf("PC1 (%.1f%%)", pca_var[1]),
      y = sprintf("PC2 (%.1f%%)", pca_var[2]),
      color = "Group",
      shape = "Analysis status",
      title = sprintf("PCA of %s complete-case proteins", scales::comma(nrow(pca_input)))
    ) +
    theme_fot()
  save_plot_both(p_pca, "QC_PCA_all6_samples", 7.2, 5.5)
}

validation <- list(
  source_rows = nrow(raw),
  source_columns = ncol(raw),
  accession_is_unique = !anyDuplicated(metadata$Protein_accession),
  explicit_decoy_or_contaminant_rows = sum(metadata$explicit_decoy_or_contaminant),
  reported_R_reconstruction_max_abs_error = relative_max_abs_error,
  fot_sum_max_abs_deviation_from_1 = max(abs(colSums(fot, na.rm = TRUE) - 1)),
  ifot_sum_max_abs_deviation_from_scale = max(abs(colSums(ifot, na.rm = TRUE) - cfg$ifot_scale)),
  all_input_intensities_positive_or_missing = all(is.na(intensity) | intensity > 0),
  excluded_sample = cfg$excluded_sample
)
jsonlite::write_json(
  validation,
  file.path(cfg$logs_dir, "preparation_validation.json"),
  pretty = TRUE,
  auto_unbox = TRUE,
  digits = NA
)

message("[01] Prepared ", nrow(raw), " proteins; FOT sums validated for all samples")
