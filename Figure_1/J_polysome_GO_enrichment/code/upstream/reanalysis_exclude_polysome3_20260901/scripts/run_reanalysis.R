#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(readxl)
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(stringr)
  library(broom)
  library(ggplot2)
  library(jsonlite)
  library(scales)
})

options(stringsAsFactors = FALSE, scipen = 999)

free_samples <- c("free_1", "free_2", "free_3")
polysome_samples <- c("polysome_1", "polysome_2")
excluded_sample <- "polysome_3"
log2fc_threshold <- 1
unique_peptide_threshold <- 5
p_threshold <- 0.05
q_threshold <- 0.05

script_args <- commandArgs(trailingOnly = FALSE)
script_arg <- grep("^--file=", script_args, value = TRUE)
script_path <- if (length(script_arg) == 1) {
  normalizePath(sub("^--file=", "", script_arg), mustWork = TRUE)
} else {
  normalizePath("run_reanalysis.R", mustWork = TRUE)
}
scripts_dir <- dirname(script_path)
analysis_dir <- dirname(scripts_dir)
output_dir <- file.path(analysis_dir, "results")
project_dir <- dirname(analysis_dir)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
source_workbook <- file.path(
  project_dir,
  "cloud_report_312553",
  "3-Functional_annotation",
  "Annotation_combine.xlsx"
)

bh_adjust_na <- function(p_values) {
  adjusted <- rep(NA_real_, length(p_values))
  valid <- which(!is.na(p_values) & is.finite(p_values))
  if (length(valid) > 0) {
    adjusted[valid] <- p.adjust(p_values[valid], method = "BH")
  }
  adjusted
}

safe_t_test <- function(free_values, polysome_values) {
  free_values <- free_values[is.finite(free_values) & free_values > 0]
  polysome_values <- polysome_values[is.finite(polysome_values) & polysome_values > 0]
  if (length(free_values) < 2 || length(polysome_values) < 2) {
    return(NA_real_)
  }
  tryCatch(
    t.test(
      log2(polysome_values),
      log2(free_values),
      alternative = "two.sided",
      var.equal = TRUE
    )$p.value,
    error = function(e) NA_real_
  )
}

calculate_statistics <- function(annotation) {
  free_matrix <- annotation |>
    select(all_of(free_samples)) |>
    as.matrix()
  polysome_matrix <- annotation |>
    select(all_of(polysome_samples)) |>
    as.matrix()
  storage.mode(free_matrix) <- "double"
  storage.mode(polysome_matrix) <- "double"
  free_matrix[!is.finite(free_matrix) | free_matrix <= 0] <- NA_real_
  polysome_matrix[!is.finite(polysome_matrix) | polysome_matrix <= 0] <- NA_real_

  n_free <- rowSums(!is.na(free_matrix))
  n_polysome <- rowSums(!is.na(polysome_matrix))
  free_mean <- rowMeans(free_matrix, na.rm = TRUE)
  polysome_mean <- rowMeans(polysome_matrix, na.rm = TRUE)
  free_mean[n_free == 0] <- NA_real_
  polysome_mean[n_polysome == 0] <- NA_real_
  ratio <- polysome_mean / free_mean
  log2fc <- log2(ratio)

  p_values <- map_dbl(seq_len(nrow(annotation)), function(index) {
    safe_t_test(free_matrix[index, ], polysome_matrix[index, ])
  })

  keep_columns <- c(
    "Protein accession", "Protein description", "Gene name", "Coverage [%]",
    "Peptides", "Razor + unique peptides", "Unique peptides",
    free_samples, polysome_samples, "Subcellular localization", "KEGG pathway",
    "Biological Process", "Cellular Component", "Molecular Function", "Pfam ID",
    "Domain description", "Reactome pathway", "WikiPathways", "Is_TF"
  )

  statistics <- annotation |>
    select(any_of(keep_columns)) |>
    mutate(
      n_free = n_free,
      n_polysome = n_polysome,
      complete_retained_5_samples = n_free == 3 & n_polysome == 2,
      ratio_excluding_polysome_3 = ratio,
      log2FC_excluding_polysome_3 = log2fc,
      p_value_excluding_polysome_3 = p_values,
      q_value_excluding_polysome_3 = bh_adjust_na(p_values),
      old_ratio_including_polysome_3 = annotation[["polysome/free Ratio"]],
      old_log2FC_including_polysome_3 = log2(annotation[["polysome/free Ratio"]]),
      old_p_value_including_polysome_3 = annotation[["polysome/free P value"]],
      old_q_value_including_polysome_3 = bh_adjust_na(annotation[["polysome/free P value"]])
    ) |>
    mutate(
      passes_effect_and_peptide_filters = coalesce(
        abs(log2FC_excluding_polysome_3) > log2fc_threshold &
          .data[["Unique peptides"]] >= unique_peptide_threshold,
        FALSE
      ),
      significant_raw_p = coalesce(
        passes_effect_and_peptide_filters &
          p_value_excluding_polysome_3 < p_threshold,
        FALSE
      ),
      significant_BH_q = coalesce(
        passes_effect_and_peptide_filters &
          q_value_excluding_polysome_3 < q_threshold,
        FALSE
      ),
      direction = case_when(
        log2FC_excluding_polysome_3 > log2fc_threshold ~ "Up",
        log2FC_excluding_polysome_3 < -log2fc_threshold ~ "Down",
        TRUE ~ "Not differential"
      ),
      old_significant_raw_p_same_filters = coalesce(
        abs(old_log2FC_including_polysome_3) > log2fc_threshold &
          .data[["Unique peptides"]] >= unique_peptide_threshold &
          old_p_value_including_polysome_3 < p_threshold,
        FALSE
      ),
      old_significant_BH_q_same_filters = coalesce(
        abs(old_log2FC_including_polysome_3) > log2fc_threshold &
          .data[["Unique peptides"]] >= unique_peptide_threshold &
          old_q_value_including_polysome_3 < q_threshold,
        FALSE
      ),
      raw_status_vs_old = case_when(
        significant_raw_p & old_significant_raw_p_same_filters ~ "Retained",
        significant_raw_p & !old_significant_raw_p_same_filters ~ "Gained after exclusion",
        !significant_raw_p & old_significant_raw_p_same_filters ~ "Lost after exclusion",
        TRUE ~ "Not significant in either"
      ),
      BH_status_vs_old = case_when(
        significant_BH_q & old_significant_BH_q_same_filters ~ "Retained",
        significant_BH_q & !old_significant_BH_q_same_filters ~ "Gained after exclusion",
        !significant_BH_q & old_significant_BH_q_same_filters ~ "Lost after exclusion",
        TRUE ~ "Not significant in either"
      )
    )

  statistics
}

split_terms <- function(values) {
  values |>
    replace_na("") |>
    str_split(";")
}

term_membership <- function(statistics, background_accessions, column_name) {
  statistics |>
    filter(.data[["Protein accession"]] %in% background_accessions) |>
    transmute(
      protein_accession = .data[["Protein accession"]],
      gene_name = .data[["Gene name"]],
      raw_term = .data[[column_name]]
    ) |>
    mutate(term = split_terms(raw_term)) |>
    unnest_longer(term) |>
    mutate(
      term = str_squish(term),
      term_id = word(term, 1),
      term_name = str_remove(term, "^[^[:space:]]+[[:space:]]*")
    ) |>
    filter(term != "", term_id != "") |>
    mutate(term_name = if_else(term_name == "", term_id, term_name)) |>
    distinct(protein_accession, gene_name, term_id, term_name)
}

calculate_enrichment <- function(statistics) {
  eligible <- !is.na(statistics$p_value_excluding_polysome_3) &
    statistics[["Unique peptides"]] >= unique_peptide_threshold
  background_accessions <- statistics[["Protein accession"]][eligible]
  background_count <- length(unique(background_accessions))

  databases <- c(
    "GO Biological Process" = "Biological Process",
    "GO Cellular Component" = "Cellular Component",
    "GO Molecular Function" = "Molecular Function",
    "KEGG" = "KEGG pathway",
    "Reactome" = "Reactome pathway",
    "Pfam" = "Domain description"
  )
  selections <- c(
    "Raw P < 0.05" = "significant_raw_p",
    "BH q < 0.05" = "significant_BH_q"
  )

  all_results <- list()
  result_index <- 1L
  for (database in names(databases)) {
    membership <- term_membership(
      statistics,
      background_accessions,
      unname(databases[[database]])
    )
    term_sizes <- membership |>
      count(term_id, term_name, name = "background_term_count")

    for (selection_method in names(selections)) {
      selection_column <- unname(selections[[selection_method]])
      for (direction_label in c("Up", "Down")) {
        foreground_accessions <- statistics |>
          filter(
            eligible,
            .data[[selection_column]],
            direction == direction_label
          ) |>
          pull(.data[["Protein accession"]]) |>
          unique()
        foreground_count <- length(foreground_accessions)
        if (foreground_count == 0) next

        overlap <- membership |>
          filter(protein_accession %in% foreground_accessions) |>
          group_by(term_id, term_name) |>
          summarise(
            mapping_count = n_distinct(protein_accession),
            genes = paste(sort(unique(na.omit(gene_name))), collapse = ";"),
            protein_accessions = paste(sort(unique(protein_accession)), collapse = ";"),
            .groups = "drop"
          ) |>
          inner_join(term_sizes, by = c("term_id", "term_name")) |>
          mutate(
            selection_method = selection_method,
            direction = direction_label,
            database = database,
            foreground_count = foreground_count,
            background_count = background_count,
            fold_enrichment =
              (mapping_count / foreground_count) /
              (background_term_count / background_count),
            p_value = phyper(
              mapping_count - 1,
              background_term_count,
              background_count - background_term_count,
              foreground_count,
              lower.tail = FALSE
            )
          ) |>
          select(
            selection_method, direction, database, term_id, term_name,
            mapping_count, foreground_count, background_term_count,
            background_count, fold_enrichment, p_value, genes,
            protein_accessions
          )
        all_results[[result_index]] <- overlap
        result_index <- result_index + 1L
      }
    }
  }

  bind_rows(all_results) |>
    group_by(selection_method, direction, database) |>
    mutate(q_value = p.adjust(p_value, method = "BH")) |>
    ungroup() |>
    mutate(
      significant_raw_p =
        p_value < p_threshold & fold_enrichment > 1.5 & mapping_count >= 2,
      significant_BH_q =
        q_value < q_threshold & fold_enrichment > 1.5 & mapping_count >= 2
    ) |>
    arrange(selection_method, direction, database, p_value, desc(fold_enrichment))
}

build_summary <- function(statistics, enrichment) {
  tested <- !is.na(statistics$p_value_excluding_polysome_3)
  peptide_eligible <- tested &
    statistics[["Unique peptides"]] >= unique_peptide_threshold
  raw <- statistics$significant_raw_p
  bh <- statistics$significant_BH_q
  old_raw <- statistics$old_significant_raw_p_same_filters
  old_bh <- statistics$old_significant_BH_q_same_filters

  list(
    analysis_date = "2026-09-01",
    comparison = "polysome_1-2 vs free_1-3",
    excluded_sample = excluded_sample,
    method = "Two-sided equal-variance t-test on log2 relative protein quantities",
    effect_filter = "abs(log2FC) > 1",
    peptide_filter = "Unique peptides >= 5",
    tested_proteins = sum(tested),
    tested_and_peptide_eligible = sum(peptide_eligible),
    complete_retained_5_samples = sum(tested & statistics$complete_retained_5_samples),
    new_raw_p_total = sum(raw, na.rm = TRUE),
    new_raw_p_up = sum(raw & statistics$direction == "Up", na.rm = TRUE),
    new_raw_p_down = sum(raw & statistics$direction == "Down", na.rm = TRUE),
    new_bh_q_total = sum(bh, na.rm = TRUE),
    new_bh_q_up = sum(bh & statistics$direction == "Up", na.rm = TRUE),
    new_bh_q_down = sum(bh & statistics$direction == "Down", na.rm = TRUE),
    new_raw_p_complete_case_total =
      sum(raw & statistics$complete_retained_5_samples, na.rm = TRUE),
    new_bh_q_complete_case_total =
      sum(bh & statistics$complete_retained_5_samples, na.rm = TRUE),
    old_raw_p_total_same_filters = sum(old_raw, na.rm = TRUE),
    old_bh_q_total_same_filters = sum(old_bh, na.rm = TRUE),
    raw_retained = sum(statistics$raw_status_vs_old == "Retained", na.rm = TRUE),
    raw_gained = sum(statistics$raw_status_vs_old == "Gained after exclusion", na.rm = TRUE),
    raw_lost = sum(statistics$raw_status_vs_old == "Lost after exclusion", na.rm = TRUE),
    bh_retained = sum(statistics$BH_status_vs_old == "Retained", na.rm = TRUE),
    bh_gained = sum(statistics$BH_status_vs_old == "Gained after exclusion", na.rm = TRUE),
    bh_lost = sum(statistics$BH_status_vs_old == "Lost after exclusion", na.rm = TRUE),
    enrichment_rows = nrow(enrichment),
    raw_significant_enrichment_terms = sum(enrichment$significant_raw_p, na.rm = TRUE),
    bh_significant_enrichment_terms = sum(enrichment$significant_BH_q, na.rm = TRUE)
  )
}

write_outputs <- function(statistics, enrichment, summary) {
  statistics <- statistics |>
    arrange(desc(significant_BH_q), q_value_excluding_polysome_3, p_value_excluding_polysome_3)

  write.csv(statistics, file.path(output_dir, "protein_statistics_all.csv"), row.names = FALSE, na = "")
  write.csv(
    filter(statistics, significant_raw_p),
    file.path(output_dir, "differential_proteins_raw_p005.csv"),
    row.names = FALSE,
    na = ""
  )
  write.csv(
    filter(statistics, significant_BH_q),
    file.path(output_dir, "differential_proteins_BH_q005.csv"),
    row.names = FALSE,
    na = ""
  )
  write.csv(
    filter(statistics, significant_raw_p, complete_retained_5_samples),
    file.path(output_dir, "differential_proteins_complete_case_raw_p005.csv"),
    row.names = FALSE,
    na = ""
  )
  write.csv(
    filter(statistics, raw_status_vs_old %in% c("Gained after exclusion", "Lost after exclusion")),
    file.path(output_dir, "sensitivity_changes_raw_p.csv"),
    row.names = FALSE,
    na = ""
  )
  write.csv(
    filter(statistics, BH_status_vs_old %in% c("Gained after exclusion", "Lost after exclusion")),
    file.path(output_dir, "sensitivity_changes_BH_q.csv"),
    row.names = FALSE,
    na = ""
  )
  write.csv(enrichment, file.path(output_dir, "enrichment_all.csv"), row.names = FALSE, na = "")
  write.csv(
    filter(enrichment, significant_raw_p),
    file.path(output_dir, "enrichment_significant_raw_p.csv"),
    row.names = FALSE,
    na = ""
  )
  write.csv(
    filter(enrichment, significant_BH_q),
    file.path(output_dir, "enrichment_significant_BH_q.csv"),
    row.names = FALSE,
    na = ""
  )
  write_json(
    summary,
    file.path(output_dir, "summary.json"),
    pretty = TRUE,
    auto_unbox = TRUE,
    na = "null"
  )

  package_names <- c(
    "readxl", "dplyr", "tidyr", "purrr", "stringr", "broom",
    "ggplot2", "jsonlite", "rmarkdown", "knitr", "scales"
  )
  package_versions <- tibble(
    package = package_names,
    version = map_chr(package_names, ~ as.character(packageVersion(.x)))
  )
  write.csv(
    package_versions,
    file.path(output_dir, "R_package_versions.csv"),
    row.names = FALSE
  )
  capture.output(sessionInfo(), file = file.path(output_dir, "sessionInfo.txt"))

  invisible(statistics)
}

make_plots <- function(statistics, enrichment) {
  volcano_data <- statistics |>
    filter(!is.na(p_value_excluding_polysome_3)) |>
    mutate(
      minus_log10_p = -log10(pmax(p_value_excluding_polysome_3, .Machine$double.xmin)),
      category = case_when(
        significant_raw_p & direction == "Up" ~ "Up",
        significant_raw_p & direction == "Down" ~ "Down",
        TRUE ~ "Not significant"
      )
    )
  volcano <- ggplot(
    volcano_data,
    aes(log2FC_excluding_polysome_3, minus_log10_p, color = category)
  ) +
    geom_point(size = 0.9, alpha = 0.65) +
    geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "grey35") +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey35") +
    scale_color_manual(values = c(Up = "#d62728", Down = "#1f77b4", `Not significant` = "#bdbdbd")) +
    labs(
      title = "Differential proteins after excluding polysome_3",
      subtitle = "|log2FC| > 1, unique peptides ≥ 5, P < 0.05",
      x = "log2 fold change (polysome / free)",
      y = "-log10 P value",
      color = NULL
    ) +
    theme_minimal(base_size = 12) +
    theme(legend.position = "top")
  ggsave(
    file.path(output_dir, "volcano_excluding_polysome3.png"),
    volcano,
    width = 8.5,
    height = 6.5,
    dpi = 300
  )

  comparison <- tibble(
    Analysis = c("Original (same filters)", "Exclude polysome_3"),
    `Raw P Up` = c(
      sum(statistics$old_significant_raw_p_same_filters & statistics$old_log2FC_including_polysome_3 > 1, na.rm = TRUE),
      sum(statistics$significant_raw_p & statistics$direction == "Up", na.rm = TRUE)
    ),
    `Raw P Down` = c(
      sum(statistics$old_significant_raw_p_same_filters & statistics$old_log2FC_including_polysome_3 < -1, na.rm = TRUE),
      sum(statistics$significant_raw_p & statistics$direction == "Down", na.rm = TRUE)
    ),
    `BH q Up` = c(
      sum(statistics$old_significant_BH_q_same_filters & statistics$old_log2FC_including_polysome_3 > 1, na.rm = TRUE),
      sum(statistics$significant_BH_q & statistics$direction == "Up", na.rm = TRUE)
    ),
    `BH q Down` = c(
      sum(statistics$old_significant_BH_q_same_filters & statistics$old_log2FC_including_polysome_3 < -1, na.rm = TRUE),
      sum(statistics$significant_BH_q & statistics$direction == "Down", na.rm = TRUE)
    )
  )
  write.csv(comparison, file.path(output_dir, "differential_count_comparison.csv"), row.names = FALSE)
  comparison_long <- comparison |>
    pivot_longer(-Analysis, names_to = "Criterion", values_to = "Protein_count")
  count_plot <- ggplot(comparison_long, aes(Analysis, Protein_count, fill = Criterion)) +
    geom_col(position = position_dodge(width = 0.8), width = 0.72) +
    scale_fill_manual(values = c(
      `Raw P Up` = "#ef8a62", `Raw P Down` = "#67a9cf",
      `BH q Up` = "#b2182b", `BH q Down` = "#2166ac"
    )) +
    labs(
      title = "Sensitivity of differential-protein counts to polysome_3 exclusion",
      x = NULL,
      y = "Protein count",
      fill = NULL
    ) +
    theme_minimal(base_size = 12) +
    theme(legend.position = "top")
  ggsave(
    file.path(output_dir, "differential_count_comparison.png"),
    count_plot,
    width = 8.5,
    height = 5.5,
    dpi = 300
  )

  top_enrichment <- enrichment |>
    filter(selection_method == "Raw P < 0.05", significant_raw_p) |>
    group_by(direction) |>
    slice_min(p_value, n = 12, with_ties = FALSE) |>
    ungroup() |>
    mutate(
      label = paste0(str_remove(database, "^GO "), ": ", term_name),
      label = reorder(label, -log10(p_value))
    )
  if (nrow(top_enrichment) > 0) {
    enrichment_plot <- ggplot(
      top_enrichment,
      aes(-log10(p_value), label, size = mapping_count, color = direction)
    ) +
      geom_point(alpha = 0.78) +
      facet_wrap(~ direction, scales = "free_y") +
      scale_color_manual(values = c(Up = "#d62728", Down = "#1f77b4")) +
      labs(
        title = "Top functional enrichments after excluding polysome_3",
        x = "-log10 enrichment P value",
        y = NULL,
        size = "Mapping count",
        color = NULL
      ) +
      theme_minimal(base_size = 10) +
      theme(legend.position = "top")
    ggsave(
      file.path(output_dir, "top_functional_enrichment.png"),
      enrichment_plot,
      width = 14,
      height = 8,
      dpi = 300
    )
  }
}

validate_results <- function(statistics, enrichment, summary) {
  raw_subset <- filter(statistics, significant_raw_p)
  bh_subset <- filter(statistics, significant_BH_q)
  rbm47 <- filter(statistics, .data[["Protein accession"]] == "A0AV96")
  checks <- tibble(
    check = c(
      "source_protein_row_count",
      "tested_protein_count",
      "raw_significant_count",
      "BH_significant_count",
      "raw_filter_integrity",
      "BH_filter_integrity",
      "raw_old_status_reconciliation",
      "raw_new_status_reconciliation",
      "BH_old_status_reconciliation",
      "BH_new_status_reconciliation",
      "vendor_t_test_spot_check_A0AV96",
      "enrichment_probability_bounds"
    ),
    passed = c(
      nrow(statistics) == 7347,
      sum(!is.na(statistics$p_value_excluding_polysome_3)) == summary$tested_proteins,
      nrow(raw_subset) == summary$new_raw_p_total,
      nrow(bh_subset) == summary$new_bh_q_total,
      all(
        abs(raw_subset$log2FC_excluding_polysome_3) > log2fc_threshold &
          raw_subset[["Unique peptides"]] >= unique_peptide_threshold &
          raw_subset$p_value_excluding_polysome_3 < p_threshold
      ),
      all(
        abs(bh_subset$log2FC_excluding_polysome_3) > log2fc_threshold &
          bh_subset[["Unique peptides"]] >= unique_peptide_threshold &
          bh_subset$q_value_excluding_polysome_3 < q_threshold
      ),
      summary$old_raw_p_total_same_filters == summary$raw_retained + summary$raw_lost,
      summary$new_raw_p_total == summary$raw_retained + summary$raw_gained,
      summary$old_bh_q_total_same_filters == summary$bh_retained + summary$bh_lost,
      summary$new_bh_q_total == summary$bh_retained + summary$bh_gained,
      nrow(rbm47) == 1 &&
        abs(
          rbm47$p_value_excluding_polysome_3 -
            rbm47$old_p_value_including_polysome_3
        ) < 1e-12,
      all(
        enrichment$p_value >= 0 & enrichment$p_value <= 1 &
          enrichment$q_value >= 0 & enrichment$q_value <= 1,
        na.rm = TRUE
      )
    )
  )
  write.csv(checks, file.path(output_dir, "validation_checks.csv"), row.names = FALSE)
  if (!all(checks$passed)) {
    stop(
      "Validation failed: ",
      paste(checks$check[!checks$passed], collapse = ", ")
    )
  }
  invisible(checks)
}

message("Reading annotation workbook...")
annotation <- read_excel(source_workbook, sheet = "Annotation_Combine") |>
  select(-starts_with("..."))
message("Calculating differential protein statistics...")
statistics <- calculate_statistics(annotation)
message("Calculating functional enrichment...")
enrichment <- calculate_enrichment(statistics)
summary <- build_summary(statistics, enrichment)
statistics <- write_outputs(statistics, enrichment, summary)
make_plots(statistics, enrichment)
validate_results(statistics, enrichment, summary)

report_template <- file.path(scripts_dir, "reanalysis_report.Rmd")
if (file.exists(report_template)) {
  message("Rendering the reproducible R Markdown report...")
  rmarkdown::render(
    input = report_template,
    output_file = "reanalysis_report.html",
    output_dir = output_dir,
    knit_root_dir = output_dir,
    quiet = TRUE,
    envir = new.env(parent = globalenv())
  )
}

cat(toJSON(summary, pretty = TRUE, auto_unbox = TRUE), "\n")
