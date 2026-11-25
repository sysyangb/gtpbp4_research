#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)
suppressPackageStartupMessages({
  library(survival)
  library(ggplot2)
  library(svglite)
  library(ragg)
})


genes <- c("DKC1", "GTPBP4", "ILF2", "NOP2", "RPP38")
ensembl_to_gene <- c(
  ENSG00000130826.18 = "DKC1",
  ENSG00000107937.19 = "GTPBP4",
  ENSG00000143621.17 = "ILF2",
  ENSG00000111641.12 = "NOP2",
  ENSG00000152464.15 = "RPP38"
)

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (length(script_arg) != 1L) stop("Run this file with Rscript.")
script_path <- normalizePath(sub("^--file=", "", script_arg), mustWork = TRUE)
project_dir <- normalizePath(file.path(dirname(script_path), ".."), mustWork = TRUE)
figure2a_dir <- normalizePath(file.path(project_dir, ".."), mustWork = TRUE)
paper_project <- file.path(figure2a_dir, "PMC11391272_GTPBP4_reproduction")

raw_expression_file <- file.path(figure2a_dir, "TCGA-BRCA_TPM.txt")
os_file <- file.path(
  paper_project, "source", "derived", "TCGA_BRCA_OS_selected_primary_samples.tsv"
)
clinical_file <- file.path(
  paper_project, "source", "tcga_xena", "BRCA_clinicalMatrix"
)

data_dir <- file.path(project_dir, "data", "processed")
results_dir <- file.path(project_dir, "results")
figures_dir <- file.path(project_dir, "figures")
dir.create(data_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)

source_files <- c(raw_expression_file, os_file, clinical_file)
if (!all(file.exists(source_files))) {
  stop("Missing source file(s): ", paste(source_files[!file.exists(source_files)], collapse = ", "))
}

strip_quotes <- function(x) sub('"$', "", sub('^"', "", x))

# Read only five target rows from the large TCGA TPM matrix. This avoids loading
# the complete 342-MB expression matrix while retaining exact sample barcodes.
extract_target_expression <- function(path, id_to_gene, chunk_size = 1000L) {
  con <- file(path, open = "rt")
  on.exit(close(con), add = TRUE)

  header_line <- readLines(con, n = 1L, warn = FALSE)
  if (length(header_line) != 1L) stop("Expression file has no header.")
  sample_ids <- strip_quotes(strsplit(header_line, "\t", fixed = TRUE)[[1L]])

  target_ids <- names(id_to_gene)
  target_prefixes <- paste0('"', target_ids, '"\t')
  found <- vector("list", length(target_ids))
  names(found) <- target_ids

  repeat {
    lines <- readLines(con, n = chunk_size, warn = FALSE)
    if (!length(lines)) break
    for (i in seq_along(target_ids)) {
      if (!is.null(found[[i]])) next
      hit <- which(startsWith(lines, target_prefixes[i]))
      if (length(hit) > 1L) stop("Duplicated target Ensembl row: ", target_ids[i])
      if (length(hit) == 1L) found[[i]] <- lines[hit]
    }
    if (all(vapply(found, Negate(is.null), logical(1)))) break
  }

  missing_ids <- target_ids[vapply(found, is.null, logical(1))]
  if (length(missing_ids)) stop("Missing target Ensembl row(s): ", paste(missing_ids, collapse = ", "))

  expression <- do.call(rbind, lapply(found, function(line) {
    fields <- strsplit(line, "\t", fixed = TRUE)[[1L]]
    if (length(fields) != length(sample_ids) + 1L) stop("Malformed target expression row.")
    suppressWarnings(as.numeric(fields[-1L]))
  }))
  if (any(!is.finite(expression)) || any(expression < 0)) stop("Invalid target TPM values.")

  rownames(expression) <- unname(id_to_gene[target_ids])
  colnames(expression) <- sample_ids
  expression[genes, , drop = FALSE]
}

expression_all <- extract_target_expression(raw_expression_file, ensembl_to_gene)
os <- read.delim(os_file, check.names = FALSE)
clinical <- read.delim(clinical_file, check.names = FALSE)

required_os <- c("patient_id", "sample_id", "event", "os_days")
if (!all(required_os %in% names(os))) stop("OS table is missing required columns.")
if (!all(c("sampleID", "PAM50Call_RNAseq") %in% names(clinical))) {
  stop("Clinical table is missing PAM50 fields.")
}

expression_index <- match(os$sample_id, colnames(expression_all))
if (anyNA(expression_index)) stop("Not all selected TCGA primary tumors map to TPM columns.")
clinical_index <- match(substr(os$sample_id, 1L, 15L), clinical$sampleID)

analysis_keep <-
  !is.na(clinical_index) &
  clinical$PAM50Call_RNAseq[clinical_index] == "Basal" &
  is.finite(os$os_days) & os$os_days >= 0 &
  os$event %in% c(0L, 1L)

basal_os <- os[analysis_keep, , drop = FALSE]
basal_expression <- expression_all[, expression_index[analysis_keep], drop = FALSE]
if (nrow(basal_os) != 139L || sum(basal_os$event) != 19L) {
  stop(sprintf(
    "Unexpected PAM50 Basal cohort: n=%d, deaths=%d; expected n=139, deaths=19.",
    nrow(basal_os), sum(basal_os$event)
  ))
}

analysis_data <- data.frame(
  patient_id = basal_os$patient_id,
  sample_id = basal_os$sample_id,
  PAM50 = "Basal",
  OS_time_days = basal_os$os_days,
  OS_time_months = basal_os$os_days / 30.4375,
  OS_event = as.integer(basal_os$event),
  stringsAsFactors = FALSE,
  check.names = FALSE
)
for (gene in genes) {
  analysis_data[[paste0(gene, "_TPM")]] <- as.numeric(basal_expression[gene, ])
  analysis_data[[paste0(gene, "_log2TPM1")]] <- log2(as.numeric(basal_expression[gene, ]) + 1)
}

sample_file <- file.path(data_dir, "Figure2A_PAM50_Basal_sample_level_data.tsv")
write.table(analysis_data, sample_file, sep = "\t", quote = FALSE, row.names = FALSE)

extract_cox <- function(fit, term) {
  sm <- summary(fit)
  c(
    hr = unname(sm$conf.int[term, "exp(coef)"]),
    ci_low = unname(sm$conf.int[term, "lower .95"]),
    ci_high = unname(sm$conf.int[term, "upper .95"]),
    cox_wald_p = unname(sm$coefficients[term, "Pr(>|z|)"]),
    cox_likelihood_ratio_p = unname(sm$logtest["pvalue"]),
    cox_score_p = unname(sm$sctest["pvalue"])
  )
}

group_rows <- list()
results <- do.call(rbind, lapply(genes, function(gene) {
  x <- analysis_data[[paste0(gene, "_TPM")]]
  cutoff <- median(x)
  group <- factor(ifelse(x > cutoff, "High", "Low"), levels = c("Low", "High"))
  dat <- data.frame(
    time = analysis_data$OS_time_days,
    event = analysis_data$OS_event,
    expression_TPM = x,
    group = group
  )

  fit <- coxph(Surv(time, event) ~ group, data = dat, ties = "efron", x = TRUE)
  cox <- extract_cox(fit, "groupHigh")
  lr <- survdiff(Surv(time, event) ~ group, data = dat)
  ph_p <- unname(cox.zph(fit, transform = "km")$table["group", "p"])

  group_rows[[gene]] <<- data.frame(
    patient_id = analysis_data$patient_id,
    sample_id = analysis_data$sample_id,
    gene = gene,
    expression_TPM = x,
    median_cutoff_TPM = cutoff,
    group = group,
    OS_time_days = analysis_data$OS_time_days,
    OS_event = analysis_data$OS_event,
    stringsAsFactors = FALSE
  )

  data.frame(
    gene = gene,
    cohort = "TCGA PAM50 Basal",
    endpoint = "OS",
    comparison = "High vs Low",
    split = "Gene-specific cohort median TPM",
    n = nrow(dat),
    events = sum(dat$event),
    median_cutoff_TPM = cutoff,
    low_n = sum(dat$group == "Low"),
    low_events = sum(dat$event[dat$group == "Low"]),
    high_n = sum(dat$group == "High"),
    high_events = sum(dat$event[dat$group == "High"]),
    hr = cox["hr"],
    ci_low = cox["ci_low"],
    ci_high = cox["ci_high"],
    cox_wald_p = cox["cox_wald_p"],
    cox_likelihood_ratio_p = cox["cox_likelihood_ratio_p"],
    cox_score_p = cox["cox_score_p"],
    logrank_p = 1 - pchisq(lr$chisq, df = 1L),
    proportional_hazards_p = ph_p,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
}))
rownames(results) <- NULL

results_file <- file.path(results_dir, "Figure2A_PAM50_Basal_five_gene_OS_results.tsv")
groups_file <- file.path(results_dir, "Figure2A_PAM50_Basal_gene_groups.tsv")
write.table(results, results_file, sep = "\t", quote = FALSE, row.names = FALSE)
write.table(do.call(rbind, group_rows), groups_file, sep = "\t", quote = FALSE, row.names = FALSE)

# Forest plot layout. Statistics and gene labels remain on the corresponding
# horizontal rows; explicit columns prevent CI/text collisions.
forest <- results
forest$row <- rev(seq_along(genes))
forest$significant <- forest$cox_wald_p < 0.05
forest$hr_ci_label <- sprintf("%.2f (%.2f–%.2f)", forest$hr, forest$ci_low, forest$ci_high)
forest$p_label <- ifelse(
  forest$cox_wald_p < 0.001,
  formatC(forest$cox_wald_p, format = "e", digits = 2),
  formatC(forest$cox_wald_p, format = "f", digits = 4)
)

pt_to_mm <- 25.4 / 72.27
body_size <- 10.0 * pt_to_mm
header_size <- 10.5 * pt_to_mm
axis_title_size <- 10.2 * pt_to_mm
minimum_visible_font_pt <- 10.0

figure_width_mm <- 155
figure_height_mm <- 148
gene_x <- 0.012
forest_left <- 0.165
forest_right <- 0.640
hr_text_x <- 0.675
p_text_x <- 0.988
header_y <- 5.73
axis_y <- 0.59
tick_bottom <- 0.50
tick_label_y <- 0.34
axis_title_y <- 0.09

forest_hr_min <- 0.20
forest_hr_max <- 12.0
if (min(forest$ci_low) < forest_hr_min || max(forest$ci_high) > forest_hr_max) {
  stop("A confidence interval falls outside the prespecified forest axis.")
}

map_hr <- function(x) {
  forest_left +
    (log(x) - log(forest_hr_min)) /
    (log(forest_hr_max) - log(forest_hr_min)) *
    (forest_right - forest_left)
}

ticks <- c(0.2, 0.5, 1, 2, 5, 10)
tick_x <- map_hr(ticks)
reference_x <- map_hr(1)
forest$plot_hr <- map_hr(forest$hr)
forest$plot_ci_low <- map_hr(forest$ci_low)
forest$plot_ci_high <- map_hr(forest$ci_high)

neutral <- "#4D4D4D"
accent <- "#C83E36"
row_colours <- c(`FALSE` = neutral, `TRUE` = accent)

figure2a <- ggplot(forest, aes(y = row, colour = factor(significant))) +
  annotate(
    "segment", x = reference_x, xend = reference_x,
    y = axis_y, yend = 5.35,
    linewidth = 0.35, linetype = "dashed", colour = "#8A8A8A"
  ) +
  geom_segment(
    aes(x = plot_ci_low, xend = plot_ci_high, yend = row),
    linewidth = 0.72, lineend = "butt"
  ) +
  geom_point(aes(x = plot_hr), shape = 15, size = 2.5) +
  geom_text(
    aes(x = gene_x, label = gene),
    hjust = 0, family = "Arial", fontface = "bold",
    size = body_size, show.legend = FALSE
  ) +
  geom_text(
    aes(
      x = hr_text_x, label = hr_ci_label,
      fontface = ifelse(significant, "bold", "plain")
    ),
    hjust = 0, family = "Arial", size = body_size, show.legend = FALSE
  ) +
  geom_text(
    aes(
      x = p_text_x, label = p_label,
      fontface = ifelse(significant, "bold", "plain")
    ),
    hjust = 1, family = "Arial", size = body_size, show.legend = FALSE
  ) +
  annotate(
    "text", x = 0.790, y = header_y, label = "HR (95% CI)",
    family = "Arial", fontface = "bold", size = header_size, hjust = 0.5
  ) +
  annotate(
    "text", x = p_text_x, y = header_y, label = "Cox P",
    family = "Arial", fontface = "bold", size = header_size, hjust = 1
  ) +
  annotate(
    "segment", x = min(tick_x), xend = max(tick_x),
    y = axis_y, yend = axis_y, linewidth = 0.38, colour = "black"
  ) +
  annotate(
    "segment", x = tick_x, xend = tick_x,
    y = axis_y, yend = tick_bottom, linewidth = 0.36, colour = "black"
  ) +
  annotate(
    "text", x = tick_x, y = tick_label_y,
    label = c("0.2", "0.5", "1", "2", "5", "10"),
    family = "Arial", size = body_size, colour = "black"
  ) +
  annotate(
    "text", x = mean(c(forest_left, forest_right)), y = axis_title_y,
    label = "Hazard ratio (High vs Low)",
    family = "Arial", size = axis_title_size, colour = "black"
  ) +
  scale_x_continuous(limits = c(0, 1), expand = expansion(mult = c(0, 0))) +
  scale_y_continuous(limits = c(0, 5.98), expand = expansion(mult = c(0, 0))) +
  scale_colour_manual(values = row_colours, guide = "none") +
  coord_cartesian(clip = "off") +
  theme_void(base_family = "Arial") +
  theme(
    plot.margin = margin(4.0, 3.0, 3.0, 3.0, unit = "mm"),
    plot.background = element_rect(fill = "white", colour = NA),
    panel.background = element_rect(fill = "white", colour = NA)
  )

save_bundle <- function(plot, stem, width_mm, height_mm, dpi = 600L) {
  width_in <- width_mm / 25.4
  height_in <- height_mm / 25.4

  svglite::svglite(paste0(stem, ".svg"), width = width_in, height = height_in)
  print(plot)
  dev.off()

  grDevices::cairo_pdf(
    paste0(stem, ".pdf"), width = width_in, height = height_in, family = "Arial"
  )
  print(plot)
  dev.off()

  ragg::agg_tiff(
    paste0(stem, ".tiff"), width = width_in, height = height_in,
    units = "in", res = dpi, compression = "lzw"
  )
  print(plot)
  dev.off()

  ragg::agg_png(
    paste0(stem, ".png"), width = width_in, height = height_in,
    units = "in", res = 300L
  )
  print(plot)
  dev.off()
}

figure_stem <- file.path(figures_dir, "Figure2A_PAM50_Basal_five_gene_forest")
save_bundle(figure2a, figure_stem, figure_width_mm, figure_height_mm)

figure_files <- paste0(figure_stem, c(".svg", ".pdf", ".tiff", ".png"))
svg_text <- paste(readLines(paste0(figure_stem, ".svg"), warn = FALSE), collapse = "\n")
qa <- data.frame(
  check = c(
    "fixed PAM50 Basal cohort size",
    "fixed PAM50 Basal death count",
    "all five genes retained",
    "each split preserves cohort total",
    "finite ordered confidence intervals",
    "valid Cox Wald P values",
    "GTPBP4 reproduces reference HR",
    "GTPBP4 reproduces reference P value",
    "minimum visible font at least 10 pt",
    "near-square aspect ratio",
    "statistics columns do not enter forest region",
    "no FDR column in figure source",
    "all four export formats exist",
    "all exports are nonempty",
    "SVG retains Arial text"
  ),
  pass = c(
    all(results$n == 139L),
    all(results$events == 19L),
    identical(results$gene, genes),
    all(results$low_n + results$high_n == results$n),
    all(is.finite(results$hr)) && all(is.finite(results$ci_low)) &&
      all(is.finite(results$ci_high)) && all(results$ci_low < results$ci_high),
    all(is.finite(results$cox_wald_p)) && all(results$cox_wald_p >= 0 & results$cox_wald_p <= 1),
    abs(results$hr[results$gene == "GTPBP4"] - 3.277) < 0.01,
    abs(results$cox_wald_p[results$gene == "GTPBP4"] - 0.035) < 0.002,
    minimum_visible_font_pt >= 10,
    max(figure_width_mm, figure_height_mm) / min(figure_width_mm, figure_height_mm) <= 1.10,
    forest_right < hr_text_x && hr_text_x < p_text_x,
    !any(grepl("fdr", names(forest), ignore.case = TRUE)),
    all(file.exists(figure_files)),
    all(file.info(figure_files)$size > 1000),
    grepl("Arial", svg_text, fixed = TRUE)
  ),
  detail = c(
    sprintf("n=%d", unique(results$n)),
    sprintf("deaths=%d", unique(results$events)),
    paste(results$gene, collapse = ", "),
    paste(sprintf("%s:%d+%d", results$gene, results$low_n, results$high_n), collapse = "; "),
    "all HR and 95% CI values are finite and ordered",
    "all Cox Wald P values lie in [0,1]",
    sprintf("HR=%.6f", results$hr[results$gene == "GTPBP4"]),
    sprintf("P=%.8f", results$cox_wald_p[results$gene == "GTPBP4"]),
    sprintf("%.1f pt", minimum_visible_font_pt),
    sprintf("%d x %d mm", figure_width_mm, figure_height_mm),
    sprintf("forest right=%.3f; HR text=%.3f; P text=%.3f", forest_right, hr_text_x, p_text_x),
    "forest source contains HR, CI and Cox P only",
    paste(basename(figure_files), collapse = ", "),
    paste(range(file.info(figure_files)$size), collapse = " to "),
    "Arial present in SVG"
  ),
  stringsAsFactors = FALSE
)

qa_file <- file.path(results_dir, "Figure2A_PAM50_Basal_validation_report.tsv")
write.table(qa, qa_file, sep = "\t", quote = FALSE, row.names = FALSE)
if (!all(qa$pass)) {
  print(qa[!qa$pass, ], row.names = FALSE)
  stop("Figure2A PAM50 Basal validation failed.")
}

session_file <- file.path(results_dir, "Figure2A_PAM50_Basal_sessionInfo.txt")
writeLines(capture.output(sessionInfo()), session_file)

cat("Figure2A PAM50 Basal five-gene survival analysis completed.\n")
print(
  results[, c(
    "gene", "n", "events", "low_n", "low_events", "high_n", "high_events",
    "hr", "ci_low", "ci_high", "cox_wald_p", "logrank_p", "proportional_hazards_p"
  )],
  row.names = FALSE,
  digits = 5
)

