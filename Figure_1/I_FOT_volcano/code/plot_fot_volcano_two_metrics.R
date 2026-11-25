#!/usr/bin/env Rscript

# The five labels come from the polysome proteomics FOT comparison:
# Polysome fraction versus Free fraction, not polysome-seq RNA counts.

args <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (length(args) != 1L) stop("Run this script with Rscript")
script_file <- normalizePath(sub("^--file=", "", args), mustWork = TRUE)
task_dir <- normalizePath(file.path(dirname(script_file), ".."), mustWork = TRUE)
polysome_dir <- normalizePath(file.path(task_dir, ".."), mustWork = TRUE)

required <- c("ggplot2", "svglite", "ragg")
missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) stop("Missing R packages: ", paste(missing, collapse = ", "))

source_file <- file.path(polysome_dir, "FOT_pipeline", "results", "tables",
                         "FOT_limma_all_tested.csv")
up_file <- file.path(polysome_dir, "FOT_pipeline", "results", "tables",
                     "FOT_TREAT_BH05_absFCgt1_UP.csv")
results_dir <- file.path(task_dir, "results")
figures_dir <- file.path(task_dir, "figures")
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)

dat <- read.csv(source_file, check.names = FALSE)
strict_up <- read.csv(up_file, check.names = FALSE)
stopifnot(nrow(dat) == 2674L, nrow(strict_up) == 218L,
          all(is.finite(dat$log2FC_Polysome_vs_Free)),
          all(is.finite(dat$TREAT_P_Value)),
          all(is.finite(dat$TREAT_BH_FDR)),
          isTRUE(all.equal(dat$TREAT_BH_FDR,
                           p.adjust(dat$TREAT_P_Value, method = "BH"),
                           check.attributes = FALSE, tolerance = 1e-12)))

targets <- c("PABPC1", "UPF1", "LARP1", "FMR1", "EIF4A3")
marked <- dat[match(targets, dat$Gene_name), , drop = FALSE]
stopifnot(identical(marked$Gene_name, targets),
          all(marked$Protein_accession %in% strict_up$Protein_accession),
          all(marked$TREAT_call == "Up"))
marked_output <- marked[, c("Gene_name", "Protein_accession", "Unique_peptides",
                            "log2FC_Polysome_vs_Free", "TREAT_P_Value",
                            "TREAT_BH_FDR", "TREAT_call")]
write.table(marked_output, file.path(results_dir, "five_labeled_proteins.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
source_data <- dat[, c("Gene_name", "Protein_accession", "Unique_peptides",
                       "log2FC_Polysome_vs_Free", "TREAT_P_Value",
                       "TREAT_BH_FDR")]
source_data$FDR_call <- ifelse(
  dat$TREAT_BH_FDR < 0.05 & dat$log2FC_Polysome_vs_Free > 1, "Up",
  ifelse(dat$TREAT_BH_FDR < 0.05 & dat$log2FC_Polysome_vs_Free < -1,
         "Down", "Not significant")
)
source_data$raw_p_call <- ifelse(
  dat$TREAT_P_Value < 0.05 & dat$log2FC_Polysome_vs_Free > 1, "Up",
  ifelse(dat$TREAT_P_Value < 0.05 & dat$log2FC_Polysome_vs_Free < -1,
         "Down", "Not significant")
)
write.table(source_data, file.path(results_dir, "volcano_source_data.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

cols <- c("Down" = "#4572A7", "Not significant" = "#BFC2C5",
          "Up" = "#AA4643")

make_plot <- function(metric = c("fdr", "pvalue")) {
  metric <- match.arg(metric)
  value_col <- if (metric == "fdr") "TREAT_BH_FDR" else "TREAT_P_Value"
  plot_data <- dat
  plot_data$y <- -log10(pmax(plot_data[[value_col]], .Machine$double.xmin))
  plot_data$call <- "Not significant"
  plot_data$call[plot_data[[value_col]] < 0.05 &
                   plot_data$log2FC_Polysome_vs_Free > 1] <- "Up"
  plot_data$call[plot_data[[value_col]] < 0.05 &
                   plot_data$log2FC_Polysome_vs_Free < -1] <- "Down"
  plot_data$call <- factor(plot_data$call,
                           levels = c("Down", "Not significant", "Up"))
  marked_plot <- plot_data[match(targets, plot_data$Gene_name), , drop = FALSE]
  # Put gene names in regions with no data points. The leader lines are
  # rendered underneath points so no data mark is hidden by annotation.
  marked_plot$label_x <- c(1.6, 0.25, 12.7, 12.7, 12.7)
  marked_plot$label_y <- c(max(plot_data$y) + 0.8,
                           marked_plot$y[2] + 0.55,
                           marked_plot$y[3],
                           marked_plot$y[4],
                           marked_plot$y[5])
  marked_plot$label_hjust <- c(0.5, 0.5, 0, 0, 0)

  y_lab <- if (metric == "fdr") expression(-log[10]("TREAT BH-FDR"))
           else expression(-log[10]("TREAT p-value"))
  p <- ggplot2::ggplot(plot_data,
                       ggplot2::aes(log2FC_Polysome_vs_Free, y)) +
    ggplot2::geom_segment(
      data = marked_plot,
      ggplot2::aes(x = log2FC_Polysome_vs_Free, y = y,
                   xend = label_x, yend = label_y),
      inherit.aes = FALSE, linewidth = 0.3, colour = "#666666",
      show.legend = FALSE
    ) +
    ggplot2::geom_point(
      data = plot_data[plot_data$call == "Not significant", , drop = FALSE],
      colour = "#BFC2C5", size = 1.0, alpha = 0.55
    ) +
    ggplot2::geom_point(
      data = plot_data[plot_data$call != "Not significant", , drop = FALSE],
      ggplot2::aes(colour = call), size = 1.2, alpha = 1
    ) +
    ggplot2::geom_hline(yintercept = -log10(0.05), linetype = "dashed",
                        linewidth = 0.4, colour = "#666666") +
    ggplot2::geom_vline(xintercept = c(-1, 1), linetype = "dashed",
                        linewidth = 0.4, colour = "#777777") +
    ggplot2::geom_vline(xintercept = 0, linewidth = 0.3,
                        colour = "#BBBBBB") +
    ggplot2::geom_point(
      data = marked_plot, ggplot2::aes(fill = call),
      shape = 21, size = 3.0, stroke = 0.8, colour = "#222222",
      show.legend = FALSE
    ) +
    ggplot2::geom_text(
      data = marked_plot,
      ggplot2::aes(x = label_x, y = label_y, label = Gene_name,
                   hjust = label_hjust),
      inherit.aes = FALSE, size = 12 / ggplot2::.pt,
      family = "Arial", fontface = "bold", colour = "#222222",
      show.legend = FALSE
    ) +
    ggplot2::scale_x_continuous(limits = c(-7.5, 16.5),
                                breaks = c(-4, 0, 4, 8, 12)) +
    ggplot2::scale_y_continuous(limits = c(0, max(plot_data$y) + 1.5)) +
    ggplot2::scale_colour_manual(values = cols, breaks = c("Up", "Down"),
                                 name = NULL, drop = FALSE) +
    ggplot2::scale_fill_manual(values = cols, guide = "none") +
    ggplot2::labs(x = "log2FC (Polysome / Free)", y = y_lab) +
    ggplot2::theme_classic(base_size = 12, base_family = "Arial") +
    ggplot2::theme(
      plot.background = ggplot2::element_blank(),
      panel.background = ggplot2::element_blank(),
      axis.line = ggplot2::element_line(linewidth = 0.5),
      axis.ticks = ggplot2::element_line(linewidth = 0.4),
      legend.position = "top", legend.direction = "horizontal",
      legend.justification = "left",
      legend.background = ggplot2::element_blank(),
      legend.key = ggplot2::element_blank(),
      legend.box.background = ggplot2::element_blank(),
      legend.margin = ggplot2::margin(0, 0, 0, 0),
      plot.margin = ggplot2::margin(10, 10, 10, 10)
    )
  list(plot = p, counts = table(plot_data$call), value_col = value_col)
}

save_figure <- function(p, suffix) {
  base <- file.path(figures_dir, paste0("FOT_volcano_five_genes_", suffix))
  svglite::svglite(paste0(base, ".svg"), width = 5.5, height = 5.5,
                   bg = "transparent")
  print(p)
  grDevices::dev.off()
  grDevices::cairo_pdf(paste0(base, ".pdf"), width = 5.5,
                       height = 5.5, family = "Arial", bg = "transparent")
  print(p)
  grDevices::dev.off()
  ragg::agg_png(paste0(base, ".png"), width = 5.5,
                height = 5.5, units = "in", res = 300,
                background = "transparent")
  print(p)
  grDevices::dev.off()
  ragg::agg_png(paste0(base, "_preview_light_grey.png"), width = 5.5,
                height = 5.5, units = "in", res = 300,
                background = "#F2F2F0")
  print(p)
  grDevices::dev.off()
  ragg::agg_tiff(paste0(base, ".tiff"), width = 5.5,
                 height = 5.5, units = "in", res = 600,
                 background = "transparent")
  print(p)
  grDevices::dev.off()
}

fdr_plot <- make_plot("fdr")
pvalue_plot <- make_plot("pvalue")
save_figure(fdr_plot$plot, "TREAT_BH_FDR")
save_figure(pvalue_plot$plot, "TREAT_raw_p")

call_counts <- data.frame(
  metric = rep(c("TREAT_BH_FDR", "TREAT_raw_p"), each = 3),
  direction = rep(c("Down", "Not significant", "Up"), times = 2),
  n = c(as.integer(fdr_plot$counts), as.integer(pvalue_plot$counts))
)
write.table(call_counts, file.path(results_dir, "volcano_call_counts.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
writeLines(c(
  paste("Source MD5:", unname(tools::md5sum(source_file))),
  paste("Strict-up MD5:", unname(tools::md5sum(up_file))),
  capture.output(sessionInfo())
), file.path(results_dir, "R_sessionInfo_and_source_md5.txt"))

print(marked_output, row.names = FALSE)
print(call_counts, row.names = FALSE)
