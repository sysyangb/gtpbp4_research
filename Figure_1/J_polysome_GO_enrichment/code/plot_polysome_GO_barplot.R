#!/usr/bin/env Rscript

# Redraw the existing ten-term enrichment panel without recomputing statistics.
# Length = -log10(raw P); fill = fold enrichment; no gene-count encoding.
suppressPackageStartupMessages({
  library(ggplot2)
  library(ragg)
  library(svglite)
})

script_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
if (length(script_arg) != 1L) stop("Run this script using Rscript.")
analysis_dir <- dirname(dirname(normalizePath(sub("^--file=", "", script_arg))))
source_path <- file.path(analysis_dir, "results", "Polysome_up_GO_BP_plus_AMFs_total_human_Top10_source_data.tsv")
figures_dir <- file.path(analysis_dir, "figures")
results_dir <- file.path(analysis_dir, "results")
qc_dir <- file.path(analysis_dir, "QC")
source <- read.delim(source_path, check.names = FALSE, stringsAsFactors = FALSE)
source <- source[order(source$top10_rank_raw_P), ]

expected_terms <- c(
  "cytoplasmic translation", "ribosome biogenesis",
  "ribonucleoprotein complex biogenesis",
  "Assembly and maturation factors (AMFs)", "translation",
  "protein biosynthetic process", "regulation of mRNA catabolic process",
  "regulation of mRNA metabolic process",
  "post-transcriptional regulation of gene expression", "rRNA processing"
)
stopifnot(
  nrow(source) == 10L,
  identical(source$term_name, expected_terms),
  all(is.finite(source$pvalue) & source$pvalue > 0),
  max(abs(source$minus_log10_pvalue + log10(source$pvalue))) < 1e-10,
  max(abs(source$fold_enrichment - source$gene_ratio / source$background_ratio)) < 1e-10
)

plot_data <- source[, c("top10_rank_raw_P", "term_id", "term_name", "source_type",
                         "pvalue", "p_adjust_BH", "minus_log10_pvalue", "fold_enrichment")]
names(plot_data)[1] <- "display_rank"
plot_data$display_name <- vapply(plot_data$term_name, function(x) {
  if (grepl("^(rRNA|mRNA)\\b", x)) return(x)
  paste0(toupper(substr(x, 1, 1)), substr(x, 2, nchar(x)))
}, character(1))
plot_data$display_label <- vapply(plot_data$display_name, function(x) {
  paste(strwrap(x, width = 36), collapse = "\n")
}, character(1))
plot_data$pathway <- factor(plot_data$display_label, levels = rev(plot_data$display_label))
write.table(plot_data[, setdiff(names(plot_data), "pathway")],
            file.path(results_dir, "FigureS1-polysome-up-GO_barplot_source_data.tsv"),
            sep = "\t", row.names = FALSE, quote = TRUE)

p <- ggplot(plot_data, aes(x = minus_log10_pvalue, y = pathway, fill = fold_enrichment)) +
  geom_col(width = 0.64, orientation = "y", colour = NA) +
  scale_x_continuous(
    name = "\u2212log\u2081\u2080(P value)",
    breaks = seq(0, 40, 10), limits = c(0, 40.5),
    expand = expansion(mult = c(0, 0))
  ) +
  scale_y_discrete(name = NULL, expand = expansion(add = c(0.6, 0.6))) +
  scale_fill_gradient(
    low = "#F8E8E5", high = "#B34D49", name = "Enrichment\nscore",
    limits = range(plot_data$fold_enrichment), breaks = 5:10,
    guide = guide_colourbar(
      title.position = "top", title.hjust = 0,
      barheight = grid::unit(44, "mm"), barwidth = grid::unit(6, "mm"),
      ticks = TRUE, frame.colour = NA
    )
  ) +
  theme_classic(base_size = 14, base_family = "Arial") +
  theme(
    text = element_text(size = 14),
    axis.line = element_line(linewidth = 0.45, colour = "#222222"),
    axis.ticks.x = element_line(linewidth = 0.45, colour = "#222222"),
    axis.ticks.y = element_blank(),
    axis.ticks.length = grid::unit(1.5, "mm"),
    axis.text.x = element_text(size = 14, colour = "#222222"),
    axis.text.y = element_text(size = 14, colour = "#222222", lineheight = 0.95,
                               margin = margin(r = 8)),
    axis.title.x = element_text(size = 15, margin = margin(t = 9)),
    legend.position = "right", legend.justification = "center",
    legend.title = element_text(size = 14.5, face = "plain", lineheight = 1.0,
                               margin = margin(b = 5)),
    legend.text = element_text(size = 14),
    legend.margin = margin(l = 3, r = 0),
    legend.box.spacing = grid::unit(3, "mm"),
    plot.margin = margin(t = 8, r = 5, b = 7, l = 5, unit = "mm"),
    panel.grid = element_blank()
  )

width_mm <- 250
height_mm <- 160
stem <- file.path(figures_dir, "FigureS1-polysome-up-GO_barplot")
ggsave(paste0(stem, ".svg"), p, device = svglite::svglite,
       width = width_mm, height = height_mm, units = "mm", bg = "white")
ggsave(paste0(stem, ".pdf"), p, device = grDevices::cairo_pdf, family = "Arial",
       width = width_mm, height = height_mm, units = "mm", bg = "white")
ggsave(paste0(stem, ".png"), p, device = ragg::agg_png,
       width = width_mm, height = height_mm, units = "mm", dpi = 600, bg = "white")
ggsave(paste0(stem, ".tiff"), p, device = ragg::agg_tiff, compression = "lzw",
       width = width_mm, height = height_mm, units = "mm", dpi = 600, bg = "white")
ggsave(paste0(stem, "_preview.png"), p, device = ragg::agg_png,
       width = width_mm, height = height_mm, units = "mm", dpi = 180, bg = "white")

# Check the actual plotted layer and vector output, not just the input table.
qa_capture <- ragg::agg_capture(width = width_mm, height = height_mm,
                              units = "mm", res = 180)
built <- ggplot_build(p)
bars <- built$data[[1]]
svg <- paste(readLines(paste0(stem, ".svg"), warn = FALSE), collapse = "\n")
collect_font_sizes <- function(grob) {
  output <- numeric()
  if (inherits(grob, "text") && !is.null(grob$gp$fontsize)) {
    output <- c(output, as.numeric(grob$gp$fontsize))
  }
  if (!is.null(grob$children)) {
    output <- c(output, unlist(lapply(grob$children, collect_font_sizes), use.names = FALSE))
  }
  if (!is.null(grob$grobs)) {
    output <- c(output, unlist(lapply(grob$grobs, collect_font_sizes), use.names = FALSE))
  }
  output
}
font_sizes <- collect_font_sizes(ggplotGrob(p))
grDevices::dev.off()
write.table(data.frame(font_size_pt = sort(unique(font_sizes))),
            file.path(qc_dir, "barplot_font_sizes.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)
checks <- c(
  all_text_at_least_14pt = length(font_sizes) > 0L && all(font_sizes >= 14),
  ten_bars = nrow(bars) == 10L,
  zero_baseline = all(bars$xmin == 0),
  exact_bar_lengths = max(abs(bars$xmax - plot_data$minus_log10_pvalue)) < 1e-10,
  exact_fold_enrichment_fill = identical(toupper(bars$fill), toupper(built$plot$scales$get_scales("fill")$map(plot_data$fold_enrichment))),
  original_term_order = identical(as.character(plot_data$term_name), expected_terms),
  no_count_aesthetic = !any(c("size", "linewidth") %in% names(p$mapping)),
  no_count_legend = !grepl("Gene count", svg, fixed = TRUE),
  editable_svg_text = grepl("<text", svg, fixed = TRUE),
  all_terms_in_svg = all(vapply(strsplit(paste(plot_data$display_label, collapse = "\n"), "\n")[[1]],
                                function(x) grepl(x, svg, fixed = TRUE), logical(1)))
)
write.table(data.frame(check = names(checks), pass = unname(checks)),
            file.path(qc_dir, "barplot_validation.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)
if (!all(checks)) stop("A figure validation check failed.")
message("Exported bar plot; all ", length(checks), " validation checks passed.")
