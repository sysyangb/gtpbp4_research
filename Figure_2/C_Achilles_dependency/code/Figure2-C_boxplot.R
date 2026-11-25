#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
  library(ragg)
  library(svglite)
})

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_path <- sub("^--file=", "", script_arg[[1]])
script_dir <- dirname(normalizePath(gsub("~\\+~", " ", script_path), mustWork = TRUE))
input_file <- file.path(script_dir, "Figure2-C_source_data.tsv.gz")
output_prefix <- file.path(script_dir, "Figure2-C")
target_label <- "polysome-AMFs"
group_levels <- c("Non cancer-genes", "Known TSGs", "Known oncogenes", target_label)
legend_levels <- c("Known TSGs", "Non cancer-genes", "Known oncogenes", target_label)
scope_levels <- c("Breast cancer cell lines", "Other cancer cell lines")

if (!file.exists(input_file)) {
  stop("Missing source data: ", input_file, ". Run Figure2-C_prepare_boxplot.py first.")
}

font_family <- "sans"
font_size_pt <- 14
text_size_mm <- font_size_pt * 25.4 / 72.27
if (requireNamespace("systemfonts", quietly = TRUE)) {
  available_fonts <- unique(systemfonts::system_fonts()$family)
  if ("Arial" %in% available_fonts) font_family <- "Arial"
}

# Exact palette used by the reference plot.
palette <- c(
  "Non cancer-genes" = "#B3B3B3",
  "Known TSGs" = "#89A54E",
  "Known oncogenes" = "#AA4643",
  "polysome-AMFs" = "#4572A7"
)

# A dedicated horizontal boxplot legend key. ggplot2's default boxplot key is
# vertical and does not inherit coord_flip(), which made the whisker, box and
# median appear offset from one another in the collected legend.
draw_key_horizontal_boxplot <- function(data, params, size) {
  key_colour <- if (is.null(data$colour) || is.na(data$colour[1])) "#333333" else data$colour[1]
  key_fill <- if (is.null(data$fill) || is.na(data$fill[1])) "#B3B3B3" else data$fill[1]
  key_alpha <- if (is.null(data$alpha) || is.na(data$alpha[1])) 1 else data$alpha[1]

  grid::grobTree(
    # Dashed whiskers, centred vertically.
    grid::segmentsGrob(
      x0 = grid::unit(0.07, "npc"), x1 = grid::unit(0.93, "npc"),
      y0 = grid::unit(0.50, "npc"), y1 = grid::unit(0.50, "npc"),
      gp = grid::gpar(col = key_colour, lwd = 1.2, lty = 2, lineend = "butt")
    ),
    # Whisker end caps.
    grid::segmentsGrob(
      x0 = grid::unit(c(0.07, 0.93), "npc"), x1 = grid::unit(c(0.07, 0.93), "npc"),
      y0 = grid::unit(c(0.31, 0.31), "npc"), y1 = grid::unit(c(0.69, 0.69), "npc"),
      gp = grid::gpar(col = key_colour, lwd = 1.2, lineend = "butt")
    ),
    # Filled box, centred on the same horizontal axis.
    grid::rectGrob(
      x = grid::unit(0.50, "npc"), y = grid::unit(0.50, "npc"),
      width = grid::unit(0.48, "npc"), height = grid::unit(0.62, "npc"),
      gp = grid::gpar(
        fill = scales::alpha(key_fill, key_alpha),
        col = key_colour, lwd = 1.2, linejoin = "mitre"
      )
    ),
    # Median line.
    grid::segmentsGrob(
      x0 = grid::unit(0.50, "npc"), x1 = grid::unit(0.50, "npc"),
      y0 = grid::unit(0.19, "npc"), y1 = grid::unit(0.81, "npc"),
      gp = grid::gpar(col = key_colour, lwd = 1.5, lineend = "butt")
    )
  )
}

format_p <- function(p) {
  if (is.na(p)) return("p=NA")
  if (p == 0) return("p<2.23e-308")
  paste0("p=", formatC(p, format = "e", digits = 2))
}

boxplot_summary <- function(values) {
  values <- values[is.finite(values)]
  if (!length(values)) stop("Cannot summarize an empty group")
  hinges <- as.numeric(quantile(values, c(0.25, 0.5, 0.75), type = 7, names = FALSE))
  iqr <- hinges[3] - hinges[1]
  inside <- values[values >= hinges[1] - 1.5 * iqr & values <= hinges[3] + 1.5 * iqr]
  data.table(
    ymin = min(inside),
    lower = hinges[1],
    middle = hinges[2],
    upper = hinges[3],
    ymax = max(inside),
    n = length(values),
    min = min(values),
    max = max(values),
    outside_display_range = sum(values < -6 | values > 6)
  )
}

if ("--render-only" %in% commandArgs(trailingOnly = TRUE)) {
  message("Using saved statistics for layout-only revision...")
  box_stats <- fread(paste0(output_prefix, "_boxplot_statistics.tsv"))
  tests <- fread(paste0(output_prefix, "_wilcoxon_tests.tsv"))
  box_stats[, Group := factor(Group, levels = group_levels)]
  box_stats[, Scope := factor(Scope, levels = scope_levels)]
  cell_audit <- fread(file.path(script_dir, "Figure2-C_cell_lines.tsv"))
  n_breast <- cell_audit[Scope == scope_levels[1], .N]
  n_other <- cell_audit[Scope == scope_levels[2], .N]
} else {
message("Reading prepared source data...")
if (requireNamespace("R.utils", quietly = TRUE)) {
  raw <- fread(
    file = input_file,
    select = c("Group", "CellLine", "Score"),
    colClasses = c(Group = "character", CellLine = "character", Score = "numeric"),
    showProgress = interactive()
  )
} else {
  raw <- fread(
    cmd = paste("gzip -dc", shQuote(input_file)),
    select = c("Group", "CellLine", "Score"),
    colClasses = c(Group = "character", CellLine = "character", Score = "numeric"),
    showProgress = interactive()
  )
}

# Normalize the source label so the legend exactly matches the reference.
raw[Group == "Non-cancer genes", Group := "Non cancer-genes"]

unexpected_groups <- setdiff(unique(raw$Group), group_levels)
if (length(unexpected_groups)) stop("Unexpected groups: ", paste(unexpected_groups, collapse = ", "))
missing_groups <- setdiff(group_levels, unique(raw$Group))
if (length(missing_groups)) stop("Missing groups: ", paste(missing_groups, collapse = ", "))

is_breast <- endsWith(raw$CellLine, "_BREAST")
n_breast <- uniqueN(raw$CellLine[is_breast])
n_other <- uniqueN(raw$CellLine[!is_breast])
cell_audit <- fread(file.path(script_dir, "Figure2-C_cell_lines.tsv"))
stopifnot(n_breast == cell_audit[Scope == scope_levels[1], .N])
stopifnot(n_other == cell_audit[Scope == scope_levels[2], .N])
message("Verified cell-line counts: breast = ", n_breast, "; other = ", n_other)

raw[, Scope := ifelse(
  endsWith(CellLine, "_BREAST"),
  "Breast cancer cell lines",
  "Other cancer cell lines"
)]
raw[, CellLine := NULL]
plot_data <- raw
rm(raw, is_breast)
gc(verbose = FALSE)

plot_data[, Scope := factor(Scope, levels = scope_levels)]
plot_data[, Group := factor(Group, levels = group_levels)]

summary_list <- list()
test_list <- list()

for (scope_name in scope_levels) {
  message("Calculating box statistics and Wilcoxon tests: ", scope_name)
  scope_data <- plot_data[Scope == scope_name]
  values_by_group <- split(scope_data$Score, scope_data$Group, drop = FALSE)

  for (group_name in group_levels) {
    group_summary <- boxplot_summary(values_by_group[[group_name]])
    group_summary[, `:=`(Scope = scope_name, Group = group_name)]
    summary_list[[length(summary_list) + 1L]] <- group_summary
  }

  reference <- values_by_group[["Non cancer-genes"]]
  for (comparison_group in c("Known oncogenes", target_label)) {
    test <- suppressWarnings(
      wilcox.test(
        values_by_group[[comparison_group]],
        reference,
        alternative = "two.sided",
        exact = FALSE,
        correct = TRUE
      )
    )
    test_list[[length(test_list) + 1L]] <- data.table(
      Scope = scope_name,
      Comparison = paste(comparison_group, "vs Non cancer-genes"),
      Group1 = comparison_group,
      Group2 = "Non cancer-genes",
      n1 = length(values_by_group[[comparison_group]]),
      n2 = length(reference),
      W = unname(test$statistic),
      p_value = test$p.value,
      p_label = format_p(test$p.value),
      test = "Two-sided Wilcoxon rank-sum test; continuity correction"
    )
  }

  rm(scope_data, values_by_group, reference)
  gc(verbose = FALSE)
}

box_stats <- rbindlist(summary_list, use.names = TRUE)
tests <- rbindlist(test_list, use.names = TRUE)
box_stats[, `:=`(
  Scope = factor(Scope, levels = scope_levels),
  Group = factor(Group, levels = group_levels)
)]

fwrite(box_stats, paste0(output_prefix, "_boxplot_statistics.tsv"), sep = "\t", quote = FALSE)
fwrite(tests, paste0(output_prefix, "_wilcoxon_tests.tsv"), sep = "\t", quote = FALSE)
}

theme_replica <- function() {
  theme_classic(base_size = font_size_pt, base_family = font_family) +
    theme(
      text = element_text(size = font_size_pt, family = font_family),
      axis.line = element_line(linewidth = 0.5, colour = "black"),
      axis.line.y = element_blank(),
      axis.ticks.x = element_line(linewidth = 0.45, colour = "black"),
      axis.ticks.y = element_blank(),
      axis.text.x = element_text(size = font_size_pt, colour = "black"),
      axis.text.y = element_blank(),
      axis.title.x = element_text(size = font_size_pt, margin = margin(t = 1)),
      legend.position = "top",
      legend.justification = "center",
      legend.title = element_blank(),
      legend.text = element_text(size = font_size_pt, colour = "black"),
      legend.key.width = grid::unit(9, "mm"),
      legend.key.height = grid::unit(5.2, "mm"),
      legend.spacing.x = grid::unit(1.5, "mm"),
      legend.key.spacing.x = grid::unit(2, "mm"),
      legend.key.spacing.y = grid::unit(0, "mm"),
      legend.margin = margin(0, 0, 0, 0),
      legend.box.margin = margin(0, 0, 0, 0),
      plot.title = element_text(size = font_size_pt, face = "plain", hjust = 0.5, margin = margin(b = 2)),
      plot.margin = margin(0.6, 3, 0.6, 3, unit = "mm")
    )
}

make_panel <- function(scope_name, title) {
  panel_data <- box_stats[Scope == scope_name]
  oncogene_p <- tests[Scope == scope_name & Group1 == "Known oncogenes", p_label]
  target_p <- tests[Scope == scope_name & Group1 == target_label, p_label]

  ggplot(panel_data, aes(x = Group, fill = Group)) +
    geom_boxplot(
      aes(ymin = ymin, lower = lower, middle = middle, upper = upper, ymax = ymax),
      stat = "identity", width = 0.66, linewidth = 0.48,
      colour = "#333333", linetype = "dashed",
      key_glyph = draw_key_horizontal_boxplot
    ) +
    geom_tile(
      aes(y = (lower + upper) / 2, height = upper - lower),
      width = 0.66, colour = "#333333", linewidth = 0.48,
      show.legend = FALSE
    ) +
    geom_errorbar(
      aes(ymin = middle, ymax = middle),
      width = 0.66, linewidth = 0.75, colour = "#333333",
      show.legend = FALSE
    ) +
    annotate("segment", x = 1, xend = 1, y = -5.35, yend = -5.58, linewidth = 0.45) +
    annotate("segment", x = 3, xend = 3, y = -5.35, yend = -5.58, linewidth = 0.45) +
    annotate("segment", x = 1, xend = 3, y = -5.58, yend = -5.58, linewidth = 0.45) +
    annotate(
      "text", x = 2, y = -5.90, label = oncogene_p, angle = 90,
      size = text_size_mm, family = font_family, hjust = 0.5
    ) +
    annotate("segment", x = 1, xend = 1, y = 3.95, yend = 4.18, linewidth = 0.45) +
    annotate("segment", x = 4, xend = 4, y = 3.95, yend = 4.18, linewidth = 0.45) +
    annotate("segment", x = 1, xend = 4, y = 4.18, yend = 4.18, linewidth = 0.45) +
    annotate(
      "text", x = 2.5, y = 4.52, label = target_p, angle = 90,
      size = text_size_mm, family = font_family, hjust = 0.5
    ) +
    scale_fill_manual(values = palette, breaks = legend_levels, drop = FALSE) +
    scale_y_continuous(breaks = seq(-6, 6, 2), expand = expansion(mult = 0)) +
    coord_flip(ylim = c(-6, 6), clip = "off") +
    labs(x = NULL, y = "shRNA log\u2082(FC)", title = title) +
    guides(fill = guide_legend(ncol = 2, byrow = TRUE)) +
    theme_replica()
}

p_breast <- make_panel(
  "Breast cancer cell lines",
  paste(n_breast, "breast cancer cell lines")
)
p_other <- make_panel(
  "Other cancer cell lines",
  paste(n_other, "other cancer cell lines")
)

combined <- guide_area() / p_breast / p_other +
  plot_layout(guides = "collect", heights = c(0.22, 1, 1)) &
  theme(
    legend.position = "top",
    legend.justification = "center",
    legend.box.margin = margin(0, 0, 0, 0)
  )

save_plot <- function(plot, stem, width_mm, height_mm, include_tiff = TRUE) {
  width_in <- width_mm / 25.4
  height_in <- height_mm / 25.4

  svglite::svglite(paste0(stem, ".svg"), width = width_in, height = height_in)
  print(plot)
  dev.off()

  grDevices::cairo_pdf(
    paste0(stem, ".pdf"), width = width_in, height = height_in, family = font_family
  )
  print(plot)
  dev.off()

  ragg::agg_png(
    paste0(stem, ".png"), width = width_in, height = height_in, units = "in", res = 600
  )
  print(plot)
  dev.off()

  if (include_tiff) {
    ragg::agg_tiff(
      paste0(stem, ".tiff"), width = width_in, height = height_in,
      units = "in", res = 600, compression = "lzw"
    )
    print(plot)
    dev.off()
  }
}

message("Exporting Figure 2C with the spreadsheet gene set...")
save_plot(combined, output_prefix, 128, 128, include_tiff = TRUE)
save_plot(p_breast + theme(legend.position = "none"), paste0(output_prefix, "_breast-cell-lines"), 128, 56, include_tiff = FALSE)
save_plot(p_other + theme(legend.position = "none"), paste0(output_prefix, "_other-cell-lines"), 128, 56, include_tiff = FALSE)

message("Done. Main output: ", output_prefix, ".pdf")
