# Figure_1B: RBP expression matrices

This directory groups analysis code, data, and statistical tables under a stable analysis identifier.

- `code/`: analysis and data-processing scripts.
- `data/`: inputs, matrices, gene lists, and saved R objects.
- `results/`: statistical tables, parameters, and computational checks.
- `FILE_MAP.tsv`: retained-file locations, checksums, and original relative paths.

Shared inputs can be stored in another analysis directory or `shared_data/`. Use `FILE_MAP.tsv` to locate them.

## Running code

Restore the original relative directory layout into a separate workspace using `tools/restore_project_layout.py` from the repository root. Run the commands below inside that restored workspace.

Working directory: `.`.

```sh
Rscript "01TCGA 筛选/Figure1A/Figure1A_plot_TCGA_BRCA_TN_HR_HER2_markers_with_NORMAL_reference_heatmap.R"
```

Software dependencies are listed in `environment/`. Runtime input paths should be checked before execution.
