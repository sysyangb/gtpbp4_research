# Figure_1J: Proteomics functional enrichment

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
Rscript 08polysome/polysome_up_GO_BP_plus_AMFs_total_human_background_20260902/scripts/plot_polysome_GO_barplot.R
```

Working directory: `.`.

```sh
Rscript 08polysome/polysome_up_GO_BP_plus_AMFs_total_human_background_20260902/scripts/run_polysome_up_GO_BP_plus_AMFs_total_human_background.R
```

Working directory: `.`.

```sh
Rscript 08polysome/reanalysis_exclude_polysome3_20260901/scripts/run_reanalysis.R
```

Software dependencies are listed in `environment/`. Runtime input paths should be checked before execution.
