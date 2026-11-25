# Figure_S2A: BRCA survival analysis

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
FIGURE2SB_MINPROP=0.10 Rscript Figure2/1A/TCGA_BRCA_full_five_gene_OS/scripts/Figure2S-B_top3_TCGA_BRCA_full_OS_KM_square.R
```

Software dependencies are listed in `environment/`. Runtime input paths should be checked before execution.
