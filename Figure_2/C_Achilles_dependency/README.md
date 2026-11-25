# Figure_2C: Achilles dependency analysis

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
Rscript "09细胞增殖/Figure2-C_AMFs36/Figure2-C_boxplot.R"
```

Working directory: `.`.

```sh
python3 "09细胞增殖/Figure2-C_AMFs36/Figure2-C_prepare_boxplot.py"
```

Software dependencies are listed in `environment/`. Runtime input paths should be checked before execution.
