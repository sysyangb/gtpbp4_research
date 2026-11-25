# Figure_1A: RBP proportions

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
Rscript outputs/Figure1A_FC2/Figure1A.R .
```

Working directory: `.`.

```sh
RBP_FC=2 Rscript outputs/TNBC_RBP_proportion_inference_20261027/scripts/run_patient_bootstrap.R
```

Rebuild QC/patient_bootstrap_model.rds with prepare_patient_bootstrap.R before execution.

Working directory: `.`.

```sh
Rscript outputs/TNBC_RBP_proportion_inference_20261027/scripts/prepare_patient_bootstrap.R
```

Rebuild QC/fresh_fitted_DESeq2.rds with run_deseq2_rbp.R before execution.

Working directory: `.`.

```sh
Rscript outputs/TNBC_DESeq2_RBP_20261026_FC1p5/scripts/run_deseq2_rbp.R .
```

Software dependencies are listed in `environment/`. Runtime input paths should be checked before execution.
