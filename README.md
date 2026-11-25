# GTPBP4 Research

This repository contains analysis code, datasets, saved statistical objects, and result tables. Analysis identifiers provide a consistent directory structure.

## Contents

- `Figure_1/`, `Figure_2/`, `Figure_S1/`, `Figure_S2/`, and `Figure_S5/`: analysis-specific code, data, and statistics.
- `shared_data/`: shared expression, clinical, proteomics, screening, and annotation inputs.
- `environment/`: package requirements and environment configuration.
- `provenance/`: retained-file mappings, checksums, and execution commands.
- `tools/`: verification and directory restoration utilities.

Use [the analysis index](FIGURE_INDEX.md) to locate a directory. Each analysis has a `README.md` and `FILE_MAP.tsv`. Shared files are referenced by their repository locations.

## Retrieve data

Large datasets use Git LFS:

```sh
git lfs install
git clone https://github.com/sysyangb/gtpbp4_research.git
cd gtpbp4_research
git lfs pull
```

## Verify files

```sh
python3 tools/verify_archive.py
python3 tools/check_public_content.py
```

`provenance/file_manifest.tsv` records repository checksums and original-source checksums. A `transform` value identifies files with filtered non-data text or workbook metadata. The optional `--source-root /path/to/original_project` flag checks original sources against their separate source checksums.

## Run analyses

Reconstruct the expected relative input paths in a new, empty directory outside this repository:

```sh
python3 tools/restore_project_layout.py --output /path/to/new_analysis_workspace
```

The tool copies retained files. It does not rewrite absolute paths unless an original-root prefix is supplied with `--original-root-prefix /path/to/original_project`. Check runtime paths and package dependencies before running code.

Execution commands and working directories are recorded in `provenance/entrypoints.json` and each analysis README. Two regenerable model checkpoints are omitted: `fresh_fitted_DESeq2.rds` and `patient_bootstrap_model.rds`. Commands that require them identify the necessary preparation step.
