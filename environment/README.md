# Software environment

The analysis code uses R and Python. Package import declarations are indexed in `R_packages_by_script.tsv`; this index is not a complete dependency lock file. Retained session and package-version tables provide additional environment information.

An R `renv.lock` is retained with its source mapping. Platform-specific installed package directories are not included. The directory restoration and checksum utilities use only the Python standard library. Other scripts may require packages such as pandas, NumPy, or openpyxl; consult their imports.

Install compatible dependencies before running an analysis. The execution commands are listed in `provenance/entrypoints.json`.
