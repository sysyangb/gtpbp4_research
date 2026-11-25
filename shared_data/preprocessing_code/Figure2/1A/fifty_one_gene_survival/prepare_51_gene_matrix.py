#!/usr/bin/env python3
"""Prepare the 51-gene AMF ∩ polysome-up TCGA-BRCA expression matrix.

The script deliberately streams the 358 MB TPM file so that only the requested
rows and the 1,093 survival-analysis samples are loaded into memory.
"""

from __future__ import annotations

import csv
from pathlib import Path


HERE = Path(__file__).resolve().parent
WORKDIR = HERE.parent
VENN_DIR = (
    WORKDIR.parent.parent
    / "08polysome"
    / "polysome_up_502_TN131_AMFs215_Venn_20260902"
)
MEMBERSHIP = VENN_DIR / "results" / "Figure1-I_gene_set_membership.tsv"
AMF_SOURCE = VENN_DIR / "source_data" / "AMFs_215_source_rows.tsv"
TPM_FILE = WORKDIR / "TCGA-BRCA_TPM.txt"
REFERENCE_EXPR = WORKDIR / "exprTrain.txt"


def read_tsv(path: Path):
    with path.open(newline="", encoding="utf-8-sig") as handle:
        yield from csv.DictReader(handle, delimiter="\t")


def main() -> None:
    HERE.mkdir(parents=True, exist_ok=True)

    members = [
        row
        for row in read_tsv(MEMBERSHIP)
        if row["AMF"] == "TRUE" and row["Polysome_upregulated"] == "TRUE"
    ]
    if len(members) != 51:
        raise RuntimeError(f"Expected 51 AMF∩polysome genes, found {len(members)}")

    symbol_to_region = {row["gene_symbol"]: row["region_code"] for row in members}
    amf_symbol_to_ensembl = {
        row["gene_symbol_clean"]: row["gene_id"].split(".")[0]
        for row in read_tsv(AMF_SOURCE)
    }
    missing_ids = sorted(set(symbol_to_region) - set(amf_symbol_to_ensembl))
    if missing_ids:
        raise RuntimeError(f"Missing Ensembl IDs for: {', '.join(missing_ids)}")

    genes = sorted(symbol_to_region)
    ensembl_to_symbol = {amf_symbol_to_ensembl[g]: g for g in genes}

    with REFERENCE_EXPR.open(newline="", encoding="utf-8-sig") as handle:
        reference_header = next(csv.reader(handle, delimiter="\t"))
    patients = reference_header[1:]
    if len(patients) != 1093 or len(set(patients)) != len(patients):
        raise RuntimeError("Unexpected reference sample header")

    with TPM_FILE.open(newline="", encoding="utf-8-sig") as handle:
        reader = csv.reader(handle, delimiter="\t")
        tpm_header = next(reader)

        # Match the legacy preparation script: only primary/metastatic tumour
        # sample codes 01-09 are accepted, and the last matching column wins.
        patient_set = set(patients)
        patient_to_col = {}
        patient_to_aliquot = {}
        for index, aliquot in enumerate(tpm_header):
            parts = aliquot.split("-")
            if len(parts) < 4:
                continue
            patient = "-".join(parts[:3])
            sample_code = parts[3][:2]
            if patient in patient_set and sample_code in {f"{x:02d}" for x in range(1, 10)}:
                patient_to_col[patient] = index
                patient_to_aliquot[patient] = aliquot

        missing_patients = [p for p in patients if p not in patient_to_col]
        if missing_patients:
            raise RuntimeError(f"Missing TPM columns for {len(missing_patients)} patients")

        expression = {}
        source_ensembl = {}
        for row in reader:
            ensembl_full = row[0]
            ensembl = ensembl_full.split(".")[0]
            if ensembl not in ensembl_to_symbol:
                continue
            symbol = ensembl_to_symbol[ensembl]
            expression[symbol] = [row[patient_to_col[p]] for p in patients]
            source_ensembl[symbol] = ensembl_full

    missing_expression = sorted(set(genes) - set(expression))
    if missing_expression:
        raise RuntimeError(f"Missing expression rows for: {', '.join(missing_expression)}")

    with (HERE / "gene_list_51.tsv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["gene_symbol", "ensembl_id", "source_ensembl_id", "region_code"])
        for gene in genes:
            writer.writerow(
                [gene, amf_symbol_to_ensembl[gene], source_ensembl[gene], symbol_to_region[gene]]
            )

    with (HERE / "expression_51_genes_raw_tpm.tsv").open(
        "w", newline="", encoding="utf-8"
    ) as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["geneSymbol", *patients])
        for gene in genes:
            writer.writerow([gene, *expression[gene]])

    with (HERE / "sample_to_tpm_aliquot.tsv").open(
        "w", newline="", encoding="utf-8"
    ) as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["patient_id", "tpm_aliquot"])
        for patient in patients:
            writer.writerow([patient, patient_to_aliquot[patient]])

    with (HERE / "preparation_qc.tsv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["check", "value", "expected", "pass"])
        checks = [
            ("requested_genes", len(genes), 51),
            ("region_011_genes", sum(symbol_to_region[g] == "011" for g in genes), 46),
            ("region_111_genes", sum(symbol_to_region[g] == "111" for g in genes), 5),
            ("expression_genes_found", len(expression), 51),
            ("survival_samples_matched", len(patients), 1093),
        ]
        for label, value, expected in checks:
            writer.writerow([label, value, expected, str(value == expected).upper()])

    print(f"Prepared {len(genes)} genes × {len(patients)} patients in {HERE}")


if __name__ == "__main__":
    main()
