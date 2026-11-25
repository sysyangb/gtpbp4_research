#!/usr/bin/env python3
"""Prepare traceable source data for the Figure 2C boxplots.

The script preserves the filtering logic of the original workflow: for each
gene/cell-line pair, all non-missing shRNA measurements are retained when at
least one measurement is <= -1 or >= 1. Gene blocks in the Achilles GCT are
processed as a stream, so the full 131 MB matrix is never held in memory.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import sys
from collections import Counter, defaultdict
from pathlib import Path

BASE_DIR = Path(__file__).resolve().parent


TARGET_LABEL = "polysome-AMFs"
NON_CANCER_LABEL = "Non-cancer genes"
ONCOGENE_LABEL = "Known oncogenes"
TSG_LABEL = "Known TSGs"


def read_gene_list(path: Path, header: str | None = None) -> set[str]:
    genes: set[str] = set()
    with path.open(encoding="utf-8-sig") as handle:
        for line_number, line in enumerate(handle, start=1):
            value = line.strip().split("\t", 1)[0].strip()
            if not value or (line_number == 1 and header and value == header):
                continue
            genes.add(value)
    if not genes:
        raise ValueError(f"No genes were read from {path}")
    return genes


def read_excel_targets(path: Path) -> list[str]:
    """Read the Gene column without changing the supplied workbook."""
    from openpyxl import load_workbook

    workbook = load_workbook(path, read_only=True, data_only=True)
    try:
        rows = workbook.worksheets[0].iter_rows(values_only=True)
        columns = [str(value).strip() if value is not None else "" for value in next(rows)]
        gene_index = columns.index("Gene")
        genes = [str(row[gene_index]).strip() for row in rows if row[gene_index] is not None]
    finally:
        workbook.close()
    if len(genes) != 36 or len(set(genes)) != 36 or any(not gene for gene in genes):
        raise ValueError("Expected 36 non-empty, unique genes in the supplied Gene column")
    return genes


def labels_for_gene(
    gene: str,
    targets: set[str],
    oncogenes: set[str],
    tumor_suppressors: set[str],
) -> tuple[str, ...]:
    """Match the original independent reference-set classification."""
    labels: list[str] = []
    if gene in targets:
        labels.append(TARGET_LABEL)
    if gene in oncogenes:
        labels.append(ONCOGENE_LABEL)
    if gene in tumor_suppressors:
        labels.append(TSG_LABEL)
    if not labels:
        labels.append(NON_CANCER_LABEL)
    return tuple(labels)


def parse_score(value: str) -> float | None:
    value = value.strip()
    if not value or value.lower() == "nan":
        return None
    return float(value)


def write_gene_block(
    gene: str,
    block: list[tuple[str, list[float | None]]],
    cell_lines: list[str],
    labels: tuple[str, ...],
    threshold: float,
    writer: csv.writer,
    counts: Counter[tuple[str, str]],
    retained_pairs: dict[tuple[str, str], set[tuple[str, str]]],
    target_metrics: dict[tuple[str, str], Counter],
) -> None:
    if not block:
        return

    for cell_index, cell_line in enumerate(cell_lines):
        observed = [
            (hairpin, scores[cell_index])
            for hairpin, scores in block
            if scores[cell_index] is not None
        ]
        if not observed:
            continue

        scopes = (
            ("Breast cancer cell lines",)
            if cell_line.endswith("_BREAST")
            else ("Other cancer cell lines",)
        )
        values = [score for _, score in observed]
        retained = min(values) <= -threshold or max(values) >= threshold
        if TARGET_LABEL in labels:
            metrics = target_metrics[(scopes[0], gene)]
            metrics["available_cell_lines"] += 1
            metrics["available_observations"] += len(observed)
            if retained:
                metrics["retained_cell_lines"] += 1
                metrics["retained_observations"] += len(observed)
        if not retained:
            continue
        for label in labels:
            for hairpin, score in observed:
                writer.writerow((label, gene, cell_line, hairpin, f"{score:.12g}"))
                for scope in scopes:
                    counts[(scope, label)] += 1
                    retained_pairs[(scope, label)].add((gene, cell_line))


def prepare(args: argparse.Namespace) -> None:
    target_order = read_excel_targets(args.target_xlsx)
    targets = set(target_order)
    with args.target_genes.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(("gene_symbol",))
        writer.writerows((gene,) for gene in target_order)
    oncogenes = read_gene_list(args.oncogenes)
    tumor_suppressors = read_gene_list(args.tumor_suppressors)

    counts: Counter[tuple[str, str]] = Counter()
    retained_pairs: dict[tuple[str, str], set[tuple[str, str]]] = defaultdict(set)
    target_metrics = defaultdict(Counter)
    hairpin_counts = Counter()
    observed_genes: set[str] = set()
    data_rows = 0

    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.gct.open(encoding="utf-8") as source, gzip.open(
        args.output, "wt", encoding="utf-8", newline="", compresslevel=3
    ) as sink:
        version = source.readline().strip()
        dimensions = source.readline().strip().split("\t")
        header = source.readline().rstrip("\n\r").split("\t")
        if not version.startswith("#") or len(dimensions) != 2 or header[:2] != ["Name", "Description"]:
            raise ValueError(f"Unexpected GCT header in {args.gct}")

        expected_rows, expected_cells = map(int, dimensions)
        cell_lines = header[2:]
        if len(cell_lines) != expected_cells:
            raise ValueError(
                f"GCT declares {expected_cells} cell lines but the header contains {len(cell_lines)}"
            )

        writer = csv.writer(sink, delimiter="\t", lineterminator="\n")
        writer.writerow(("Group", "Gene", "CellLine", "Hairpin", "Score"))

        current_gene: str | None = None
        current_block: list[tuple[str, list[float | None]]] = []

        for line_number, line in enumerate(source, start=4):
            fields = line.rstrip("\n\r").split("\t")
            if len(fields) != expected_cells + 2:
                raise ValueError(
                    f"Line {line_number} has {len(fields)} columns; expected {expected_cells + 2}"
                )
            hairpin, gene = fields[:2]
            scores = [parse_score(value) for value in fields[2:]]
            data_rows += 1
            observed_genes.add(gene)
            if gene in targets:
                hairpin_counts[gene] += 1

            if current_gene is None:
                current_gene = gene
            elif gene != current_gene:
                write_gene_block(
                    current_gene,
                    current_block,
                    cell_lines,
                    labels_for_gene(current_gene, targets, oncogenes, tumor_suppressors),
                    args.threshold,
                    writer,
                    counts,
                    retained_pairs,
                    target_metrics,
                )
                current_gene = gene
                current_block = []
            current_block.append((hairpin, scores))

        if current_gene is not None:
            write_gene_block(
                current_gene,
                current_block,
                cell_lines,
                labels_for_gene(current_gene, targets, oncogenes, tumor_suppressors),
                args.threshold,
                writer,
                counts,
                retained_pairs,
                target_metrics,
            )

    if data_rows != expected_rows:
        raise ValueError(f"GCT declares {expected_rows} rows but {data_rows} were read")

    missing_targets = sorted(targets - observed_genes)
    with args.gene_audit.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow((
            "Input_order", "Gene", "Exact_match_in_GCT", "GCT_hairpin_rows",
            "Breast_available_cell_lines", "Breast_available_observations",
            "Breast_retained_cell_lines", "Breast_retained_observations",
            "Other_available_cell_lines", "Other_available_observations",
            "Other_retained_cell_lines", "Other_retained_observations", "Note"
        ))
        for rank, gene in enumerate(target_order, 1):
            note = ""
            if gene == "PPAN" and gene not in observed_genes and "PPAN-P2RY11" in observed_genes:
                note = "PPAN-P2RY11 is a different GCT symbol; not substituted for PPAN"
            metrics = [
                target_metrics[(scope, gene)][key]
                for scope in ("Breast cancer cell lines", "Other cancer cell lines")
                for key in ("available_cell_lines", "available_observations", "retained_cell_lines", "retained_observations")
            ]
            writer.writerow((rank, gene, gene in observed_genes, hairpin_counts[gene], *metrics, note))

    with args.cell_line_audit.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(("CellLine", "Scope"))
        writer.writerows((cell, "Breast cancer cell lines" if cell.endswith("_BREAST") else "Other cancer cell lines") for cell in cell_lines)
    if missing_targets:
        print(
            "WARNING: target genes absent from the Achilles GCT and omitted from the plot: "
            + ", ".join(missing_targets),
            file=sys.stderr,
        )

    group_order = (NON_CANCER_LABEL, TSG_LABEL, ONCOGENE_LABEL, TARGET_LABEL)
    with args.audit.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(
            (
                "Scope",
                "Group",
                "Target_gene_count",
                "Matched_target_gene_count",
                "Cell_line_count",
                "Retained_gene_count",
                "Retained_gene_cell_pairs",
                "Observations",
                "Missing_target_genes",
            )
        )
        for scope in ("Breast cancer cell lines", "Other cancer cell lines"):
            for label in group_order:
                pairs = retained_pairs[(scope, label)]
                writer.writerow(
                    (
                        scope,
                        label,
                        len(targets) if label == TARGET_LABEL else "",
                        len(targets & observed_genes) if label == TARGET_LABEL else "",
                        sum(cell.endswith("_BREAST") == (scope == "Breast cancer cell lines") for cell in cell_lines),
                        len({gene for gene, _ in pairs}),
                        len(pairs),
                        counts[(scope, label)],
                        ",".join(missing_targets) if label == TARGET_LABEL else "",
                    )
                )

    print(f"Target genes: {len(targets)}")
    print(f"Exact matches in GCT: {len(targets & observed_genes)}/{len(targets)}")
    for scope in ("Breast cancer cell lines", "Other cancer cell lines"):
        retained = [gene for gene in target_order if target_metrics[(scope, gene)]["retained_observations"]]
        print(f"Target genes retained in {scope}: {len(retained)}")
    print(f"GCT rows: {data_rows}; cell lines: {expected_cells}")
    print(f"Source data: {args.output}")
    print(f"Audit table: {args.audit}")
    for scope in ("Breast cancer cell lines", "Other cancer cell lines"):
        summary = ", ".join(f"{label}={counts[(scope, label)]:,}" for label in group_order)
        print(f"{scope}: {summary}")


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--gct", type=Path, default=BASE_DIR.parent / "Achilles_QC_v2.4.3.rnai.gct")
    parser.add_argument("--target-xlsx", type=Path, default=BASE_DIR / "Figure2-C_AMFs36_input.xlsx")
    parser.add_argument(
        "--target-genes",
        type=Path,
        default=BASE_DIR / "Figure2-C_polysome-AMFs_genes.tsv",
    )
    parser.add_argument("--oncogenes", type=Path, default=BASE_DIR.parent / "oncogene.txt")
    parser.add_argument("--tumor-suppressors", type=Path, default=BASE_DIR.parent / "tumor_sp.txt")
    parser.add_argument("--threshold", type=float, default=1.0)
    parser.add_argument(
        "--output", type=Path, default=BASE_DIR / "Figure2-C_source_data.tsv.gz"
    )
    parser.add_argument("--audit", type=Path, default=BASE_DIR / "Figure2-C_data_audit.tsv")
    parser.add_argument("--gene-audit", type=Path, default=BASE_DIR / "Figure2-C_gene_coverage.tsv")
    parser.add_argument("--cell-line-audit", type=Path, default=BASE_DIR / "Figure2-C_cell_lines.tsv")
    return parser


if __name__ == "__main__":
    prepare(build_parser().parse_args())
