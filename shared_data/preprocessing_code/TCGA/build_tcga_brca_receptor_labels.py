#!/usr/bin/env python3
"""Build TCGA-BRCA TN / non-TN / normal labels with source traceability.

The script matches the exact column order in TCGA-BRCA_counts.txt and
TCGA-BRCA_TPM.txt. It uses:

1. GDC TCGA-BRCA Clinical Supplement BCR Biotab for ER/PR IHC and
   HER2 IHC/FISH fields.
2. The TCGA Nature 2012 curated ER/PR/HER2 final-status fields distributed
   by cBioPortal and mirrored by UCSC Xena.
3. The GDC Nature 2012 Final Full BRCA Sample Summary direct TN marker as
   an audit field.

No expression-derived receptor calls are made.
"""

from __future__ import annotations

import csv
import hashlib
from collections import Counter
from pathlib import Path


ROOT = Path(__file__).resolve().parent
SOURCES = ROOT / "external_receptor_sources"

COUNTS = ROOT / "TCGA-BRCA_counts.txt"
TPM = ROOT / "TCGA-BRCA_TPM.txt"
GDC_BIOTAB = SOURCES / "nationwidechildrens.org_clinical_patient_brca.txt"
NATURE_CBIO = SOURCES / "cbioportal_brca_tcga_pub_data_clinical_sample.txt"
NATURE_XENA = SOURCES / "UCSC_Xena_TCGA_BRCA_clinicalMatrix.tsv"
NATURE_DIRECT_TN = SOURCES / "GDC_BRCA_2012_Final_Full_Sample_Summary.tsv"

SAMPLE_OUTPUT = ROOT / "TCGA-BRCA_receptor_subtype_by_sample.tsv"
PATIENT_OUTPUT = ROOT / "TCGA-BRCA_receptor_subtype_by_patient.tsv"
CONDITION_OUTPUT = ROOT / "TCGA-BRCA_TN_NON_TN_NORMAL.condition.tsv"
QC_OUTPUT = ROOT / "TCGA-BRCA_receptor_label_QC.txt"

CONCLUSIVE = {"Positive", "Negative"}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def normalize_binary(value: str | None) -> str:
    value = (value or "").strip().lower()
    if value == "positive":
        return "Positive"
    if value == "negative":
        return "Negative"
    return "Unknown"


def receptor_class(er: str, pr: str, her2: str) -> str:
    if er == pr == her2 == "Negative":
        return "TN"
    if "Positive" in (er, pr, her2):
        return "NON_TN"
    return "UNKNOWN"


def read_expression_samples(path: Path) -> list[str]:
    with path.open(newline="") as handle:
        return [value.strip('"') for value in next(csv.reader(handle, delimiter="\t"))]


def read_gdc_biotab(path: Path) -> dict[str, dict[str, str]]:
    with path.open(newline="", encoding="utf-8-sig") as handle:
        reader = csv.reader(handle, delimiter="\t")
        header = next(reader)
        next(reader)  # XML element names
        next(reader)  # CDE identifiers
        records = [dict(zip(header, row)) for row in reader]
    return {row["bcr_patient_barcode"]: row for row in records}


def read_cbio_nature(path: Path) -> dict[str, dict[str, str]]:
    with path.open(newline="") as handle:
        reader = csv.DictReader(
            (line for line in handle if not line.startswith("#")), delimiter="\t"
        )
        records = list(reader)
    return {row["PATIENT_ID"]: row for row in records}


def read_xena_nature(path: Path) -> dict[str, dict[str, str]]:
    with path.open(newline="") as handle:
        records = list(csv.DictReader(handle, delimiter="\t"))
    return {row["sampleID"]: row for row in records}


def read_direct_tn(path: Path) -> dict[str, str]:
    with path.open(newline="") as handle:
        records = list(csv.DictReader(handle, delimiter="\t"))
    return {row["bcr_patient_barcode"]: row["Triple Negative"] for row in records}


def resolve_gdc_her2(row: dict[str, str] | None) -> tuple[str, bool]:
    if row is None:
        return "Unknown", False
    ihc = normalize_binary(row.get("her2_status_by_ihc"))
    fish = normalize_binary(row.get("her2_fish_status"))
    conclusive = {value for value in (ihc, fish) if value in CONCLUSIVE}
    assay_conflict = len(conclusive) > 1

    # A positive result in either submitted assay prevents a TN call. If there
    # is no positive result, a conclusive negative IHC or FISH result supports
    # HER2-negative status. The original fields and conflict flag are retained.
    if "Positive" in conclusive:
        return "Positive", assay_conflict
    if "Negative" in conclusive:
        return "Negative", assay_conflict
    return "Unknown", assay_conflict


def source_consensus(gdc_value: str, nature_value: str) -> tuple[str, bool]:
    values = {v for v in (gdc_value, nature_value) if v in CONCLUSIVE}
    if len(values) > 1:
        return "Conflict", True
    if len(values) == 1:
        return next(iter(values)), False
    return "Unknown", False


def nature_first(gdc_value: str, nature_value: str) -> str:
    if nature_value in CONCLUSIVE:
        return nature_value
    if gdc_value in CONCLUSIVE:
        return gdc_value
    return "Unknown"


def validate_xena_mirror(
    cbio: dict[str, dict[str, str]], xena: dict[str, dict[str, str]]
) -> tuple[int, int]:
    mapping = {
        "ER_STATUS": "ER_Status_nature2012",
        "PR_STATUS": "PR_Status_nature2012",
        "HER2_STATUS": "HER2_Final_Status_nature2012",
    }
    compared = 0
    discrepancies = 0
    for patient_id, cbio_row in cbio.items():
        xena_row = xena.get(f"{patient_id}-01")
        if xena_row is None:
            continue
        for cbio_field, xena_field in mapping.items():
            a = normalize_binary(cbio_row.get(cbio_field))
            b = normalize_binary(xena_row.get(xena_field))
            if a in CONCLUSIVE and b in CONCLUSIVE:
                compared += 1
                if a != b:
                    discrepancies += 1
    if discrepancies:
        raise ValueError(
            f"cBioPortal and UCSC Xena Nature 2012 fields disagree in "
            f"{discrepancies}/{compared} conclusive comparisons"
        )
    return compared, discrepancies


def sample_type(sample_id: str) -> tuple[str, str]:
    code = sample_id[13:15]
    labels = {
        "01": "PRIMARY_TUMOR",
        "06": "METASTATIC",
        "11": "SOLID_TISSUE_NORMAL",
    }
    return code, labels.get(code, "OTHER")


def build_tumor_record(
    patient_id: str,
    gdc: dict[str, dict[str, str]],
    nature: dict[str, dict[str, str]],
    direct_tn: dict[str, str],
) -> dict[str, str]:
    gdc_row = gdc.get(patient_id)
    nature_row = nature.get(patient_id)

    gdc_er = normalize_binary(gdc_row.get("er_status_by_ihc") if gdc_row else None)
    gdc_pr = normalize_binary(gdc_row.get("pr_status_by_ihc") if gdc_row else None)
    gdc_her2, gdc_her2_assay_conflict = resolve_gdc_her2(gdc_row)
    gdc_class = receptor_class(gdc_er, gdc_pr, gdc_her2)

    nature_er = normalize_binary(nature_row.get("ER_STATUS") if nature_row else None)
    nature_pr = normalize_binary(nature_row.get("PR_STATUS") if nature_row else None)
    nature_her2 = normalize_binary(
        nature_row.get("HER2_STATUS") if nature_row else None
    )
    nature_class = receptor_class(nature_er, nature_pr, nature_her2)

    consensus_er, er_conflict = source_consensus(gdc_er, nature_er)
    consensus_pr, pr_conflict = source_consensus(gdc_pr, nature_pr)
    consensus_her2, her2_conflict = source_consensus(gdc_her2, nature_her2)

    selected_er = nature_first(gdc_er, nature_er)
    selected_pr = nature_first(gdc_pr, nature_pr)
    selected_her2 = nature_first(gdc_her2, nature_her2)
    best_label = receptor_class(selected_er, selected_pr, selected_her2)

    class_conflict = (
        gdc_class in {"TN", "NON_TN"}
        and nature_class in {"TN", "NON_TN"}
        and gdc_class != nature_class
    )
    recommended_label = "UNKNOWN" if class_conflict else best_label

    used_nature = any(v in CONCLUSIVE for v in (nature_er, nature_pr, nature_her2))
    used_gdc_fill = any(
        n not in CONCLUSIVE and g in CONCLUSIVE
        for g, n in (
            (gdc_er, nature_er),
            (gdc_pr, nature_pr),
            (gdc_her2, nature_her2),
        )
    )
    if used_nature and used_gdc_fill:
        label_source = "NATURE2012_CURATED+GDC_FILL"
    elif used_nature:
        label_source = "NATURE2012_CURATED"
    elif gdc_row:
        label_source = "GDC_CLINICAL_SUPPLEMENT"
    else:
        label_source = "NO_RECEPTOR_SOURCE"

    flags: list[str] = []
    conflicts = [
        name
        for name, is_conflict in (
            ("ER", er_conflict),
            ("PR", pr_conflict),
            ("HER2", her2_conflict),
        )
        if is_conflict
    ]
    if class_conflict:
        flags.append("SOURCE_CLASS_CONFLICT")
    elif conflicts:
        flags.append("RECEPTOR_CONFLICT_CLASS_STABLE")
    if gdc_her2_assay_conflict:
        flags.append("GDC_HER2_IHC_FISH_CONFLICT")
    if best_label == "UNKNOWN":
        flags.append("INCOMPLETE_RECEPTOR_STATUS")
    if gdc_row is None:
        flags.append("MISSING_GDC_PATIENT_ROW")
    if not flags:
        flags.append("OK")

    if class_conflict:
        confidence = "CONFLICT"
    elif best_label == "UNKNOWN":
        confidence = "UNKNOWN"
    elif gdc_class == nature_class and gdc_class in {"TN", "NON_TN"}:
        confidence = "HIGH"
    elif direct_tn.get(patient_id) == "TN" and best_label == "TN":
        confidence = "HIGH"
    else:
        confidence = "MEDIUM"

    return {
        "patient_id": patient_id,
        "ER_status_selected": selected_er,
        "PR_status_selected": selected_pr,
        "HER2_status_selected": selected_her2,
        "ER_status_consensus": consensus_er,
        "PR_status_consensus": consensus_pr,
        "HER2_status_consensus": consensus_her2,
        "label_recommended": recommended_label,
        "label_best_available": best_label,
        "label_source": label_source,
        "confidence": confidence,
        "quality_flag": "|".join(flags),
        "receptor_source_conflicts": ",".join(conflicts),
        "gdc_er_status_by_ihc": gdc_row.get("er_status_by_ihc", "") if gdc_row else "",
        "gdc_pr_status_by_ihc": gdc_row.get("pr_status_by_ihc", "") if gdc_row else "",
        "gdc_her2_status_by_ihc": gdc_row.get("her2_status_by_ihc", "") if gdc_row else "",
        "gdc_her2_ihc_score": gdc_row.get("her2_ihc_score", "") if gdc_row else "",
        "gdc_her2_fish_status": gdc_row.get("her2_fish_status", "") if gdc_row else "",
        "gdc_receptor_class": gdc_class,
        "nature2012_er_status": nature_row.get("ER_STATUS", "") if nature_row else "",
        "nature2012_pr_status": nature_row.get("PR_STATUS", "") if nature_row else "",
        "nature2012_her2_final_status": nature_row.get("HER2_STATUS", "") if nature_row else "",
        "nature2012_receptor_class": nature_class,
        "nature2012_direct_tn_marker": direct_tn.get(patient_id, ""),
    }


def write_tsv(path: Path, rows: list[dict[str, str]], fields: list[str]) -> None:
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    counts_samples = read_expression_samples(COUNTS)
    tpm_samples = read_expression_samples(TPM)
    if counts_samples != tpm_samples:
        raise ValueError("Counts and TPM sample columns differ or are in different order")
    if len(counts_samples) != len(set(counts_samples)):
        raise ValueError("Expression matrix contains duplicate sample columns")

    gdc = read_gdc_biotab(GDC_BIOTAB)
    nature = read_cbio_nature(NATURE_CBIO)
    xena = read_xena_nature(NATURE_XENA)
    direct_tn = read_direct_tn(NATURE_DIRECT_TN)
    xena_compared, xena_discrepancies = validate_xena_mirror(nature, xena)

    patients = sorted({sample[:12] for sample in counts_samples})
    patient_rows = [build_tumor_record(p, gdc, nature, direct_tn) for p in patients]
    patient_map = {row["patient_id"]: row for row in patient_rows}

    sample_rows: list[dict[str, str]] = []
    condition_rows: list[dict[str, str]] = []
    for sample in counts_samples:
        patient_id = sample[:12]
        code, kind = sample_type(sample)
        row = dict(patient_map[patient_id])
        row = {"sample_id": sample, "sample_type_code": code, "sample_type": kind, **row}
        if code == "11":
            row.update(
                {
                    "ER_status_selected": "Not_Applicable",
                    "PR_status_selected": "Not_Applicable",
                    "HER2_status_selected": "Not_Applicable",
                    "ER_status_consensus": "Not_Applicable",
                    "PR_status_consensus": "Not_Applicable",
                    "HER2_status_consensus": "Not_Applicable",
                    "label_recommended": "NORMAL",
                    "label_best_available": "NORMAL",
                    "label_source": "TCGA_BARCODE_SAMPLE_TYPE_11",
                    "confidence": "HIGH",
                    "quality_flag": "NORMAL_NOT_RECEPTOR_CLASSIFIED",
                    "receptor_source_conflicts": "",
                }
            )
        elif code == "06":
            # TCGA receptor annotations are patient/diagnosis level and may not
            # have been measured on the metastatic aliquot itself.
            if row["quality_flag"] == "OK":
                row["quality_flag"] = "METASTATIC_USES_PATIENT_LEVEL_RECEPTORS"
            else:
                row["quality_flag"] += "|METASTATIC_USES_PATIENT_LEVEL_RECEPTORS"
            if row["confidence"] == "HIGH":
                row["confidence"] = "MEDIUM"
        elif code not in {"01", "06"}:
            row.update(
                {
                    "label_recommended": "UNKNOWN",
                    "label_best_available": "UNKNOWN",
                    "confidence": "UNKNOWN",
                    "quality_flag": "UNSUPPORTED_SAMPLE_TYPE",
                }
            )
        sample_rows.append(row)
        condition_rows.append(
            {
                "sample_id": sample,
                "patient_id": patient_id,
                "condition": row["label_recommended"],
                "condition_best_available": row["label_best_available"],
                "quality_flag": row["quality_flag"],
            }
        )

    patient_fields = list(patient_rows[0])
    sample_fields = ["sample_id", "sample_type_code", "sample_type"] + patient_fields
    write_tsv(SAMPLE_OUTPUT, sample_rows, sample_fields)
    write_tsv(PATIENT_OUTPUT, patient_rows, patient_fields)
    write_tsv(
        CONDITION_OUTPUT,
        condition_rows,
        ["sample_id", "patient_id", "condition", "condition_best_available", "quality_flag"],
    )

    recommended_counts = Counter(row["label_recommended"] for row in sample_rows)
    best_counts = Counter(row["label_best_available"] for row in sample_rows)
    type_counts = Counter(row["sample_type"] for row in sample_rows)
    patient_recommended_counts = Counter(row["label_recommended"] for row in patient_rows)
    quality_counts = Counter(row["quality_flag"] for row in sample_rows)

    qc_lines = [
        "TCGA-BRCA receptor label QC",
        "",
        f"Expression samples: {len(counts_samples)}",
        f"Expression patients: {len(patients)}",
        f"Counts/TPM headers identical: YES",
        f"Sample types: {dict(sorted(type_counts.items()))}",
        f"Recommended sample labels: {dict(sorted(recommended_counts.items()))}",
        f"Best-available sample labels: {dict(sorted(best_counts.items()))}",
        f"Recommended patient labels: {dict(sorted(patient_recommended_counts.items()))}",
        f"GDC Clinical Supplement patients: {len(gdc)}",
        f"Nature 2012 curated patients: {len(nature)}",
        f"Nature 2012 Xena conclusive comparisons: {xena_compared}",
        f"Nature 2012 Xena conclusive discrepancies: {xena_discrepancies}",
        "",
        "Quality flag counts:",
    ]
    qc_lines.extend(f"  {key}: {value}" for key, value in sorted(quality_counts.items()))
    qc_lines.extend(
        [
            "",
            "Input SHA-256:",
            f"  {COUNTS.name}: {sha256(COUNTS)}",
            f"  {TPM.name}: {sha256(TPM)}",
            f"  {GDC_BIOTAB.name}: {sha256(GDC_BIOTAB)}",
            f"  {NATURE_CBIO.name}: {sha256(NATURE_CBIO)}",
            f"  {NATURE_XENA.name}: {sha256(NATURE_XENA)}",
            f"  {NATURE_DIRECT_TN.name}: {sha256(NATURE_DIRECT_TN)}",
            "",
            "Rule: TN requires ER=Negative, PR=Negative, and HER2=Negative.",
            "Rule: NON_TN requires at least one receptor=Positive.",
            "Rule: missing/equivocal/indeterminate/conflicting evidence is not forced into NON_TN.",
            "Rule: sample type code 11 is NORMAL; codes 01 and 06 are tumor samples.",
        ]
    )
    QC_OUTPUT.write_text("\n".join(qc_lines) + "\n")

    print(f"Wrote {SAMPLE_OUTPUT.name}: {len(sample_rows)} rows")
    print(f"Wrote {PATIENT_OUTPUT.name}: {len(patient_rows)} rows")
    print(f"Wrote {CONDITION_OUTPUT.name}: {len(condition_rows)} rows")
    print(f"Recommended sample labels: {dict(sorted(recommended_counts.items()))}")
    print(f"Best-available sample labels: {dict(sorted(best_counts.items()))}")


if __name__ == "__main__":
    main()
