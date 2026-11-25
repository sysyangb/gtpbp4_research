#!/usr/bin/env python3
"""Independently validate gene coverage and target box statistics."""

import csv
import math
from collections import defaultdict
from pathlib import Path
from statistics import median

from openpyxl import load_workbook


BASE = Path(__file__).resolve().parent


def read_tsv(name):
    with (BASE / name).open(encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def quantile(values, probability):
    position = (len(values) - 1) * probability
    low, high = math.floor(position), math.ceil(position)
    return values[low] + (values[high] - values[low]) * (position - low)


workbook = load_workbook(BASE / "Figure2-C_AMFs36_input.xlsx", read_only=True, data_only=True)
genes = [row[0] for row in workbook.active.iter_rows(min_row=2, min_col=2, max_col=2, values_only=True) if row[0]]
workbook.close()
assert len(genes) == len(set(genes)) == 36

blocks = defaultdict(list)
with (BASE.parent / "Achilles_QC_v2.4.3.rnai.gct").open() as handle:
    next(handle)
    rows, n_cells = map(int, next(handle).strip().split("\t"))
    cells = next(handle).rstrip().split("\t")[2:]
    observed_rows = 0
    for line in handle:
        observed_rows += 1
        fields = line.rstrip().split("\t")
        if fields[1] in genes:
            blocks[fields[1]].append(fields[2:])
assert observed_rows == rows and len(cells) == n_cells == 216
assert sum(cell.endswith("_BREAST") for cell in cells) == 13
assert len(blocks) == 20

coverage = {row["Gene"]: row for row in read_tsv("Figure2-C_gene_coverage.tsv")}
target_values = defaultdict(list)
for gene in genes:
    row = coverage[gene]
    assert (row["Exact_match_in_GCT"] == "True") == (gene in blocks)
    assert int(row["GCT_hairpin_rows"]) == len(blocks[gene])
    for is_breast, name, scope in (
        (True, "Breast", "Breast cancer cell lines"),
        (False, "Other", "Other cancer cell lines"),
    ):
        available_cells = available_n = retained_cells = retained_n = 0
        for index, cell in enumerate(cells):
            if cell.endswith("_BREAST") != is_breast:
                continue
            scores = [float(hairpin[index]) for hairpin in blocks[gene] if hairpin[index].strip().lower() not in ("", "nan")]
            if not scores:
                continue
            available_cells += 1
            available_n += len(scores)
            if any(abs(score) >= 1 for score in scores):
                retained_cells += 1
                retained_n += len(scores)
                target_values[scope].extend(scores)
        for key, value in (
            ("available_cell_lines", available_cells), ("available_observations", available_n),
            ("retained_cell_lines", retained_cells), ("retained_observations", retained_n),
        ):
            assert int(row[f"{name}_{key}"]) == value, (gene, name, key)

summaries = read_tsv("Figure2-C_boxplot_statistics.tsv")
for scope, values in target_values.items():
    values.sort()
    summary = next(row for row in summaries if row["Scope"] == scope and row["Group"] == "polysome-AMFs")
    q1, q2, q3 = (quantile(values, p) for p in (0.25, 0.5, 0.75))
    inside = [value for value in values if q1 - 1.5 * (q3 - q1) <= value <= q3 + 1.5 * (q3 - q1)]
    for column, expected in (("lower", q1), ("middle", q2), ("upper", q3), ("ymin", min(inside)), ("ymax", max(inside))):
        assert math.isclose(float(summary[column]), expected, abs_tol=1e-10), (scope, column)
    assert int(summary["n"]) == len(values)

tests = read_tsv("Figure2-C_wilcoxon_tests.tsv")
assert len(tests) == 4 and all(0 <= float(row["p_value"]) <= 1 for row in tests)
svg = (BASE / "Figure2-C.svg").read_text()
for label in ("13 breast cancer cell lines", "203 other cancer cell lines", "polysome-AMFs"):
    assert label in svg
for extension in ("svg", "pdf", "png", "tiff"):
    assert (BASE / f"Figure2-C.{extension}").stat().st_size > 1000

matched = [gene for gene in genes if gene in blocks]
missing = [gene for gene in genes if gene not in blocks]
notes = [
    "# Figure 2C：36 基因替换版",
    "",
    f"输入包含 36 个不重复基因，其中 {len(matched)} 个（{len(matched)/36:.1%}）在当前 Achilles shRNA 矩阵中精确匹配。",
    "乳腺癌与其他癌细胞系面板均实际纳入这 20 个基因。",
    "",
    "## 纳入基因",
    ", ".join(matched),
    "",
    "## 未匹配基因",
    ", ".join(missing),
    "",
    "PPAN-P2RY11 是矩阵中的不同符号，本次不替代 PPAN；未作未验证的别名合并。",
    "",
    "## 数据和统计定义",
    "当前矩阵共 216 个细胞系，按列名组织后缀归类为 13 个乳腺癌与 203 个其他癌细胞系。原图的 16/200 标签已更正。",
    "沿用原筛选规则：每个基因和细胞系至少有一个 shRNA 测量值 ≤−1 或 ≥1 时，保留该组合全部非缺失测量值。",
    "箱体表示第 25–75 百分位数，中线表示中位数，须线延伸至 1.5×IQR 内的极值；离群点不显示，横轴范围为 −6 至 6。",
    "P 值沿用原方法，为各参照组及目标组与非癌基因组之间的双侧 Wilcoxon 秩和检验，带连续性校正，未作多重比较校正。",
    "统计单位为汇总的 shRNA 测量值，不等同于独立基因数或独立生物学重复数。",
    "",
    "| 面板 | 实际纳入基因 | 基因×细胞系组合 | shRNA 测量值 | 目标组中位数 | 目标组 vs 非癌基因 P 值 |",
    "| --- | ---: | ---: | ---: | ---: | --- |",
]
audit = read_tsv("Figure2-C_data_audit.tsv")
for scope, label in (("Breast cancer cell lines", "乳腺癌（13 个细胞系）"), ("Other cancer cell lines", "其他癌（203 个细胞系）")):
    row = next(row for row in audit if row["Scope"] == scope and row["Group"] == "polysome-AMFs")
    test = next(row for row in tests if row["Scope"] == scope and row["Group1"] == "polysome-AMFs")
    notes.append(f"| {label} | {row['Retained_gene_count']} | {row['Retained_gene_cell_pairs']} | {row['Observations']} | {median(target_values[scope]):.4f} | {test['p_label']} |")
notes += ["", "## 图形与复现", "共享图例置顶居中并分两列，两个面板上下紧凑排列；原有配色、箱体顺序和比较括号保留。", "所有可见文字均为 14 pt。R 生成全部图像，128 × 128 mm，PNG/TIFF 为 600 dpi，SVG/PDF 保留文本。", "原四基因版结果保留在上级目录。", "", "在本目录运行：", "```sh", "python3 Figure2-C_prepare_boxplot.py", "Rscript Figure2-C_boxplot.R", "python3 Figure2-C_validate.py", "```", "Python 需要 openpyxl；R 需要 data.table、ggplot2、patchwork、ragg 和 svglite。", ""]
(BASE / "Figure2-C_report.md").write_text("\n".join(notes), encoding="utf-8")
with (BASE / "Figure2-C_validation.tsv").open("w", encoding="utf-8", newline="") as handle:
    writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
    writer.writerow(("Check", "Pass"))
    for check in ("36 unique input genes", "20 exact GCT matches", "13 breast and 203 other cell lines", "independent per-gene availability and filtering counts", "independent target boxplot quantiles and whiskers", "four finite valid P values", "correct labels in editable SVG", "four exported figure formats"):
        writer.writerow((check, True))
print("Validated: 36 input genes, 20 matched and plotted in each panel; 13/203 cell lines.")
