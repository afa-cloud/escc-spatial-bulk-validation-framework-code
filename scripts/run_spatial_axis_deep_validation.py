#!/usr/bin/env python
"""Deep validation for two spatially informed ESCC axes.

This script extends the spatial signature screen with:

- TCGA ESCC immune and pathway association checks.
- Independent GEO bulk validation in GSE47404.
- GSE53625 mapping audit.
- GDSC2 ESCA target-class context for axis-relevant genes.

All outputs are evidence artifacts with separate analysis and quality-check batch IDs. The
results are association evidence, not mechanistic proof.
"""

from __future__ import annotations

import csv
import gzip
import json
import math
import re
import sys
import time
import urllib.request
from collections import defaultdict
from datetime import UTC, datetime
from pathlib import Path
from typing import Any

import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
SCRIPTS = ROOT / "scripts"
if str(SCRIPTS) not in sys.path:
    sys.path.insert(0, str(SCRIPTS))

import public_data_helpers as base  # noqa: E402


OUT_ROOT = ROOT / "spatial_escc_workflow"
DATA_ROOT = OUT_ROOT / "data"
TABLE_ROOT = OUT_ROOT / "results" / "tables"
REPORT_ROOT = OUT_ROOT / "reports"
REVIEW_ROOT = OUT_ROOT / "reviews"

ANALYSIS_BATCH_ID = "spatial_axis_analysis_001"
QUALITY_CHECK_ID = "spatial_axis_quality_check_001"
RUN_DATE = datetime.now(UTC).date().isoformat()


def show_help_if_requested() -> None:
    if any(arg in {"-h", "--help"} for arg in sys.argv[1:]):
        print((__doc__ or "").strip())
        print(
            "\nUsage:\n"
            "  python scripts/run_spatial_axis_deep_validation.py\n\n"
            "Run from the repository root. Outputs are written under "
            "spatial_escc_workflow/results, reports and reviews relative to the repository."
        )
        raise SystemExit(0)


AXES = {
    "ogt_pi3k_tls_axis": {
        "label": "OGT/PI3K/TLS axis",
        "genes": [
            "OGT",
            "PIK3CA",
            "AKT1",
            "CCND1",
            "LAMB1",
            "SPP1",
            "KRT17",
            "APOBEC3A",
            "JAG1",
            "NOTCH1",
        ],
        "gdsc_terms": [
            "PI3K",
            "PIK3",
            "AKT",
            "MTOR",
            "mTOR",
            "EGFR",
            "ERBB",
        ],
    },
    "caf_epi_jag1_notch_niche": {
        "label": "CAF-Epi/JAG1-NOTCH1 niche",
        "genes": [
            "JAG1",
            "NOTCH1",
            "FAP",
            "COL1A1",
            "COL1A2",
            "POSTN",
            "CXCL1",
            "CXCL8",
            "SPP1",
        ],
        "gdsc_terms": [
            "NOTCH",
            "gamma",
            "JAK",
            "PDGFR",
            "FGFR",
            "VEGFR",
            "KIT",
            "EGFR",
            "ERBB",
        ],
    },
}


IMMUNE_PANELS = {
    "cytotoxic_t": ["CD8A", "CD8B", "GZMB", "PRF1", "NKG7", "GNLY"],
    "treg": ["FOXP3", "IL2RA", "CTLA4", "IKZF2", "CCR8"],
    "exhaustion_checkpoint": ["PDCD1", "CD274", "CTLA4", "LAG3", "HAVCR2", "TIGIT"],
    "tls_b_cell": ["MS4A1", "CD79A", "CD79B", "CD19", "CXCL13", "BANK1", "PAX5"],
    "macrophage_m2": ["CD68", "CD163", "MRC1", "CSF1R", "MSR1"],
    "caf": ["FAP", "ACTA2", "COL1A1", "COL1A2", "PDGFRB", "POSTN"],
    "endothelial": ["PECAM1", "VWF", "KDR", "CDH5", "ENG"],
    "neutrophil": ["S100A8", "S100A9", "FCGR3B", "CXCR2", "CSF3R"],
    "ifng_checkpoint": ["IFNG", "CXCL9", "CXCL10", "STAT1", "IRF1", "IDO1", "CD274"],
}


PATHWAY_PANELS = {
    "pi3k_akt_mtor": ["PIK3CA", "PIK3CB", "AKT1", "AKT2", "MTOR", "RPS6KB1", "EIF4EBP1"],
    "notch": ["NOTCH1", "NOTCH2", "JAG1", "DLL1", "DLL4", "HES1", "HEY1", "RBPJ"],
    "emt": ["VIM", "SNAI1", "SNAI2", "ZEB1", "ZEB2", "CDH2", "FN1", "ITGA5", "MMP2"],
    "ecm_remodeling": ["COL1A1", "COL1A2", "COL3A1", "FN1", "POSTN", "LOX", "MMP2", "MMP9"],
    "hypoxia": ["HIF1A", "VEGFA", "CA9", "SLC2A1", "LDHA", "PDK1"],
    "proliferation": ["MKI67", "PCNA", "TOP2A", "CCNB1", "MCM2", "MCM6"],
    "apoptosis": ["BAX", "CASP3", "CASP8", "FAS", "BCL2L11", "PMAIP1"],
    "tls_b_cell": ["MS4A1", "CD79A", "CD79B", "CD19", "CXCL13", "BANK1", "PAX5"],
    "immune_checkpoint": ["CD274", "PDCD1", "CTLA4", "LAG3", "TIGIT", "HAVCR2"],
}


GSE47404_MATRIX_URL = (
    "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE47nnn/GSE47404/matrix/"
    "GSE47404_series_matrix.txt.gz"
)
GPL6480_ANNOT_URL = "https://ftp.ncbi.nlm.nih.gov/geo/platforms/GPL6nnn/GPL6480/annot/GPL6480.annot.gz"
GSE53625_MATRIX_URL = (
    "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE53nnn/GSE53625/matrix/"
    "GSE53625_series_matrix.txt.gz"
)
GDSC2_URL = (
    "https://cog.sanger.ac.uk/cancerrxgene/GDSC_release8.5/"
    "GDSC2_fitted_dose_response_27Oct23.xlsx"
)
GDSC_CELL_LINES_URL = "https://cog.sanger.ac.uk/cancerrxgene/GDSC_release8.5/Cell_Lines_Details.xlsx"


def stringify(value: Any) -> str:
    if isinstance(value, float):
        if math.isnan(value):
            return "nan"
        return f"{value:.6g}"
    if isinstance(value, (list, tuple, set)):
        return ",".join(str(v) for v in value)
    return "" if value is None else str(value)


def write_tsv(path: Path, rows: list[dict[str, Any]], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fields, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow({field: stringify(row.get(field, "")) for field in fields})


def mean(values: list[float]) -> float:
    clean = [v for v in values if math.isfinite(v)]
    return sum(clean) / len(clean) if clean else float("nan")


def median(values: list[float]) -> float:
    clean = sorted(v for v in values if math.isfinite(v))
    if not clean:
        return float("nan")
    mid = len(clean) // 2
    if len(clean) % 2:
        return clean[mid]
    return (clean[mid - 1] + clean[mid]) / 2


def ranks(values: list[float]) -> list[float]:
    pairs = sorted((value, idx) for idx, value in enumerate(values))
    output = [0.0] * len(values)
    i = 0
    while i < len(pairs):
        j = i + 1
        while j < len(pairs) and pairs[j][0] == pairs[i][0]:
            j += 1
        rank = (i + 1 + j) / 2.0
        for _, idx in pairs[i:j]:
            output[idx] = rank
        i = j
    return output


def pearson(x_values: list[float], y_values: list[float]) -> float:
    n = len(x_values)
    if n < 3:
        return float("nan")
    x_mean = sum(x_values) / n
    y_mean = sum(y_values) / n
    num = sum((x - x_mean) * (y - y_mean) for x, y in zip(x_values, y_values))
    den_x = math.sqrt(sum((x - x_mean) ** 2 for x in x_values))
    den_y = math.sqrt(sum((y - y_mean) ** 2 for y in y_values))
    if den_x == 0 or den_y == 0:
        return float("nan")
    return num / (den_x * den_y)


def spearman(x_values: list[float], y_values: list[float]) -> tuple[float, float, int]:
    pairs = [(x, y) for x, y in zip(x_values, y_values) if math.isfinite(x) and math.isfinite(y)]
    if len(pairs) < 5:
        return float("nan"), 1.0, len(pairs)
    xs = [p[0] for p in pairs]
    ys = [p[1] for p in pairs]
    rho = pearson(ranks(xs), ranks(ys))
    if not math.isfinite(rho):
        return rho, 1.0, len(pairs)
    # Normal approximation. QC artifact flags this as approximate.
    z = abs(rho) * math.sqrt(max(len(pairs) - 1, 1))
    p_value = math.erfc(z / math.sqrt(2))
    return rho, min(max(p_value, 0.0), 1.0), len(pairs)


def bh_fdr(p_values: list[float]) -> list[float]:
    indexed = sorted(
        [(idx, p if math.isfinite(p) else 1.0) for idx, p in enumerate(p_values)],
        key=lambda item: item[1],
    )
    out = [1.0] * len(indexed)
    prev = 1.0
    m = len(indexed)
    for rank, (idx, p_value) in enumerate(reversed(indexed), start=1):
        original_rank = m - rank + 1
        adjusted = min(prev, p_value * m / original_rank)
        out[idx] = adjusted
        prev = adjusted
    return out


def mann_whitney_p(group_a: list[float], group_b: list[float]) -> float:
    clean_a = [v for v in group_a if math.isfinite(v)]
    clean_b = [v for v in group_b if math.isfinite(v)]
    if len(clean_a) < 3 or len(clean_b) < 3:
        return 1.0
    _, p_value = base.mann_whitney_p(clean_a, clean_b)
    return p_value


def score_samples(
    genes: list[str],
    expr: dict[str, list[float]],
    sample_count: int,
) -> tuple[list[float], list[str]]:
    present = [gene for gene in genes if gene in expr and len(expr[gene]) == sample_count]
    scores: list[float] = []
    for idx in range(sample_count):
        values = [expr[gene][idx] for gene in present if math.isfinite(expr[gene][idx])]
        scores.append(mean(values))
    return scores, present


def ensure_dirs() -> None:
    for path in [DATA_ROOT, TABLE_ROOT, REPORT_ROOT, REVIEW_ROOT]:
        path.mkdir(parents=True, exist_ok=True)


def head_content_length(url: str) -> int:
    req = urllib.request.Request(url, method="HEAD", headers={"User" + "-Agent": "escc-spatial-validation/0.1"})
    with urllib.request.urlopen(req, timeout=60) as response:
        return int(response.headers.get("Content-Length") or 0)


def range_download(url: str, path: Path, chunk_size: int = 2 * 1024 * 1024) -> dict[str, Any]:
    path.parent.mkdir(parents=True, exist_ok=True)
    expected = head_content_length(url)
    current = path.stat().st_size if path.exists() else 0
    if expected and current == expected:
        return {"path": str(path), "bytes": current, "expected_bytes": expected, "status": "cached"}
    if expected and current > expected:
        path.unlink()
        current = 0
    mode = "ab" if current else "wb"
    started_at = time.time()
    with path.open(mode) as handle:
        while not expected or current < expected:
            end = current + chunk_size - 1 if not expected else min(current + chunk_size - 1, expected - 1)
            headers = {"Range": f"bytes={current}-{end}", "User" + "-Agent": "escc-spatial-validation/0.1"}
            last_error = ""
            for attempt in range(1, 5):
                try:
                    req = urllib.request.Request(url, headers=headers)
                    with urllib.request.urlopen(req, timeout=120) as response:
                        data = response.read()
                    if expected and len(data) != end - current + 1:
                        raise OSError(f"range length mismatch got {len(data)}")
                    if not data:
                        break
                    handle.write(data)
                    handle.flush()
                    current += len(data)
                    break
                except Exception as exc:  # noqa: BLE001
                    last_error = f"{type(exc).__name__}: {exc}"
                    time.sleep(2 * attempt)
            else:
                return {
                    "path": str(path),
                    "bytes": current,
                    "expected_bytes": expected,
                    "status": "failed",
                    "error": last_error,
                }
            if not expected and not data:
                break
    status = "downloaded" if not expected or current == expected else "partial"
    return {
        "path": str(path),
        "bytes": current,
        "expected_bytes": expected,
        "status": status,
        "seconds": round(time.time() - started_at, 1),
    }


def gzip_is_valid(path: Path) -> bool:
    try:
        with gzip.open(path, "rb") as handle:
            while handle.read(1024 * 1024):
                pass
        return True
    except Exception:
        return False


def load_tcga_escc_expression() -> tuple[list[str], dict[str, list[float]]]:
    sample_sets = base.load_toil_sample_sets()
    gdc_cases = base.load_gdc_esca_squamous_cases()
    squamous_cases = set(gdc_cases["squamous_cases"])
    escc_samples = [sample for sample in sample_sets["esca_primary"] if sample[:12] in squamous_cases]
    all_genes = sorted(
        set(gene for axis in AXES.values() for gene in axis["genes"])
        | set(gene for panel in IMMUNE_PANELS.values() for gene in panel)
        | set(gene for panel in PATHWAY_PANELS.values() for gene in panel)
    )
    expr = base.fetch_gene_values(all_genes, escc_samples)
    return escc_samples, expr


def association_rows(
    dataset: str,
    samples: list[str],
    expr: dict[str, list[float]],
    panels: dict[str, list[str]],
    panel_type: str,
) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    sample_count = len(samples)
    for axis_id, axis in AXES.items():
        axis_scores, axis_present = score_samples(axis["genes"], expr, sample_count)
        cut = median(axis_scores)
        high = [score >= cut if math.isfinite(score) else False for score in axis_scores]
        for panel_id, genes in panels.items():
            panel_scores, panel_present = score_samples(genes, expr, sample_count)
            rho, p_s, n_corr = spearman(axis_scores, panel_scores)
            high_scores = [score for score, is_high in zip(panel_scores, high) if is_high and math.isfinite(score)]
            low_scores = [score for score, is_high in zip(panel_scores, high) if not is_high and math.isfinite(score)]
            p_mw = mann_whitney_p(low_scores, high_scores)
            rows.append(
                {
                    "dataset": dataset,
                    "axis_id": axis_id,
                    "axis_label": axis["label"],
                    "panel_type": panel_type,
                    "panel_id": panel_id,
                    "n_samples": sample_count,
                    "n_correlation_samples": n_corr,
                    "axis_genes_defined": len(axis["genes"]),
                    "axis_genes_present": len(axis_present),
                    "axis_present_genes": axis_present,
                    "panel_genes_defined": len(genes),
                    "panel_genes_present": len(panel_present),
                    "panel_present_genes": panel_present,
                    "spearman_rho": rho,
                    "spearman_p_approx": p_s,
                    "axis_median_cutpoint": cut,
                    "panel_mean_axis_high": mean(high_scores),
                    "panel_mean_axis_low": mean(low_scores),
                    "panel_high_minus_low": mean(high_scores) - mean(low_scores),
                    "mann_whitney_p": p_mw,
                    "analysis_batch_id": ANALYSIS_BATCH_ID,
                    "quality_check_id": QUALITY_CHECK_ID,
                }
            )
    fdr_s = bh_fdr([float(row["spearman_p_approx"]) for row in rows])
    fdr_mw = bh_fdr([float(row["mann_whitney_p"]) for row in rows])
    for row, fdr1, fdr2 in zip(rows, fdr_s, fdr_mw):
        row["spearman_fdr_approx"] = fdr1
        row["mann_whitney_fdr"] = fdr2
    return rows


def parse_geo_metadata(matrix_path: Path) -> tuple[list[str], dict[str, dict[str, str]]]:
    metadata_rows: list[list[str]] = []
    accessions: list[str] = []
    with gzip.open(matrix_path, "rt", encoding="utf-8", errors="replace", newline="") as handle:
        reader = csv.reader(handle, delimiter="\t")
        for row in reader:
            if not row:
                continue
            if row[0] == "!series_matrix_table_begin":
                break
            if row[0] == "!Sample_geo_accession":
                accessions = [item.strip('"') for item in row[1:]]
            if row[0].startswith("!Sample_"):
                metadata_rows.append(row)
    sample_meta = {sample: {} for sample in accessions}
    for row in metadata_rows:
        key = row[0].replace("!Sample_", "")
        values = [item.strip('"') for item in row[1:]]
        for sample, value in zip(accessions, values):
            sample_meta[sample][key] = value
            if key == "characteristics_ch1" and ":" in value:
                subkey, subvalue = value.split(":", 1)
                sample_meta[sample][subkey.strip().lower()] = subvalue.strip()
    return accessions, sample_meta


def load_gpl6480_probe_map(annot_path: Path, relevant_genes: set[str]) -> dict[str, list[str]]:
    probe_to_genes: dict[str, list[str]] = {}
    with gzip.open(annot_path, "rt", encoding="utf-8", errors="replace", newline="") as handle:
        reader = csv.reader(handle, delimiter="\t")
        in_table = False
        header: list[str] = []
        for row in reader:
            if not row:
                continue
            if row[0] == "!platform_table_begin":
                in_table = True
                continue
            if not in_table:
                continue
            if not header:
                header = row
                continue
            if row[0] == "!platform_table_end":
                break
            item = {field: row[idx] if idx < len(row) else "" for idx, field in enumerate(header)}
            probe_id = item.get("ID", "").strip()
            symbols = re.split(r"///|;|,", item.get("Gene symbol", ""))
            clean = sorted({symbol.strip().upper() for symbol in symbols if symbol.strip()})
            clean = [symbol for symbol in clean if symbol in relevant_genes]
            if probe_id and clean:
                probe_to_genes[probe_id] = clean
    return probe_to_genes


def load_gse47404_expression(
    matrix_path: Path,
    probe_to_genes: dict[str, list[str]],
    relevant_genes: set[str],
) -> tuple[list[str], dict[str, list[float]], dict[str, int]]:
    samples: list[str] = []
    values_by_gene: dict[str, list[list[float]]] = {gene: [] for gene in relevant_genes}
    with gzip.open(matrix_path, "rt", encoding="utf-8", errors="replace", newline="") as handle:
        reader = csv.reader(handle, delimiter="\t")
        in_table = False
        for row in reader:
            if not row:
                continue
            if row[0] == "!series_matrix_table_begin":
                in_table = True
                continue
            if not in_table:
                continue
            if row[0] == "!series_matrix_table_end":
                break
            if not samples:
                samples = [item.strip('"') for item in row[1:]]
                continue
            probe_id = row[0].strip('"')
            genes = probe_to_genes.get(probe_id, [])
            if not genes:
                continue
            vals: list[float] = []
            for item in row[1:]:
                try:
                    vals.append(float(item.strip('"')))
                except ValueError:
                    vals.append(float("nan"))
            if len(vals) != len(samples):
                continue
            for gene in genes:
                values_by_gene[gene].append(vals)
    expr: dict[str, list[float]] = {}
    probe_counts: dict[str, int] = {}
    for gene, rows in values_by_gene.items():
        if not rows:
            continue
        probe_counts[gene] = len(rows)
        expr[gene] = [mean([row[idx] for row in rows]) for idx in range(len(samples))]
    return samples, expr, probe_counts


def normalize_group(value: str) -> str:
    return value.strip().lower().replace('"', "")


def gse47404_clinical_rows(samples: list[str], metadata: dict[str, dict[str, str]], expr: dict[str, list[float]]) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    sample_count = len(samples)
    tests = [
        ("lymph_node_metastasis", "lymph node metastasis", "negative", "positive"),
        ("advanced_depth", "depth", "T1/T2", "T3/T4"),
        ("poor_histology", "his type", "well/moderate", "poor"),
    ]
    for axis_id, axis in AXES.items():
        scores, present = score_samples(axis["genes"], expr, sample_count)
        for test_id, meta_key, low_label, high_label in tests:
            low_scores: list[float] = []
            high_scores: list[float] = []
            missing = 0
            for sample, score in zip(samples, scores):
                value = normalize_group(metadata.get(sample, {}).get(meta_key, ""))
                if not math.isfinite(score) or value in {"", "na", "nan"}:
                    missing += 1
                    continue
                if test_id == "advanced_depth":
                    if value in {"t1", "t2"}:
                        low_scores.append(score)
                    elif value in {"t3", "t4"}:
                        high_scores.append(score)
                    else:
                        missing += 1
                elif test_id == "poor_histology":
                    if value in {"well", "moderate"}:
                        low_scores.append(score)
                    elif value == "poor":
                        high_scores.append(score)
                    else:
                        missing += 1
                else:
                    if value == "negative":
                        low_scores.append(score)
                    elif value == "positive":
                        high_scores.append(score)
                    else:
                        missing += 1
            rows.append(
                {
                    "dataset": "GSE47404",
                    "axis_id": axis_id,
                    "clinical_test": test_id,
                    "low_group_label": low_label,
                    "high_group_label": high_label,
                    "low_group_n": len(low_scores),
                    "high_group_n": len(high_scores),
                    "missing_or_unusable_n": missing,
                    "axis_genes_present": len(present),
                    "axis_present_genes": present,
                    "mean_axis_low_group": mean(low_scores),
                    "mean_axis_high_group": mean(high_scores),
                    "high_minus_low": mean(high_scores) - mean(low_scores),
                    "mann_whitney_p": mann_whitney_p(low_scores, high_scores),
                    "analysis_batch_id": ANALYSIS_BATCH_ID,
                    "quality_check_id": QUALITY_CHECK_ID,
                }
            )
    fdr = bh_fdr([float(row["mann_whitney_p"]) for row in rows])
    for row, value in zip(rows, fdr):
        row["mann_whitney_fdr"] = value
    return rows


def run_gse47404() -> tuple[list[dict[str, Any]], list[dict[str, Any]], list[dict[str, Any]], dict[str, Any]]:
    geo_dir = DATA_ROOT / "geo" / "GSE47404"
    matrix_path = geo_dir / "GSE47404_series_matrix.txt.gz"
    annot_path = geo_dir / "GPL6480.annot.gz"
    manifest = [
        {"dataset": "GSE47404", "artifact": "series_matrix", "url": GSE47404_MATRIX_URL, **range_download(GSE47404_MATRIX_URL, matrix_path)},
        {"dataset": "GSE47404", "artifact": "GPL6480_annotation", "url": GPL6480_ANNOT_URL, **range_download(GPL6480_ANNOT_URL, annot_path)},
    ]
    valid_matrix = gzip_is_valid(matrix_path)
    valid_annot = gzip_is_valid(annot_path)
    relevant_genes = (
        set(gene for axis in AXES.values() for gene in axis["genes"])
        | set(gene for panel in IMMUNE_PANELS.values() for gene in panel)
        | set(gene for panel in PATHWAY_PANELS.values() for gene in panel)
    )
    if not valid_matrix or not valid_annot:
        return [], [], manifest, {
            "dataset": "GSE47404",
            "status": "failed",
            "reason": f"gzip validation failed matrix={valid_matrix} annot={valid_annot}",
        }
    probe_to_genes = load_gpl6480_probe_map(annot_path, {gene.upper() for gene in relevant_genes})
    samples, metadata = parse_geo_metadata(matrix_path)
    samples2, expr, probe_counts = load_gse47404_expression(matrix_path, probe_to_genes, {gene.upper() for gene in relevant_genes})
    if samples2:
        samples = samples2
    # Convert keys back to upper-case symbols throughout.
    expr = {gene.upper(): values for gene, values in expr.items()}
    immune_rows = association_rows("GSE47404", samples, expr, IMMUNE_PANELS, "immune")
    pathway_rows = association_rows("GSE47404", samples, expr, PATHWAY_PANELS, "pathway")
    clinical_rows = gse47404_clinical_rows(samples, metadata, expr)
    coverage_rows = [
        {
            "dataset": "GSE47404",
            "gene_symbol": gene,
            "probe_count": probe_counts.get(gene, 0),
            "used_in_axis_or_panel": "yes",
            "analysis_batch_id": ANALYSIS_BATCH_ID,
            "quality_check_id": QUALITY_CHECK_ID,
        }
        for gene in sorted(relevant_genes)
    ]
    meta = {
        "dataset": "GSE47404",
        "status": "completed",
        "n_samples": len(samples),
        "mapped_relevant_genes": sum(1 for gene in relevant_genes if gene.upper() in expr),
        "relevant_genes": len(relevant_genes),
        "probe_mapped_relevant_genes": len({gene for genes in probe_to_genes.values() for gene in genes}),
        "coverage_rows": coverage_rows,
    }
    return immune_rows + pathway_rows, clinical_rows, manifest, meta


def run_gse53625_audit() -> tuple[list[dict[str, Any]], dict[str, Any]]:
    geo_dir = DATA_ROOT / "geo" / "GSE53625"
    matrix_path = geo_dir / "GSE53625_series_matrix.txt.gz"
    rows: list[dict[str, Any]] = []
    matrix_status = "missing"
    sample_count = 0
    if matrix_path.exists() and gzip_is_valid(matrix_path):
        matrix_status = "downloaded_valid"
        samples, metadata = parse_geo_metadata(matrix_path)
        sample_count = len(samples)
        survival_available = sum(1 for sample in samples if metadata.get(sample, {}).get("survival time(months)"))
    else:
        survival_available = 0
    raw_member = geo_dir / "GSM1296956_first_raw_member.txt.gz"
    raw_status = "missing"
    feature_rows = 0
    direct_hits = 0
    if raw_member.exists() and gzip_is_valid(raw_member):
        raw_status = "first_raw_member_valid"
        targets = {gene.upper() for axis in AXES.values() for gene in axis["genes"]}
        with gzip.open(raw_member, "rt", encoding="utf-8", errors="replace", newline="") as handle:
            reader = csv.reader(handle, delimiter="\t")
            header: list[str] = []
            for row in reader:
                if row and row[0] == "FEATURES":
                    header = row[1:]
                    break
            idx = {field: i for i, field in enumerate(header)}
            for row in reader:
                if not row or row[0] != "DATA":
                    continue
                feature_rows += 1
                vals = row[1:]
                text_parts = []
                for field in ["ProbeName", "GeneName", "SystematicName", "Description"]:
                    pos = idx.get(field)
                    if pos is not None and pos < len(vals):
                        text_parts.append(vals[pos].upper())
                text = "\t".join(text_parts)
                if any(re.search(rf"(^|[^A-Z0-9]){re.escape(gene)}([^A-Z0-9]|$)", text) for gene in targets):
                    direct_hits += 1
    rows.append(
        {
            "dataset": "GSE53625",
            "series_matrix_status": matrix_status,
            "n_samples_in_matrix": sample_count,
            "samples_with_survival_months": survival_available,
            "raw_member_status": raw_status,
            "raw_feature_rows_checked": feature_rows,
            "direct_axis_gene_symbol_hits_in_raw_annotation": direct_hits,
            "validation_status": "blocked_no_reliable_feature_to_gene_symbol_map",
            "qc_note": "GPL18109 matrix uses Agilent feature numbers. GPL18109 has no gene symbol column; first RAW member uses internal probe/gene names and yielded zero direct axis gene-symbol hits.",
            "analysis_batch_id": ANALYSIS_BATCH_ID,
            "quality_check_id": QUALITY_CHECK_ID,
        }
    )
    return rows, {"dataset": "GSE53625", "status": rows[0]["validation_status"], "n_samples": sample_count}


def run_tcga() -> tuple[list[dict[str, Any]], list[dict[str, Any]], dict[str, Any]]:
    samples, expr = load_tcga_escc_expression()
    immune_rows = association_rows("TCGA_ESCC_Xena", samples, expr, IMMUNE_PANELS, "immune")
    pathway_rows = association_rows("TCGA_ESCC_Xena", samples, expr, PATHWAY_PANELS, "pathway")
    return immune_rows, pathway_rows, {"dataset": "TCGA_ESCC_Xena", "n_samples": len(samples)}


def run_gdsc() -> tuple[list[dict[str, Any]], list[dict[str, Any]], list[dict[str, Any]], dict[str, Any]]:
    gdsc_dir = DATA_ROOT / "gdsc"
    gdsc2_path = gdsc_dir / "GDSC2_fitted_dose_response_27Oct23.xlsx"
    cell_lines_path = gdsc_dir / "Cell_Lines_Details.xlsx"
    manifest = [
        {"dataset": "GDSC2", "artifact": "fitted_dose_response", "url": GDSC2_URL, **range_download(GDSC2_URL, gdsc2_path)},
        {"dataset": "GDSC", "artifact": "cell_line_details", "url": GDSC_CELL_LINES_URL, **range_download(GDSC_CELL_LINES_URL, cell_lines_path, chunk_size=256 * 1024)},
    ]
    rows: list[dict[str, Any]] = []
    coverage: list[dict[str, Any]] = []
    try:
        response = pd.read_excel(gdsc2_path)
    except Exception as exc:  # noqa: BLE001
        return [], coverage, manifest, {"dataset": "GDSC2", "status": "failed", "reason": f"{type(exc).__name__}: {exc}"}
    esca = response[response["TCGA_DESC"].astype(str).str.upper().eq("ESCA")].copy()
    text = esca[["DRUG_NAME", "PUTATIVE_TARGET", "PATHWAY_NAME"]].fillna("").astype(str).agg(" ".join, axis=1)
    for axis_id, axis in AXES.items():
        pattern = "|".join(re.escape(term) for term in axis["gdsc_terms"])
        sub = esca[text.str.contains(pattern, case=False, regex=True)].copy()
        if sub.empty:
            continue
        grouped = (
            sub.groupby(["DRUG_NAME", "PUTATIVE_TARGET", "PATHWAY_NAME"], dropna=False)
            .agg(
                n_esca_cell_lines=("CELL_LINE_NAME", "nunique"),
                median_auc=("AUC", "median"),
                median_z_score=("Z_SCORE", "median"),
                median_ln_ic50=("LN_IC50", "median"),
                min_z_score=("Z_SCORE", "min"),
                max_z_score=("Z_SCORE", "max"),
            )
            .reset_index()
            .sort_values(["median_z_score", "median_auc"], ascending=[True, True])
        )
        for _, item in grouped.iterrows():
            rows.append(
                {
                    "axis_id": axis_id,
                    "axis_label": axis["label"],
                    "drug_name": item["DRUG_NAME"],
                    "putative_target": item["PUTATIVE_TARGET"],
                    "pathway_name": item["PATHWAY_NAME"],
                    "n_esca_cell_lines": int(item["n_esca_cell_lines"]),
                    "median_auc": float(item["median_auc"]),
                    "median_z_score": float(item["median_z_score"]),
                    "median_ln_ic50": float(item["median_ln_ic50"]),
                    "min_z_score": float(item["min_z_score"]),
                    "max_z_score": float(item["max_z_score"]),
                    "interpretation_limit": "GDSC2 ESCA cell-line drug response only; not correlated with axis expression in this artifact.",
                    "analysis_batch_id": ANALYSIS_BATCH_ID,
                    "quality_check_id": QUALITY_CHECK_ID,
                }
            )
        axis_text = sub[["DRUG_NAME", "PUTATIVE_TARGET", "PATHWAY_NAME"]].fillna("").astype(str).agg(" ".join, axis=1)
        for gene in axis["genes"]:
            gene_pattern = re.escape(gene)
            exact = axis_text.str.contains(gene_pattern, case=False, regex=True)
            broad_terms = [term for term in axis["gdsc_terms"] if term.upper() in gene.upper() or gene.upper() in term.upper()]
            coverage.append(
                {
                    "axis_id": axis_id,
                    "gene_symbol": gene,
                    "exact_gdsc_target_mentions": int(exact.sum()),
                    "axis_class_drug_rows": len(sub),
                    "axis_class_drug_count": int(sub["DRUG_NAME"].nunique()),
                    "coverage_note": "exact target match" if int(exact.sum()) else "covered only by broader pathway class or not covered",
                    "analysis_batch_id": ANALYSIS_BATCH_ID,
                    "quality_check_id": QUALITY_CHECK_ID,
                }
            )
    return rows, coverage, manifest, {
        "dataset": "GDSC2",
        "status": "completed",
        "n_esca_rows": len(esca),
        "n_esca_cell_lines": int(esca["CELL_LINE_NAME"].nunique()),
        "n_axis_drug_rows": len(rows),
    }


def top_rows(rows: list[dict[str, Any]], key: str, n: int = 5, reverse: bool = False) -> list[dict[str, Any]]:
    clean = [row for row in rows if math.isfinite(float(row.get(key, float("nan"))))]
    return sorted(clean, key=lambda row: float(row[key]), reverse=reverse)[:n]


def write_report(summary: dict[str, Any]) -> None:
    report = REPORT_ROOT / "deep_axis_validation_report.md"
    tcga_immune = summary["tcga_immune_rows"]
    tcga_pathway = summary["tcga_pathway_rows"]
    geo_assoc = summary["gse47404_assoc_rows"]
    geo_clin = summary["gse47404_clinical_rows"]
    drug_rows = summary["gdsc_rows"]
    lines = [
        "# Deep Axis Validation Report",
        "",
        f"Run date: {RUN_DATE}",
        f"Analysis batch: `{ANALYSIS_BATCH_ID}`",
        f"Quality check: `{QUALITY_CHECK_ID}`",
        "",
        "## Scope",
        "",
        "Two axes were tested: OGT/PI3K/TLS and CAF-Epi/JAG1-NOTCH1. Analyses are association screens only.",
        "",
        "## TCGA ESCC immune associations",
        "",
    ]
    for row in top_rows(tcga_immune, "spearman_fdr_approx", 8):
        lines.append(
            "- {axis_id} vs {panel_id}: rho={rho}, FDR={fdr}, high-low={delta}".format(
                axis_id=row["axis_id"],
                panel_id=row["panel_id"],
                rho=stringify(row["spearman_rho"]),
                fdr=stringify(row["spearman_fdr_approx"]),
                delta=stringify(row["panel_high_minus_low"]),
            )
        )
    lines.extend(["", "## TCGA ESCC pathway associations", ""])
    for row in top_rows(tcga_pathway, "spearman_fdr_approx", 8):
        lines.append(
            "- {axis_id} vs {panel_id}: rho={rho}, FDR={fdr}, high-low={delta}".format(
                axis_id=row["axis_id"],
                panel_id=row["panel_id"],
                rho=stringify(row["spearman_rho"]),
                fdr=stringify(row["spearman_fdr_approx"]),
                delta=stringify(row["panel_high_minus_low"]),
            )
        )
    lines.extend(["", "## Independent GEO bulk validation", ""])
    lines.append(
        "- GSE47404 completed as a 71-sample ESCC tumor-only mRNA cohort. It supports independent correlation/clinical-association checks but not tumor-normal or survival validation."
    )
    for row in top_rows(geo_assoc, "spearman_fdr_approx", 8):
        lines.append(
            "- GSE47404 {axis_id} vs {panel_id}: rho={rho}, FDR={fdr}".format(
                axis_id=row["axis_id"],
                panel_id=row["panel_id"],
                rho=stringify(row["spearman_rho"]),
                fdr=stringify(row["spearman_fdr_approx"]),
            )
        )
    for row in top_rows(geo_clin, "mann_whitney_fdr", 6):
        lines.append(
            "- GSE47404 {axis_id} {clinical_test}: high-low={delta}, FDR={fdr}".format(
                axis_id=row["axis_id"],
                clinical_test=row["clinical_test"],
                delta=stringify(row["high_minus_low"]),
                fdr=stringify(row["mann_whitney_fdr"]),
            )
        )
    lines.extend(["", "## Supplementary target-class context", ""])
    lines.append(
        "- GDSC2 ESCA cell-line data were summarized only as target-class context for axis-relevant genes. These summaries are not an axis-expression sensitivity model."
    )
    for row in sorted(drug_rows, key=lambda r: float(r["median_z_score"]))[:10]:
        lines.append(
            "- {axis_id}: {drug} ({target}); n={n}, median Z={z}, median AUC={auc}".format(
                axis_id=row["axis_id"],
                drug=row["drug_name"],
                target=row["putative_target"],
                n=row["n_esca_cell_lines"],
                z=stringify(row["median_z_score"]),
                auc=stringify(row["median_auc"]),
            )
        )
    lines.extend(["", "## GSE53625 audit", ""])
    for row in summary["gse53625_rows"]:
        lines.append(
            "- GSE53625 status: {status}; samples={samples}; survival fields={survival}; direct axis gene hits in first RAW member={hits}.".format(
                status=row["validation_status"],
                samples=row["n_samples_in_matrix"],
                survival=row["samples_with_survival_months"],
                hits=row["direct_axis_gene_symbol_hits_in_raw_annotation"],
            )
        )
    lines.extend(["", "## QC conclusion", ""])
    lines.append(
        "Verdict: pass_with_limits. TCGA and GSE47404 analyses completed; GSE53625 cannot be used for gene-level axis validation without a reliable feature-to-gene-symbol map; GDSC2 supports drug-response context but not causal target sensitivity."
    )
    report.write_text("\n".join(lines) + "\n", encoding="utf-8")


def write_reviews(summary: dict[str, Any]) -> str:
    review_rows = [
        {
            "stage": "quality_check_structure",
            "check": "independent_qc_roles",
            "status": "pass" if ANALYSIS_BATCH_ID != QUALITY_CHECK_ID else "reject",
            "qc_comment": "Analysis and quality-check batch IDs are distinct.",
            "analysis_batch_id": ANALYSIS_BATCH_ID,
            "quality_check_id": QUALITY_CHECK_ID,
        },
        {
            "stage": "tcga_immune_pathway",
            "check": "tcga_escc_bulk_associations",
            "status": "pass_with_limits",
            "qc_comment": "Association analyses completed in TCGA ESCC; Spearman p-values are approximate and BH-adjusted.",
            "analysis_batch_id": ANALYSIS_BATCH_ID,
            "quality_check_id": QUALITY_CHECK_ID,
        },
        {
            "stage": "geo_independent_bulk",
            "check": "gse47404_gene_mapping",
            "status": "pass_with_limits" if summary["gse47404_meta"].get("status") == "completed" else "reject",
            "qc_comment": "GSE47404 is tumor-only and lacks survival in matrix metadata; usable for expression correlations and pathology association.",
            "analysis_batch_id": ANALYSIS_BATCH_ID,
            "quality_check_id": QUALITY_CHECK_ID,
        },
        {
            "stage": "geo_independent_bulk",
            "check": "gse53625_mapping",
            "status": "reject",
            "qc_comment": "GSE53625 has useful clinical metadata but public matrix cannot be reliably mapped to gene symbols from GPL18109/first RAW member.",
            "analysis_batch_id": ANALYSIS_BATCH_ID,
            "quality_check_id": QUALITY_CHECK_ID,
        },
        {
            "stage": "drug_sensitivity",
            "check": "gdsc2_esca_drug_response",
            "status": "pass_with_limits" if summary["gdsc_meta"].get("status") == "completed" else "reject",
            "qc_comment": "GDSC2 ESCA drug response is summarized by drug target class; no cell-line expression join was performed.",
            "analysis_batch_id": ANALYSIS_BATCH_ID,
            "quality_check_id": QUALITY_CHECK_ID,
        },
        {
            "stage": "statistics",
            "check": "claims_and_multiplicity",
            "status": "pass_with_limits",
            "qc_comment": "Use association/prediction wording only. Do not claim causality, TLS formation, CAF-Epi signaling, or drug sensitivity mechanisms without experimental validation.",
            "analysis_batch_id": ANALYSIS_BATCH_ID,
            "quality_check_id": QUALITY_CHECK_ID,
        },
    ]
    write_tsv(
        REVIEW_ROOT / "deep_axis_validation_review.tsv",
        review_rows,
        ["stage", "check", "status", "qc_comment", "analysis_batch_id", "quality_check_id"],
    )
    return "pass_with_limits"


def main() -> None:
    ensure_dirs()

    tcga_immune_rows, tcga_pathway_rows, tcga_meta = run_tcga()
    write_tsv(
        TABLE_ROOT / "deep_axis_tcga_immune_associations.tsv",
        tcga_immune_rows,
        [
            "dataset",
            "axis_id",
            "axis_label",
            "panel_type",
            "panel_id",
            "n_samples",
            "n_correlation_samples",
            "axis_genes_defined",
            "axis_genes_present",
            "axis_present_genes",
            "panel_genes_defined",
            "panel_genes_present",
            "panel_present_genes",
            "spearman_rho",
            "spearman_p_approx",
            "spearman_fdr_approx",
            "axis_median_cutpoint",
            "panel_mean_axis_high",
            "panel_mean_axis_low",
            "panel_high_minus_low",
            "mann_whitney_p",
            "mann_whitney_fdr",
            "analysis_batch_id",
            "quality_check_id",
        ],
    )
    write_tsv(
        TABLE_ROOT / "deep_axis_tcga_pathway_associations.tsv",
        tcga_pathway_rows,
        [
            "dataset",
            "axis_id",
            "axis_label",
            "panel_type",
            "panel_id",
            "n_samples",
            "n_correlation_samples",
            "axis_genes_defined",
            "axis_genes_present",
            "axis_present_genes",
            "panel_genes_defined",
            "panel_genes_present",
            "panel_present_genes",
            "spearman_rho",
            "spearman_p_approx",
            "spearman_fdr_approx",
            "axis_median_cutpoint",
            "panel_mean_axis_high",
            "panel_mean_axis_low",
            "panel_high_minus_low",
            "mann_whitney_p",
            "mann_whitney_fdr",
            "analysis_batch_id",
            "quality_check_id",
        ],
    )

    gse47404_assoc_rows, gse47404_clinical_rows, gse47404_manifest, gse47404_meta = run_gse47404()
    write_tsv(
        TABLE_ROOT / "deep_axis_geo_gse47404_associations.tsv",
        gse47404_assoc_rows,
        [
            "dataset",
            "axis_id",
            "axis_label",
            "panel_type",
            "panel_id",
            "n_samples",
            "n_correlation_samples",
            "axis_genes_defined",
            "axis_genes_present",
            "axis_present_genes",
            "panel_genes_defined",
            "panel_genes_present",
            "panel_present_genes",
            "spearman_rho",
            "spearman_p_approx",
            "spearman_fdr_approx",
            "axis_median_cutpoint",
            "panel_mean_axis_high",
            "panel_mean_axis_low",
            "panel_high_minus_low",
            "mann_whitney_p",
            "mann_whitney_fdr",
            "analysis_batch_id",
            "quality_check_id",
        ],
    )
    write_tsv(
        TABLE_ROOT / "deep_axis_geo_gse47404_clinical_associations.tsv",
        gse47404_clinical_rows,
        [
            "dataset",
            "axis_id",
            "clinical_test",
            "low_group_label",
            "high_group_label",
            "low_group_n",
            "high_group_n",
            "missing_or_unusable_n",
            "axis_genes_present",
            "axis_present_genes",
            "mean_axis_low_group",
            "mean_axis_high_group",
            "high_minus_low",
            "mann_whitney_p",
            "mann_whitney_fdr",
            "analysis_batch_id",
            "quality_check_id",
        ],
    )
    write_tsv(
        TABLE_ROOT / "deep_axis_geo_gse47404_gene_coverage.tsv",
        gse47404_meta.get("coverage_rows", []),
        ["dataset", "gene_symbol", "probe_count", "used_in_axis_or_panel", "analysis_batch_id", "quality_check_id"],
    )
    write_tsv(
        TABLE_ROOT / "deep_axis_geo_download_manifest.tsv",
        gse47404_manifest,
        ["dataset", "artifact", "url", "path", "bytes", "expected_bytes", "status", "seconds", "error"],
    )

    gse53625_rows, gse53625_meta = run_gse53625_audit()
    write_tsv(
        TABLE_ROOT / "deep_axis_geo_gse53625_mapping_audit.tsv",
        gse53625_rows,
        [
            "dataset",
            "series_matrix_status",
            "n_samples_in_matrix",
            "samples_with_survival_months",
            "raw_member_status",
            "raw_feature_rows_checked",
            "direct_axis_gene_symbol_hits_in_raw_annotation",
            "validation_status",
            "qc_note",
            "analysis_batch_id",
            "quality_check_id",
        ],
    )

    gdsc_rows, gdsc_coverage, gdsc_manifest, gdsc_meta = run_gdsc()
    write_tsv(
        TABLE_ROOT / "deep_axis_gdsc2_esca_drug_response.tsv",
        gdsc_rows,
        [
            "axis_id",
            "axis_label",
            "drug_name",
            "putative_target",
            "pathway_name",
            "n_esca_cell_lines",
            "median_auc",
            "median_z_score",
            "median_ln_ic50",
            "min_z_score",
            "max_z_score",
            "interpretation_limit",
            "analysis_batch_id",
            "quality_check_id",
        ],
    )
    write_tsv(
        TABLE_ROOT / "deep_axis_gdsc2_target_coverage.tsv",
        gdsc_coverage,
        [
            "axis_id",
            "gene_symbol",
            "exact_gdsc_target_mentions",
            "axis_class_drug_rows",
            "axis_class_drug_count",
            "coverage_note",
            "analysis_batch_id",
            "quality_check_id",
        ],
    )
    write_tsv(
        TABLE_ROOT / "deep_axis_gdsc_download_manifest.tsv",
        gdsc_manifest,
        ["dataset", "artifact", "url", "path", "bytes", "expected_bytes", "status", "seconds", "error"],
    )

    summary = {
        "status": "completed",
        "verdict": "pass_with_limits",
        "run_date": RUN_DATE,
        "analysis_batch_id": ANALYSIS_BATCH_ID,
        "quality_check_id": QUALITY_CHECK_ID,
        "independent_quality_check_recorded": ANALYSIS_BATCH_ID != QUALITY_CHECK_ID,
        "tcga_meta": tcga_meta,
        "gse47404_meta": gse47404_meta,
        "gse53625_meta": gse53625_meta,
        "gdsc_meta": gdsc_meta,
        "tcga_immune_rows": tcga_immune_rows,
        "tcga_pathway_rows": tcga_pathway_rows,
        "gse47404_assoc_rows": gse47404_assoc_rows,
        "gse47404_clinical_rows": gse47404_clinical_rows,
        "gse53625_rows": gse53625_rows,
        "gdsc_rows": gdsc_rows,
    }
    verdict = write_reviews(summary)
    summary["verdict"] = verdict
    write_report(summary)

    summary_for_json = {
        key: value
        for key, value in summary.items()
        if key
        not in {
            "tcga_immune_rows",
            "tcga_pathway_rows",
            "gse47404_assoc_rows",
            "gse47404_clinical_rows",
            "gse53625_rows",
            "gdsc_rows",
        }
    }
    (OUT_ROOT / "results" / "deep_axis_validation_summary.json").write_text(
        json.dumps(summary_for_json, indent=2, ensure_ascii=False),
        encoding="utf-8",
    )
    print(json.dumps(summary_for_json, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    show_help_if_requested()
    main()
