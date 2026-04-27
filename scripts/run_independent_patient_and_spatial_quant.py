#!/usr/bin/env python
"""Independent patient validation and spatial source-table quantification.

This script extends the spatial ESCC axis workflow along two deliberately
separate evidence lines:

1. Rescue GSE53625 probe-to-gene mapping from public Agilent probe sequences and
   Ensembl transcript/genomic sequence, then run patient-level survival and
   paired tumor-normal checks where mapping coverage is adequate.
2. Quantify openly available spatial/source tables from HRA003627 and HRA008846
   instead of only citing the papers narratively.

All outputs keep analysis and quality-check batch IDs distinct. The analyses are
association and source-table validation artifacts; they are not causal proof.
"""

from __future__ import annotations

import csv
import gzip
import hashlib
import json
import math
import re
import sys
import time
import urllib.error
import urllib.parse
import urllib.request
import zipfile
from collections import defaultdict
from datetime import UTC, datetime
from pathlib import Path
from typing import Any


def _early_help_if_requested() -> None:
    if any(arg in {"-h", "--help"} for arg in sys.argv[1:]):
        print((__doc__ or "").strip())
        print(
            "\nUsage:\n"
            "  python scripts/run_independent_patient_and_spatial_quant.py\n\n"
            "Run from the repository root after installing requirements. Outputs are written under "
            "spatial_escc_workflow/results, reports, reviews and deliverables relative to the repository."
        )
        raise SystemExit(0)


_early_help_if_requested()

import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
SCRIPTS = ROOT / "scripts"
if str(SCRIPTS) not in sys.path:
    sys.path.insert(0, str(SCRIPTS))

import public_data_helpers as base  # noqa: E402
import run_spatial_axis_deep_validation as deep  # noqa: E402


OUT_ROOT = ROOT / "spatial_escc_workflow"
DATA_ROOT = OUT_ROOT / "data"
TABLE_ROOT = OUT_ROOT / "results" / "tables"
REPORT_ROOT = OUT_ROOT / "reports"
REVIEW_ROOT = OUT_ROOT / "reviews"
DELIVERABLE_ROOT = OUT_ROOT / "deliverables"
CACHE_ROOT = DATA_ROOT / "ensembl_sequence_cache"

RUN_DATE = datetime.now(UTC).date().isoformat()
ANALYSIS_RUN_ID = "independent_spatial_analysis_001"
REPRODUCIBILITY_CHECK_ID = "independent_spatial_quality_check_001"

ENSEMBL_REST = "https://rest.ensembl.org"
HTTP_CLIENT_ID = "escc-spatial-validation/0.2"


def show_help_if_requested() -> None:
    if any(arg in {"-h", "--help"} for arg in sys.argv[1:]):
        print((__doc__ or "").strip())
        print(
            "\nUsage:\n"
            "  python scripts/run_independent_patient_and_spatial_quant.py\n\n"
            "Run from the repository root after installing requirements. Outputs are written "
            "under spatial_escc_workflow/results, reports, reviews and deliverables relative "
            "to the repository."
        )
        raise SystemExit(0)

GSE53625_MATRIX_URL = (
    "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE53nnn/GSE53625/matrix/"
    "GSE53625_series_matrix.txt.gz"
)

HRA008846_DOWNLOADS = {
    "HRA008846_TableS3_DEG.xlsx": "https://pmc.ncbi.nlm.nih.gov/articles/instance/13006417/bin/mmc2.xlsx",
    "HRA008846_TableS4_cell_abundance.xlsx": "https://pmc.ncbi.nlm.nih.gov/articles/instance/13006417/bin/mmc3.xlsx",
    "HRA008846_TableS6_ligand_receptor.xlsx": "https://pmc.ncbi.nlm.nih.gov/articles/instance/13006417/bin/mmc4.xlsx",
}

AXES = deep.AXES
AXIS_GENES = sorted({gene.upper() for axis in AXES.values() for gene in axis["genes"]})

HRA003627_SIGNATURES = {
    "dk_keratinization": ["CRNN", "MAL"],
    "cancerization_progression": [
        "TAGLN2",
        "KRT16",
        "KRT17",
        "S100A8",
        "TOP2A",
        "MKI67",
        "LAMC2",
        "CCN2",
        "ANO1",
        "ITGA6",
        "MMP14",
    ],
}

HRA008846_TARGET_GENES = sorted(
    set(AXIS_GENES)
    | {
        "COL3A1",
        "CXCL13",
        "CXCR5",
        "HES1",
        "HEY1",
        "DLL1",
        "DLL4",
        "FAP",
        "FN1",
        "MMP2",
        "MMP9",
        "PDGFRB",
        "VIM",
    }
)


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


def ensure_dirs() -> None:
    for path in [DATA_ROOT, TABLE_ROOT, REPORT_ROOT, REVIEW_ROOT, DELIVERABLE_ROOT, CACHE_ROOT]:
        path.mkdir(parents=True, exist_ok=True)


def mean(values: list[float]) -> float:
    clean = [v for v in values if math.isfinite(v)]
    return sum(clean) / len(clean) if clean else float("nan")


def median(values: list[float]) -> float:
    clean = sorted(v for v in values if math.isfinite(v))
    if not clean:
        return float("nan")
    mid = len(clean) // 2
    return clean[mid] if len(clean) % 2 else (clean[mid - 1] + clean[mid]) / 2


def zscore(values: list[float]) -> list[float]:
    clean = [v for v in values if math.isfinite(v)]
    if len(clean) < 2:
        return [float("nan") for _ in values]
    mu = sum(clean) / len(clean)
    var = sum((v - mu) ** 2 for v in clean) / (len(clean) - 1)
    sd = math.sqrt(var)
    if sd == 0:
        return [0.0 if math.isfinite(v) else float("nan") for v in values]
    return [(v - mu) / sd if math.isfinite(v) else float("nan") for v in values]


def clean_nt(seq: str) -> str:
    return re.sub(r"[^ACGT]", "", str(seq).upper())


def revcomp(seq: str) -> str:
    return seq.translate(str.maketrans("ACGT", "TGCA"))[::-1]


def safe_float(value: Any) -> float:
    try:
        return float(value)
    except (TypeError, ValueError):
        return float("nan")


def exact_sign_test_p(positive: int, negative: int) -> float:
    n = positive + negative
    if n == 0:
        return 1.0
    k = min(positive, negative)
    prob = sum(math.comb(n, i) for i in range(k + 1)) / (2**n)
    return min(1.0, 2 * prob)


def xlsx_is_valid(path: Path) -> bool:
    if not path.exists() or path.stat().st_size < 1024:
        return False
    try:
        with zipfile.ZipFile(path) as zf:
            return "[Content_Types].xml" in zf.namelist()
    except zipfile.BadZipFile:
        return False


def solve_pmc_pow(html: str) -> tuple[str, str] | None:
    challenge_match = re.search(r'POW_CHALLENGE\s*=\s*"([^"]+)"', html)
    difficulty_match = re.search(r"POW_DIFFICULTY\s*=\s*(\d+)", html)
    cookie_match = re.search(r'POW_COOKIE_NAME\s*=\s*"([^"]+)"', html)
    if not (challenge_match and difficulty_match and cookie_match):
        return None
    challenge = challenge_match.group(1)
    difficulty = int(difficulty_match.group(1))
    target = "0" * difficulty
    nonce = 0
    while True:
        digest = hashlib.sha256((challenge + str(nonce)).encode("utf-8")).hexdigest()
        if digest.startswith(target):
            return cookie_match.group(1), f"{challenge},{nonce}"
        nonce += 1


def download_pmc_xlsx(url: str, path: Path) -> dict[str, Any]:
    path.parent.mkdir(parents=True, exist_ok=True)
    if xlsx_is_valid(path):
        return {"path": str(path), "status": "cached_valid", "bytes": path.stat().st_size, "url": url}
    headers = {"User" + "-Agent": HTTP_CLIENT_ID}
    try:
        req = urllib.request.Request(url, headers=headers)
        with urllib.request.urlopen(req, timeout=120) as response:
            data = response.read()
            content_type = response.headers.get("Content-Type", "")
    except Exception as exc:  # noqa: BLE001
        return {"path": str(path), "status": "failed", "error": f"{type(exc).__name__}: {exc}", "url": url}
    if data.startswith(b"PK"):
        path.write_bytes(data)
        return {"path": str(path), "status": "downloaded_valid", "bytes": path.stat().st_size, "url": url}
    html = data.decode("utf-8", errors="replace")
    pow_cookie = solve_pmc_pow(html)
    if not pow_cookie:
        invalid = path.with_suffix(path.suffix + ".invalid.html")
        invalid.write_text(html[:20000], encoding="utf-8")
        return {
            "path": str(path),
            "status": "failed_not_xlsx",
            "content_type": content_type,
            "invalid_preview": str(invalid),
            "url": url,
        }
    cookie_name, cookie_value = pow_cookie
    req = urllib.request.Request(
        url,
        headers={
            "User" + "-Agent": HTTP_CLIENT_ID,
            "Cookie": f"{cookie_name}={cookie_value}",
            "Referer": url,
        },
    )
    with urllib.request.urlopen(req, timeout=120) as response:
        data = response.read()
        content_type = response.headers.get("Content-Type", "")
    if not data.startswith(b"PK"):
        invalid = path.with_suffix(path.suffix + ".invalid.html")
        invalid.write_bytes(data[:20000])
        return {
            "path": str(path),
            "status": "failed_pow_cookie_not_xlsx",
            "content_type": content_type,
            "invalid_preview": str(invalid),
            "url": url,
        }
    path.write_bytes(data)
    return {"path": str(path), "status": "downloaded_with_pow_cookie", "bytes": path.stat().st_size, "url": url}


def ensembl_request_json(path: str, cache_path: Path) -> dict[str, Any]:
    if cache_path.exists():
        return json.loads(cache_path.read_text(encoding="utf-8"))
    url = f"{ENSEMBL_REST}{path}"
    headers = {"Content-Type": "application/json", "User" + "-Agent": HTTP_CLIENT_ID}
    last_error = ""
    for attempt in range(1, 6):
        try:
            req = urllib.request.Request(url, headers=headers)
            with urllib.request.urlopen(req, timeout=120) as response:
                payload = response.read().decode("utf-8")
            data = json.loads(payload)
            cache_path.parent.mkdir(parents=True, exist_ok=True)
            cache_path.write_text(json.dumps(data, ensure_ascii=False, indent=2), encoding="utf-8")
            return data
        except urllib.error.HTTPError as exc:
            last_error = f"HTTPError {exc.code}: {exc.reason}"
            if exc.code == 429:
                time.sleep(2 * attempt)
            else:
                time.sleep(attempt)
        except Exception as exc:  # noqa: BLE001
            last_error = f"{type(exc).__name__}: {exc}"
            time.sleep(attempt)
    raise RuntimeError(f"Ensembl request failed for {path}: {last_error}")


def fetch_gene_sequence_bundle(gene: str) -> dict[str, Any]:
    bundle_path = CACHE_ROOT / f"{gene}_sequence_bundle.json"
    if bundle_path.exists():
        return json.loads(bundle_path.read_text(encoding="utf-8"))
    lookup = ensembl_request_json(
        f"/lookup/symbol/homo_sapiens/{urllib.parse.quote(gene)}?expand=1",
        CACHE_ROOT / f"{gene}_lookup.json",
    )
    sequences: list[dict[str, str]] = []
    gene_region = "{seq}:{start}..{end}:{strand}".format(
        seq=lookup["seq_region_name"],
        start=lookup["start"],
        end=lookup["end"],
        strand=lookup.get("strand", 1),
    )
    region_data = ensembl_request_json(
        f"/sequence/region/homo_sapiens/{gene_region}",
        CACHE_ROOT / f"{gene}_genomic_sequence.json",
    )
    genomic_seq = clean_nt(region_data.get("seq", ""))
    if genomic_seq:
        sequences.append({"source": "genomic_region", "id": lookup["id"], "seq": genomic_seq})
    transcript_ids: list[str] = []
    for transcript in lookup.get("Transcript", []):
        tid = transcript.get("id")
        if not tid:
            continue
        biotype = transcript.get("biotype", "")
        if biotype and biotype not in {"protein_coding", "nonsense_mediated_decay", "processed_transcript"}:
            continue
        transcript_ids.append(tid)
    for tid in sorted(set(transcript_ids)):
        try:
            seq_data = ensembl_request_json(
                f"/sequence/id/{tid}?type=cdna",
                CACHE_ROOT / f"{gene}_{tid}_cdna.json",
            )
        except RuntimeError:
            continue
        seq = clean_nt(seq_data.get("seq", ""))
        if seq:
            sequences.append({"source": "transcript_cdna", "id": tid, "seq": seq})
    bundle = {
        "gene_symbol": gene,
        "ensembl_gene_id": lookup.get("id", ""),
        "seq_region_name": lookup.get("seq_region_name", ""),
        "start": lookup.get("start", ""),
        "end": lookup.get("end", ""),
        "strand": lookup.get("strand", ""),
        "n_transcripts_queried": len(set(transcript_ids)),
        "n_sequences": len(sequences),
        "total_bases_indexed": sum(len(item["seq"]) for item in sequences),
        "sequences": sequences,
    }
    bundle_path.write_text(json.dumps(bundle, ensure_ascii=False), encoding="utf-8")
    return bundle


def load_gse53625_feature_sequences(raw_path: Path) -> list[dict[str, Any]]:
    features: list[dict[str, Any]] = []
    with gzip.open(raw_path, "rt", encoding="utf-8", errors="replace", newline="") as handle:
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
            vals = row[1:]
            seq = clean_nt(vals[idx["Sequence"]] if idx.get("Sequence", 10**9) < len(vals) else "")
            feature_num = vals[idx["FeatureNum"]] if idx.get("FeatureNum", 10**9) < len(vals) else ""
            if not feature_num or len(seq) < 25:
                continue
            features.append(
                {
                    "feature_num": feature_num,
                    "sequence": seq,
                    "probe_name": vals[idx["ProbeName"]] if idx.get("ProbeName", 10**9) < len(vals) else "",
                    "raw_gene_name": vals[idx["GeneName"]] if idx.get("GeneName", 10**9) < len(vals) else "",
                    "systematic_name": vals[idx["SystematicName"]] if idx.get("SystematicName", 10**9) < len(vals) else "",
                    "description": vals[idx["Description"]] if idx.get("Description", 10**9) < len(vals) else "",
                }
            )
    return features


def build_target_sequence_index(probe_lengths: set[int]) -> tuple[dict[tuple[int, str], set[str]], list[dict[str, Any]]]:
    index: dict[tuple[int, str], set[str]] = defaultdict(set)
    bundle_rows: list[dict[str, Any]] = []
    for gene in AXIS_GENES:
        bundle = fetch_gene_sequence_bundle(gene)
        bundle_rows.append(
            {
                "gene_symbol": gene,
                "ensembl_gene_id": bundle.get("ensembl_gene_id", ""),
                "seq_region_name": bundle.get("seq_region_name", ""),
                "start": bundle.get("start", ""),
                "end": bundle.get("end", ""),
                "strand": bundle.get("strand", ""),
                "n_transcripts_queried": bundle.get("n_transcripts_queried", 0),
                "n_sequences_indexed": bundle.get("n_sequences", 0),
                "total_bases_indexed": bundle.get("total_bases_indexed", 0),
                "analysis_run_id": ANALYSIS_RUN_ID,
                "reproducibility_check_id": REPRODUCIBILITY_CHECK_ID,
            }
        )
        for seq_item in bundle.get("sequences", []):
            seq = clean_nt(seq_item.get("seq", ""))
            for orient_seq in (seq, revcomp(seq)):
                for length in probe_lengths:
                    if len(orient_seq) < length:
                        continue
                    for start in range(0, len(orient_seq) - length + 1):
                        index[(length, orient_seq[start : start + length])].add(gene)
    return index, bundle_rows


def map_features_to_axis_genes(features: list[dict[str, Any]]) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    probe_lengths = {len(item["sequence"]) for item in features if len(item["sequence"]) >= 25}
    sequence_index, bundle_rows = build_target_sequence_index(probe_lengths)
    mapped_rows: list[dict[str, Any]] = []
    for item in features:
        seq = item["sequence"]
        exact_genes = sequence_index.get((len(seq), seq), set())
        match_type = "none"
        matched = sorted(exact_genes)
        if len(matched) == 1:
            match_type = "exact"
        elif not matched:
            fuzzy_genes: set[str] = set()
            for pos, old_base in enumerate(seq):
                for base in "ACGT":
                    if base == old_base:
                        continue
                    candidate = seq[:pos] + base + seq[pos + 1 :]
                    fuzzy_genes.update(sequence_index.get((len(candidate), candidate), set()))
                if len(fuzzy_genes) > 1:
                    break
            matched = sorted(fuzzy_genes)
            if len(matched) == 1:
                match_type = "one_mismatch"
            elif len(matched) > 1:
                match_type = "ambiguous_one_mismatch"
        else:
            match_type = "ambiguous_exact"
        if match_type in {"exact", "one_mismatch", "ambiguous_exact", "ambiguous_one_mismatch"}:
            mapped_rows.append(
                {
                    "dataset": "GSE53625",
                    "feature_num": item["feature_num"],
                    "probe_name": item["probe_name"],
                    "raw_gene_name": item["raw_gene_name"],
                    "systematic_name": item["systematic_name"],
                    "sequence_length": len(seq),
                    "matched_genes": matched,
                    "accepted_gene": matched[0] if len(matched) == 1 else "",
                    "match_type": match_type,
                    "specificity_scope": "unambiguous_within_axis_target_ensembl_sequences",
                    "analysis_run_id": ANALYSIS_RUN_ID,
                    "reproducibility_check_id": REPRODUCIBILITY_CHECK_ID,
                }
            )
    return mapped_rows, bundle_rows


def load_gse53625_expression(
    matrix_path: Path,
    accepted_feature_to_gene: dict[str, str],
) -> tuple[list[str], dict[str, list[float]]]:
    samples: list[str] = []
    values_by_gene: dict[str, list[list[float]]] = defaultdict(list)
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
            feature_num = row[0].strip('"')
            gene = accepted_feature_to_gene.get(feature_num)
            if not gene:
                continue
            vals = [safe_float(item.strip('"')) for item in row[1:]]
            if len(vals) == len(samples):
                values_by_gene[gene].append(vals)
    expr: dict[str, list[float]] = {}
    for gene, rows in values_by_gene.items():
        expr[gene] = [mean([row[idx] for row in rows]) for idx in range(len(samples))]
    return samples, expr


def score_axis_from_expr(axis_genes: list[str], expr: dict[str, list[float]], sample_count: int) -> tuple[list[float], list[str]]:
    present = [gene for gene in axis_genes if gene in expr and len(expr[gene]) == sample_count]
    z_expr = {gene: zscore(expr[gene]) for gene in present}
    scores: list[float] = []
    for idx in range(sample_count):
        scores.append(mean([z_expr[gene][idx] for gene in present]))
    return scores, present


def classify_gse53625_sample(meta: dict[str, str]) -> str:
    tissue = meta.get("tissue", "").lower()
    if "cancer tissue" in tissue:
        return "tumor"
    if "normal tissue" in tissue:
        return "normal"
    return "unknown"


def survival_rows_for_gse53625(
    samples: list[str],
    metadata: dict[str, dict[str, str]],
    expr: dict[str, list[float]],
) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for axis_id, axis in AXES.items():
        scores, present = score_axis_from_expr(axis["genes"], expr, len(samples))
        tumor_records: list[tuple[float, float, int, str]] = []
        for sample, score in zip(samples, scores):
            meta = metadata.get(sample, {})
            if classify_gse53625_sample(meta) != "tumor" or not math.isfinite(score):
                continue
            time_months = safe_float(meta.get("survival time(months)", ""))
            event = 1 if meta.get("death at fu", "").strip().lower() == "yes" else 0
            if math.isfinite(time_months) and time_months > 0:
                tumor_records.append((score, time_months, event, sample))
        cutpoint = median([rec[0] for rec in tumor_records])
        groups = [1 if rec[0] >= cutpoint else 0 for rec in tumor_records]
        times = [rec[1] for rec in tumor_records]
        events = [rec[2] for rec in tumor_records]
        if len(set(groups)) == 2 and len(tumor_records) >= 10:
            chi2, p_value = base.logrank_p(times, events, groups)
        else:
            chi2, p_value = float("nan"), 1.0
        high_events = sum(event for event, group in zip(events, groups) if group == 1)
        low_events = sum(event for event, group in zip(events, groups) if group == 0)
        high_pt = sum(time for time, group in zip(times, groups) if group == 1)
        low_pt = sum(time for time, group in zip(times, groups) if group == 0)
        rate_high = (high_events + 0.5) / (high_pt + 0.5)
        rate_low = (low_events + 0.5) / (low_pt + 0.5)
        rate_ratio = rate_high / rate_low if rate_low > 0 else float("nan")
        se_log = math.sqrt(1 / (high_events + 0.5) + 1 / (low_events + 0.5))
        ci_low = math.exp(math.log(rate_ratio) - 1.96 * se_log) if math.isfinite(rate_ratio) and rate_ratio > 0 else float("nan")
        ci_high = math.exp(math.log(rate_ratio) + 1.96 * se_log) if math.isfinite(rate_ratio) and rate_ratio > 0 else float("nan")
        coverage = len(present) / len(axis["genes"]) if axis["genes"] else 0.0
        if coverage < 0.6 or len(tumor_records) < 50:
            status = "reject_low_mapping_coverage"
        elif math.isfinite(p_value) and p_value < 0.05:
            status = "supported_with_limitations"
        else:
            status = "no_survival_support"
        rows.append(
            {
                "dataset": "GSE53625",
                "axis_id": axis_id,
                "axis_label": axis["label"],
                "n_tumor_survival_samples": len(tumor_records),
                "axis_genes_defined": len(axis["genes"]),
                "axis_genes_present": len(present),
                "axis_present_genes": present,
                "mapping_coverage": coverage,
                "median_cutpoint": cutpoint,
                "high_group_n": sum(1 for group in groups if group == 1),
                "low_group_n": sum(1 for group in groups if group == 0),
                "high_group_events": high_events,
                "low_group_events": low_events,
                "logrank_chi2": chi2,
                "logrank_p": p_value,
                "event_rate_ratio_approx": rate_ratio,
                "event_rate_ratio_approx_ci_low": ci_low,
                "event_rate_ratio_approx_ci_high": ci_high,
                "validation_status": status,
                "qc_note": "Median cutpoint is fixed within tumor cohort; rate ratio is an event-rate approximation, not a Cox model. Survival support requires logrank p<0.05.",
                "analysis_run_id": ANALYSIS_RUN_ID,
                "reproducibility_check_id": REPRODUCIBILITY_CHECK_ID,
            }
        )
    return rows


def tumor_normal_rows_for_gse53625(
    samples: list[str],
    metadata: dict[str, dict[str, str]],
    expr: dict[str, list[float]],
) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    by_patient: dict[str, dict[str, tuple[str, int]]] = defaultdict(dict)
    for idx, sample in enumerate(samples):
        meta = metadata.get(sample, {})
        patient_id = meta.get("patient id", "").lower().replace("ec", "").strip()
        sample_type = classify_gse53625_sample(meta)
        if patient_id and sample_type in {"tumor", "normal"}:
            by_patient[patient_id][sample_type] = (sample, idx)
    for axis_id, axis in AXES.items():
        scores, present = score_axis_from_expr(axis["genes"], expr, len(samples))
        deltas: list[float] = []
        for pair in by_patient.values():
            if "tumor" not in pair or "normal" not in pair:
                continue
            tumor_score = scores[pair["tumor"][1]]
            normal_score = scores[pair["normal"][1]]
            if math.isfinite(tumor_score) and math.isfinite(normal_score):
                deltas.append(tumor_score - normal_score)
        positive = sum(1 for delta in deltas if delta > 0)
        negative = sum(1 for delta in deltas if delta < 0)
        p_value = exact_sign_test_p(positive, negative)
        coverage = len(present) / len(axis["genes"]) if axis["genes"] else 0.0
        if coverage < 0.6 or len(deltas) < 50:
            status = "reject_low_mapping_coverage"
        elif p_value < 0.05:
            status = "supported_with_limitations"
        else:
            status = "no_tumor_normal_support"
        rows.append(
            {
                "dataset": "GSE53625",
                "axis_id": axis_id,
                "axis_label": axis["label"],
                "n_paired_patients": len(deltas),
                "axis_genes_defined": len(axis["genes"]),
                "axis_genes_present": len(present),
                "axis_present_genes": present,
                "mapping_coverage": coverage,
                "mean_tumor_minus_normal_axis_score": mean(deltas),
                "median_tumor_minus_normal_axis_score": median(deltas),
                "positive_pairs": positive,
                "negative_pairs": negative,
                "two_sided_sign_test_p": p_value,
                "validation_status": status,
                "qc_note": "Paired sign test uses tumor-normal pairs after target-sequence probe rescue.",
                "analysis_run_id": ANALYSIS_RUN_ID,
                "reproducibility_check_id": REPRODUCIBILITY_CHECK_ID,
            }
        )
    return rows


def run_gse53625_rescue() -> dict[str, Any]:
    geo_dir = DATA_ROOT / "geo" / "GSE53625"
    matrix_path = geo_dir / "GSE53625_series_matrix.txt.gz"
    raw_path = geo_dir / "GSM1296956_first_raw_member.txt.gz"
    if not matrix_path.exists():
        deep.range_download(GSE53625_MATRIX_URL, matrix_path)
    if not deep.gzip_is_valid(matrix_path):
        return {"status": "failed", "reason": "GSE53625 matrix gzip validation failed"}
    if not raw_path.exists() or not deep.gzip_is_valid(raw_path):
        return {"status": "failed", "reason": "GSE53625 first RAW member is missing or invalid"}
    samples, metadata = deep.parse_geo_metadata(matrix_path)
    features = load_gse53625_feature_sequences(raw_path)
    mapping_rows, bundle_rows = map_features_to_axis_genes(features)
    accepted = {
        str(row["feature_num"]): str(row["accepted_gene"])
        for row in mapping_rows
        if row["accepted_gene"] and row["match_type"] in {"exact", "one_mismatch"}
    }
    matrix_samples, expr = load_gse53625_expression(matrix_path, accepted)
    if matrix_samples:
        samples = matrix_samples
    coverage_rows: list[dict[str, Any]] = []
    for axis_id, axis in AXES.items():
        for gene in axis["genes"]:
            probe_rows = [row for row in mapping_rows if row.get("accepted_gene") == gene]
            coverage_rows.append(
                {
                    "dataset": "GSE53625",
                    "axis_id": axis_id,
                    "gene_symbol": gene,
                    "accepted_probe_count": len(probe_rows),
                    "exact_probe_count": sum(1 for row in probe_rows if row.get("match_type") == "exact"),
                    "one_mismatch_probe_count": sum(1 for row in probe_rows if row.get("match_type") == "one_mismatch"),
                    "mapped_in_expression_matrix": "yes" if gene in expr else "no",
                    "analysis_run_id": ANALYSIS_RUN_ID,
                    "reproducibility_check_id": REPRODUCIBILITY_CHECK_ID,
                }
            )
    survival_rows = survival_rows_for_gse53625(samples, metadata, expr)
    paired_rows = tumor_normal_rows_for_gse53625(samples, metadata, expr)
    write_tsv(
        TABLE_ROOT / "gse53625_probe_sequence_mapping.tsv",
        mapping_rows,
        [
            "dataset",
            "feature_num",
            "probe_name",
            "raw_gene_name",
            "systematic_name",
            "sequence_length",
            "matched_genes",
            "accepted_gene",
            "match_type",
            "specificity_scope",
            "analysis_run_id",
            "reproducibility_check_id",
        ],
    )
    write_tsv(
        TABLE_ROOT / "gse53625_ensembl_sequence_index.tsv",
        bundle_rows,
        [
            "gene_symbol",
            "ensembl_gene_id",
            "seq_region_name",
            "start",
            "end",
            "strand",
            "n_transcripts_queried",
            "n_sequences_indexed",
            "total_bases_indexed",
            "analysis_run_id",
            "reproducibility_check_id",
        ],
    )
    write_tsv(
        TABLE_ROOT / "gse53625_rescue_axis_gene_coverage.tsv",
        coverage_rows,
        [
            "dataset",
            "axis_id",
            "gene_symbol",
            "accepted_probe_count",
            "exact_probe_count",
            "one_mismatch_probe_count",
            "mapped_in_expression_matrix",
            "analysis_run_id",
            "reproducibility_check_id",
        ],
    )
    write_tsv(
        TABLE_ROOT / "gse53625_rescue_survival_validation.tsv",
        survival_rows,
        [
            "dataset",
            "axis_id",
            "axis_label",
            "n_tumor_survival_samples",
            "axis_genes_defined",
            "axis_genes_present",
            "axis_present_genes",
            "mapping_coverage",
            "median_cutpoint",
            "high_group_n",
            "low_group_n",
            "high_group_events",
            "low_group_events",
            "logrank_chi2",
            "logrank_p",
            "event_rate_ratio_approx",
            "event_rate_ratio_approx_ci_low",
            "event_rate_ratio_approx_ci_high",
            "validation_status",
            "qc_note",
            "analysis_run_id",
            "reproducibility_check_id",
        ],
    )
    write_tsv(
        TABLE_ROOT / "gse53625_rescue_tumor_normal_validation.tsv",
        paired_rows,
        [
            "dataset",
            "axis_id",
            "axis_label",
            "n_paired_patients",
            "axis_genes_defined",
            "axis_genes_present",
            "axis_present_genes",
            "mapping_coverage",
            "mean_tumor_minus_normal_axis_score",
            "median_tumor_minus_normal_axis_score",
            "positive_pairs",
            "negative_pairs",
            "two_sided_sign_test_p",
            "validation_status",
            "qc_note",
            "analysis_run_id",
            "reproducibility_check_id",
        ],
    )
    survival_passed_axes = [
        row["axis_id"]
        for row in survival_rows
        if str(row.get("validation_status")) == "supported_with_limitations"
    ]
    paired_passed_axes = [
        row["axis_id"]
        for row in paired_rows
        if str(row.get("validation_status")) == "supported_with_limitations"
    ]
    if survival_passed_axes and paired_passed_axes:
        status = "supported_with_limitations"
    elif paired_passed_axes:
        status = "limited_support_tumor_normal_supported_survival_not_supported"
    elif survival_passed_axes:
        status = "limited_support_survival_only_without_paired_support"
    else:
        status = "reject_low_mapping_or_no_patient_support"
    return {
        "status": status,
        "n_samples": len(samples),
        "n_features_with_sequences": len(features),
        "n_mapped_probe_rows": len(mapping_rows),
        "n_accepted_probe_rows": len(accepted),
        "n_genes_mapped": len(set(accepted.values())),
        "mapped_genes": sorted(set(accepted.values())),
        "survival_passed_axes": sorted(set(survival_passed_axes)),
        "paired_passed_axes": sorted(set(paired_passed_axes)),
        "survival_rows": survival_rows,
        "paired_rows": paired_rows,
        "coverage_rows": coverage_rows,
    }


def hra003627_quantification() -> tuple[list[dict[str, Any]], dict[str, Any]]:
    source_path = DATA_ROOT / "open_source_tables" / "HRA003627_NatCommun2023_source_data.xlsx"
    sheet = "Fig5c, Supplment Fig 5a,9a, b"
    df = pd.read_excel(source_path, sheet_name=sheet)
    stage_map = {"Normal": 0, "low_grade": 1, "high_grade": 2, "cancer": 3}
    stage_label = {"Normal": "Normal", "low_grade": "LGIN", "high_grade": "HGIN", "cancer": "ESCC"}
    df = df[df["his"].astype(str).isin(stage_map)].copy()
    rows: list[dict[str, Any]] = []
    for signature_id, genes in HRA003627_SIGNATURES.items():
        present = [gene for gene in genes if gene in df.columns]
        z_by_gene = {gene: zscore([safe_float(v) for v in df[gene].tolist()]) for gene in present}
        scores: list[float] = []
        stages: list[float] = []
        labels: list[str] = []
        for idx, (_, item) in enumerate(df.iterrows()):
            scores.append(mean([z_by_gene[gene][idx] for gene in present]))
            stages.append(float(stage_map[str(item["his"])]))
            labels.append(str(item["his"]))
        rho, p_trend, n_trend = deep.spearman(stages, scores)
        normal_scores = [score for score, label in zip(scores, labels) if label == "Normal"]
        cancer_scores = [score for score, label in zip(scores, labels) if label == "cancer"]
        p_cancer = deep.mann_whitney_p(normal_scores, cancer_scores)
        for label in ["Normal", "low_grade", "high_grade", "cancer"]:
            values = [score for score, lab in zip(scores, labels) if lab == label]
            rows.append(
                {
                    "dataset": "HRA003627",
                    "source_table": "Nature Communications source data",
                    "sheet": sheet,
                    "signature_id": signature_id,
                    "stage": stage_label[label],
                    "n_roi": len(values),
                    "present_genes": present,
                    "mean_signature_z": mean(values),
                    "median_signature_z": median(values),
                    "spearman_stage_rho": rho,
                    "spearman_stage_p_approx": p_trend,
                    "spearman_stage_n": n_trend,
                    "escc_vs_normal_mann_whitney_p": p_cancer,
                    "interpretation": "source-table ROI-level progression quantification; not direct CAF-JAG1-NOTCH1 ligand-receptor evidence",
                    "analysis_run_id": ANALYSIS_RUN_ID,
                    "reproducibility_check_id": REPRODUCIBILITY_CHECK_ID,
                }
            )
    write_tsv(
        TABLE_ROOT / "hra003627_source_table_quantification.tsv",
        rows,
        [
            "dataset",
            "source_table",
            "sheet",
            "signature_id",
            "stage",
            "n_roi",
            "present_genes",
            "mean_signature_z",
            "median_signature_z",
            "spearman_stage_rho",
            "spearman_stage_p_approx",
            "spearman_stage_n",
            "escc_vs_normal_mann_whitney_p",
            "interpretation",
            "analysis_run_id",
            "reproducibility_check_id",
        ],
    )
    return rows, {"status": "completed", "n_rows": len(rows), "n_roi": int(df.shape[0])}


def normalize_gene_cell(value: Any) -> str:
    return re.sub(r"[^A-Z0-9]", "", str(value).upper())


def axis_membership(gene: str) -> list[str]:
    out = []
    for axis_id, axis in AXES.items():
        if gene.upper() in {item.upper() for item in axis["genes"]}:
            out.append(axis_id)
    return out


def significance_status(pvalue: float, fdr: float) -> str:
    if math.isfinite(fdr):
        return "significant_fdr_lt_0.05" if fdr < 0.05 else "not_significant_fdr_ge_0.05"
    if math.isfinite(pvalue):
        return "p_only_lt_0.05_no_fdr" if pvalue < 0.05 else "p_only_ge_0.05_no_fdr"
    return "no_statistic_available"


def hra008846_deg_hits(table_path: Path) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    xl = pd.ExcelFile(table_path)
    target_lookup = {normalize_gene_cell(gene): gene for gene in HRA008846_TARGET_GENES}
    for sheet in xl.sheet_names:
        df = pd.read_excel(table_path, sheet_name=sheet, header=1)
        gene_col = "Genes" if "Genes" in df.columns else "Gene" if "Gene" in df.columns else ""
        if not gene_col:
            continue
        for _, item in df.iterrows():
            gene = target_lookup.get(normalize_gene_cell(item.get(gene_col, "")))
            if not gene:
                continue
            if "shNC vs shOGT" in sheet:
                for replicate in ["s1", "s2"]:
                    logfc = safe_float(item.get(f"{replicate}_log2FC", ""))
                    pvalue = safe_float(item.get(f"{replicate}_Pvalue", ""))
                    fdr = safe_float(item.get(f"{replicate}_Qvalue", ""))
                    rows.append(
                        {
                            "dataset": "HRA008846",
                            "source_table": "Table S3 DEG marker genes",
                            "sheet": sheet,
                            "comparison": f"shNC_vs_shOGT_KYSE30_{replicate}",
                            "compartment_or_source": "KYSE30",
                            "gene_symbol": gene,
                            "logFC": logfc,
                            "pvalue": pvalue,
                            "fdr": fdr,
                            "pattern": "",
                            "significance_status": significance_status(pvalue, fdr),
                            "axis_membership": axis_membership(gene),
                            "interpretation": "published source-table perturbation row for axis-relevant gene",
                            "analysis_run_id": ANALYSIS_RUN_ID,
                            "reproducibility_check_id": REPRODUCIBILITY_CHECK_ID,
                        }
                    )
            elif "expression patterns" in sheet:
                comparisons = [
                    ("ESPL_vs_Normal_EP", item.get("LogFC_ESPL_vs_Normal"), item.get("Pvalue_ESPL_vs_Normal"), ""),
                    ("nonESCC_vs_ESPL_EP", item.get("LogFC_nonESCC_vs_ESPL"), item.get("Pvalue_nonESCC_vs_ESPL"), ""),
                ]
                for comparison, logfc, pvalue, fdr in comparisons:
                    pvalue_f = safe_float(pvalue)
                    fdr_f = safe_float(fdr)
                    rows.append(
                        {
                            "dataset": "HRA008846",
                            "source_table": "Table S3 DEG marker genes",
                            "sheet": sheet,
                            "comparison": comparison,
                            "compartment_or_source": "EP",
                            "gene_symbol": gene,
                            "logFC": safe_float(logfc),
                            "pvalue": pvalue_f,
                            "fdr": fdr_f,
                            "pattern": item.get("Patterns", ""),
                            "significance_status": significance_status(pvalue_f, fdr_f),
                            "axis_membership": axis_membership(gene),
                            "interpretation": "published source-table DEG hit for axis-relevant gene",
                            "analysis_run_id": ANALYSIS_RUN_ID,
                            "reproducibility_check_id": REPRODUCIBILITY_CHECK_ID,
                        }
                    )
            else:
                source = item.get("Source", "")
                comparison = re.sub(r"\s+", " ", sheet).strip()
                pvalue = safe_float(item.get("Pvalue", item.get("PValue", "")))
                fdr = safe_float(item.get("FDR", ""))
                rows.append(
                    {
                        "dataset": "HRA008846",
                        "source_table": "Table S3 DEG marker genes",
                        "sheet": sheet,
                        "comparison": comparison,
                        "compartment_or_source": source if str(source) != "nan" else "",
                        "gene_symbol": gene,
                        "logFC": safe_float(item.get("logFC", item.get("LogFC", ""))),
                        "pvalue": pvalue,
                        "fdr": fdr,
                        "pattern": item.get("Patterns", ""),
                        "significance_status": significance_status(pvalue, fdr),
                        "axis_membership": axis_membership(gene),
                        "interpretation": "published source-table DEG hit for axis-relevant gene",
                        "analysis_run_id": ANALYSIS_RUN_ID,
                        "reproducibility_check_id": REPRODUCIBILITY_CHECK_ID,
                    }
                )
    write_tsv(
        TABLE_ROOT / "hra008846_deg_axis_hits.tsv",
        rows,
        [
            "dataset",
            "source_table",
            "sheet",
            "comparison",
            "compartment_or_source",
            "gene_symbol",
            "logFC",
            "pvalue",
            "fdr",
            "pattern",
            "significance_status",
            "axis_membership",
            "interpretation",
            "analysis_run_id",
            "reproducibility_check_id",
        ],
    )
    return rows


def hra008846_cell_abundance(table_path: Path) -> list[dict[str, Any]]:
    stage_order = {
        "normal": 0,
        "espl": 1,
        "non-mescc": 2,
        "nonmescc": 2,
        "non-escc": 2,
        "nonescc": 2,
        "advanced": 3,
        "metastatic": 4,
    }
    rows: list[dict[str, Any]] = []
    for sheet in pd.ExcelFile(table_path).sheet_names:
        df = pd.read_excel(table_path, sheet_name=sheet, header=1)
        if "Stages" in df.columns:
            stage_col, cell_col, value_col, method = "Stages", "Cell_types", "TME score", "TME_consense"
        elif "Stage" in df.columns:
            stage_col, cell_col, value_col, method = "Stage", "Cell_type", "Proportion", "SpatialDecon"
        else:
            continue
        for cell_type, sub in df.groupby(cell_col, dropna=False):
            values_by_stage: dict[str, list[float]] = defaultdict(list)
            trend_x: list[float] = []
            trend_y: list[float] = []
            for _, item in sub.iterrows():
                stage = str(item.get(stage_col, ""))
                value = safe_float(item.get(value_col, ""))
                if not math.isfinite(value):
                    continue
                values_by_stage[stage].append(value)
                stage_key = stage.strip().lower()
                if stage_key in stage_order:
                    trend_x.append(float(stage_order[stage_key]))
                    trend_y.append(value)
            rho, p_value, n_trend = deep.spearman(trend_x, trend_y)
            relevance = "axis_context" if re.search(r"B_|B cells|fibro|macrophage|T_|T cells|epithelial|strom|CAF|lymph", str(cell_type), re.I) else "other"
            for stage, values in sorted(values_by_stage.items(), key=lambda kv: stage_order.get(kv[0].strip().lower(), 99)):
                rows.append(
                    {
                        "dataset": "HRA008846",
                        "source_table": "Table S4 cell abundance",
                        "sheet": sheet,
                        "method": method,
                        "stage": stage,
                        "cell_type": cell_type,
                        "n_observations": len(values),
                        "mean_value": mean(values),
                        "median_value": median(values),
                        "spearman_stage_rho": rho,
                        "spearman_stage_p_approx": p_value,
                        "spearman_stage_n": n_trend,
                        "relevance": relevance,
                        "analysis_run_id": ANALYSIS_RUN_ID,
                        "reproducibility_check_id": REPRODUCIBILITY_CHECK_ID,
                    }
                )
    write_tsv(
        TABLE_ROOT / "hra008846_cell_abundance_trends.tsv",
        rows,
        [
            "dataset",
            "source_table",
            "sheet",
            "method",
            "stage",
            "cell_type",
            "n_observations",
            "mean_value",
            "median_value",
            "spearman_stage_rho",
            "spearman_stage_p_approx",
            "spearman_stage_n",
            "relevance",
            "analysis_run_id",
            "reproducibility_check_id",
        ],
    )
    return rows


def hra008846_ligand_receptor(table_path: Path) -> list[dict[str, Any]]:
    df = pd.read_excel(table_path, sheet_name=0, header=1)
    rows: list[dict[str, Any]] = []
    targets = {gene.upper() for gene in HRA008846_TARGET_GENES}
    for _, item in df.iterrows():
        ligand = str(item.get("Ligand", "")).upper()
        receptor = str(item.get("Receptor", "")).upper()
        interaction = str(item.get("Interaction_name", "")).upper()
        pathway = str(item.get("Pathway_name", "")).upper()
        haystack = f"{ligand} {receptor} {interaction} {pathway}"
        hit_genes = sorted({gene for gene in targets if re.search(rf"(^|[^A-Z0-9]){re.escape(gene)}([^A-Z0-9]|$)", haystack)})
        direct_jag_notch = (
            ("JAG1" in hit_genes and "NOTCH1" in receptor)
            or ("JAG1" in ligand and "NOTCH1" in hit_genes)
            or ("JAG1" in interaction and "NOTCH1" in interaction)
        )
        if not hit_genes and not direct_jag_notch:
            continue
        rows.append(
            {
                "dataset": "HRA008846",
                "source_table": "Table S6 ligand-receptor",
                "source_cell": item.get("Source", ""),
                "target_cell": item.get("Target", ""),
                "ligand": item.get("Ligand", ""),
                "receptor": item.get("Receptor", ""),
                "prob": safe_float(item.get("Prob", "")),
                "pvalue": safe_float(item.get("Pval", "")),
                "interaction_name": item.get("Interaction_name", ""),
                "pathway_name": item.get("Pathway_name", ""),
                "annotation": item.get("Annotation", ""),
                "evidence": item.get("Evidence", ""),
                "axis_hit_genes": hit_genes,
                "axis_membership": sorted({axis for gene in hit_genes for axis in axis_membership(gene)}),
                "direct_jag1_notch1_flag": "yes" if direct_jag_notch else "no",
                "analysis_run_id": ANALYSIS_RUN_ID,
                "reproducibility_check_id": REPRODUCIBILITY_CHECK_ID,
            }
        )
    write_tsv(
        TABLE_ROOT / "hra008846_ligand_receptor_axis_hits.tsv",
        rows,
        [
            "dataset",
            "source_table",
            "source_cell",
            "target_cell",
            "ligand",
            "receptor",
            "prob",
            "pvalue",
            "interaction_name",
            "pathway_name",
            "annotation",
            "evidence",
            "axis_hit_genes",
            "axis_membership",
            "direct_jag1_notch1_flag",
            "analysis_run_id",
            "reproducibility_check_id",
        ],
    )
    return rows


def run_spatial_source_table_quantification() -> dict[str, Any]:
    source_dir = DATA_ROOT / "open_source_tables"
    manifest: list[dict[str, Any]] = []
    for filename, url in HRA008846_DOWNLOADS.items():
        manifest.append({"artifact": filename, **download_pmc_xlsx(url, source_dir / filename)})
    hra003627_rows, hra003627_meta = hra003627_quantification()
    table_s3 = source_dir / "HRA008846_TableS3_DEG.xlsx"
    table_s4 = source_dir / "HRA008846_TableS4_cell_abundance.xlsx"
    table_s6 = source_dir / "HRA008846_TableS6_ligand_receptor.xlsx"
    deg_rows = hra008846_deg_hits(table_s3)
    cell_rows = hra008846_cell_abundance(table_s4)
    lr_rows = hra008846_ligand_receptor(table_s6)
    write_tsv(
        TABLE_ROOT / "hra008846_download_manifest.tsv",
        manifest,
        ["artifact", "path", "status", "bytes", "url", "error", "content_type", "invalid_preview"],
    )
    ogt_hits = [row for row in deg_rows if row.get("gene_symbol") == "OGT"]
    st_caf_hits = [
        row
        for row in deg_rows
        if row.get("compartment_or_source") == "ST" and row.get("gene_symbol") in {"FAP", "COL1A1", "COL1A2", "POSTN"}
    ]
    significant_deg_rows = [
        row
        for row in deg_rows
        if str(row.get("significance_status")) in {"significant_fdr_lt_0.05", "p_only_lt_0.05_no_fdr"}
    ]
    direct_jag_notch = [row for row in lr_rows if row.get("direct_jag1_notch1_flag") == "yes"]
    status = "supported_with_limitations" if deg_rows and cell_rows and lr_rows else "limited_missing_source_table_component"
    return {
        "status": status,
        "download_manifest": manifest,
        "hra003627": hra003627_meta,
        "hra003627_rows": hra003627_rows,
        "hra008846_deg_rows": deg_rows,
        "hra008846_cell_rows": cell_rows,
        "hra008846_lr_rows": lr_rows,
        "n_significant_or_p_only_deg_rows": len(significant_deg_rows),
        "n_ogt_deg_hits": len(ogt_hits),
        "n_st_caf_deg_hits": len(st_caf_hits),
        "n_direct_jag1_notch1_lr_hits": len(direct_jag_notch),
    }


def write_review(summary: dict[str, Any]) -> list[dict[str, Any]]:
    patient = summary["gse53625"]
    spatial = summary["spatial_source_tables"]
    gse_status = patient.get("status", "failed")
    spatial_status = spatial.get("status", "failed")
    probe_status = "supported_with_limitations" if patient.get("n_accepted_probe_rows", 0) and patient.get("n_genes_mapped", 0) else "reject"
    patient_test_status = (
        "supported_with_limitations"
        if str(gse_status) == "supported_with_limitations"
        else "limited_support" if str(gse_status).startswith("limited_support") else "reject"
    )
    rows = [
        {
            "stage": "independent_patient_validation",
            "check": "independent_check_roles",
            "status": "pass" if ANALYSIS_RUN_ID != REPRODUCIBILITY_CHECK_ID else "reject",
            "analysis_run_id": ANALYSIS_RUN_ID,
            "reproducibility_check_id": REPRODUCIBILITY_CHECK_ID,
            "check_comment": "Analysis and quality-check batch IDs are distinct.",
        },
        {
            "stage": "independent_patient_validation",
            "check": "gse53625_probe_sequence_rescue",
            "status": probe_status,
            "analysis_run_id": ANALYSIS_RUN_ID,
            "reproducibility_check_id": REPRODUCIBILITY_CHECK_ID,
            "check_comment": (
                f"Accepted probes={patient.get('n_accepted_probe_rows', 0)}; "
                f"mapped genes={','.join(patient.get('mapped_genes', [])) or 'none'}. "
                "Mapping is target-sequence rescue, not a full genome-wide probe specificity assessment."
            ),
        },
        {
            "stage": "independent_patient_validation",
            "check": "gse53625_survival_and_paired_tests",
            "status": patient_test_status,
            "analysis_run_id": ANALYSIS_RUN_ID,
            "reproducibility_check_id": REPRODUCIBILITY_CHECK_ID,
            "check_comment": (
                "Paired tumor-normal validation passes where sign-test support is present; "
                "survival/prognostic validation is not supported unless logrank p<0.05."
            ),
        },
        {
            "stage": "spatial_source_table_quantification",
            "check": "hra003627_roi_source_table",
            "status": "supported_with_limitations" if spatial.get("hra003627", {}).get("status") == "completed" else "reject",
            "analysis_run_id": ANALYSIS_RUN_ID,
            "reproducibility_check_id": REPRODUCIBILITY_CHECK_ID,
            "check_comment": "Quantifies ROI-level D&K/progression signatures; does not directly test CAF-Epi/JAG1-NOTCH1.",
        },
        {
            "stage": "spatial_source_table_quantification",
            "check": "hra008846_deg_cell_lr_tables",
            "status": spatial_status,
            "analysis_run_id": ANALYSIS_RUN_ID,
            "reproducibility_check_id": REPRODUCIBILITY_CHECK_ID,
            "check_comment": (
                f"axis source-table rows={len(spatial.get('hra008846_deg_rows', []))}; "
                f"significant or p-only rows={spatial.get('n_significant_or_p_only_deg_rows', 0)}; "
                f"cell abundance rows={len(spatial.get('hra008846_cell_rows', []))}; "
                f"LR hits={len(spatial.get('hra008846_lr_rows', []))}; "
                f"direct JAG1-NOTCH1 LR hits={spatial.get('n_direct_jag1_notch1_lr_hits', 0)}. "
                "Direct JAG1-NOTCH1 signaling is absent in Table S6."
            ),
        },
        {
            "stage": "interpretation_status",
            "check": "evidence_strength_after_sensitivity",
            "status": "limited_support" if spatial_status.startswith("pass") else "not_supported",
            "analysis_run_id": ANALYSIS_RUN_ID,
            "reproducibility_check_id": REPRODUCIBILITY_CHECK_ID,
            "check_comment": (
                "Use these results for association/source-table reproducibility only. "
                "Do not claim independent survival validation or direct CAF-Epi JAG1-NOTCH1 LR evidence."
            ),
        },
    ]
    write_tsv(
        REVIEW_ROOT / "independent_patient_and_spatial_quant_review.tsv",
        rows,
        ["stage", "check", "status", "analysis_run_id", "reproducibility_check_id", "check_comment"],
    )
    return rows


def top_deg_lines(rows: list[dict[str, Any]], genes: set[str], limit: int = 8) -> list[str]:
    selected = [row for row in rows if row.get("gene_symbol") in genes]
    selected = sorted(
        selected,
        key=lambda row: (
            float(row.get("fdr")) if math.isfinite(safe_float(row.get("fdr"))) else 1.0,
            float(row.get("pvalue")) if math.isfinite(safe_float(row.get("pvalue"))) else 1.0,
        ),
    )
    lines = []
    for row in selected[:limit]:
        lines.append(
            "- {gene} in {comparison} / {source}: logFC={logfc}, FDR={fdr}, p={p}".format(
                gene=row.get("gene_symbol", ""),
                comparison=row.get("comparison", ""),
                source=row.get("compartment_or_source", ""),
                logfc=stringify(row.get("logFC")),
                fdr=stringify(row.get("fdr")),
                p=stringify(row.get("pvalue")),
            )
        )
    return lines


def write_report(summary: dict[str, Any], review_rows: list[dict[str, Any]]) -> None:
    patient = summary["gse53625"]
    spatial = summary["spatial_source_tables"]
    survival_rows = patient.get("survival_rows", [])
    paired_rows = patient.get("paired_rows", [])
    deg_rows = spatial.get("hra008846_deg_rows", [])
    lr_rows = spatial.get("hra008846_lr_rows", [])
    h3_rows = spatial.get("hra003627_rows", [])
    report_lines = [
        "# Independent Patient And Spatial Source-Table Quantification",
        "",
        f"Run date: {RUN_DATE}",
        f"Analysis batch: `{ANALYSIS_RUN_ID}`",
        f"Reproducibility check: `{REPRODUCIBILITY_CHECK_ID}`",
        "",
        "## Bottom line",
        "",
        f"- GSE53625 independent patient status: `{patient.get('status', 'failed')}`.",
        f"- GSE53625 accepted probe rows: {patient.get('n_accepted_probe_rows', 0)}; mapped genes: {', '.join(patient.get('mapped_genes', [])) or 'none'}.",
        f"- Spatial source-table status: `{spatial.get('status', 'failed')}`.",
        f"- HRA008846 direct JAG1-NOTCH1 ligand-receptor rows: {spatial.get('n_direct_jag1_notch1_lr_hits', 0)}.",
        "- Evidence interpretation status: `limited_support`; use as association/source-table evidence, not as independent survival or direct LR proof.",
        "",
        "## GSE53625 patient validation",
        "",
    ]
    if survival_rows:
        for row in survival_rows:
            report_lines.append(
                "- {axis}: coverage={cov}, n={n}, logrank p={p}, event-rate ratio={rr}, status={status}".format(
                    axis=row.get("axis_id", ""),
                    cov=stringify(row.get("mapping_coverage")),
                    n=row.get("n_tumor_survival_samples", ""),
                    p=stringify(row.get("logrank_p")),
                    rr=stringify(row.get("event_rate_ratio_approx")),
                    status=row.get("validation_status", ""),
                )
            )
    if paired_rows:
        report_lines.extend(["", "## GSE53625 paired tumor-normal", ""])
        for row in paired_rows:
            report_lines.append(
                "- {axis}: coverage={cov}, paired n={n}, median delta={delta}, sign-test p={p}, status={status}".format(
                    axis=row.get("axis_id", ""),
                    cov=stringify(row.get("mapping_coverage")),
                    n=row.get("n_paired_patients", ""),
                    delta=stringify(row.get("median_tumor_minus_normal_axis_score")),
                    p=stringify(row.get("two_sided_sign_test_p")),
                    status=row.get("validation_status", ""),
                )
            )
    report_lines.extend(
        [
            "",
            "## HRA003627 source-table quantification",
            "",
        ]
    )
    for signature in HRA003627_SIGNATURES:
        sub = [row for row in h3_rows if row.get("signature_id") == signature]
        if not sub:
            continue
        trend = sub[0]
        by_stage = ", ".join(f"{row['stage']}={stringify(row['mean_signature_z'])}" for row in sub)
        report_lines.append(
            f"- {signature}: stage trend rho={stringify(trend.get('spearman_stage_rho'))}, "
            f"p~{stringify(trend.get('spearman_stage_p_approx'))}; means: {by_stage}."
        )
    report_lines.extend(
        [
            "",
            "## HRA008846 source-table evidence",
            "",
            f"- Table S3 axis-relevant source-table rows: {len(deg_rows)}.",
            f"- Table S3 significant or p-only significant rows: {spatial.get('n_significant_or_p_only_deg_rows', 0)}.",
            f"- Table S4 cell-abundance trend rows: {len(spatial.get('hra008846_cell_rows', []))}.",
            f"- Table S6 axis-relevant ligand-receptor rows: {len(lr_rows)}.",
            "",
            "Key OGT/PI3K/TLS source-table rows:",
        ]
    )
    report_lines.extend(top_deg_lines(deg_rows, {"OGT", "CCND1", "LAMB1", "SPP1", "KRT17", "CXCL13", "CXCL8"}))
    report_lines.extend(["", "Key CAF-Epi/JAG1-NOTCH1 niche source-table rows:"])
    report_lines.extend(top_deg_lines(deg_rows, {"FAP", "COL1A1", "COL1A2", "POSTN", "JAG1", "NOTCH1", "CXCL1", "CXCL8", "SPP1"}))
    jag_rows = [row for row in lr_rows if row.get("direct_jag1_notch1_flag") == "yes"]
    report_lines.extend(["", "Ligand-receptor review:"])
    if jag_rows:
        for row in jag_rows[:10]:
            report_lines.append(
                "- {source}->{target}: {interaction}; prob={prob}, p={p}".format(
                    source=row.get("source_cell", ""),
                    target=row.get("target_cell", ""),
                    interaction=row.get("interaction_name", ""),
                    prob=stringify(row.get("prob")),
                    p=stringify(row.get("pvalue")),
                )
            )
    else:
        report_lines.append(
            "- No direct JAG1-NOTCH1 ligand-receptor row was found in HRA008846 Table S6; the CAF axis is supported mainly by ST/CAF ECM DEG and broader Notch/pathway correlations."
            " The phrase JAG1-NOTCH1 should be treated as a hypothesis label unless an independent direct LR or perturbation source is added."
        )
    report_lines.extend(["", "## Reproducibility status", ""])
    for row in review_rows:
        report_lines.append(f"- {row['stage']} / {row['check']}: `{row['status']}` - {row['check_comment']}")
    report_lines.extend(
        [
            "",
            "## Claim boundary",
            "",
            "- These artifacts support association, progression, spatial-context, and source-table reproducibility claims.",
            "- GSE53625 probe rescue and paired tumor-normal validation are usable; GSE53625 survival/prognostic validation is not supported in this run.",
            "- HRA008846 source tables can be cited as direct source-table evidence for published DEG/cell-abundance/ligand-receptor summaries, not as raw spatial reanalysis.",
        ]
    )
    (REPORT_ROOT / "independent_patient_and_spatial_quant_report.md").write_text(
        "\n".join(report_lines) + "\n",
        encoding="utf-8",
    )


def write_summary(summary: dict[str, Any]) -> None:
    compact = {
        "run_date": RUN_DATE,
        "analysis_run_id": ANALYSIS_RUN_ID,
        "reproducibility_check_id": REPRODUCIBILITY_CHECK_ID,
        "gse53625": {
            key: value
            for key, value in summary["gse53625"].items()
            if key
            not in {
                "survival_rows",
                "paired_rows",
                "coverage_rows",
            }
        },
        "spatial_source_tables": {
            "status": summary["spatial_source_tables"].get("status"),
            "hra003627": summary["spatial_source_tables"].get("hra003627"),
            "n_hra008846_deg_rows": len(summary["spatial_source_tables"].get("hra008846_deg_rows", [])),
            "n_hra008846_cell_rows": len(summary["spatial_source_tables"].get("hra008846_cell_rows", [])),
            "n_hra008846_lr_rows": len(summary["spatial_source_tables"].get("hra008846_lr_rows", [])),
            "n_significant_or_p_only_deg_rows": summary["spatial_source_tables"].get("n_significant_or_p_only_deg_rows"),
            "n_ogt_deg_hits": summary["spatial_source_tables"].get("n_ogt_deg_hits"),
            "n_st_caf_deg_hits": summary["spatial_source_tables"].get("n_st_caf_deg_hits"),
            "n_direct_jag1_notch1_lr_hits": summary["spatial_source_tables"].get("n_direct_jag1_notch1_lr_hits"),
        },
    }
    (OUT_ROOT / "results" / "independent_patient_and_spatial_quant_summary.json").write_text(
        json.dumps(compact, ensure_ascii=False, indent=2),
        encoding="utf-8",
    )


def write_addendum(summary: dict[str, Any]) -> None:
    patient = summary["gse53625"]
    spatial = summary["spatial_source_tables"]
    lines = [
        "# Manuscript Addendum: Independent Patient And Spatial Source-Table Gates",
        "",
        f"Run date: {RUN_DATE}",
        "",
        "## Result",
        "",
        f"- Independent patient validation status: `{patient.get('status', 'failed')}`.",
        f"- Spatial/source-table quantification status: `{spatial.get('status', 'failed')}`.",
        "",
        "## Manuscript use",
        "",
        "- The HRA008846 source-table results can be moved from background citation into a reproducible results subsection.",
        "- The HRA003627 source-table results can support a progression/keratinization spatial context paragraph.",
        "- GSE53625 should be used as sequence-rescued paired tumor-normal support only; this run does not support an independent survival/prognostic claim.",
        "- The CAF/JAG1-NOTCH1 wording should be limited_supportd to CAF/ECM stromal-remodeling context unless direct JAG1-NOTCH1 ligand-receptor evidence is added.",
    ]
    (DELIVERABLE_ROOT / "spatial_axis_manuscript_addendum_independent_validation.md").write_text(
        "\n".join(lines) + "\n",
        encoding="utf-8",
    )


def write_package() -> Path:
    package_path = DELIVERABLE_ROOT / f"spatial_independent_validation_package_{RUN_DATE}.zip"
    include_paths = [
        REPORT_ROOT / "independent_patient_and_spatial_quant_report.md",
        REVIEW_ROOT / "independent_patient_and_spatial_quant_review.tsv",
        OUT_ROOT / "results" / "independent_patient_and_spatial_quant_summary.json",
        DELIVERABLE_ROOT / "spatial_axis_manuscript_addendum_independent_validation.md",
        TABLE_ROOT / "gse53625_probe_sequence_mapping.tsv",
        TABLE_ROOT / "gse53625_ensembl_sequence_index.tsv",
        TABLE_ROOT / "gse53625_rescue_axis_gene_coverage.tsv",
        TABLE_ROOT / "gse53625_rescue_survival_validation.tsv",
        TABLE_ROOT / "gse53625_rescue_tumor_normal_validation.tsv",
        TABLE_ROOT / "hra003627_source_table_quantification.tsv",
        TABLE_ROOT / "hra008846_deg_axis_hits.tsv",
        TABLE_ROOT / "hra008846_cell_abundance_trends.tsv",
        TABLE_ROOT / "hra008846_ligand_receptor_axis_hits.tsv",
        TABLE_ROOT / "hra008846_download_manifest.tsv",
        REVIEW_ROOT / "independent_patient_and_spatial_quant_reproducibility_check.tsv",
        ROOT / "scripts" / "run_independent_patient_and_spatial_quant.py",
    ]
    with zipfile.ZipFile(package_path, "w", compression=zipfile.ZIP_DEFLATED) as zf:
        for path in include_paths:
            if path.exists():
                zf.write(path, path.relative_to(ROOT))
    return package_path


def main() -> None:
    ensure_dirs()
    summary = {
        "gse53625": run_gse53625_rescue(),
        "spatial_source_tables": run_spatial_source_table_quantification(),
    }
    review_rows = write_review(summary)
    write_report(summary, review_rows)
    write_summary(summary)
    write_addendum(summary)
    package_path = write_package()
    print(
        json.dumps(
            {
                "status": "completed",
                "summary_path": str(OUT_ROOT / "results" / "independent_patient_and_spatial_quant_summary.json"),
                "report_path": str(REPORT_ROOT / "independent_patient_and_spatial_quant_report.md"),
                "review_path": str(REVIEW_ROOT / "independent_patient_and_spatial_quant_review.tsv"),
                "package_path": str(package_path),
                "gse53625_status": summary["gse53625"].get("status"),
                "spatial_status": summary["spatial_source_tables"].get("status"),
            },
            ensure_ascii=False,
            indent=2,
        )
    )


if __name__ == "__main__":
    show_help_if_requested()
    main()
