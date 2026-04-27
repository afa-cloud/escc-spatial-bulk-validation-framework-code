#!/usr/bin/env python
"""Shared public-data helpers for the ESCC spatial-to-bulk workflow.

This module contains only data-access and statistical helper functions used by
the submitted analysis scripts. It intentionally does not generate manuscripts,
cover letters, DOCX files, or submission packages.
"""

from __future__ import annotations

import json
import math
import urllib.parse
import urllib.request
from collections import defaultdict
from pathlib import Path
from typing import Any

import xenaPython as xena

ROOT = Path(__file__).resolve().parents[1]

TOIL_HUB = "https://toil.xenahubs.net"
GENE_DATASET = "TcgaTargetGtex_rsem_gene_tpm"
PHENOTYPE_DATASET = "TcgaTargetGTEX_phenotype.txt"


def http_json(url: str, params: dict[str, Any] | None = None, timeout: int = 120) -> dict[str, Any]:
    if params:
        url = url + "?" + urllib.parse.urlencode(params)
    request = urllib.request.Request(url, headers={"User-Agent": "escc-public-data-workflow/1.0"})
    with urllib.request.urlopen(request, timeout=timeout) as response:
        return json.loads(response.read().decode("utf-8"))


def decode_xena_codes(codes: list[dict[str, str]]) -> dict[str, list[str]]:
    decoded: dict[str, list[str]] = {}
    for item in codes:
        code = item.get("code")
        if code:
            decoded[item["name"]] = code.split("\t")
    return decoded


def categorical_decode(value: Any, labels: list[str] | None) -> str:
    if value == "NaN" or value is None:
        return ""
    if labels is None:
        return str(value)
    try:
        idx = int(value)
    except (TypeError, ValueError):
        return str(value)
    if 0 <= idx < len(labels):
        return labels[idx]
    return str(value)


def load_toil_sample_sets() -> dict[str, Any]:
    samples = xena.dataset_samples(TOIL_HUB, PHENOTYPE_DATASET, None)
    fields = ["_study", "_primary_site", "_sample_type", "primary disease or tissue", "detailed_category"]
    _, value_sets = xena.dataset_probe_values(TOIL_HUB, PHENOTYPE_DATASET, samples, fields)
    codes = decode_xena_codes(xena.field_codes(TOIL_HUB, PHENOTYPE_DATASET, fields))

    records = []
    for idx, sample in enumerate(samples):
        item = {"sample": sample}
        for field, values in zip(fields, value_sets):
            item[field] = categorical_decode(values[idx], codes.get(field))
        records.append(item)

    esca_primary = [
        record["sample"]
        for record in records
        if record["_study"] == "TCGA"
        and record["_sample_type"] == "Primary Tumor"
        and record["primary disease or tissue"] == "Esophageal Carcinoma"
    ]
    gtex_esophagus = [
        record["sample"]
        for record in records
        if record["_study"] == "GTEX"
        and record["_sample_type"] == "Normal Tissue"
        and record["_primary_site"] == "Esophagus"
    ]
    return {"records": records, "esca_primary": esca_primary, "gtex_esophagus": gtex_esophagus}


def load_gdc_esca_squamous_cases() -> dict[str, Any]:
    filters = {"op": "in", "content": {"field": "project.project_id", "value": ["TCGA-ESCA"]}}
    data = http_json(
        "https://api.gdc.cancer.gov/cases",
        {
            "filters": json.dumps(filters),
            "fields": (
                "case_id,submitter_id,diagnoses.primary_diagnosis,diagnoses.morphology,"
                "diagnoses.vital_status,diagnoses.days_to_death,diagnoses.days_to_last_follow_up"
            ),
            "format": "JSON",
            "size": 1000,
        },
    )
    hits = data.get("data", {}).get("hits", [])
    squamous = set()
    adenocarcinoma = set()
    for hit in hits:
        submitter = hit.get("submitter_id", "")
        diagnoses = hit.get("diagnoses", [])
        diag_text = " ".join(
            str(diag.get("primary_diagnosis", "")) + " " + str(diag.get("morphology", ""))
            for diag in diagnoses
        ).lower()
        if "squamous" in diag_text or "8070/3" in diag_text or "8071/3" in diag_text:
            squamous.add(submitter)
        if "adenocarcinoma" in diag_text or "8140/3" in diag_text:
            adenocarcinoma.add(submitter)
    return {
        "hits": hits,
        "squamous_cases": sorted(squamous),
        "adenocarcinoma_cases": sorted(adenocarcinoma),
    }


def log2_to_tpm(value: Any) -> float:
    try:
        numeric = float(value)
    except (TypeError, ValueError):
        return float("nan")
    if math.isnan(numeric):
        return float("nan")
    return max(0.0, (2.0**numeric) - 0.001)


def fetch_gene_values(gene_symbols: list[str], samples: list[str]) -> dict[str, list[float]]:
    if not gene_symbols or not samples:
        return {}
    rows = xena.dataset_gene_probe_avg(TOIL_HUB, GENE_DATASET, samples, gene_symbols)
    output: dict[str, list[float]] = {}
    for row in rows:
        scores = row.get("scores") or []
        values = scores[0] if scores else []
        output[row["gene"]] = [log2_to_tpm(value) for value in values]
    return output


def logrank_p(times: list[float], events: list[int], groups: list[int]) -> tuple[float, float]:
    unique_times = sorted({time for time, event in zip(times, events) if event == 1})
    observed = 0.0
    expected = 0.0
    variance = 0.0
    for time in unique_times:
        at_risk_1 = sum(1 for item_time, group in zip(times, groups) if item_time >= time and group == 1)
        at_risk_0 = sum(1 for item_time, group in zip(times, groups) if item_time >= time and group == 0)
        events_1 = sum(
            1
            for item_time, event, group in zip(times, events, groups)
            if item_time == time and event == 1 and group == 1
        )
        events_total = sum(1 for item_time, event in zip(times, events) if item_time == time and event == 1)
        at_risk_total = at_risk_1 + at_risk_0
        if at_risk_total <= 1:
            continue
        observed += events_1
        expected += events_total * (at_risk_1 / at_risk_total)
        variance += (
            at_risk_1
            * at_risk_0
            * events_total
            * (at_risk_total - events_total)
            / ((at_risk_total**2) * (at_risk_total - 1))
        )
    if variance <= 0:
        return 0.0, 1.0
    chi_square = ((observed - expected) ** 2) / variance
    p_value = math.erfc(math.sqrt(chi_square / 2.0))
    return chi_square, p_value


def mann_whitney_p(group_a: list[float], group_b: list[float]) -> tuple[float, float]:
    a = [value for value in group_a if not math.isnan(value)]
    b = [value for value in group_b if not math.isnan(value)]
    n1, n2 = len(a), len(b)
    if n1 == 0 or n2 == 0:
        return 0.0, 1.0
    ranked = sorted([(value, 0) for value in a] + [(value, 1) for value in b], key=lambda item: item[0])
    ranks = [0.0] * len(ranked)
    i = 0
    while i < len(ranked):
        j = i
        while j < len(ranked) and ranked[j][0] == ranked[i][0]:
            j += 1
        avg_rank = (i + 1 + j) / 2.0
        for k in range(i, j):
            ranks[k] = avg_rank
        i = j
    rank_sum_a = sum(rank for rank, item in zip(ranks, ranked) if item[1] == 0)
    u1 = rank_sum_a - (n1 * (n1 + 1) / 2.0)
    mean_u = n1 * n2 / 2.0
    sd_u = math.sqrt(n1 * n2 * (n1 + n2 + 1) / 12.0)
    if sd_u == 0:
        return u1, 1.0
    z = (u1 - mean_u) / sd_u
    p_value = math.erfc(abs(z) / math.sqrt(2.0))
    return u1, p_value
