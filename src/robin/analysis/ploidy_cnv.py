"""
Ploidy / sex karyotype estimation from sample bin-level CNV profiles (r_cnv).

Uses the same bin-level coverage signal produced by ``cnv_from_bam`` as the main
CNV pipeline (``r_cnv``), not normalized sample-minus-reference tracks (CNV3),
which are centered for CNV calling and are a poor match for X/Y vs autosome
depth ratios.

When ``r2_cnv`` (reference pass, persisted as ``CNV2.npy``) is provided, autosomal
pooled values and X/Y medians use **only bins where the reference bin value is
non-zero**, so Y (and other chromosomes) are not dominated by unmappable / zero
reference bins.
"""

from __future__ import annotations

import logging
from typing import Any, Dict, Iterable, List, Optional, Tuple

import numpy as np

# Skewness thresholds for classifying WGS vs WES from autosomal bin distribution
SKEWNESS_WGS_MAX = 0.2
SKEWNESS_WES_MIN = 0.6

# (name, x_min, x_max, y_min, y_max) — checked in order; first match wins.
_SEX_KARYOTYPE_RULES: Tuple[Tuple[str, float, float, float, float], ...] = (
    ("XXX", 1.25, 1.75, 0.00, 0.25),
    ("XXXY", 1.25, 1.75, 0.25, 0.75),
    ("XXY", 0.75, 1.25, 0.25, 0.75),
    ("XYY", 0.25, 0.75, 0.75, 1.25),
    ("XX", 0.75, 1.25, 0.00, 0.25),
    ("XY", 0.25, 0.75, 0.25, 0.75),
    ("X0", 0.25, 0.75, 0.00, 0.25),
)


def _autosome_keys(keys: Iterable[str]) -> List[str]:
    out: List[str] = []
    for k in keys:
        if not k.startswith("chr"):
            continue
        rest = k[3:]
        if rest.isdigit():
            n = int(rest)
            if 1 <= n <= 22:
                out.append(k)
    return sorted(out, key=lambda x: int(x[3:]))


def _pooled_autosome_bins(r_cnv: Dict[str, np.ndarray]) -> np.ndarray:
    parts: List[np.ndarray] = []
    for c in _autosome_keys(r_cnv.keys()):
        if c not in r_cnv:
            continue
        arr = np.asarray(r_cnv[c], dtype=float).ravel()
        if arr.size:
            parts.append(arr)
    if not parts:
        return np.array([], dtype=float)
    return np.concatenate(parts)


def _align_pair(
    samp: np.ndarray, ref: np.ndarray
) -> Tuple[np.ndarray, np.ndarray]:
    s = np.asarray(samp, dtype=float).ravel()
    r = np.asarray(ref, dtype=float).ravel()
    n = min(len(s), len(r))
    return s[:n], r[:n]


def _pooled_autosome_bins_ref_masked(
    r_cnv: Dict[str, np.ndarray],
    r2_cnv: Dict[str, np.ndarray],
) -> np.ndarray:
    """Concatenate sample autosome bins where **reference** bin > 0."""
    parts: List[np.ndarray] = []
    for c in _autosome_keys(r_cnv.keys()):
        if c not in r_cnv:
            continue
        s = np.asarray(r_cnv[c], dtype=float).ravel()
        if s.size == 0:
            continue
        if c not in r2_cnv:
            parts.append(s)
            continue
        s2, r2 = _align_pair(s, r2_cnv[c])
        mask = r2 > 0
        picked = s2[mask]
        if picked.size:
            parts.append(picked)
        else:
            parts.append(s)
    if not parts:
        return np.array([], dtype=float)
    return np.concatenate(parts)


def _median_chr(r_cnv: Dict[str, np.ndarray], chrom: str) -> float:
    if chrom not in r_cnv:
        return 0.0
    arr = np.asarray(r_cnv[chrom], dtype=float).ravel()
    if arr.size == 0:
        return 0.0
    return float(np.median(arr))


def _median_chr_ref_masked(
    r_cnv: Dict[str, np.ndarray],
    r2_cnv: Optional[Dict[str, np.ndarray]],
    chrom: str,
) -> float:
    """Median of sample bins where reference has non-zero bins; else unmasked median."""
    if chrom not in r_cnv:
        return 0.0
    s = np.asarray(r_cnv[chrom], dtype=float).ravel()
    if s.size == 0:
        return 0.0
    if not r2_cnv or chrom not in r2_cnv:
        return float(np.median(s))
    s2, r2 = _align_pair(s, r2_cnv[chrom])
    mask = r2 > 0
    sub = s2[mask]
    if sub.size == 0:
        return float(np.median(s))
    return float(np.median(sub))


def _seq_type_from_skewness(skewness: float) -> str:
    if np.isnan(skewness):
        return "unknown"
    if skewness <= SKEWNESS_WGS_MAX:
        return "WGS"
    if skewness >= SKEWNESS_WES_MIN:
        return "WES"
    return "unknown"


def _match_karyotype(x_ratio: float, y_ratio: float) -> Optional[str]:
    for name, xmin, xmax, ymin, ymax in _SEX_KARYOTYPE_RULES:
        if xmin <= x_ratio <= xmax and ymin <= y_ratio <= ymax:
            return name
    return None


def estimate_ploidy_from_sample_cnv(
    r_cnv: Dict[str, np.ndarray],
    bin_width: int,
    logger: Optional[logging.Logger] = None,
    r2_cnv: Optional[Dict[str, np.ndarray]] = None,
) -> Dict[str, Any]:
    """
    Estimate sequencing type (WGS vs WES heuristic), coverage proxies, and sex karyotype
    from sample bin-level CNV data (``r_cnv``).

    Ratios follow the common convention: median(X) / median(autosomes) and
    median(Y) / median(autosomes), using pooled autosomal bin values for the
    autosome denominator.

    When ``r2_cnv`` is provided (reference pass bins, ``CNV2.npy``), only bins with
    **reference value > 0** are used for autosome pooling and for X/Y medians.
    If a chromosome has no reference non-zero bins, that chromosome falls back to
    unmasked sample medians / pooling for that chromosome.

    When skewness falls between ``SKEWNESS_WGS_MAX`` and ``SKEWNESS_WES_MIN``,
    sequencing type is unknown and reported median coverage proxies are zeroed,
    matching the reference behaviour you described.

    Absolute "2x genome coverage" is not available from uncalibrated bin counts;
    this function requires a positive pooled autosomal median and a resolved WGS/WES
    class before assigning a karyotype.

    Parameters
    ----------
    r_cnv
        Per-chromosome bin arrays from the sample pass of ``cnv_from_bam`` (same
        object persisted as CNV.npy / ``result.cnv``).
    bin_width
        Bin width in bases (for reporting; does not change ratios).
    logger
        Optional logger for debug lines.
    r2_cnv
        Optional reference pass bin arrays (``CNV2.npy``). When set, defines which
        bins are treated as informative in the reference.

    Returns
    -------
    dict
        Fields include ``seq_type``, ``skewness``, ``median_autosome``,
        ``median_chrX``, ``median_chrY``, ``x_ratio``, ``y_ratio``,
        ``sex_karyotype`` (or None), ``wes_median_exome_proxy``, and
        ``reported_median_autosome`` / ``reported_median_exome`` (zeroed when
        seq_type is unknown), and ``reference_guided`` (True when ``r2_cnv``
        was used for masking).
    """
    log = logger or logging.getLogger("robin.ploidy_cnv")

    if r2_cnv:
        pooled = _pooled_autosome_bins_ref_masked(r_cnv, r2_cnv)
        reference_guided = True
    else:
        pooled = _pooled_autosome_bins(r_cnv)
        reference_guided = False
    if pooled.size == 0:
        out = _empty_result(bin_width, reason="no_autosomal_bins")
        log.debug("ploidy_cnv: no autosomal bins")
        return out

    autosome_mean = float(np.mean(pooled))
    autosome_median = float(np.median(pooled))
    if autosome_mean <= 0:
        out = _empty_result(bin_width, reason="zero_autosome_mean")
        log.debug("ploidy_cnv: zero autosome mean")
        return out

    skewness = abs(autosome_mean - autosome_median) / autosome_mean
    seq_type = _seq_type_from_skewness(skewness)

    if reference_guided:
        median_x = _median_chr_ref_masked(r_cnv, r2_cnv, "chrX")
        median_y = _median_chr_ref_masked(r_cnv, r2_cnv, "chrY")
    else:
        median_x = _median_chr(r_cnv, "chrX")
        median_y = _median_chr(r_cnv, "chrY")

    reported_auto = 0.0
    reported_exome = 0.0
    wes_exome_proxy = 0.0

    if seq_type != "unknown":
        reported_auto = autosome_median
        if seq_type == "WES":
            p99_per_chr: List[float] = []
            for c in _autosome_keys(r_cnv.keys()):
                if c not in r_cnv:
                    continue
                s = np.asarray(r_cnv[c], dtype=float).ravel()
                if s.size == 0:
                    continue
                if reference_guided and r2_cnv is not None and c in r2_cnv:
                    s2, r2 = _align_pair(s, r2_cnv[c])
                    mask = r2 > 0
                    arr = s2[mask]
                    if arr.size == 0:
                        arr = s
                else:
                    arr = s
                p99_per_chr.append(float(np.percentile(arr, 99)))
            if p99_per_chr:
                wes_exome_proxy = float(np.median(p99_per_chr))
                reported_exome = wes_exome_proxy

    x_ratio = (median_x / autosome_median) if autosome_median > 0 else float("nan")
    y_ratio = (median_y / autosome_median) if autosome_median > 0 else float("nan")

    sex_karyotype: Optional[str] = None
    coverage_ok = autosome_median > 0 and seq_type in ("WGS", "WES")
    if coverage_ok and not np.isnan(x_ratio) and not np.isnan(y_ratio):
        sex_karyotype = _match_karyotype(x_ratio, y_ratio)

    result: Dict[str, Any] = {
        "seq_type": seq_type,
        "skewness": float(skewness),
        "bin_width": int(bin_width),
        "median_autosome": autosome_median,
        "median_chrX": median_x,
        "median_chrY": median_y,
        "x_ratio": float(x_ratio) if not np.isnan(x_ratio) else None,
        "y_ratio": float(y_ratio) if not np.isnan(y_ratio) else None,
        "sex_karyotype": sex_karyotype,
        "reported_median_autosome": reported_auto,
        "reported_median_exome": reported_exome,
        "wes_median_exome_proxy": wes_exome_proxy,
        "autosome_mean": autosome_mean,
        "reference_guided": reference_guided,
    }
    log.debug(
        "ploidy_cnv: seq_type=%s skew=%.4f x_ratio=%s y_ratio=%s karyotype=%s ref_guided=%s",
        seq_type,
        skewness,
        result["x_ratio"],
        result["y_ratio"],
        sex_karyotype,
        reference_guided,
    )
    return result


def _empty_result(bin_width: int, reason: str) -> Dict[str, Any]:
    return {
        "seq_type": "unknown",
        "skewness": float("nan"),
        "bin_width": int(bin_width),
        "median_autosome": 0.0,
        "median_chrX": 0.0,
        "median_chrY": 0.0,
        "x_ratio": None,
        "y_ratio": None,
        "sex_karyotype": None,
        "reported_median_autosome": 0.0,
        "reported_median_exome": 0.0,
        "wes_median_exome_proxy": 0.0,
        "autosome_mean": 0.0,
        "reference_guided": False,
        "failure_reason": reason,
    }
