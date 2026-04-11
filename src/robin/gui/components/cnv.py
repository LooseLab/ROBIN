from __future__ import annotations

from typing import Any, Dict, List, Optional, Tuple
from pathlib import Path

import asyncio
import natsort
import numpy as np
from functools import lru_cache
import logging
import pickle
import time
import json
import importlib.resources as importlib_resources
import pandas as pd

try:
    from nicegui import ui
except ImportError:  # pragma: no cover
    ui = None

from robin.gui.theme import styled_table, register_theme_sync_callback
from robin.analysis.cnv_classification import detect_cnv_events, get_cnv_summary, CNVEvent
from robin.analysis.cnv_analysis import moving_average, pad_arrays, persist_ichor_py_outputs
from robin.analysis.ploidy_cnv import estimate_ploidy_from_sample_cnv
from robin.classification_config import get_cnv_thresholds

# Same chromosome set as reporting (plotting.py): chr0–chr22, chrX, chrY only
CNV_PLOT_CONTIGS = frozenset(
    ["chr" + str(i) for i in range(0, 23)] + ["chrX", "chrY"]
)

# Integer-clonality metrics use this genomic bin width (not the GUI plot bin width).
INTEGER_CLONAL_TARGET_BIN_BP = 500_000


def _cnv_contig_ok(contig: str) -> bool:
    """True if contig should be included in CNV plots (matches report behaviour)."""
    return contig in CNV_PLOT_CONTIGS


def _cnv_sample_bin_map_from_state(state: Dict[str, Any]) -> Optional[Dict[str, np.ndarray]]:
    """Resolve per-chromosome bin arrays from GUI state (matches ``CNV.npy`` layout)."""
    raw = state.get("cnv")
    if not raw:
        return None
    if isinstance(raw, dict) and "cnv" in raw:
        raw = raw["cnv"]
    return raw if isinstance(raw, dict) else None


def _cnv_ref_bin_map_from_state(state: Dict[str, Any]) -> Optional[Dict[str, np.ndarray]]:
    """Reference pass bins (``CNV2.npy``), same layout as sample."""
    raw = state.get("cnv2")
    if not raw:
        return None
    if isinstance(raw, dict) and "cnv" in raw:
        raw = raw["cnv"]
    return raw if isinstance(raw, dict) else None


def _format_ploidy_gui_line(est: Dict[str, Any]) -> str:
    """Single-line summary for the CNV panel."""
    fr = est.get("failure_reason")
    if fr:
        return f"Ploidy (CNV bins): unavailable ({fr})"
    st = est.get("seq_type", "unknown")
    sk = est.get("sex_karyotype")
    rg = est.get("reference_guided")
    xr, yr = est.get("x_ratio"), est.get("y_ratio")
    skew = est.get("skewness")
    parts = [f"class: {st}"]
    if rg:
        parts.append("ref-masked bins")
    if sk:
        parts.append(f"karyotype: {sk}")
    if xr is not None and yr is not None:
        parts.append(f"X ratio {xr:.3f}, Y ratio {yr:.3f}")
    if isinstance(skew, (int, float)) and not (isinstance(skew, float) and np.isnan(skew)):
        parts.append(f"skew {skew:.3f}")
    return "Ploidy (CNV bins): " + " · ".join(parts)


def _apply_ploidy_rescale(vals: np.ndarray, assumed_ploidy: float) -> np.ndarray:
    """Rescale saved CNV/CNV3 arrays as if ``cnv_from_bam`` had used this ploidy (default lib = 2).

    Library output is proportional to ``(bin_sum / median) * ploidy_run`` with ``ploidy_run=2``.
    Multiplying stored values by ``assumed_ploidy / 2`` matches changing only that factor.
    """
    p = float(assumed_ploidy)
    if p <= 0:
        p = 2.0
    return np.asarray(vals, dtype=float) * (p / 2.0)


def _downsample_cnv_bins_for_plot(
    values_1d: np.ndarray,
    analysis_bin_width: int,
    plot_bin_width: int,
) -> Tuple[np.ndarray, np.ndarray]:
    """Downsample CNV values when plot_bin_width > analysis_bin_width (same as CNV plots).

    Returns ``(x_positions_bp, values)``. If plot_bin_width <= analysis_bin_width,
    returns original positions and values unchanged.
    """
    if plot_bin_width <= analysis_bin_width or plot_bin_width <= 0:
        x_bp = np.arange(len(values_1d), dtype=float) * analysis_bin_width
        return x_bp, np.asarray(values_1d, dtype=float)
    group_size = int(plot_bin_width / analysis_bin_width)
    if group_size < 1:
        x_bp = np.arange(len(values_1d), dtype=float) * analysis_bin_width
        return x_bp, np.asarray(values_1d, dtype=float)
    n = len(values_1d)
    n_trim = (n // group_size) * group_size
    if n_trim == 0:
        x_bp = np.arange(n, dtype=float) * analysis_bin_width
        return x_bp, np.asarray(values_1d, dtype=float)
    trimmed = np.asarray(values_1d[:n_trim], dtype=float)
    grouped = trimmed.reshape(-1, group_size)
    values_out = np.mean(grouped, axis=1)
    x_bp = (np.arange(len(values_out)) + 0.5) * plot_bin_width
    return x_bp, values_out


def _cnv_autosome_keys(cnv_map: Dict[str, np.ndarray]) -> List[str]:
    """chr1–chr22 only (excludes chr0 and sex chromosomes)."""
    out: List[str] = []
    for k in natsort.natsorted(cnv_map.keys()):
        if k.startswith("chr") and k[3:].isdigit():
            n = int(k[3:])
            if 1 <= n <= 22:
                out.append(k)
    return out


def _cnv_unwrap_inner_map(raw: Any) -> Optional[Dict[str, np.ndarray]]:
    if not raw or not isinstance(raw, dict):
        return None
    if "cnv" in raw:
        raw = raw["cnv"]
    return raw if isinstance(raw, dict) else None


def _cnv_autosomes_have_bins(cnv_dict: Optional[Dict[str, np.ndarray]]) -> bool:
    if not cnv_dict:
        return False
    for k in _cnv_autosome_keys(cnv_dict):
        if k in cnv_dict and np.asarray(cnv_dict[k]).size > 0:
            return True
    return False


def _masked_contiguous_runs(mask: np.ndarray) -> List[Tuple[int, int]]:
    """Return ``[start, end)`` index spans where ``mask`` is True."""
    m = np.asarray(mask, dtype=bool)
    if m.size == 0 or not np.any(m):
        return []
    idx = np.where(m)[0]
    runs: List[Tuple[int, int]] = []
    start = int(idx[0])
    prev = int(idx[0])
    for i in idx[1:]:
        ii = int(i)
        if ii == prev + 1:
            prev = ii
        else:
            runs.append((start, prev + 1))
            start = ii
            prev = ii
    runs.append((start, prev + 1))
    return runs


def _chunk_1d(arr: np.ndarray, chunk_bins: int) -> List[np.ndarray]:
    if chunk_bins < 1:
        chunk_bins = 1
    a = np.asarray(arr, dtype=float).ravel()
    if a.size == 0:
        return []
    out: List[np.ndarray] = []
    for i in range(0, len(a), chunk_bins):
        out.append(a[i : i + chunk_bins])
    return out


def _merge_adjacent_chunks(
    chunks: List[np.ndarray], merge_tol: float
) -> List[np.ndarray]:
    """Merge consecutive chunks when their medians differ by less than ``merge_tol``."""
    if not chunks:
        return []
    med = [float(np.median(c)) for c in chunks]
    groups: List[List[np.ndarray]] = [[chunks[0]]]
    for i in range(1, len(chunks)):
        if abs(med[i] - med[i - 1]) < merge_tol:
            groups[-1].append(chunks[i])
        else:
            groups.append([chunks[i]])
    return [np.concatenate(g) for g in groups]


def _segment_medians_from_run(
    r_run: np.ndarray, chunk_bins: int, merge_tol: float
) -> List[float]:
    """Split a residual run into coarse chunks, merge neighbours, return one median per segment."""
    chunks = _chunk_1d(r_run, chunk_bins)
    if not chunks:
        return []
    merged = _merge_adjacent_chunks(chunks, merge_tol)
    return [float(np.median(block)) for block in merged]


def _segment_spans_from_series(
    x: np.ndarray,
    y: np.ndarray,
    *,
    chunk_bins: int = 5,
    merge_tol: float = 0.12,
) -> List[Tuple[float, float, float]]:
    """Return horizontal segment spans ``(x_start, x_end, y_median)`` from series."""
    xx = np.asarray(x, dtype=float).ravel()
    yy = np.asarray(y, dtype=float).ravel()
    n = min(len(xx), len(yy))
    if n == 0:
        return []
    xx = xx[:n]
    yy = yy[:n]
    if n < 2:
        return [(float(xx[0]), float(xx[0]), float(yy[0]))]
    chunks_x = _chunk_1d(xx, max(1, int(chunk_bins)))
    chunks_y = _chunk_1d(yy, max(1, int(chunk_bins)))
    if not chunks_x or not chunks_y:
        return []
    med = [float(np.median(c)) for c in chunks_y]
    groups: List[List[int]] = [[0]]
    for i in range(1, len(chunks_y)):
        if abs(med[i] - med[i - 1]) < float(merge_tol):
            groups[-1].append(i)
        else:
            groups.append([i])
    spans: List[Tuple[float, float, float]] = []
    for g in groups:
        i0, i1 = g[0], g[-1]
        x0 = float(chunks_x[i0][0])
        x1 = float(chunks_x[i1][-1])
        ymed = float(np.median(np.concatenate([chunks_y[i] for i in g])))
        spans.append((x0, x1, ymed))
    return spans


def _resample_residuals_to_target_bin_width(
    r: np.ndarray,
    mask: np.ndarray,
    analysis_bin_width_bp: int,
    target_bin_width_bp: int,
) -> Tuple[np.ndarray, np.ndarray]:
    """Resample per-bin residuals to ~``target_bin_width_bp`` windows (mean + ref mask).

    If the analysis bin is **wider** than the target (e.g. 1 Mb vs 500 kb), each value is
    repeated across ``analysis/target`` target bins (constant within the analysis bin).
    If the analysis bin is **narrower**, consecutive bins are averaged until each block
    covers ~target width; all-masked blocks are dropped.
    """
    r = np.asarray(r, dtype=float).ravel()
    mask = np.asarray(mask, dtype=bool).ravel()
    n = min(len(r), len(mask))
    if n == 0:
        return np.array([], dtype=float), np.array([], dtype=bool)
    r, mask = r[:n], mask[:n]
    aw = max(1, int(analysis_bin_width_bp))
    tw = max(1, int(target_bin_width_bp))
    if aw == tw:
        return r, mask
    if tw > aw:
        k = max(1, int(round(float(tw) / float(aw))))
        n_out = n // k
        if n_out == 0:
            return np.array([], dtype=float), np.array([], dtype=bool)
        rs: List[float] = []
        ms: List[bool] = []
        for i in range(0, n_out * k, k):
            blk_r = r[i : i + k]
            blk_m = mask[i : i + k]
            if not np.any(blk_m):
                continue
            rs.append(float(np.mean(blk_r[blk_m])))
            ms.append(True)
        return np.asarray(rs, dtype=float), np.asarray(ms, dtype=bool)
    rep = max(1, int(round(float(aw) / float(tw))))
    return np.repeat(r, rep), np.repeat(mask, rep)


def _compute_data_driven_integer_band(
    r_all: np.ndarray,
    mask_all: np.ndarray,
    *,
    variable_threshold: float,
    k_sigma: float = 2.5,
    band_floor: float = 0.04,
    band_ceiling: float = 0.22,
    min_neutral_bins: int = 40,
) -> Tuple[float, float, int]:
    """Return ``(band, sigma, neutral_bin_count)`` for classifying non-integer residuals.

    ``sigma`` is the standard deviation of residuals in **neutral** bins (|r| below
    ``variable_threshold`` when enough bins exist; otherwise a widened percentile band).
    ``band = clip(k_sigma * sigma, band_floor, band_ceiling)``.
    """
    r = np.asarray(r_all, dtype=float).ravel()
    m = np.asarray(mask_all, dtype=bool).ravel()
    n = min(len(r), len(m))
    if n == 0:
        return float(band_floor), 0.0, 0
    r, m = r[:n], m[:n]
    if not np.any(m):
        return float(band_floor), 0.0, 0
    ra = r[m]
    if ra.size == 0:
        return float(band_floor), 0.0, 0
    vt = float(variable_threshold)
    neutral = np.abs(ra) < vt
    nn = int(np.sum(neutral))
    if nn < min_neutral_bins and ra.size >= 10:
        pct = float(
            np.percentile(
                np.abs(ra),
                min(55.0, max(25.0, 100.0 * min_neutral_bins / max(ra.size, 1))),
            )
        )
        neutral = np.abs(ra) <= pct
        nn = int(np.sum(neutral))
    if nn >= 10:
        sigma = float(np.std(ra[neutral]))
    else:
        sigma = float(np.std(ra))
    band = float(np.clip(k_sigma * sigma, band_floor, band_ceiling))
    return band, sigma, nn


def _build_autosome_residuals(
    cnv_map: Dict[str, np.ndarray],
    *,
    assumed_ploidy: float,
    cnv3_map: Optional[Dict[str, np.ndarray]],
    cnv2_map: Optional[Dict[str, np.ndarray]],
) -> Tuple[str, List[Tuple[str, np.ndarray, np.ndarray]]]:
    """Return ``(signal_source, list of (chrom, residual_r, mask))`` for autosomes.

    Preference order matches the batch pipeline: use persisted **CNV3** (smoothed
    sample−reference) when present; else **moving_average(sample) − moving_average(ref)**
    when **CNV2** exists; else absolute **CNV** minus assumed ploidy.

    **CNV2** bins with value ``<= 0`` are excluded (same unmappable / zero-control idea
    as ploidy estimation).
    """
    p = float(assumed_ploidy)
    if p <= 0:
        p = 2.0

    use_cnv3 = _cnv_autosomes_have_bins(cnv3_map)
    use_ref = cnv2_map is not None and isinstance(cnv2_map, dict)

    rows: List[Tuple[str, np.ndarray, np.ndarray]] = []

    if use_cnv3 and cnv3_map is not None:
        for contig in _cnv_autosome_keys(cnv3_map):
            if contig not in cnv3_map:
                continue
            raw3 = np.asarray(cnv3_map[contig], dtype=float).ravel()
            if raw3.size == 0:
                continue
            v = _apply_ploidy_rescale(raw3, p)
            r = v
            mask = np.ones(len(r), dtype=bool)
            if use_ref and contig in cnv2_map:
                ref = np.asarray(cnv2_map[contig], dtype=float).ravel()
                n = min(len(r), len(ref))
                r = r[:n]
                mask = ref[:n] > 0
            rows.append((contig, r, mask))
        if rows:
            return ("cnv3", rows)

    if use_ref:
        for contig in _cnv_autosome_keys(cnv_map):
            if contig not in cnv_map or contig not in cnv2_map:
                continue
            s = np.asarray(cnv_map[contig], dtype=float).ravel()
            ref = np.asarray(cnv2_map[contig], dtype=float).ravel()
            if s.size == 0:
                continue
            s_ma = moving_average(s)
            r_ma = moving_average(ref)
            s_ma, r_ma = pad_arrays(s_ma, r_ma)
            diff = s_ma - r_ma
            v = _apply_ploidy_rescale(diff, p)
            r = v
            mask = r_ma > 0
            rows.append((contig, np.asarray(r, dtype=float), mask))
        if rows:
            return ("ref_ma_diff", rows)

    for contig in _cnv_autosome_keys(cnv_map):
        raw = np.asarray(cnv_map[contig], dtype=float).ravel()
        if raw.size == 0:
            continue
        v = _apply_ploidy_rescale(raw, p)
        r = v - p
        mask = np.ones(len(r), dtype=bool)
        if use_ref and contig in cnv2_map:
            ref = np.asarray(cnv2_map[contig], dtype=float).ravel()
            n = min(len(r), len(ref))
            r = r[:n]
            mask = ref[:n] > 0
        rows.append((contig, r, mask))
    return ("absolute", rows)


def estimate_cnv_integer_clonal_diagnostics(
    cnv_map: Dict[str, np.ndarray],
    *,
    bin_width_analysis: int,
    plot_bin_width: int,
    assumed_ploidy: float,
    cnv3_map: Optional[Dict[str, np.ndarray]] = None,
    cnv2_map: Optional[Dict[str, np.ndarray]] = None,
    band: Optional[float] = None,
    band_k_sigma: float = 2.5,
    band_floor: float = 0.04,
    band_ceiling: float = 0.22,
    variable_threshold: float = 0.25,
    segment_chunk_bins: int = 5,
    segment_merge_tol: Optional[float] = None,
    target_bin_width_bp: int = INTEGER_CLONAL_TARGET_BIN_BP,
) -> Optional[Dict[str, Any]]:
    """Summarise integer clonality using **control-aware** residuals and **segment medians**.

    When **CNV3** (sample−reference, as in the batch pipeline) is available, residuals are
    ``rescale(CNV3)`` with expected **0** — reference pass **CNV2** removes shared GC /
    mappability bias. Bins with **reference ≤ 0** are skipped (unmappable control bins).

    If **CNV3** is absent but **CNV2** exists, we use ``MA(sample) − MA(reference)``
    like ``save_cnv_files`` / ``process_single_bam``.

    Otherwise we fall back to absolute ``rescale(CNV) − assumed_ploidy`` (optionally
    **ref-masked**).

    **Autosomes only** (chr1–chr22). Residuals are **resampled to** ``target_bin_width_bp``
    (default **500 kb**) before the non-integer band and segmentation — **not** the GUI plot
    bin width. The **band** for ``|r − round(r)|`` is **data-driven**:
    ``clip(band_k_sigma * σ, band_floor, band_ceiling)`` where **σ** is the standard deviation
    of residuals in **neutral** bins (``|r| < variable_threshold``, with a percentile fallback
    if too few). Pass ``band`` to override with a fixed width.

    Segments: on the **500 kb** (or target) series, chunks of ``segment_chunk_bins`` bins,
    merged when adjacent chunk medians differ by less than ``segment_merge_tol`` (default
    tied to the computed band).

    ``plot_bin_width`` is unused here (kept for API compatibility).

    Returns keys including ``signal_source``, ``band``, ``band_sigma``, ``frac_outside_band``,
    ``median_frac_dev_variable``, ``purity_proxy``, ``n_bins_all``, ``target_bin_width_bp``.
    """
    _ = plot_bin_width
    if not cnv_map:
        return None
    aw = int(bin_width_analysis)
    if aw <= 0:
        return None
    tw = max(1, int(target_bin_width_bp))

    p = float(assumed_ploidy)
    if p <= 0:
        p = 2.0

    source, per_chrom = _build_autosome_residuals(
        cnv_map,
        assumed_ploidy=p,
        cnv3_map=cnv3_map,
        cnv2_map=cnv2_map,
    )
    if not per_chrom:
        return None

    resampled: List[Tuple[np.ndarray, np.ndarray]] = []
    r_parts: List[np.ndarray] = []
    m_parts: List[np.ndarray] = []
    n_bins = 0
    for _contig, r, mask in per_chrom:
        r = np.asarray(r, dtype=float).ravel()
        mask = np.asarray(mask, dtype=bool).ravel()
        n = min(len(r), len(mask))
        if n == 0:
            continue
        r = r[:n]
        mask = mask[:n]
        r5, m5 = _resample_residuals_to_target_bin_width(r, mask, aw, tw)
        if r5.size == 0:
            continue
        resampled.append((r5, m5))
        r_parts.append(r5)
        m_parts.append(m5)
        n_bins += int(np.sum(m5))

    if not r_parts:
        return None

    r_cat = np.concatenate(r_parts)
    m_cat = np.concatenate(m_parts)
    band_auto, sigma, neutral_n = _compute_data_driven_integer_band(
        r_cat,
        m_cat,
        variable_threshold=variable_threshold,
        k_sigma=float(band_k_sigma),
        band_floor=float(band_floor),
        band_ceiling=float(band_ceiling),
    )
    band_eff = float(band) if band is not None else band_auto
    merge_tol = (
        float(segment_merge_tol)
        if segment_merge_tol is not None
        else max(0.08, 0.95 * band_eff)
    )
    chunk_bins = max(1, int(segment_chunk_bins))

    seg_medians: List[float] = []
    for r5, m5 in resampled:
        for start, end in _masked_contiguous_runs(m5):
            run = r5[start:end]
            if run.size == 0:
                continue
            seg_medians.extend(
                _segment_medians_from_run(run, chunk_bins, merge_tol)
            )

    if not seg_medians:
        return None

    sm = np.asarray(seg_medians, dtype=float)
    frac_dev = np.abs(sm - np.round(sm))
    n_seg = int(sm.size)
    mean_all = float(np.mean(frac_dev))
    frac_out = float(np.mean(frac_dev > band_eff))

    var_mask = np.abs(sm) > variable_threshold
    n_var_seg = int(np.sum(var_mask))
    base_meta = {
        "signal_source": source,
        "band": band_eff,
        "band_sigma": float(sigma),
        "band_auto": float(band_auto),
        "band_neutral_bins": float(neutral_n),
        "band_override": band is not None,
        "target_bin_width_bp": float(tw),
        "variable_threshold": variable_threshold,
        "n_bins_all": float(n_bins),
        "n_segments_all": float(n_seg),
        "mean_frac_dev_all": mean_all,
        "frac_outside_band": frac_out,
    }
    if n_var_seg == 0:
        return {
            **base_meta,
            "n_segments_variable": 0.0,
            "median_frac_dev_variable": float("nan"),
            "purity_proxy": float("nan"),
        }

    fv = frac_dev[var_mask]
    med_v = float(np.median(fv))
    purity = max(0.0, min(1.0, 1.0 - 2.0 * med_v))
    return {
        **base_meta,
        "n_segments_variable": float(n_var_seg),
        "median_frac_dev_variable": med_v,
        "purity_proxy": purity,
    }


def _format_integer_clonal_line(d: Optional[Dict[str, Any]]) -> str:
    """One-line GUI summary for :func:`estimate_cnv_integer_clonal_diagnostics`."""
    if not d:
        return "Integer clonality (heuristic): --"
    src = str(d.get("signal_source", ""))
    src_lbl = {
        "cnv3": "CNV3+ref-mask",
        "ref_ma_diff": "MA−ref",
        "absolute": "absolute CNV",
    }.get(src, src or "?")
    tw = int(d.get("target_bin_width_bp", INTEGER_CLONAL_TARGET_BIN_BP))
    tw_kb = max(1, tw // 1000)
    n_bins = int(d.get("n_bins_all", 0))
    n_seg = int(d.get("n_segments_all", 0))
    n_var = int(d.get("n_segments_variable", 0))
    band = float(d.get("band", 0.12))
    frac_out = d.get("frac_outside_band")
    med_v = d.get("median_frac_dev_variable")
    pur = d.get("purity_proxy")
    sig = d.get("band_sigma")
    parts = [f"{src_lbl}", f"{tw_kb} kb", f"{n_bins} bins", f"{n_seg} seg."]
    if d.get("band_override"):
        parts.append("band=fixed")
    elif sig is not None and not (isinstance(sig, float) and np.isnan(sig)):
        parts.append(f"band≈{band:.3f} (σ={float(sig):.3f})")
    else:
        parts.append(f"band≈{band:.3f}")
    if frac_out is not None and not (isinstance(frac_out, float) and np.isnan(frac_out)):
        parts.append(f"non-int. seg. Δ >{band:.2f}: {frac_out:.1%}")
    if n_var == 0:
        vt = float(d.get("variable_threshold", 0.25))
        parts.append(f"no seg. |r|>{vt:.2f} (proxy N/A)")
        return "Integer clonality (heuristic): " + " · ".join(parts)
    if med_v is not None and not (isinstance(med_v, float) and np.isnan(med_v)):
        parts.append(f"median |r−round(r)| (var. seg.): {med_v:.3f}")
    if pur is not None and not (isinstance(pur, float) and np.isnan(pur)):
        parts.append(f"proxy purity ≈ {pur:.0%}")
    return "Integer clonality (heuristic): " + " · ".join(parts)


def _cnv_load_binary_payload(
    sample_dir: Path,
    *,
    cnv_dict_npy_changed: bool,
    cnv_npy_changed: bool,
    cnv2_npy_changed: bool,
    cnv3_npy_changed: bool,
    data_array_reload: bool,
    xy_pkl_changed: bool,
) -> Dict[str, Any]:
    """Load CNV numpy/pickle data from disk (no UI). Safe for ``asyncio.to_thread``."""
    out: Dict[str, Any] = {}
    cnv_dict_npy = sample_dir / "CNV_dict.npy"
    cnv_npy = sample_dir / "CNV.npy"
    cnv2_npy = sample_dir / "CNV2.npy"
    cnv3_npy = sample_dir / "CNV3.npy"
    data_array_npy = sample_dir / "cnv_data_array.npy"
    xy_pkl = sample_dir / "XYestimate.pkl"

    if cnv_dict_npy_changed and cnv_dict_npy.exists():
        out["cnv_dict"] = np.load(cnv_dict_npy, allow_pickle=True).item()

    if xy_pkl_changed and xy_pkl.exists():
        try:
            with xy_pkl.open("rb") as f:
                out["xy"] = pickle.load(f)
        except Exception:
            pass

    if cnv_npy_changed and cnv_npy.exists():
        try:
            out["cnv"] = np.load(cnv_npy, allow_pickle=True).item()
        except Exception:
            out["cnv"] = None

    if cnv2_npy_changed and cnv2_npy.exists():
        try:
            out["cnv2"] = np.load(cnv2_npy, allow_pickle=True).item()
        except Exception:
            out["cnv2"] = None

    if cnv3_npy_changed and cnv3_npy.exists():
        try:
            out["cnv3"] = np.load(cnv3_npy, allow_pickle=True).item()
        except Exception:
            out["cnv3"] = None

    if data_array_reload and data_array_npy.exists():
        try:
            out["bp_array"] = np.load(data_array_npy, allow_pickle=True)
        except Exception:
            pass

    return out


def _cnv3_map_and_binwidth_for_ichor(
    state: Dict[str, Any], sample_dir: Path
) -> Tuple[Optional[Dict[str, np.ndarray]], int]:
    """Resolve CNV3 per-chromosome arrays and analysis bin width from GUI state or disk."""
    cnv3 = state.get("cnv3")
    if not cnv3 and (sample_dir / "CNV3.npy").exists():
        try:
            cnv3 = np.load(sample_dir / "CNV3.npy", allow_pickle=True).item()
        except Exception:
            cnv3 = None
    raw: Any = cnv3
    if isinstance(raw, dict) and "cnv" in raw:
        raw = raw["cnv"]
    cnv3_map = raw if isinstance(raw, dict) else None

    bw = 0
    cd = state.get("cnv_dict")
    if isinstance(cd, dict) and isinstance(cd.get("bin_width"), (int, float)):
        bw = int(cd["bin_width"])
    if bw <= 0 and (sample_dir / "CNV_dict.npy").exists():
        try:
            d = np.load(sample_dir / "CNV_dict.npy", allow_pickle=True).item()
            if isinstance(d, dict) and isinstance(d.get("bin_width"), (int, float)):
                bw = int(d["bin_width"])
        except Exception:
            pass
    return cnv3_map, bw


def _is_dark_mode() -> bool:
    """Quasar ``body--dark`` via app storage (see theme.frame)."""
    try:
        from nicegui import app

        return bool(app.storage.user.get("dark_mode"))
    except Exception:
        return False


def _cnv_chromosome_scatter_palette(dark: bool) -> List[str]:
    """Distinct scatter colours per chromosome.

    ECharts defaults include very dark greys that disappear on midnight backgrounds;
    dark mode uses lighter, saturated hues (design.md §5–§6).
    """
    if dark:
        return [
            "#34d399",
            "#38bdf8",
            "#fbbf24",
            "#fb7185",
            "#a78bfa",
            "#2dd4bf",
            "#f472b6",
            "#facc15",
            "#4ade80",
            "#60a5fa",
            "#f97316",
            "#e879f9",
            "#c084fc",
            "#22d3ee",
            "#fde047",
            "#93c5fd",
            "#f87171",
            "#bef264",
            "#5eead4",
            "#fcd34d",
            "#7dd3fc",
            "#fda4af",
            "#86efac",
            "#d8b4fe",
            "#eab308",
            "#67e8f9",
        ]
    return [
        "#059669",
        "#0284c7",
        "#b45309",
        "#dc2626",
        "#7c3aed",
        "#0d9488",
        "#db2777",
        "#ca8a04",
        "#16a34a",
        "#2563eb",
        "#ea580c",
        "#c026d3",
        "#9333ea",
        "#0891b2",
        "#ca8a04",
        "#3b82f6",
        "#ef4444",
        "#65a30d",
        "#14b8a6",
        "#eab308",
        "#0ea5e9",
        "#ec4899",
        "#4ade80",
        "#8b5cf6",
        "#ca8a04",
        "#06b6d4",
    ]


def _cnv_value_mode_colors(dark: bool) -> Tuple[str, str, str]:
    """High / low / in-range scatter colours (semantic green / rose / slate)."""
    if dark:
        return ("#34d399", "#fb7185", "#94a3b8")
    return ("#007AFF", "#FF3B30", "#8E8E93")


def _cnv_echart_palette(dark: bool) -> Dict[str, str]:
    """Axis, title, and tooltip colours for CNV scatter plots (design.md §5, §6)."""
    if dark:
        return {
            "text": "#e2e8f0",
            "muted": "#94a3b8",
            "axis_line": "#64748b",
            "split": "rgba(148, 163, 184, 0.52)",
            "tooltip_bg": "rgba(15, 23, 42, 0.96)",
            "tooltip_border": "#334155",
        }
    return {
        "text": "#0f172a",
        "muted": "#475569",
        "axis_line": "#94a3b8",
        "split": "rgba(71, 85, 105, 0.33)",
        "tooltip_bg": "rgba(255, 255, 255, 0.98)",
        "tooltip_border": "#e2e8f0",
    }


def _apply_cnv_echart_chrome(echart: Any, dark: bool) -> None:
    """Apply light/dark readable chrome without touching series data."""
    p = _cnv_echart_palette(dark)
    try:
        o = echart.options
        if not isinstance(o, dict):
            return
        o["backgroundColor"] = "transparent"
        o["textStyle"] = {"color": p["text"]}
        title = o.get("title")
        if isinstance(title, dict):
            title["textStyle"] = {"color": p["text"], "fontSize": 14}
        o["tooltip"] = {
            **(o.get("tooltip") or {}),
            "backgroundColor": p["tooltip_bg"],
            "borderColor": p["tooltip_border"],
            "textStyle": {"color": p["text"]},
        }
        xa = o.get("xAxis")
        if isinstance(xa, dict):
            xa["axisLine"] = {"lineStyle": {"color": p["axis_line"]}}
            xa["axisLabel"] = {**(xa.get("axisLabel") or {}), "color": p["muted"]}
            xa["splitLine"] = {"lineStyle": {"color": p["split"]}}
        ya_list = o.get("yAxis")
        if isinstance(ya_list, list):
            for ya in ya_list:
                if not isinstance(ya, dict):
                    continue
                ya["axisLine"] = {"lineStyle": {"color": p["axis_line"]}}
                ya["axisLabel"] = {**(ya.get("axisLabel") or {}), "color": p["muted"]}
                ya["nameTextStyle"] = {"color": p["muted"]}
                ya["splitLine"] = {"lineStyle": {"color": p["split"]}}
        leg = o.get("legend")
        if isinstance(leg, dict):
            leg["textStyle"] = {"color": p["muted"]}
        try:
            dz_list = o.get("dataZoom")
            if isinstance(dz_list, list):
                for dz in dz_list:
                    if not isinstance(dz, dict):
                        continue
                    dz["borderColor"] = p["axis_line"]
                    dz["fillerColor"] = (
                        "rgba(51, 65, 85, 0.35)"
                        if dark
                        else "rgba(148, 163, 184, 0.2)"
                    )
                    dz["handleStyle"] = {
                        "color": p["text"],
                        "borderColor": p["axis_line"],
                    }
                    dz["moveHandleStyle"] = {"color": p["muted"]}
                    dz["emphasis"] = {
                        "handleStyle": {"borderColor": p["text"]},
                    }
                    dbb = dz.get("dataBackground") or {}
                    if isinstance(dbb, dict):
                        dbb["lineStyle"] = {
                            **(dbb.get("lineStyle") or {}),
                            "color": p["muted"],
                        }
                        dbb["areaStyle"] = {
                            **(dbb.get("areaStyle") or {}),
                            "color": (
                                "rgba(148, 163, 184, 0.12)"
                                if dark
                                else "rgba(71, 85, 105, 0.08)"
                            ),
                        }
                        dz["dataBackground"] = dbb
        except Exception:
            pass
    except Exception:
        pass


def add_cnv_section(launcher: Any, sample_dir: Path) -> None:
    """Build the CNV UI section and attach refresh timers.

    Uses `launcher._cnv_state` for per-sample cache/state.
    Expects CNV.npy, CNV3.npy, CNV_dict.npy, XYestimate.pkl, cnv_data_array.npy in sample folder.

    The "Ploidy (CNV bins)" line is computed in the GUI from sample bin data (``CNV.npy`` +
    ``bin_width`` in ``CNV_dict.npy``), and when present uses ``CNV2.npy`` (reference bins)
    to restrict medians to bins where the reference bin is non-zero; it is not written during
    batch CNV analysis.

    **Assumed ploidy** (default 2) rescales loaded ``CNV.npy`` and ``CNV3.npy`` values by
    ``assumed / 2`` for plotting only, matching a different ``ploidy`` argument to
    ``cnv_from_bam`` without re-running the BAM.

    **Integer clonality** is a heuristic: prefer **CNV3** (sample−reference) with **CNV2**
    ref-masking to reduce GC/mappability noise; otherwise ``MA(sample)−MA(ref)`` or absolute
    CNV. Residuals are resampled to **500 kb** genomic bins (not the plot bin width); the
    non-integer band is **data-driven** from neutral-bin variance. **Segment medians**
    (chunk + merge) tighten noise. It is not a substitute for BAF-based tumour purity.

    Controls now trigger immediate refresh instead of waiting for timer updates.
    """
    with ui.element("div").classes("w-full min-w-0").props("id=analysis-detail-cnv"):
        with ui.element("div").classes("classification-insight-shell w-full min-w-0"):
            ui.label("Copy number (CNV)").classes(
                "classification-insight-heading text-headline-small"
            )
            with ui.element("div").classes(
                "classification-insight-card w-full min-w-0"
            ):
                with ui.column().classes("w-full min-w-0 gap-2 p-2 md:p-3"):
                    with ui.row().classes("items-center gap-2 min-w-0"):
                        ui.icon("person").classes("classification-insight-icon")
                        ui.label("Genome-wide profile").classes(
                            "classification-insight-model flex-1 min-w-0"
                        )
                    cnv_status = ui.label("Status: Awaiting Data").classes(
                        "classification-insight-result w-full"
                    )
                    cnv_xy = ui.label("Genetic sex: --").classes(
                        "classification-insight-meta w-full"
                    )
                    cnv_ploidy = ui.label("Ploidy (CNV bins): --").classes(
                        "classification-insight-meta w-full"
                    )
                    cnv_integer_clonal = ui.label(
                        "Integer clonality (heuristic): --"
                    ).classes("classification-insight-meta w-full")
                    with ui.row().classes("w-full items-center gap-2 flex-wrap"):
                        cnv_ichor_py = ui.label(
                            "Tumour fraction/ploidy (ichor_py): --"
                        ).classes("classification-insight-meta flex-grow min-w-0")
                        cnv_run_ichor_py = (
                            ui.button("Run ichor_py", icon="calculate")
                            .props("flat dense no-caps")
                            .tooltip(
                                "Compute tumour fraction and ploidy from CNV3 in memory "
                                "or from CNV3.npy on disk — no BAM re-run required."
                            )
                        )
                    with ui.row().classes(
                        "w-full justify-end gap-4 flex-wrap items-baseline"
                    ):
                        cnv_bin = ui.label("Bin width: --").classes(
                            "classification-insight-meta"
                        )
                        cnv_var = ui.label("Variance: --").classes(
                            "classification-insight-meta"
                        )
            with ui.row().classes(
                "w-full gap-3 items-center mb-2 flex-wrap mt-2"
            ):
                ui.label("Chromosome").classes("classification-insight-meta")
                cnv_chrom_select = ui.select(options={"All": "All"}, value="All").style(
                    "width: 160px"
                )
                ui.label("Gene").classes("classification-insight-meta ml-2")
                cnv_gene_select = ui.select(options={"All": "All"}, value="All").style(
                    "width: 200px"
                )
                ui.label("Color by").classes("classification-insight-meta ml-2")
                cnv_color = ui.toggle(
                    options={"chromosome": "Chromosome", "value": "Up/Down"},
                    value="chromosome",
                ).classes("mt-1")
                ui.label("Y-axis").classes("classification-insight-meta ml-2")
                cnv_scale = ui.toggle(
                    options={"linear": "Linear", "log": "Log"}, value="linear"
                ).classes("mt-1")
                ui.label("Plot bin").classes("classification-insight-meta ml-2")
                # NiceGUI select with dict uses keys as option values; map key -> bp (None = use data)
                _PLOT_BIN_KEY_DEFAULT = "Data default"
                _PLOT_BIN_OPTIONS = {
                    _PLOT_BIN_KEY_DEFAULT: "Data default",
                    "500 kb": "500 kb",
                    "1 Mb": "1 Mb",
                    "2 Mb": "2 Mb",
                    "5 Mb": "5 Mb",
                    "10 Mb": "10 Mb",
                }
                _PLOT_BIN_KEY_TO_BP = {
                    _PLOT_BIN_KEY_DEFAULT: None,
                    "500 kb": 500_000,
                    "1 Mb": 1_000_000,
                    "2 Mb": 2_000_000,
                    "5 Mb": 5_000_000,
                    "10 Mb": 10_000_000,
                }
                cnv_plot_bin = ui.select(
                    options=_PLOT_BIN_OPTIONS,
                    value=_PLOT_BIN_KEY_DEFAULT,
                ).style("width: 120px")
                ui.label("Segments").classes("classification-insight-meta ml-2")
                cnv_segments = ui.toggle(
                    options={"hide": "Hide", "show": "Show"}, value="show"
                ).classes("mt-1")
                ui.label("Segment source").classes("classification-insight-meta ml-2")
                cnv_segment_source = ui.toggle(
                    options={"heuristic": "Heuristic", "ichor_py": "IchorPy"},
                    value="heuristic",
                ).classes("mt-1")
                cnv_bp_label = ui.label("Breakpoints").classes(
                    "classification-insight-meta ml-2"
                ).style("display: none")
                cnv_bp = ui.toggle(
                    options={"hide": "Hide", "show": "Show"}, value="show"
                ).classes("mt-1").style("display: none")
            with ui.row().classes(
                "w-full gap-3 items-center mb-1 flex-wrap"
            ):
                ui.label("Assumed ploidy (display)").classes(
                    "classification-insight-meta"
                )
                cnv_assumed_ploidy = ui.number(
                    value=2.0, min=0.5, max=16.0, step=0.5, format="%.2f"
                ).classes("max-w-[140px]").props("dense")
                try:
                    _st_cnv = launcher._cnv_state.setdefault(str(sample_dir), {})
                    if "plot_assumed_ploidy" not in _st_cnv and "plot_abs_baseline_copy" in _st_cnv:
                        try:
                            b = float(_st_cnv.get("plot_abs_baseline_copy", 2.0))
                            if b > 0:
                                _st_cnv["plot_assumed_ploidy"] = 4.0 / b
                        except Exception:
                            pass
                    _ap = float(_st_cnv.get("plot_assumed_ploidy", 2.0))
                    if _ap > 0:
                        cnv_assumed_ploidy.value = _ap
                except Exception:
                    pass
            with ui.element("div").classes("w-full target-coverage-panel__plot-wrap"):
                cnv_abs = ui.echart(
                    {
                        "backgroundColor": "transparent",
                        "title": {"text": "CNV scatter plot", "left": "center", "top": 10},
                        "grid": {
                            "left": "5%",
                            "right": "5%",
                            "bottom": "10%",
                            "top": "20%",
                            "containLabel": True,
                        },
                        "tooltip": {"trigger": "axis"},
                        "xAxis": {"type": "value", "max": "dataMax"},
                        "yAxis": [
                            {"type": "value", "name": "Ploidy"},
                            {
                                "type": "value",
                                "name": "Breakpoint density",
                                "position": "right",
                            },
                        ],
                        "dataZoom": [
                            {"type": "slider", "xAxisIndex": [0]},
                            {
                                "type": "slider",
                                "yAxisIndex": [0, 1],
                                "right": 20,
                                "startValue": 0,
                                "endValue": 6,
                            },
                        ],
                        "series": [
                            {"type": "scatter", "name": "CNV", "symbolSize": 3, "data": []},
                            {
                                "type": "scatter",
                                "name": "centromeres_highlight",
                                "data": [],
                                "symbolSize": 3,
                                "markArea": {
                                    "itemStyle": {"color": "rgba(135, 206, 250, 0.4)"},
                                    "data": [],
                                },
                            },
                            {
                                "type": "scatter",
                                "name": "cytobands_highlight",
                                "data": [],
                                "symbolSize": 3,
                                "markArea": {
                                    "itemStyle": {"color": "rgba(200, 200, 200, 0.4)"},
                                    "data": [],
                                },
                                "markLine": {"symbol": "none", "data": []},
                            },
                        ],
                    }
                ).classes("w-full h-72")
            with ui.element("div").classes("w-full target-coverage-panel__plot-wrap mt-2"):
                cnv_diff = ui.echart(
                    {
                        "backgroundColor": "transparent",
                        "title": {"text": "Difference plot", "left": "center", "top": 10},
                        "grid": {
                            "left": "5%",
                            "right": "5%",
                            "bottom": "10%",
                            "top": "20%",
                            "containLabel": True,
                        },
                        "tooltip": {"trigger": "axis"},
                        "xAxis": {"type": "value", "max": "dataMax"},
                        "yAxis": [
                            {"type": "value", "name": "Relative"},
                            {
                                "type": "value",
                                "name": "Breakpoint density",
                                "position": "right",
                            },
                        ],
                        "dataZoom": [
                            {"type": "slider", "xAxisIndex": [0]},
                            {
                                "type": "slider",
                                "yAxisIndex": [0, 1],
                                "right": 20,
                                "startValue": -4,
                                "endValue": 4,
                            },
                        ],
                        "series": [
                            {
                                "type": "scatter",
                                "name": "CNV Δ",
                                "symbolSize": 3,
                                "data": [],
                            },
                            {
                                "type": "scatter",
                                "name": "centromeres_highlight",
                                "data": [],
                                "symbolSize": 3,
                                "markArea": {
                                    "itemStyle": {"color": "rgba(135, 206, 250, 0.4)"},
                                    "data": [],
                                },
                            },
                            {
                                "type": "scatter",
                                "name": "cytobands_highlight",
                                "data": [],
                                "symbolSize": 3,
                                "markArea": {
                                    "itemStyle": {"color": "rgba(200, 200, 200, 0.4)"},
                                    "data": [],
                                },
                                "markLine": {"symbol": "none", "data": []},
                            },
                        ],
                    }
                ).classes("w-full h-72")

            ui.separator().classes("mgmt-detail-separator")
            ui.label("CNV events").classes("target-coverage-panel__meta-label mt-2 mb-1")
            cnv_events_summary = ui.label("No CNV events detected").classes(
                "classification-insight-meta mb-2"
            )
        
        # CNV Events Table
        cnv_events_columns = [
            {"name": "chromosome", "label": "Chr", "field": "chromosome", "sortable": True},
            {"name": "event_type", "label": "Event Type", "field": "event_type", "sortable": True},
            {"name": "arm", "label": "Arm", "field": "arm", "sortable": True},
            {"name": "start_mb", "label": "Start (Mb)", "field": "start_mb", "sortable": True, "align": "right"},
            {"name": "end_mb", "label": "End (Mb)", "field": "end_mb", "sortable": True, "align": "right"},
            {"name": "length_mb", "label": "Length (Mb)", "field": "length_mb", "sortable": True, "align": "right"},
            {"name": "mean_cnv_str", "label": "Mean CNV", "field": "mean_cnv_str", "sortable": True, "align": "right"},
            {"name": "confidence", "label": "Confidence", "field": "confidence", "sortable": True, "align": "center"},
            {"name": "proportion_affected", "label": "% Affected", "field": "proportion_affected", "sortable": True, "align": "right"},
            {"name": "genes_str", "label": "Genes", "field": "genes_str"},
        ]
        _, cnv_events_table = styled_table(
            columns=cnv_events_columns, rows=[], pagination=20, class_size="table-xs"
        )
        try:
            cnv_events_table.props('multi-sort rows-per-page-options="[10,20,50,0]"')
        except Exception:
            pass
        ui.separator().classes("mgmt-detail-separator")
        ui.label("Cytoband summary").classes(
            "target-coverage-panel__meta-label mt-2 mb-1"
        )
        cyto_columns = [
            {"name": "chrom", "label": "Chr", "field": "chrom", "sortable": True},
            {"name": "region", "label": "Region", "field": "region", "sortable": True},
            {
                "name": "start_mb",
                "label": "Start (Mb)",
                "field": "start_mb",
                "sortable": True,
                "align": "right",
            },
            {
                "name": "end_mb",
                "label": "End (Mb)",
                "field": "end_mb",
                "sortable": True,
                "align": "right",
            },
            {
                "name": "length_mb",
                "label": "Length (Mb)",
                "field": "length_mb",
                "sortable": True,
                "align": "right",
            },
            {
                "name": "mean_cnv",
                "label": "Mean CNV",
                "field": "mean_cnv",
                "sortable": True,
                "align": "right",
            },
            {
                "name": "state",
                "label": "State",
                "field": "state",
                "sortable": True,
                "align": "center",
            },
            {"name": "genes", "label": "Genes", "field": "genes"},
        ]
        _, cyto_table = styled_table(
            columns=cyto_columns, rows=[], pagination=20, class_size="table-xs"
        )
        try:
            cyto_table.props('multi-sort rows-per-page-options="[10,20,50,0]"')
        except Exception:
            pass

    # Adaptive thinning helpers
    MAX_POINTS_PER_CHART = 10000

    def _get_visible_range(chart, series_list):
        try:
            # Compute overall data span
            min_x = None
            max_x = None
            for s in series_list:
                data = s.get("data") or []
                if not data:
                    continue
                sx = data[0][0] if isinstance(data[0], (list, tuple)) else None
                ex = data[-1][0] if isinstance(data[-1], (list, tuple)) else None
                if sx is None or ex is None:
                    # Fallback compute min/max
                    for p in data:
                        try:
                            x = float(p[0])
                        except Exception:
                            continue
                        min_x = x if min_x is None else min(min_x, x)
                        max_x = x if max_x is None else max(max_x, x)
                else:
                    vmin = min(float(sx), float(ex))
                    vmax = max(float(sx), float(ex))
                    min_x = vmin if min_x is None else min(min_x, vmin)
                    max_x = vmax if max_x is None else max(max_x, vmax)
            dz = None
            try:
                dzo = chart.options.get("dataZoom")
                if isinstance(dzo, list) and dzo:
                    dz = dzo[0]
            except Exception:
                dz = None
            if not dz:
                return (min_x, max_x)
            # Prefer explicit values
            sv = dz.get("startValue") if isinstance(dz, dict) else None
            ev = dz.get("endValue") if isinstance(dz, dict) else None
            if sv is not None or ev is not None:
                left = float(sv) if sv is not None else min_x
                right = float(ev) if ev is not None else max_x
                return (left, right)
            # Fallback to percentage range
            sp = dz.get("start") if isinstance(dz, dict) else None
            ep = dz.get("end") if isinstance(dz, dict) else None
            if (
                (sp is not None or ep is not None)
                and min_x is not None
                and max_x is not None
            ):
                width = (
                    max_x - min_x if max_x is not None and min_x is not None else None
                )
                if width and width > 0:
                    left = min_x + (float(sp or 0) / 100.0) * width
                    right = min_x + (float(ep or 100) / 100.0) * width
                    return (left, right)
            return (min_x, max_x)
        except Exception:
            return (None, None)

    def _evenly_sample(seq, k):
        try:
            n = len(seq)
            if k <= 0 or n <= k:
                return list(seq)
            if k == 1:
                return [seq[n // 2]]
            # Choose k indices evenly across [0, n-1]
            return [seq[int(round(i * (n - 1) / (k - 1)))] for i in range(k)]
        except Exception:
            return list(seq)[:k]

    def _thin_chart_series(chart, max_points: int = MAX_POINTS_PER_CHART) -> None:
        try:
            series = chart.options.get("series", [])
            # Identify data series to thin (exclude overlays)
            data_idx = []
            data_series = []
            for idx, s in enumerate(series):
                name = s.get("name", "")
                if s.get("type") == "scatter" and name not in (
                    "centromeres_highlight",
                    "cytobands_highlight",
                ):
                    data = s.get("data") or []
                    if isinstance(data, list) and data:
                        data_idx.append(idx)
                        data_series.append(s)
            if not data_series:
                return
            x_range = _get_visible_range(chart, data_series)
            # Gather visible counts and data within range
            vis_data = []
            total = 0
            left, right = x_range
            for s in data_series:
                pts = s.get("data") or []
                if left is not None and right is not None:
                    sub = [
                        p
                        for p in pts
                        if isinstance(p, (list, tuple)) and left <= float(p[0]) <= right
                    ]
                else:
                    sub = list(pts)
                vis_data.append(sub)
                total += len(sub)
            if total <= max_points:
                return
            # Allocate budgets proportional to visible counts with a small floor
            budgets = []
            
            for sub in vis_data:
                share = int(max(1, round((len(sub) / total) * max_points)))
                budgets.append(share)
            # Normalize budgets to exactly max_points
            adj = sum(budgets) - max_points
            i = 0
            while adj != 0 and budgets:
                if adj > 0 and budgets[i] > 1:
                    budgets[i] -= 1
                    adj -= 1
                elif adj < 0:
                    budgets[i] += 1
                    adj += 1
                i = (i + 1) % len(budgets)
            # Apply sampling and replace data (preserve points outside range sparsely)
            for (idx, s), sub, k in zip(zip(data_idx, data_series), vis_data, budgets):
                original = s.get("data") or []
                # Keep outside-range points sparsely so context remains when zoomed out/in
                if left is not None and right is not None:
                    outside = [
                        p
                        for p in original
                        if isinstance(p, (list, tuple))
                        and not (left <= float(p[0]) <= right)
                    ]
                    outside_keep = _evenly_sample(
                        outside, max(0, k // 10)
                    )  # at most 10% of budget
                else:
                    outside_keep = []
                inside_keep = _evenly_sample(sub, max(1, k - len(outside_keep)))
                new_data = inside_keep + outside_keep
                chart.options["series"][idx]["data"] = new_data
        except Exception:
            pass

    def _clear_segment_series(chart: Any) -> None:
        """Remove previously rendered segment overlay line series."""
        try:
            base = chart.options.get("series", [])
            chart.options["series"] = [
                s
                for s in base
                if not str(s.get("name", "")).startswith(("seg_abs", "seg_diff"))
            ]
        except Exception:
            pass

    @lru_cache(maxsize=1)
    def _load_centromere_regions() -> Dict[str, List[Tuple[int, int, str]]]:
        """Load centromere/satellite regions from packaged resources.
        Returns mapping: chrom -> list of (start_bp, end_bp, name).
        """
        regions: Dict[str, List[Tuple[int, int, str]]] = {}
        try:
            res_path = (
                importlib_resources.files("robin.resources") / "cenSatRegions.bed"
            )
            with res_path.open("r") as fh:
                for line in fh:
                    parts = line.strip().split("\t")
                    if len(parts) < 4:
                        continue
                    chrom, start, end, name = (
                        parts[0],
                        int(parts[1]),
                        int(parts[2]),
                        parts[3],
                    )
                    regions.setdefault(chrom, []).append((start, end, name))
        except Exception:
            pass
        return regions

    @lru_cache(maxsize=1)
    def _load_cytobands_df() -> pd.DataFrame:
        try:
            res_path = importlib_resources.files("robin.resources") / "cytoBand.txt"
            df = pd.read_csv(
                res_path,
                sep="\t",
                header=None,
                names=["chrom", "start_pos", "end_pos", "name", "stain"],
            )
            return df
        except Exception:
            return pd.DataFrame(columns=["chrom", "start_pos", "end_pos", "name", "stain"])

    def _load_gene_bed(sample_dir: Path = None) -> pd.DataFrame:
        """Load gene BED file based on the analysis panel used for the sample"""
        try:
            # Determine which panel to use
            panel = ""  # No default fallback
            if sample_dir:
                try:
                    master_csv_path = sample_dir / "master.csv"
                    if master_csv_path.exists():
                        import pandas as pd
                        df = pd.read_csv(master_csv_path)
                        if not df.empty and "analysis_panel" in df.columns:
                            panel_val = df.iloc[0]["analysis_panel"]
                            if panel_val and str(panel_val).strip() != "":
                                panel = str(panel_val).strip()
                except Exception:
                    pass
            
            # Map panel to BED filename
            bed_filename = None
            if not panel:
                # No panel found - return empty DataFrame
                return pd.DataFrame(columns=["chrom", "start_pos", "end_pos", "gene"])
            elif panel == "rCNS2":
                bed_filename = "rCNS2_panel_name_uniq.bed"
            elif panel == "AML":
                bed_filename = "AML_panel_name_uniq.bed"
            elif panel == "PanCan":
                bed_filename = "PanCan_panel_name_uniq.bed"
            else:
                # Check for custom panel
                bed_filename = f"{panel}_panel_name_uniq.bed"
            
            # Try to load the panel-specific BED file
            try:
                res_path = importlib_resources.files("robin.resources") / bed_filename
                if res_path.exists():
                    return pd.read_csv(
                        res_path,
                        sep="\t",
                        header=None,
                        names=["chrom", "start_pos", "end_pos", "gene"],
                    )
            except Exception:
                pass
            
            # Fallback to unique_genes.bed if panel-specific file not found
            try:
                res_path = importlib_resources.files("robin.resources") / "unique_genes.bed"
                if res_path.exists():
                    return pd.read_csv(
                        res_path,
                        sep="\t",
                        header=None,
                        names=["chrom", "start_pos", "end_pos", "gene"],
                    )
            except Exception:
                pass
                
        except Exception:
            pass
            
        return pd.DataFrame(columns=["chrom", "start_pos", "end_pos", "gene"])

    def _sex_label(xy_val: Any) -> str:
        try:
            s = str(xy_val).strip().upper()
            if s in ("MALE", "XY"):
                return "Male"
            if s in ("FEMALE", "XX"):
                return "Female"
        except Exception:
            pass
        return "Unknown"

    def _update_ploidy_label(state: Dict[str, Any]) -> None:
        """Compute ploidy / karyotype from sample CNV bins (GUI-only; not written by analysis)."""
        try:
            cnv_map = _cnv_sample_bin_map_from_state(state)
            if not cnv_map:
                cnv_ploidy.set_text("Ploidy (CNV bins): --")
                return
            bw = int(state.get("cnv_dict", {}).get("bin_width", 1_000_000))
            r2_map = _cnv_ref_bin_map_from_state(state)
            est = estimate_ploidy_from_sample_cnv(
                cnv_map,
                bw,
                logging.getLogger("robin.gui.cnv.ploidy"),
                r2_cnv=r2_map,
            )
            cnv_ploidy.set_text(_format_ploidy_gui_line(est))
        except Exception as exc:
            logging.debug("CNV ploidy GUI update failed: %s", exc)
            cnv_ploidy.set_text("Ploidy (CNV bins): unavailable")

    def _update_integer_clonal_label(state: Dict[str, Any]) -> None:
        """Heuristic tumour purity proxy from deviation of residuals from integers (autosomes)."""
        try:
            cnv_map = _cnv_sample_bin_map_from_state(state)
            if not cnv_map:
                cnv_integer_clonal.set_text("Integer clonality (heuristic): --")
                return
            st = launcher._cnv_state.setdefault(str(sample_dir), {})
            if "plot_assumed_ploidy" not in st and "plot_abs_baseline_copy" in st:
                try:
                    b = float(st.get("plot_abs_baseline_copy", 2.0))
                    if b > 0:
                        st["plot_assumed_ploidy"] = 4.0 / b
                except Exception:
                    pass
            ploidy_assumed = float(st.get("plot_assumed_ploidy", 2.0))
            if ploidy_assumed <= 0:
                ploidy_assumed = 2.0
            try:
                ploidy_assumed = float(cnv_assumed_ploidy.value)
                if ploidy_assumed > 0:
                    st["plot_assumed_ploidy"] = ploidy_assumed
            except Exception:
                pass
            binw = int(state.get("cnv_dict", {}).get("bin_width", 1_000_000))
            plot_bw = state.get("plot_bin_width") or binw
            if plot_bw is not None and int(plot_bw) < binw:
                plot_bw = binw
            cnv3_map = _cnv_unwrap_inner_map(state.get("cnv3"))
            cnv2_map = _cnv_ref_bin_map_from_state(state)
            diag = estimate_cnv_integer_clonal_diagnostics(
                cnv_map,
                bin_width_analysis=binw,
                plot_bin_width=int(plot_bw) if plot_bw else binw,
                assumed_ploidy=ploidy_assumed,
                cnv3_map=cnv3_map,
                cnv2_map=cnv2_map,
            )
            cnv_integer_clonal.set_text(_format_integer_clonal_line(diag))
        except Exception as exc:
            logging.debug("CNV integer clonality GUI update failed: %s", exc)
            cnv_integer_clonal.set_text("Integer clonality (heuristic): unavailable")

    def _update_ichor_py_label() -> None:
        """Display persisted ichor_py purity/ploidy estimates if available."""
        try:
            p = sample_dir / "ichor_py_params.json"
            if not p.exists():
                cnv_ichor_py.set_text("Tumour fraction/ploidy (ichor_py): --")
                return
            with p.open("r", encoding="utf-8") as f:
                obj = json.load(f)
            tf = obj.get("tumor_fraction")
            pl = obj.get("tumor_ploidy")
            ag = obj.get("altered_genome_fraction")
            ns = obj.get("n_segments")
            if isinstance(tf, (int, float)) and isinstance(pl, (int, float)):
                text = (
                    f"Tumour fraction/ploidy (ichor_py): {float(tf):.1%}, {float(pl):.2f}"
                )
                pl_hmm = obj.get("tumor_ploidy_hmm_only")
                if isinstance(pl_hmm, (int, float)):
                    text += f" · HMM-only ploidy {float(pl_hmm):.2f}"
                if isinstance(ag, (int, float)):
                    text += f" · altered genome {float(ag):.1%}"
                if isinstance(ns, (int, float)):
                    text += f" · segments {int(ns)}"
                cnv_ichor_py.set_text(text)
            else:
                cnv_ichor_py.set_text("Tumour fraction/ploidy (ichor_py): unavailable")
        except Exception as exc:
            logging.debug("ichor_py GUI label update failed: %s", exc)
            cnv_ichor_py.set_text("Tumour fraction/ploidy (ichor_py): unavailable")

    def _load_ichor_segments_cached(state: Dict[str, Any]) -> Optional[pd.DataFrame]:
        """Load ichor_py segments TSV with mtime cache in state."""
        try:
            seg_path = sample_dir / "ichor_py_segments.tsv"
            if not seg_path.exists():
                state["ichor_segments_df"] = None
                state["ichor_segments_m"] = 0
                return None
            m = seg_path.stat().st_mtime
            if state.get("ichor_segments_m") == m and isinstance(
                state.get("ichor_segments_df"), pd.DataFrame
            ):
                return state.get("ichor_segments_df")
            df = pd.read_csv(seg_path, sep="\t")
            state["ichor_segments_df"] = df
            state["ichor_segments_m"] = m
            return df
        except Exception:
            state["ichor_segments_df"] = None
            return None

    def _downsample_cnv_for_plot(
        values_1d: np.ndarray,
        analysis_bin_width: int,
        plot_bin_width: int,
    ) -> Tuple[np.ndarray, np.ndarray]:
        """Downsample CNV values for display when plot_bin_width > analysis_bin_width."""
        return _downsample_cnv_bins_for_plot(
            values_1d, analysis_bin_width, plot_bin_width
        )

    def _analyze_cytoband_cnv(
        cnv_data: Dict[str, np.ndarray],
        chromosome: str,
        bin_width: int,
        sex_estimate: str,
    ) -> pd.DataFrame:
        try:
            import numpy as _np
            import pandas as _pd

            if not cnv_data or chromosome not in cnv_data:
                return _pd.DataFrame()
            if bin_width > 10_000_000:
                return _pd.DataFrame()

            cyto_df = _load_cytobands_df()
            chromosome_cytobands = cyto_df[cyto_df["chrom"] == chromosome].copy()
            if chromosome_cytobands.empty:
                return _pd.DataFrame()

            # Centromere mask
            centro = _load_centromere_regions()
            mask = _np.ones(len(cnv_data[chromosome]), dtype=bool)
            cent_regions = centro.get(chromosome, [])
            if cent_regions:
                # Combine all annotated centromere/satellite spans
                for s_bp, e_bp, _nm in cent_regions:
                    s_bin = max(0, int(s_bp // bin_width))
                    e_bin = min(len(mask), int(e_bp // bin_width))
                    if e_bin > s_bin:
                        mask[s_bin:e_bin] = False
            chr_cnv = _np.asarray(cnv_data[chromosome])
            chr_cnv = chr_cnv[mask] if mask.any() else _np.asarray(cnv_data[chromosome])

            # Chromosome-wide stats
            chr_mean = float(_np.mean(chr_cnv)) if chr_cnv.size else 0.0
            chr_std = float(_np.std(chr_cnv)) if chr_cnv.size else 1.0

            # Autosomal distribution for whole-chrom thresholds
            chrom_means: List[float] = []
            for ck, arr in cnv_data.items():
                if ck.startswith("chr") and ck[3:].isdigit():
                    a = _np.asarray(arr)
                    if a.size:
                        chrom_means.append(float(_np.mean(a)))
            means_mean = float(_np.mean(chrom_means)) if chrom_means else 0.0
            means_std = float(_np.std(chrom_means)) if chrom_means else 1.0

            # Use the same thresholds as centralized CNV detection for consistency
            from robin.classification_config import get_cnv_thresholds
            gain_threshold, loss_threshold = get_cnv_thresholds(chromosome, sex_estimate)
            
            # Use the same thresholds for cytoband-level analysis
            cyto_gain_th = gain_threshold
            cyto_loss_th = loss_threshold

            # Whole chromosome event detection
            bins_above_gain = float(
                (_np.asarray(cnv_data[chromosome]) > gain_threshold).sum()
            ) / max(1, len(cnv_data[chromosome]))
            bins_below_loss = float(
                (_np.asarray(cnv_data[chromosome]) < loss_threshold).sum()
            ) / max(1, len(cnv_data[chromosome]))
            min_prop = 0.7
            whole_chr_event = False
            whole_chr_state = "NORMAL"
            if bins_above_gain > min_prop:
                whole_chr_event = True
                whole_chr_state = "GAIN"
            elif bins_below_loss > min_prop:
                whole_chr_event = True
                whole_chr_state = "LOSS"

            merged_rows: List[dict] = []
            if whole_chr_event:
                gene_df = _load_gene_bed(sample_dir)
                genes_in_chr = (
                    gene_df[gene_df["chrom"] == chromosome]["gene"].astype(str).tolist()
                )
                merged_rows.append(
                    {
                        "chrom": chromosome,
                        "start_pos": int(chromosome_cytobands["start_pos"].min()),
                        "end_pos": int(chromosome_cytobands["end_pos"].max()),
                        "name": f"{chromosome} WHOLE CHROMOSOME {whole_chr_state}",
                        "mean_cnv": chr_mean,
                        "cnv_state": whole_chr_state,
                        "length": int(chromosome_cytobands["end_pos"].max())
                        - int(chromosome_cytobands["start_pos"].min()),
                        "genes": genes_in_chr,
                    }
                )

            # Group contiguous cytobands by state
            current_group = None
            vals = _np.asarray(cnv_data[chromosome])
            for _, band in chromosome_cytobands.iterrows():
                s_bp = int(band["start_pos"])
                e_bp = int(band["end_pos"])
                s_bin = max(0, s_bp // bin_width)
                e_bin = min(len(vals) - 1, max(0, e_bp // bin_width))
                region = (
                    vals[s_bin : e_bin + 1]
                    if len(vals) and e_bin >= s_bin
                    else _np.array([])
                )
                mean_val = float(_np.mean(region)) if region.size else 0.0
                # Determine cytoband state - always use standard thresholds for regional analysis
                # This ensures we capture all significant regional variations regardless of whole chromosome events
                state = (
                    "GAIN"
                    if mean_val > cyto_gain_th
                    else ("LOSS" if mean_val < cyto_loss_th else "NORMAL")
                )

                if current_group is None:
                    current_group = {
                        "chrom": chromosome,
                        "start_pos": s_bp,
                        "end_pos": e_bp,
                        "bands": [str(band["name"])],
                        "mean_vals": [mean_val],
                        "cnv_state": state,
                    }
                elif state == current_group["cnv_state"]:
                    current_group["end_pos"] = e_bp
                    current_group["bands"].append(str(band["name"]))
                    current_group["mean_vals"].append(mean_val)
                else:
                    # finalize
                    mean_cnv = (
                        float(_np.mean(current_group["mean_vals"]))
                        if current_group["mean_vals"]
                        else 0.0
                    )
                    row = {
                        "chrom": current_group["chrom"],
                        "start_pos": int(current_group["start_pos"]),
                        "end_pos": int(current_group["end_pos"]),
                        "name": f"{current_group['chrom']} {current_group['bands'][0]}-{current_group['bands'][-1]}",
                        "mean_cnv": mean_cnv,
                        "cnv_state": current_group["cnv_state"],
                        "length": int(current_group["end_pos"])
                        - int(current_group["start_pos"]),
                    }
                    if row["cnv_state"] in ("GAIN", "LOSS", "HIGH_GAIN", "DEEP_LOSS"):
                        gene_df = _load_gene_bed(sample_dir)
                        genes = (
                            gene_df[
                                (gene_df["chrom"] == row["chrom"])
                                & (gene_df["start_pos"] <= row["end_pos"])
                                & (gene_df["end_pos"] >= row["start_pos"])
                            ]["gene"]
                            .astype(str)
                            .tolist()
                        )
                        row["genes"] = genes
                    merged_rows.append(row)
                    # start new group
                    current_group = {
                        "chrom": chromosome,
                        "start_pos": s_bp,
                        "end_pos": e_bp,
                        "bands": [str(band["name"])],
                        "mean_vals": [mean_val],
                        "cnv_state": state,
                    }
            # finalize last
            if current_group is not None:
                mean_cnv = (
                    float(_np.mean(current_group["mean_vals"]))
                    if current_group["mean_vals"]
                    else 0.0
                )
                row = {
                    "chrom": current_group["chrom"],
                    "start_pos": int(current_group["start_pos"]),
                    "end_pos": int(current_group["end_pos"]),
                    "name": f"{current_group['chrom']} {current_group['bands'][0]}-{current_group['bands'][-1]}",
                    "mean_cnv": mean_cnv,
                    "cnv_state": current_group["cnv_state"],
                    "length": int(current_group["end_pos"])
                    - int(current_group["start_pos"]),
                }
                if row["cnv_state"] in ("GAIN", "LOSS", "HIGH_GAIN", "DEEP_LOSS"):
                    gene_df = _load_gene_bed(sample_dir)
                    genes = (
                        gene_df[
                            (gene_df["chrom"] == row["chrom"])
                            & (gene_df["start_pos"] <= row["end_pos"])
                            & (gene_df["end_pos"] >= row["start_pos"])
                        ]["gene"]
                        .astype(str)
                        .tolist()
                    )
                    row["genes"] = genes
                merged_rows.append(row)

            df = _pd.DataFrame(merged_rows)
            if not df.empty:
                df = df.sort_values("start_pos")
            return df
        except Exception:
            return pd.DataFrame()

    def _get_cytoband_cnv_summary(
        cnv_data: Dict[str, np.ndarray],
        chromosome: str,
        bin_width: int,
        sex_estimate: str,
    ) -> str:
        try:
            df = _analyze_cytoband_cnv(cnv_data, chromosome, bin_width, sex_estimate)
            if df.empty:
                return "No significant CNV changes detected"
            gains = df[df["cnv_state"] == "GAIN"]
            losses = df[df["cnv_state"] == "LOSS"]
            parts: List[str] = []
            if not gains.empty:
                parts.append(
                    "Gains: "
                    + ", ".join(
                        f"{r['name']} ({r['mean_cnv']:.2f})"
                        for _, r in gains.iterrows()
                    )
                )
            if not losses.empty:
                parts.append(
                    "Losses: "
                    + ", ".join(
                        f"{r['name']} ({r['mean_cnv']:.2f})"
                        for _, r in losses.iterrows()
                    )
                )
            return "\n".join(parts) if parts else "No significant CNV changes detected"
        except Exception:
            return "No CNV data available"

    def _compute_all_cytoband_df(
        cnv_data: Dict[str, np.ndarray], bin_width: int, sex_estimate: str
    ) -> pd.DataFrame:
        try:
            frames: List[pd.DataFrame] = []
            for chrom in natsort.natsorted(cnv_data.keys()):
                if not _cnv_contig_ok(chrom):
                    continue
                df = _analyze_cytoband_cnv(cnv_data, chrom, bin_width, sex_estimate)
                if not df.empty:
                    frames.append(df)
            if frames:
                out = pd.concat(frames, ignore_index=True)
                if not out.empty:
                    # Natural chromosome sort: chr1..chr22, chrX, chrY
                    def _rank(label: Any) -> int:
                        try:
                            s = str(label)
                            if s.startswith("chr"):
                                s = s[3:]
                            mapping = {"X": 23, "Y": 24, "M": 25}
                            return int(s) if s.isdigit() else mapping.get(s, 1000)
                        except Exception:
                            return 1000

                    out["_chrom_rank"] = out["chrom"].map(_rank)
                    out = out.sort_values(["_chrom_rank", "start_pos"]).drop(
                        columns=["_chrom_rank"]
                    )
                return out
            return pd.DataFrame()
        except Exception:
            return pd.DataFrame()

    def _build_cyto_rows(df: pd.DataFrame) -> List[Dict[str, Any]]:
        rows: List[Dict[str, Any]] = []
        try:
            if df is None or df.empty:
                return rows
            for _, r in df.iterrows():
                if r.get("cnv_state") in ("GAIN", "LOSS", "HIGH_GAIN", "DEEP_LOSS"):
                    rows.append(
                        {
                            "chrom": str(r["chrom"]).replace("chr", ""),
                            "region": str(r["name"]).replace(f"{r['chrom']} ", ""),
                            "start_mb": f"{float(r['start_pos'])/1e6:.2f}",
                            "end_mb": f"{float(r['end_pos'])/1e6:.2f}",
                            "length_mb": f"{float(r['length'])/1e6:.2f}",
                            "mean_cnv": f"{float(r['mean_cnv']):.3f}",
                            "state": str(r["cnv_state"]),
                            "genes": ", ".join(
                                r.get("genes", [])
                                if isinstance(r.get("genes"), list)
                                else []
                            ),
                        }
                    )
        except Exception:
            return rows
        return rows

    def _update_cnv_events_analysis(state: Dict[str, Any]) -> None:
        """Update CNV events analysis using centralized classification rules."""
        try:
            cnv_map = state.get("cnv")
            cnv3_map = state.get("cnv3")
            if isinstance(cnv_map, dict) and "cnv" in cnv_map:
                cnv_map = cnv_map["cnv"]
            if isinstance(cnv3_map, dict) and "cnv" in cnv3_map:
                cnv3_map = cnv3_map["cnv"]
            
            if not cnv_map:
                cnv_events_table.rows = []
                cnv_events_summary.set_text("No CNV data available")
                return
            
            binw = state.get("cnv_dict", {}).get("bin_width", 1000000)
            sex_lbl = _sex_label(state.get("xy"))
            
            # Use difference map (CNV3) for calling if available, otherwise absolute
            data = cnv3_map if isinstance(cnv3_map, dict) else cnv_map
            
            if data and binw:
                # Load cytobands and genes
                cyto_df = _load_cytobands_df()
                gene_df = _load_gene_bed(sample_dir)
                
                # Detect CNV events using centralized rules
                events = detect_cnv_events(
                    cnv_data=data,
                    bin_width=int(binw),
                    sex_estimate=sex_lbl,
                    cytobands_df=cyto_df,
                    gene_df=gene_df
                )
                
                # Update events table
                events_rows = []
                for event in events:
                    event_dict = event.to_dict()
                    # Format proportion as percentage
                    event_dict["proportion_affected"] = f"{event.proportion_affected:.1%}"
                    events_rows.append(event_dict)
                
                cnv_events_table.rows = events_rows
                try:
                    cnv_events_table.update()
                except Exception:
                    pass
                
                # Update summary
                summary = get_cnv_summary(events)
                if summary["total_events"] > 0:
                    summary_text = f"Detected {summary['total_events']} CNV events: "
                    parts = []
                    if summary["whole_chromosome_events"]:
                        parts.append(f"{len(summary['whole_chromosome_events'])} whole chromosome")
                    if summary["arm_events"]:
                        parts.append(f"{len(summary['arm_events'])} arm-specific")
                    if summary["gene_containing_events"]:
                        parts.append(f"{summary['total_genes_affected']} genes affected")
                    summary_text += ", ".join(parts)
                    cnv_events_summary.set_text(summary_text)
                else:
                    cnv_events_summary.set_text("No threshold triggered CNV events detected")
            else:
                cnv_events_table.rows = []
                cnv_events_summary.set_text("CNV data not available")
        except Exception as e:
            logging.error(f"Error updating CNV events analysis: {e}")
            cnv_events_table.rows = []
            cnv_events_summary.set_text("Error analyzing CNV events")

    def _render_cnv_from_state(state: Dict[str, Any]) -> None:
        try:
            _clear_segment_series(cnv_abs)
            _clear_segment_series(cnv_diff)
            cnv_map = state.get("cnv")
            cnv3_map = state.get("cnv3")
            if isinstance(cnv_map, dict) and "cnv" in cnv_map:
                cnv_map = cnv_map["cnv"]
            if isinstance(cnv3_map, dict) and "cnv" in cnv3_map:
                cnv3_map = cnv3_map["cnv"]
            if not cnv_map:
                return
            _st_plot = launcher._cnv_state.setdefault(str(sample_dir), {})
            if "plot_assumed_ploidy" not in _st_plot and "plot_abs_baseline_copy" in _st_plot:
                try:
                    b = float(_st_plot.get("plot_abs_baseline_copy", 2.0))
                    if b > 0:
                        _st_plot["plot_assumed_ploidy"] = 4.0 / b
                except Exception:
                    pass
            ploidy_assumed = float(_st_plot.get("plot_assumed_ploidy", 2.0))
            if ploidy_assumed <= 0:
                ploidy_assumed = 2.0
            dark_ui = _is_dark_mode()
            chrom_palette = _cnv_chromosome_scatter_palette(dark_ui)
            col_high, col_low, col_norm = _cnv_value_mode_colors(dark_ui)
            seg_col_abs = "#e5e7eb" if dark_ui else "#111827"
            seg_col_diff = "#cbd5e1" if dark_ui else "#6b7280"
            seg_lbl_col = "#f8fafc" if dark_ui else "#0f172a"
            seg_lbl_bg = "rgba(15, 23, 42, 0.75)" if dark_ui else "rgba(255, 255, 255, 0.9)"
            chrom_divider = _cnv_echart_palette(dark_ui)["muted"]
            binw_analysis = state.get("cnv_dict", {}).get("bin_width", 1_000_000)
            plot_bin_width = state.get("plot_bin_width") or binw_analysis
            if plot_bin_width < binw_analysis:
                plot_bin_width = binw_analysis
            binw = plot_bin_width  # used for x positions and padding in the plot
            # Keep plot bin width dropdown in sync with state
            try:
                pb = state.get("plot_bin_width")
                if pb is None:
                    cnv_plot_bin.value = _PLOT_BIN_KEY_DEFAULT
                else:
                    key = next(
                        (k for k, v in _PLOT_BIN_KEY_TO_BP.items() if v == pb),
                        _PLOT_BIN_KEY_DEFAULT,
                    )
                    cnv_plot_bin.value = key
                cnv_plot_bin.update()
            except Exception:
                pass
            try:
                cnv_segments.value = "show" if state.get("show_segments", True) else "hide"
                cnv_segments.update()
            except Exception:
                pass
            try:
                cnv_segment_source.value = state.get("segment_source", "heuristic")
                cnv_segment_source.update()
            except Exception:
                pass
            selected = state.get("selected_chrom", "All")
            use_log = state.get("y_scale", "linear") == "log"
            show_segments = bool(state.get("show_segments", True))
            segment_source = str(state.get("segment_source", "heuristic")).strip().lower()
            if segment_source not in ("heuristic", "ichor_py"):
                segment_source = "heuristic"
            raw_color_mode = state.get("color_mode", "chromosome")
            # normalize color mode to expected keys
            lval = str(raw_color_mode).strip().lower()
            if lval in ("chromosome", "chromosomes"):
                color_mode = "chromosome"
            elif lval in (
                "value",
                "up/down",
                "updown",
                "up_down",
                "updown ",
                "up down",
            ):
                color_mode = "value"
            else:
                color_mode = "chromosome"
            
            # Show/hide breakpoint density y-axis based on selection
            should_show_breakpoint_density = selected != "All"
            for chart in (cnv_abs, cnv_diff):
                if len(chart.options["yAxis"]) > 1:
                    # Index 1 is the "Breakpoint density" axis
                    chart.options["yAxis"][1]["show"] = should_show_breakpoint_density
            
            cnv_abs.options["yAxis"][0]["type"] = "log" if use_log else "value"
            cnv_abs.options["yAxis"][0]["logBase"] = 10 if use_log else None
            # X-axis is always in genomic base pairs; use actual genome/chromosome length
            # so the scale does not change when plot bin width changes (dataMax would shrink
            # with fewer downsampled points).
            if selected == "All":
                x_axis_max = sum(
                    len(cnv_map[c]) * binw_analysis
                    for c in natsort.natsorted(cnv_map.keys())
                    if _cnv_contig_ok(c)
                )
            else:
                chr_cnv = cnv_map.get(selected)
                x_axis_max = len(chr_cnv) * binw_analysis if chr_cnv is not None else 0
            cnv_abs.options["xAxis"]["min"] = 0
            cnv_abs.options["xAxis"]["max"] = x_axis_max
            cnv_diff.options["xAxis"]["min"] = 0
            cnv_diff.options["xAxis"]["max"] = x_axis_max
            # Clear any previous zoom constraints when viewing All
            if selected == "All":
                try:
                    if (
                        isinstance(cnv_abs.options.get("dataZoom"), list)
                        and cnv_abs.options["dataZoom"]
                    ):
                        dz = cnv_abs.options["dataZoom"][0]
                        dz.pop("startValue", None)
                        dz.pop("endValue", None)
                        dz.update({"start": 0, "end": 100})
                    if (
                        isinstance(cnv_diff.options.get("dataZoom"), list)
                        and cnv_diff.options["dataZoom"]
                    ):
                        dz2 = cnv_diff.options["dataZoom"][0]
                        dz2.pop("startValue", None)
                        dz2.pop("endValue", None)
                        dz2.update({"start": 0, "end": 100})
                except Exception:
                    pass
            logging.debug(
                f"CNV render: selected={selected}, y_scale={state.get('y_scale')}, color_mode={state.get('color_mode')}"
            )

            def _segment_spans_for_plot(
                x_coords: np.ndarray,
                y_vals: np.ndarray,
                *,
                contig: Optional[str] = None,
                contig_offset_bp: float = 0.0,
            ) -> List[Tuple[float, float, float]]:
                if not len(x_coords) or not len(y_vals):
                    return []
                if segment_source == "heuristic":
                    return _segment_spans_from_series(
                        x_coords, y_vals, chunk_bins=5, merge_tol=0.12
                    )
                seg_df = _load_ichor_segments_cached(state)
                if seg_df is None or seg_df.empty or not contig:
                    return []
                try:
                    seg_chr = seg_df[seg_df["chrom"].astype(str) == str(contig)]
                except Exception:
                    return []
                if seg_chr.empty:
                    return []
                x_local = np.asarray(x_coords, dtype=float) - float(contig_offset_bp)
                y_arr = np.asarray(y_vals, dtype=float)
                spans: List[Tuple[float, float, float]] = []
                for _, row in seg_chr.iterrows():
                    try:
                        s0 = float(row["start"])
                        s1 = float(row["end"])
                    except Exception:
                        continue
                    if s1 <= s0:
                        continue
                    mask = (x_local >= s0) & (x_local <= s1)
                    if not np.any(mask):
                        continue
                    sy = float(np.nanmedian(y_arr[mask]))
                    spans.append((s0 + contig_offset_bp, s1 + contig_offset_bp, sy))
                return spans

            # Absolute plot
            series_abs = []
            segment_abs = []
            # Prepare chromosome partitions for labels/areas when viewing All
            chrom_bounds = []  # list of (name, start_bp, end_bp)
            chrom_offsets: Dict[str, float] = {}
            if selected == "All":
                # X-axis is in genomic base pairs; chromosome bounds use actual lengths.
                offset_bp = 0
                for contig, cnv in natsort.natsorted(cnv_map.items()):
                    if not _cnv_contig_ok(contig):
                        continue
                    x_local, vals = _downsample_cnv_for_plot(
                        np.asarray(cnv), binw_analysis, int(plot_bin_width)
                    )
                    vals = _apply_ploidy_rescale(vals, ploidy_assumed)
                    x_global = offset_bp + x_local
                    pts = list(zip(x_global.tolist(), [float(v) for v in vals]))
                    if show_segments and len(vals):
                        spans = _segment_spans_for_plot(
                            x_global,
                            vals,
                            contig=contig,
                            contig_offset_bp=offset_bp,
                        )
                        for si, (sx0, sx1, sy) in enumerate(spans):
                            segment_abs.append(
                                {
                                    "type": "line",
                                    "name": f"seg_abs_{contig}",
                                    "showSymbol": False,
                                    "symbol": "none",
                                    "lineStyle": {"color": seg_col_abs, "width": 2, "opacity": 0.9},
                                    "label": {
                                        "show": False,
                                        "formatter": f"{sy:.2f}",
                                        "color": seg_lbl_col,
                                        "fontSize": 11,
                                        "fontWeight": "bold",
                                        "backgroundColor": seg_lbl_bg,
                                        "padding": [2, 4],
                                        "borderRadius": 3,
                                    },
                                    "silent": True,
                                    "data": [[sx0, sy], [sx1, sy]],
                                    "z": 5,
                                }
                            )
                    start_bp = offset_bp
                    end_bp = offset_bp + len(cnv) * binw_analysis
                    chrom_offsets[contig] = start_bp
                    chrom_bounds.append((contig, start_bp, end_bp))
                    offset_bp = end_bp
                    if color_mode == "chromosome":
                        ci = len(
                            [s for s in series_abs if s.get("type") == "scatter"]
                        )
                        series_abs.append(
                            {
                                "type": "scatter",
                                "name": contig,
                                "symbolSize": 3,
                                "itemStyle": {
                                    "color": chrom_palette[
                                        ci % len(chrom_palette)
                                    ]
                                },
                                "data": pts,
                            }
                        )
                    else:
                        try:
                            autosome_vals = [
                                float(v) * (ploidy_assumed / 2.0)
                                for k, arr in cnv_map.items()
                                if k.startswith("chr") and k[3:].isdigit()
                                for v in arr
                            ]
                            mean_val = (
                                float(np.mean(autosome_vals)) if autosome_vals else 2.0
                            )
                            std_val = (
                                float(np.std(autosome_vals)) if autosome_vals else 1.0
                            )
                        except Exception:
                            mean_val, std_val = 2.0, 1.0
                        high, low, norm = [], [], []
                        for xi, vi in pts:
                            z = (vi - mean_val) / std_val if std_val > 0 else 0.0
                            (high if z > 0.5 else low if z < -0.5 else norm).append(
                                [xi, vi]
                            )
                        if high:
                            series_abs.append(
                                {
                                    "type": "scatter",
                                    "name": f"High {contig}",
                                    "symbolSize": 4,
                                    "itemStyle": {"color": col_high},
                                    "data": high,
                                }
                            )
                        if low:
                            series_abs.append(
                                {
                                    "type": "scatter",
                                    "name": f"Low {contig}",
                                    "symbolSize": 4,
                                    "itemStyle": {"color": col_low},
                                    "data": low,
                                }
                            )
                        if norm:
                            series_abs.append(
                                {
                                    "type": "scatter",
                                    "name": f"Normal {contig}",
                                    "symbolSize": 2,
                                    "itemStyle": {"color": col_norm},
                                    "data": norm,
                                }
                            )
            else:
                cnv = cnv_map.get(selected)
                if cnv is not None:
                    x_local, vals = _downsample_cnv_for_plot(
                        np.asarray(cnv), binw_analysis, int(plot_bin_width)
                    )
                    vals = _apply_ploidy_rescale(vals, ploidy_assumed)
                    pts = list(zip(x_local.tolist(), [float(v) for v in vals]))
                    if show_segments and len(vals):
                        spans = _segment_spans_for_plot(x_local, vals, contig=selected)
                        for si, (sx0, sx1, sy) in enumerate(spans):
                            show_label = (si % 2 == 0)
                            segment_abs.append(
                                {
                                    "type": "line",
                                    "name": "seg_abs",
                                    "showSymbol": False,
                                    "symbol": "none",
                                    "lineStyle": {"color": seg_col_abs, "width": 2, "opacity": 0.9},
                                    "label": {
                                        "show": show_label,
                                        "position": "end",
                                        "formatter": f"{sy:.2f}",
                                        "color": seg_lbl_col,
                                        "fontSize": 11,
                                        "fontWeight": "bold",
                                        "backgroundColor": seg_lbl_bg,
                                        "padding": [2, 4],
                                        "borderRadius": 3,
                                    },
                                    "silent": True,
                                    "data": [[sx0, sy], [sx1, sy]],
                                    "z": 5,
                                }
                            )
                    if color_mode == "chromosome":
                        series_abs.append(
                            {
                                "type": "scatter",
                                "name": selected,
                                "symbolSize": 3,
                                "itemStyle": {"color": chrom_palette[0]},
                                "data": pts,
                            }
                        )
                    else:
                        expected = float(ploidy_assumed)
                        try:
                            if selected in ("chrX", "chrY") and str(
                                state.get("xy", "")
                            ).upper().startswith("MALE"):
                                expected = float(ploidy_assumed) / 2.0
                        except Exception:
                            pass
                        vals = [v for _, v in pts]
                        std_val = float(np.std(vals)) if vals else 1.0
                        high, low, norm = [], [], []
                        for xi, vi in pts:
                            z = (vi - expected) / std_val if std_val > 0 else 0
                            (high if z > 0.5 else low if z < -0.5 else norm).append(
                                [xi, vi]
                            )
                        if high:
                            series_abs.append(
                                {
                                    "type": "scatter",
                                    "name": f"High {selected}",
                                    "symbolSize": 4,
                                    "itemStyle": {"color": col_high},
                                    "data": high,
                                }
                            )
                        if low:
                            series_abs.append(
                                {
                                    "type": "scatter",
                                    "name": f"Low {selected}",
                                    "symbolSize": 4,
                                    "itemStyle": {"color": col_low},
                                    "data": low,
                                }
                            )
                        if norm:
                            series_abs.append(
                                {
                                    "type": "scatter",
                                    "name": f"Normal {selected}",
                                    "symbolSize": 2,
                                    "itemStyle": {"color": col_norm},
                                    "data": norm,
                                }
                            )
            # Preserve highlight series (centromeres, cytobands) and replace data series only
            keep = [
                s
                for s in cnv_abs.options["series"]
                if s.get("name") in ("centromeres_highlight", "cytobands_highlight")
            ]
            cnv_abs.options["series"] = series_abs + segment_abs + keep
            # Build background chromosome areas and vertical labels when showing All
            try:
                if selected == "All" and chrom_bounds:
                    # Alternating shaded bands per chromosome for readability
                    areas_data = []
                    lines_data = []
                    for contig, start_bp, end_bp in chrom_bounds:
                        areas_data.append(
                            [{"xAxis": float(start_bp)}, {"xAxis": float(end_bp)}]
                        )
                        # Label at the center of chromosome region
                        center_bp = (start_bp + end_bp) / 2
                        lines_data.append(
                            {
                                "xAxis": float(center_bp),
                                "lineStyle": {
                                    "type": "dashed",
                                    "color": chrom_divider,
                                },
                                "label": {
                                    "show": True,
                                    "formatter": contig,
                                    "color": chrom_divider,
                                },
                            }
                        )

                    # Helper to set overlays by series name regardless of index
                    def _apply_overlays(
                        chart, band_areas, band_lines, centro_areas_list
                    ):
                        try:
                            idx_cyto = next(
                                (
                                    i
                                    for i, s in enumerate(chart.options["series"])
                                    if s.get("name") == "cytobands_highlight"
                                ),
                                None,
                            )
                            idx_centro = next(
                                (
                                    i
                                    for i, s in enumerate(chart.options["series"])
                                    if s.get("name") == "centromeres_highlight"
                                ),
                                None,
                            )
                            if idx_cyto is not None:
                                # Remove background shading per request; keep dashed vertical lines only
                                chart.options["series"][idx_cyto]["markArea"][
                                    "data"
                                ] = []
                                chart.options["series"][idx_cyto].setdefault(
                                    "markLine", {}
                                )
                                chart.options["series"][idx_cyto]["markLine"][
                                    "data"
                                ] = band_lines
                                # Disable animation for marker lines so they appear instantly
                                chart.options["series"][idx_cyto]["markLine"][
                                    "animation"
                                ] = False
                            # Do not show centromeres in All view
                            if idx_centro is not None:
                                chart.options["series"][idx_centro]["markArea"][
                                    "data"
                                ] = []
                        except Exception:
                            pass

                    _apply_overlays(cnv_abs, areas_data, lines_data, [])
                    _apply_overlays(cnv_diff, areas_data, lines_data, [])
                else:
                    # Clear overlays when focusing on a single chromosome
                    def _clear_overlays(chart):
                        try:
                            for s in chart.options["series"]:
                                if s.get("name") in (
                                    "cytobands_highlight",
                                    "centromeres_highlight",
                                ):
                                    if "markArea" in s and "data" in s["markArea"]:
                                        s["markArea"]["data"] = []
                                    if "markLine" in s and "data" in s["markLine"]:
                                        s["markLine"]["data"] = []
                        except Exception:
                            pass

                    _clear_overlays(cnv_abs)
                    _clear_overlays(cnv_diff)
                    # In single-chromosome view, draw centromeres, cytobands (by state), genes, and breakpoint candidates
                    if selected != "All":
                        try:
                            centro = _load_centromere_regions()
                            idx_centro_abs = next(
                                (
                                    i
                                    for i, s in enumerate(cnv_abs.options["series"])
                                    if s.get("name") == "centromeres_highlight"
                                ),
                                None,
                            )
                            idx_centro_diff = next(
                                (
                                    i
                                    for i, s in enumerate(cnv_diff.options["series"])
                                    if s.get("name") == "centromeres_highlight"
                                ),
                                None,
                            )
                            areas = []
                            for s, e, _n in centro.get(selected, []):
                                areas.append([{"xAxis": float(s)}, {"xAxis": float(e)}])
                            if idx_centro_abs is not None:
                                cnv_abs.options["series"][idx_centro_abs]["markArea"][
                                    "data"
                                ] = areas
                            if idx_centro_diff is not None:
                                cnv_diff.options["series"][idx_centro_diff]["markArea"][
                                    "data"
                                ] = areas
                        except Exception:
                            pass
                        # Cytobands colored by CNV state (from CNV3 values)
                        try:
                            idx_cyto_abs = next(
                                (
                                    i
                                    for i, s in enumerate(cnv_abs.options["series"])
                                    if s.get("name") == "cytobands_highlight"
                                ),
                                None,
                            )
                            if (
                                idx_cyto_abs is not None
                                and cnv3_map
                                and selected in cnv3_map
                            ):
                                cyto_df = _load_cytobands_df()
                                bands = cyto_df[cyto_df["chrom"] == selected]
                                vals = (
                                    np.array(cnv3_map[selected])
                                    if isinstance(
                                        cnv3_map[selected], (list, np.ndarray)
                                    )
                                    else np.array([])
                                )
                                band_areas = []
                                
                                # Get CNV events for this chromosome to highlight significant events
                                events = []
                                try:
                                    sex_lbl = _sex_label(state.get("xy"))
                                    gene_df = _load_gene_bed(sample_dir)
                                    events = detect_cnv_events(
                                        cnv_data={selected: vals},
                                        bin_width=int(binw_analysis),
                                        sex_estimate=sex_lbl,
                                        cytobands_df=cyto_df,
                                        gene_df=gene_df
                                    )
                                except Exception:
                                    pass
                                
                                # Create event lookup for highlighting
                                event_regions = {}
                                for event in events:
                                    key = f"{event.start_pos}-{event.end_pos}"
                                    event_regions[key] = event
                                
                                for _, row in bands.iterrows():
                                    s_bp, e_bp = int(row["start_pos"]), int(row["end_pos"])
                                    s_bin = max(0, s_bp // binw_analysis)
                                    e_bin = min(len(vals) - 1, max(0, e_bp // binw_analysis))
                                    if len(vals) > 0 and e_bin >= s_bin:
                                        mean_val = float(
                                            np.mean(vals[s_bin : e_bin + 1])
                                        )
                                    else:
                                        mean_val = 0.0
                                    
                                    # Check if this region has a significant CNV event
                                    region_key = f"{s_bp}-{e_bp}"
                                    event = event_regions.get(region_key)
                                    
                                    fill_neutral = (
                                        "rgba(255, 255, 255, 0.07)"
                                        if dark_ui
                                        else "rgba(0, 0, 0, 0.03)"
                                    )
                                    if event:
                                        # Highlight significant events with stronger colors
                                        if event.event_type in ("GAIN", "WHOLE_CHR_GAIN"):
                                            color = "rgba(52, 199, 89, 0.3)"  # gains
                                        elif event.event_type in ("LOSS", "WHOLE_CHR_LOSS"):
                                            color = "rgba(255, 45, 85, 0.3)"  # losses
                                        else:
                                            color = fill_neutral
                                    else:
                                        # Standard cytoband coloring
                                        if mean_val > 0.5:
                                            color = "rgba(52, 199, 89, 0.12)"
                                        elif mean_val < -0.5:
                                            color = "rgba(255, 45, 85, 0.12)"
                                        else:
                                            color = fill_neutral

                                    band_areas.append(
                                        [
                                            {
                                                "name": str(row["name"]),
                                                "xAxis": float(s_bp),
                                                "itemStyle": {"color": color},
                                                "label": {
                                                    "show": True,
                                                    "position": "insideTop",
                                                    "color": (
                                                        "#cbd5e1"
                                                        if dark_ui
                                                        else "#555"
                                                    ),
                                                    "fontSize": 11,
                                                },
                                            },
                                            {"xAxis": float(e_bp)},
                                        ]
                                    )
                                cnv_abs.options["series"][idx_cyto_abs]["markArea"][
                                    "data"
                                ] = band_areas
                        except Exception:
                            pass
                        # Genes of interest and gene selector options
                        try:
                            gene_df = _load_gene_bed(sample_dir)
                            gchr = gene_df[gene_df["chrom"] == selected]
                            # limit to reduce clutter; still add many labels
                            gene_opts = {"All": "All"}
                            for _, gr in gchr.iterrows():
                                gene_opts[str(gr["gene"])] = str(gr["gene"])
                            try:
                                cnv_gene_select.set_options(gene_opts)
                            except Exception as e:
                                pass
                            
                            # annotate genes on main series
                            if series_abs:
                                main = series_abs[0]
                                # attach gene regions via markArea on main series after replacement
                                mark = []
                                for _, gr in gchr.iterrows():
                                    mark.append(
                                        [
                                            {
                                                "name": str(gr["gene"]),
                                                "xAxis": float(gr["start_pos"]),
                                                "label": {
                                                    "position": "insideTop",
                                                    "color": (
                                                        "#e2e8f0"
                                                        if dark_ui
                                                        else "#000"
                                                    ),
                                                    "fontSize": 11,
                                                }
                                            },
                                            {"xAxis": float(gr["end_pos"])},
                                        ]
                                    )
                                main.setdefault("markArea", {"data": []})
                                main["markArea"]["data"] = (
                                    main["markArea"]["data"] or []
                                ) + mark
                                # ensure series_abs[0] updated
                                series_abs[0] = main
                            
                        except Exception as e:
                            pass
                        
                        # Breakpoint candidates as dashed vertical lines
                        try:
                            idx_cyto_abs = next(
                                (
                                    i
                                    for i, s in enumerate(cnv_abs.options["series"])
                                    if s.get("name") == "cytobands_highlight"
                                ),
                                None,
                            )
                            if (
                                idx_cyto_abs is not None
                                and state.get("bp_array") is not None
                            ):
                                cnv_abs.options["series"][idx_cyto_abs].setdefault(
                                    "markLine", {"symbol": "none", "data": []}
                                )
                                if state.get("show_bp", True):
                                    arr = state["bp_array"]
                                    pos = [
                                        int(r["end_pos"])
                                        for r in arr
                                        if r["name"] == selected
                                    ]
                                    lines = [
                                        {
                                            "xAxis": float(p),
                                            "lineStyle": {
                                                "type": "dashed",
                                                "color": "#E0162B",
                                            },
                                        }
                                        for p in pos
                                    ]
                                    cnv_abs.options["series"][idx_cyto_abs]["markLine"][
                                        "data"
                                    ] = lines
                                else:
                                    cnv_abs.options["series"][idx_cyto_abs]["markLine"][
                                        "data"
                                    ] = []
                        except Exception:
                            pass
                # debug label removed
            except Exception:
                pass
            # Adaptive thinning based on current zoom and cap total points
            _thin_chart_series(cnv_abs, MAX_POINTS_PER_CHART)
            
            # Apply gene zoom before updating chart
            try:
                sel_gene = launcher._cnv_state.setdefault(
                    str(sample_dir), {}
                ).get("selected_gene", "All")
                
                if sel_gene and sel_gene != "All":
                    # Load gene data for zoom if not already loaded
                    gene_df = _load_gene_bed(sample_dir)
                    gchr = gene_df[gene_df["chrom"] == selected]
                    
                    row = gchr[gchr["gene"] == sel_gene]
                    if not row.empty:
                        s_bp = int(row.iloc[0]["start_pos"])
                        e_bp = int(row.iloc[0]["end_pos"])
                        # Use 10x analysis bin width for padding (in bp)
                        pad = 10 * binw_analysis
                        zoom_start = max(0, s_bp - pad)
                        zoom_end = e_bp + pad
                        try:
                            cnv_abs.options["dataZoom"][0].update(
                                {
                                    "startValue": zoom_start,
                                    "endValue": zoom_end,
                                    "start": None,  # Remove percentage-based zoom
                                    "end": None,    # Remove percentage-based zoom
                                }
                            )
                        except Exception as e:
                            pass
                else:
                    # Reset zoom when "All" is selected
                    try:
                        if isinstance(cnv_abs.options.get("dataZoom"), list) and cnv_abs.options["dataZoom"]:
                            dz = cnv_abs.options["dataZoom"][0]
                            dz.pop("startValue", None)
                            dz.pop("endValue", None)
                            dz.update({"start": 0, "end": 100})
                    except Exception as e:
                        pass
            except Exception as e:
                pass
            
            try:
                abs_title = "CNV scatter plot"
                if abs(ploidy_assumed - 2.0) > 1e-6:
                    abs_title = (
                        f"CNV scatter (rescaled ×{ploidy_assumed:.2f}/2)"
                    )
                cnv_abs.options["title"]["text"] = abs_title
            except Exception:
                pass

            _apply_cnv_echart_chrome(cnv_abs, _is_dark_mode())
            cnv_abs.update()
            # Difference plot
            if cnv3_map:
                series_diff = []
                segment_diff = []
                if selected == "All":
                    offset_bp = 0
                    dj = 0
                    for contig, cnv in natsort.natsorted(cnv3_map.items()):
                        if not _cnv_contig_ok(contig):
                            continue
                        x_local, vals = _downsample_cnv_for_plot(
                            np.asarray(cnv), binw_analysis, int(plot_bin_width)
                        )
                        vals = _apply_ploidy_rescale(vals, ploidy_assumed)
                        x_global = offset_bp + x_local
                        pts = list(zip(x_global.tolist(), [float(v) for v in vals]))
                        if show_segments and len(vals):
                            spans = _segment_spans_for_plot(
                                x_global,
                                vals,
                                contig=contig,
                                contig_offset_bp=offset_bp,
                            )
                            for si, (sx0, sx1, sy) in enumerate(spans):
                                segment_diff.append(
                                    {
                                        "type": "line",
                                        "name": f"seg_diff_{contig}",
                                        "showSymbol": False,
                                        "symbol": "none",
                                        "lineStyle": {"color": seg_col_diff, "width": 2, "opacity": 0.9},
                                        "label": {
                                            "show": False,
                                            "formatter": f"{sy:.2f}",
                                            "color": seg_lbl_col,
                                            "fontSize": 11,
                                            "fontWeight": "bold",
                                            "backgroundColor": seg_lbl_bg,
                                            "padding": [2, 4],
                                            "borderRadius": 3,
                                        },
                                        "silent": True,
                                        "data": [[sx0, sy], [sx1, sy]],
                                        "z": 5,
                                    }
                                )
                        offset_bp += len(cnv) * binw_analysis
                        series_diff.append(
                            {
                                "type": "scatter",
                                "name": contig,
                                "symbolSize": 3,
                                "itemStyle": {
                                    "color": chrom_palette[
                                        dj % len(chrom_palette)
                                    ]
                                },
                                "data": pts,
                            }
                        )
                        dj += 1
                else:
                    cnv = cnv3_map.get(selected)
                    if cnv is not None:
                        x_local, vals = _downsample_cnv_for_plot(
                            np.asarray(cnv), binw_analysis, int(plot_bin_width)
                        )
                        vals = _apply_ploidy_rescale(vals, ploidy_assumed)
                        pts = list(zip(x_local.tolist(), [float(v) for v in vals]))
                        if show_segments and len(vals):
                            spans = _segment_spans_for_plot(x_local, vals, contig=selected)
                            for si, (sx0, sx1, sy) in enumerate(spans):
                                show_label = (si % 2 == 0)
                                segment_diff.append(
                                    {
                                        "type": "line",
                                        "name": "seg_diff",
                                        "showSymbol": False,
                                        "symbol": "none",
                                        "lineStyle": {"color": seg_col_diff, "width": 2, "opacity": 0.9},
                                        "label": {
                                            "show": show_label,
                                            "position": "end",
                                            "formatter": f"{sy:.2f}",
                                            "color": seg_lbl_col,
                                            "fontSize": 11,
                                            "fontWeight": "bold",
                                            "backgroundColor": seg_lbl_bg,
                                            "padding": [2, 4],
                                            "borderRadius": 3,
                                        },
                                        "silent": True,
                                        "data": [[sx0, sy], [sx1, sy]],
                                        "z": 5,
                                    }
                                )
                        series_diff.append(
                            {
                                "type": "scatter",
                                "name": selected,
                                "symbolSize": 3,
                                "itemStyle": {"color": chrom_palette[0]},
                                "data": pts,
                            }
                        )
                # Preserve highlight series by name
                try:
                    base_series = cnv_diff.options["series"]
                    keep = [
                        s
                        for s in base_series
                        if s.get("name")
                        in ("centromeres_highlight", "cytobands_highlight")
                    ]
                except Exception:
                    keep = []
                cnv_diff.options["series"] = series_diff + segment_diff + keep
                _thin_chart_series(cnv_diff, MAX_POINTS_PER_CHART)
                try:
                    diff_title = "Difference plot"
                    if abs(ploidy_assumed - 2.0) > 1e-6:
                        diff_title = (
                            f"Difference plot (rescaled ×{ploidy_assumed:.2f}/2)"
                        )
                    cnv_diff.options["title"]["text"] = diff_title
                except Exception:
                    pass
                
                # Apply gene zoom to difference chart before updating
                try:
                    sel_gene = launcher._cnv_state.setdefault(
                        str(sample_dir), {}
                    ).get("selected_gene", "All")
                    
                    if sel_gene and sel_gene != "All":
                        gene_df = _load_gene_bed(sample_dir)
                        gchr = gene_df[gene_df["chrom"] == selected]
                        
                        row = gchr[gchr["gene"] == sel_gene]
                        if not row.empty:
                            s_bp = int(row.iloc[0]["start_pos"])
                            e_bp = int(row.iloc[0]["end_pos"])
                            # Use 10x analysis bin width for padding (in bp)
                            pad = 10 * binw_analysis
                            zoom_start = max(0, s_bp - pad)
                            zoom_end = e_bp + pad
                            try:
                                cnv_diff.options["dataZoom"][0].update(
                                    {
                                        "startValue": zoom_start,
                                        "endValue": zoom_end,
                                        "start": None,  # Remove percentage-based zoom
                                        "end": None,    # Remove percentage-based zoom
                                    }
                                )
                            except Exception as e:
                                pass
                    else:
                        # Reset zoom when "All" is selected
                        try:
                            if isinstance(cnv_diff.options.get("dataZoom"), list) and cnv_diff.options["dataZoom"]:
                                dz = cnv_diff.options["dataZoom"][0]
                                dz.pop("startValue", None)
                                dz.pop("endValue", None)
                                dz.update({"start": 0, "end": 100})
                        except Exception as e:
                            pass
                except Exception as e:
                    pass
                
                _apply_cnv_echart_chrome(cnv_diff, _is_dark_mode())
                cnv_diff.update()
            else:
                try:
                    keep = [
                        s
                        for s in (cnv_diff.options.get("series") or [])
                        if s.get("name") in ("centromeres_highlight", "cytobands_highlight")
                    ]
                    cnv_diff.options["series"] = keep
                    cnv_diff.options["title"]["text"] = "Difference plot"
                except Exception:
                    pass
                _apply_cnv_echart_chrome(cnv_diff, _is_dark_mode())
                cnv_diff.update()

            # Cytoband CNV table update (whole-genome table with per-chromosome subsetting)
            try:
                selected = state.get("selected_chrom", "All")
                binw = state.get("cnv_dict", {}).get("bin_width", 1_000_000)
                sex_lbl = _sex_label(state.get("xy"))
                # Prefer difference map (CNV3) for calling; fall back to absolute
                if isinstance(cnv3_map, dict):
                    data = cnv3_map
                    source = "cnv3"
                else:
                    data = cnv_map if isinstance(cnv_map, dict) else None
                    source = "cnv"
                if data and binw:
                    cache_key = (
                        f"{source}:{state.get(source+'_m')}:{int(binw)}:{sex_lbl}"
                    )
                    if state.get("cyto_cache_key") != cache_key:
                        df_all = _compute_all_cytoband_df(data, int(binw), sex_lbl)
                        state["cyto_df_all"] = df_all
                        state["cyto_cache_key"] = cache_key
                    df_all = state.get("cyto_df_all")
                    if isinstance(df_all, pd.DataFrame) and not df_all.empty:
                        if selected and selected != "All":
                            df_show = df_all[df_all["chrom"] == selected]
                        else:
                            df_show = df_all
                        cyto_table.rows = _build_cyto_rows(df_show)
                        try:
                            cyto_table.update()
                        except Exception:
                            pass
                        if selected and selected != "All":
                            cyto_summary.set_text(
                                _get_cytoband_cnv_summary(
                                    data, selected, int(binw), sex_lbl
                                )
                            )
                        else:
                            cyto_summary.set_text(
                                f"Whole genome cytoband events: {len(cyto_table.rows)}"
                            )
                    else:
                        cyto_table.rows = []
                        try:
                            cyto_table.update()
                        except Exception:
                            pass
                        cyto_summary.set_text("No significant CNV changes detected")
                else:
                    cyto_summary.set_text("CNV data not available")
            except Exception:
                pass
        except Exception:
            pass

    def _prepare_cnv_refresh() -> Optional[Dict[str, Any]]:
        """Main-thread debounce, UI sync, mtime checks. Returns None if no work needed."""
        key = str(sample_dir)
        state = launcher._cnv_state.get(key, {})

        current_time = time.time()
        last_refresh = state.get("_last_refresh", 0)
        if current_time - last_refresh < 0.1:
            return None
        state["_last_refresh"] = current_time

        ui_changed = False
        try:
            ui_sel = getattr(cnv_chrom_select, "value", None)
            if ui_sel and ui_sel != state.get("selected_chrom"):
                state["selected_chrom"] = ui_sel
                ui_changed = True
            ui_scale = getattr(cnv_scale, "value", None)
            if ui_scale and ui_scale != state.get("y_scale"):
                state["y_scale"] = ui_scale
                ui_changed = True
            ui_bp = getattr(cnv_bp, "value", None)
            if ui_bp is not None:
                vbp = str(ui_bp).strip().lower()
                desired = vbp in ("show", "true", "1")
                if desired != state.get("show_bp", True):
                    state["show_bp"] = desired
                    ui_changed = True
            ui_color = getattr(cnv_color, "value", None)
            current_color_mode = state.get("color_mode", "chromosome")
            if ui_color and ui_color != current_color_mode:
                state["color_mode"] = ui_color
                ui_changed = True
            ui_plot_bin = getattr(cnv_plot_bin, "value", None)
            if ui_plot_bin is not None:
                want_bin = _PLOT_BIN_KEY_TO_BP.get(
                    ui_plot_bin, _PLOT_BIN_KEY_TO_BP[_PLOT_BIN_KEY_DEFAULT]
                )
                if want_bin != state.get("plot_bin_width"):
                    state["plot_bin_width"] = want_bin
                    ui_changed = True
            ui_segments = getattr(cnv_segments, "value", None)
            if ui_segments is not None:
                show_seg = ui_segments == "show"
                if show_seg != state.get("show_segments", True):
                    state["show_segments"] = show_seg
                    ui_changed = True
            ui_seg_source = getattr(cnv_segment_source, "value", None)
            if ui_seg_source is not None:
                ssrc = str(ui_seg_source).strip().lower()
                if ssrc not in ("heuristic", "ichor_py"):
                    ssrc = "heuristic"
                if ssrc != state.get("segment_source", "heuristic"):
                    state["segment_source"] = ssrc
                    ui_changed = True
        except Exception:
            pass

        is_fresh_visit = "last_visit_time" not in state
        if is_fresh_visit:
            state["last_visit_time"] = time.time()

        cnv_npy = sample_dir / "CNV.npy"
        cnv2_npy = sample_dir / "CNV2.npy"
        cnv3_npy = sample_dir / "CNV3.npy"
        cnv_dict_npy = sample_dir / "CNV_dict.npy"
        data_array_npy = sample_dir / "cnv_data_array.npy"
        xy_pkl = sample_dir / "XYestimate.pkl"
        ichor_params_json = sample_dir / "ichor_py_params.json"
        ichor_seg_tsv = sample_dir / "ichor_py_segments.tsv"

        cnv_npy_mtime = cnv_npy.stat().st_mtime if cnv_npy.exists() else 0
        cnv2_npy_mtime = cnv2_npy.stat().st_mtime if cnv2_npy.exists() else 0
        cnv3_npy_mtime = cnv3_npy.stat().st_mtime if cnv3_npy.exists() else 0
        cnv_dict_npy_mtime = cnv_dict_npy.stat().st_mtime if cnv_dict_npy.exists() else 0
        data_array_npy_mtime = data_array_npy.stat().st_mtime if data_array_npy.exists() else 0
        xy_pkl_mtime = xy_pkl.stat().st_mtime if xy_pkl.exists() else 0
        ichor_params_mtime = (
            ichor_params_json.stat().st_mtime if ichor_params_json.exists() else 0
        )
        ichor_seg_mtime = ichor_seg_tsv.stat().st_mtime if ichor_seg_tsv.exists() else 0

        prev_cnv_npy_mtime = state.get("cnv_m", 0)
        prev_cnv2_npy_mtime = state.get("cnv2_m", 0)
        prev_cnv3_npy_mtime = state.get("cnv3_m", 0)
        prev_cnv_dict_npy_mtime = state.get("dict_m", 0)
        prev_data_array_npy_mtime = state.get("bp_array_mtime", 0)
        prev_xy_pkl_mtime = state.get("xy_m", 0)
        prev_ichor_params_mtime = state.get("ichor_params_m", 0)
        prev_ichor_seg_mtime = state.get("ichor_seg_m", 0)

        cnv_npy_changed = prev_cnv_npy_mtime != cnv_npy_mtime
        cnv2_npy_changed = prev_cnv2_npy_mtime != cnv2_npy_mtime
        cnv3_npy_changed = prev_cnv3_npy_mtime != cnv3_npy_mtime
        cnv_dict_npy_changed = prev_cnv_dict_npy_mtime != cnv_dict_npy_mtime
        data_array_npy_changed = prev_data_array_npy_mtime != data_array_npy_mtime
        xy_pkl_changed = prev_xy_pkl_mtime != xy_pkl_mtime
        ichor_params_changed = prev_ichor_params_mtime != ichor_params_mtime
        ichor_seg_changed = prev_ichor_seg_mtime != ichor_seg_mtime
        ichor_outputs_changed = ichor_params_changed or ichor_seg_changed

        force_color_refresh = state.get("_force_color_refresh", False)
        force_gene_refresh = state.get("_force_gene_refresh", False)
        force_chrom_refresh = state.get("_force_chrom_refresh", False)

        files_changed = (
            cnv_npy_changed
            or cnv2_npy_changed
            or cnv3_npy_changed
            or cnv_dict_npy_changed
            or data_array_npy_changed
            or xy_pkl_changed
            or ichor_outputs_changed
        )

        needs_update = (
            is_fresh_visit
            or files_changed
            or ui_changed
            or force_color_refresh
            or force_gene_refresh
            or force_chrom_refresh
            or not state.get("_rendered_once")
        )

        if not needs_update:
            logging.debug("[CNV] ⏭ Skipping CNV update - no changes detected")
            launcher._cnv_state[key] = state
            return None

        reasons = []
        if is_fresh_visit:
            reasons.append("fresh_visit")
        if cnv_npy_changed:
            reasons.append("CNV.npy")
        if cnv2_npy_changed:
            reasons.append("CNV2.npy")
        if cnv3_npy_changed:
            reasons.append("CNV3.npy")
        if cnv_dict_npy_changed:
            reasons.append("CNV_dict.npy")
        if data_array_npy_changed:
            reasons.append("cnv_data_array.npy")
        if xy_pkl_changed:
            reasons.append("XYestimate.pkl")
        if ichor_params_changed:
            reasons.append("ichor_py_params.json")
        if ichor_seg_changed:
            reasons.append("ichor_py_segments.tsv")
        if ui_changed:
            reasons.append("ui_changed")
        if force_color_refresh:
            reasons.append("force_color_refresh")
        if force_gene_refresh:
            reasons.append("force_gene_refresh")
        if force_chrom_refresh:
            reasons.append("force_chrom_refresh")
        if not state.get("_rendered_once"):
            reasons.append("first_render")
        logging.debug(f"[CNV] Update needed. Reasons: {', '.join(reasons)}")

        data_array_reload = data_array_npy_changed or is_fresh_visit
        need_load = (
            cnv_dict_npy_changed
            or cnv_npy_changed
            or cnv2_npy_changed
            or cnv3_npy_changed
            or data_array_reload
            or xy_pkl_changed
        )

        return {
            "state": state,
            "key": key,
            "ui_changed": ui_changed,
            "is_fresh_visit": is_fresh_visit,
            "cnv_npy": cnv_npy,
            "cnv2_npy": cnv2_npy,
            "cnv3_npy": cnv3_npy,
            "cnv_dict_npy": cnv_dict_npy,
            "data_array_npy": data_array_npy,
            "xy_pkl": xy_pkl,
            "cnv_npy_mtime": cnv_npy_mtime,
            "cnv2_npy_mtime": cnv2_npy_mtime,
            "cnv3_npy_mtime": cnv3_npy_mtime,
            "cnv_dict_npy_mtime": cnv_dict_npy_mtime,
            "data_array_npy_mtime": data_array_npy_mtime,
            "xy_pkl_mtime": xy_pkl_mtime,
            "cnv_dict_npy_changed": cnv_dict_npy_changed,
            "cnv_npy_changed": cnv_npy_changed,
            "cnv2_npy_changed": cnv2_npy_changed,
            "cnv3_npy_changed": cnv3_npy_changed,
            "data_array_npy_changed": data_array_npy_changed,
            "xy_pkl_changed": xy_pkl_changed,
            "data_array_reload": data_array_reload,
            "need_load": need_load,
            "ichor_params_mtime": ichor_params_mtime,
            "ichor_seg_mtime": ichor_seg_mtime,
            "ichor_outputs_changed": ichor_outputs_changed,
        }

    def _apply_cnv_refresh_after_load(
        plan: Dict[str, Any], payload: Dict[str, Any]
    ) -> None:
        """Merge disk payload into state and update charts (main thread only)."""
        p = plan
        state = p["state"]
        key = p["key"]
        ui_changed = p["ui_changed"]
        is_fresh_visit = p["is_fresh_visit"]
        cnv_dict_npy = p["cnv_dict_npy"]
        cnv_npy = p["cnv_npy"]
        cnv2_npy = p["cnv2_npy"]
        cnv3_npy = p["cnv3_npy"]
        data_array_npy = p["data_array_npy"]
        xy_pkl = p["xy_pkl"]
        cnv_npy_mtime = p["cnv_npy_mtime"]
        cnv2_npy_mtime = p["cnv2_npy_mtime"]
        cnv3_npy_mtime = p["cnv3_npy_mtime"]
        cnv_dict_npy_mtime = p["cnv_dict_npy_mtime"]
        data_array_npy_mtime = p["data_array_npy_mtime"]
        xy_pkl_mtime = p["xy_pkl_mtime"]
        cnv_dict_npy_changed = p["cnv_dict_npy_changed"]
        cnv_npy_changed = p["cnv_npy_changed"]
        cnv2_npy_changed = p["cnv2_npy_changed"]
        cnv3_npy_changed = p["cnv3_npy_changed"]
        data_array_npy_changed = p["data_array_npy_changed"]
        data_array_reload = p["data_array_reload"]
        xy_pkl_changed = p["xy_pkl_changed"]
        ichor_outputs_changed = p.get("ichor_outputs_changed", False)
        ichor_params_mtime = p.get("ichor_params_mtime", 0)
        ichor_seg_mtime = p.get("ichor_seg_mtime", 0)

        changed = ("cnv" in payload or "cnv2" in payload or "cnv3" in payload)

        if cnv_dict_npy.exists():
            m = cnv_dict_npy_mtime
            if cnv_dict_npy_changed and "cnv_dict" in payload:
                state["cnv_dict"] = payload["cnv_dict"]
                cnv_bin.set_text(
                    f"Bin width: {state['cnv_dict'].get('bin_width', '--'):,}"
                )
                cnv_var.set_text(
                    f"Variance: {state['cnv_dict'].get('variance','--'):.3f}"
                    if isinstance(state["cnv_dict"].get("variance"), (int, float))
                    else "Variance: --"
                )
                state["dict_m"] = m
                if state.get("cnv"):
                    _update_ploidy_label(state)
                    _update_integer_clonal_label(state)
        if xy_pkl.exists():
            m = xy_pkl_mtime
            if xy_pkl_changed:
                try:
                    if "xy" in payload:
                        xy = payload["xy"]
                        cnv_xy.set_text(f"Genetic sex: {xy}")
                        state["xy"] = xy
                except Exception:
                    pass
                state["xy_m"] = m

        if "cnv" in payload:
            state["cnv"] = payload["cnv"]
            state["cnv_m"] = cnv_npy_mtime
        if "cnv2" in payload:
            state["cnv2"] = payload["cnv2"]
            state["cnv2_m"] = cnv2_npy_mtime
        if "cnv3" in payload:
            state["cnv3"] = payload["cnv3"]
            state["cnv3_m"] = cnv3_npy_mtime

        if state.get("cnv"):
            if changed:
                cnv_status.set_text("Status: CNV data loaded")
            if changed or not state.get("_chrom_opts_set"):
                try:
                    cnv_map = state["cnv"].get("cnv", state["cnv"])
                    chrom_opts = {"All": "All"}
                    for contig in natsort.natsorted(cnv_map.keys()):
                        if _cnv_contig_ok(contig):
                            chrom_opts[contig] = contig
                    cnv_chrom_select.set_options(chrom_opts)
                    chrom_opts = {"All": "All"}
                    for contig in natsort.natsorted(cnv_map.keys()):
                        if _cnv_contig_ok(contig):
                            chrom_opts[contig] = contig
                    cnv_chrom_select.set_options(chrom_opts)
                    sel = launcher._cnv_state.setdefault(str(sample_dir), {}).get(
                        "selected_chrom", "All"
                    )
                    if sel not in chrom_opts:
                        sel = "All"
                    cnv_chrom_select.value = sel
                    try:
                        cnv_chrom_select.update()
                    except Exception:
                        pass
                    state["_chrom_opts_set"] = True
                except Exception:
                    pass
            force_color_refresh = state.get("_force_color_refresh", False)
            force_gene_refresh = state.get("_force_gene_refresh", False)
            force_chrom_refresh = state.get("_force_chrom_refresh", False)
            if (
                changed
                or ui_changed
                or not state.get("_rendered_once")
                or force_color_refresh
                or force_gene_refresh
                or force_chrom_refresh
                or ichor_outputs_changed
            ):
                _render_cnv_from_state(state)
                state["_rendered_once"] = True
                if force_color_refresh:
                    state["_force_color_refresh"] = False
                if force_gene_refresh:
                    state["_force_gene_refresh"] = False
                if force_chrom_refresh:
                    state["_force_chrom_refresh"] = False

            _update_cnv_events_analysis(state)
            _update_ploidy_label(state)
            _update_integer_clonal_label(state)
            _update_ichor_py_label()
        else:
            try:
                cnv_ploidy.set_text("Ploidy (CNV bins): --")
                cnv_integer_clonal.set_text("Integer clonality (heuristic): --")
                cnv_ichor_py.set_text("Tumour fraction/ploidy (ichor_py): --")
            except Exception:
                pass

        if data_array_npy.exists() and (data_array_npy_changed or is_fresh_visit):
            try:
                arr = payload.get("bp_array")
                if arr is None:
                    arr = np.load(data_array_npy, allow_pickle=True)
                if hasattr(arr, "dtype") and "name" in arr.dtype.names:
                    state["bp_array"] = arr
                    selected = launcher._cnv_state.setdefault(
                        str(sample_dir), {}
                    ).get("selected_chrom", "All")
                    breakpoint_lines = []
                    for r in arr:
                        if selected == "All" or r["name"] == selected:
                            start_pos = int(r["start"])
                            end_pos = int(r["end"])
                            midpoint = (start_pos + end_pos) // 2
                            breakpoint_lines.append(midpoint)
                    current_series = [
                        s
                        for s in cnv_diff.options["series"]
                        if not s.get("name", "").startswith("Breakpoint")
                    ]
                    if selected != "All" and state.get("show_bp", True):
                        if current_series:
                            main_series = current_series[0]
                            markLine_data = []
                            for bp_pos in breakpoint_lines:
                                markLine_data.append(
                                    {
                                        "xAxis": bp_pos,
                                        "lineStyle": {
                                            "type": "dashed",
                                            "color": "#ff6b6b",
                                            "width": 3,
                                        },
                                    }
                                )
                            main_series["markLine"] = {
                                "data": markLine_data,
                                "symbol": "none",
                                "lineStyle": {"type": "dashed", "color": "#ff6b6b", "width": 3},
                            }
                    else:
                        if current_series:
                            current_series[0].pop("markLine", None)
                    cnv_diff.options["series"] = current_series
                    _apply_cnv_echart_chrome(cnv_diff, _is_dark_mode())
                    cnv_diff.update()
                state["bp_array_mtime"] = data_array_npy_mtime
            except Exception:
                pass
        elif data_array_npy.exists() and not data_array_npy_changed:
            if state.get("bp_array") is not None:
                try:
                    selected = launcher._cnv_state.setdefault(
                        str(sample_dir), {}
                    ).get("selected_chrom", "All")
                    current_series = [
                        s
                        for s in cnv_diff.options["series"]
                        if not s.get("name", "").startswith("Breakpoint")
                    ]
                    if selected != "All" and state.get("show_bp", True):
                        arr = state["bp_array"]
                        breakpoint_lines = []
                        for r in arr:
                            if selected == "All" or r["name"] == selected:
                                start_pos = int(r["start"])
                                end_pos = int(r["end"])
                                midpoint = (start_pos + end_pos) // 2
                                breakpoint_lines.append(midpoint)
                        if current_series:
                            main_series = current_series[0]
                            markLine_data = []
                            for bp_pos in breakpoint_lines:
                                markLine_data.append(
                                    {
                                        "xAxis": bp_pos,
                                        "lineStyle": {
                                            "type": "dashed",
                                            "color": "#ff6b6b",
                                            "width": 3,
                                        },
                                    }
                                )
                            main_series["markLine"] = {
                                "data": markLine_data,
                                "symbol": "none",
                                "lineStyle": {"type": "dashed", "color": "#ff6b6b", "width": 3},
                            }
                    else:
                        if current_series:
                            current_series[0].pop("markLine", None)
                    cnv_diff.options["series"] = current_series
                    _apply_cnv_echart_chrome(cnv_diff, _is_dark_mode())
                    cnv_diff.update()
                except Exception:
                    pass

        state["cnv_m"] = cnv_npy_mtime
        state["cnv2_m"] = cnv2_npy_mtime
        state["cnv3_m"] = cnv3_npy_mtime
        state["dict_m"] = cnv_dict_npy_mtime
        state["xy_m"] = xy_pkl_mtime
        state["ichor_params_m"] = ichor_params_mtime
        state["ichor_seg_m"] = ichor_seg_mtime
        state["last_visit_time"] = state.get("last_visit_time", time.time())
        state["cnv_plot_theme_dark"] = _is_dark_mode()

        launcher._cnv_state[key] = state
        _update_breakpoints_visibility()

    async def _refresh_cnv_async() -> None:
        try:
            plan = _prepare_cnv_refresh()
            if plan is None:
                return
            if plan["need_load"]:
                payload = await asyncio.to_thread(
                    _cnv_load_binary_payload,
                    sample_dir,
                    cnv_dict_npy_changed=plan["cnv_dict_npy_changed"],
                    cnv_npy_changed=plan["cnv_npy_changed"],
                    cnv2_npy_changed=plan["cnv2_npy_changed"],
                    cnv3_npy_changed=plan["cnv3_npy_changed"],
                    data_array_reload=plan["data_array_reload"],
                    xy_pkl_changed=plan["xy_pkl_changed"],
                )
            else:
                payload = {}
            _apply_cnv_refresh_after_load(plan, payload)
        except Exception:
            pass

    def _refresh_cnv_sync(sample_dir: Path, launcher: Any) -> None:
        """Synchronous CNV refresh (no event loop)."""
        try:
            plan = _prepare_cnv_refresh()
            if plan is None:
                return
            if plan["need_load"]:
                payload = _cnv_load_binary_payload(
                    sample_dir,
                    cnv_dict_npy_changed=plan["cnv_dict_npy_changed"],
                    cnv_npy_changed=plan["cnv_npy_changed"],
                    cnv2_npy_changed=plan["cnv2_npy_changed"],
                    cnv3_npy_changed=plan["cnv3_npy_changed"],
                    data_array_reload=plan["data_array_reload"],
                    xy_pkl_changed=plan["xy_pkl_changed"],
                )
            else:
                payload = {}
            _apply_cnv_refresh_after_load(plan, payload)
        except Exception:
            pass

    def _refresh_cnv() -> None:
        """Refresh CNV data; offload binary loads when a loop is running."""
        try:
            if not sample_dir or not sample_dir.exists():
                logging.warning(f"[CNV] Sample directory not found: {sample_dir}")
                return
            try:
                asyncio.get_running_loop()
            except RuntimeError:
                _refresh_cnv_sync(sample_dir, launcher)
                return
            asyncio.create_task(_refresh_cnv_async())
        except Exception as e:
            logging.exception(f"[CNV] Refresh failed: {e}")

    def _update_breakpoints_visibility() -> None:
        """Show/hide breakpoints controls based on chromosome selection."""
        try:
            key = str(sample_dir)
            state = launcher._cnv_state.get(key, {})
            selected = state.get("selected_chrom", "All")
            
            # Show breakpoints controls only when viewing individual chromosomes
            should_show = selected != "All"
            
            try:
                display_value = "block" if should_show else "none"
                cnv_bp_label.style(f"display: {display_value}")
                cnv_bp.style(f"display: {display_value}")
            except Exception:
                pass
        except Exception:
            pass

    # Bind control events
    try:

        def _val(ev, default=None):
            # Handle toggle events which have args like [index, {'value': X, 'label': 'Y'}]
            if hasattr(ev, "args") and ev.args and isinstance(ev.args, list) and len(ev.args) >= 2:
                if isinstance(ev.args[1], dict):
                    # Extract the label from the toggle event structure
                    return ev.args[1].get("label", default)
            # Handle direct value objects like {'value': 2, 'label': 'GNB1'}
            if hasattr(ev, "value") and isinstance(ev.value, dict):
                return ev.value.get("label", default)
            # Handle args that are directly a dictionary with label
            if hasattr(ev, "args") and isinstance(ev.args, dict) and "label" in ev.args:
                return ev.args.get("label", default)
            # Fallback to standard value extraction
            return (
                getattr(ev, "value", None)
                if hasattr(ev, "value")
                else (getattr(ev, "args", None) or default)
            )

        def _on_chrom(ev):
            st = launcher._cnv_state.setdefault(str(sample_dir), {})
            st["selected_chrom"] = _val(ev, "All") or "All"
            st["_force_chrom_refresh"] = True  # Force refresh for chromosome selection
            logging.debug(f"CNV select changed -> {st['selected_chrom']}")
            # Update breakpoints visibility
            _update_breakpoints_visibility()
            # reset x zoom when switching scope
            try:
                for chart in (cnv_abs, cnv_diff):
                    if (
                        isinstance(chart.options.get("dataZoom"), list)
                        and chart.options["dataZoom"]
                    ):
                        chart.options["dataZoom"][0].pop("startValue", None)
                        chart.options["dataZoom"][0].pop("endValue", None)
            except Exception:
                pass
            # Trigger immediate refresh to update all UI elements
            ui.timer(0.1, _refresh_cnv, once=True)

        def _on_scale(ev):
            st = launcher._cnv_state.setdefault(str(sample_dir), {})
            st["y_scale"] = _val(ev, "linear") or "linear"
            # Trigger immediate refresh to update all UI elements
            ui.timer(0.1, _refresh_cnv, once=True)

        def _on_bp(ev):
            st = launcher._cnv_state.setdefault(str(sample_dir), {})
            raw = _val(ev, "show")
            v = str(raw).strip().lower()
            show_bp = v in ("show", "true", "1")
            st["show_bp"] = show_bp
            st["_force_chrom_refresh"] = True
            st["_last_refresh"] = 0  # bypass debounce for immediate repaint
            # Trigger immediate refresh to update all UI elements
            ui.timer(0.12, _refresh_cnv, once=True)

        def _on_segments(ev):
            st = launcher._cnv_state.setdefault(str(sample_dir), {})
            sval = str(_val(ev, "show")).strip().lower()
            st["show_segments"] = sval in ("show", "true", "1")
            st["_force_color_refresh"] = True
            ui.timer(0.1, _refresh_cnv, once=True)

        def _on_segment_source(ev):
            st = launcher._cnv_state.setdefault(str(sample_dir), {})
            sval = str(_val(ev, "heuristic")).strip().lower()
            if sval not in ("heuristic", "ichor_py"):
                sval = "heuristic"
            st["segment_source"] = sval
            st["_force_color_refresh"] = True
            st["_last_refresh"] = 0
            ui.timer(0.1, _refresh_cnv, once=True)

        def _on_color(ev):
            st = launcher._cnv_state.setdefault(str(sample_dir), {})
            val = _val(ev, "chromosome") or "chromosome"
            # Accept either keys or labels from the toggle
            vlow = str(val).strip().lower()
            if vlow in ("chromosome", "chromosomes"):
                st["color_mode"] = "chromosome"
            elif vlow in ("value", "up/down", "updown", "up_down", "up down"):
                st["color_mode"] = "value"
            else:
                st["color_mode"] = "chromosome"
            # Force a refresh by setting a flag that bypasses the state sync logic
            st["_force_color_refresh"] = True
            
            # Trigger immediate refresh to update all UI elements
            ui.timer(0.1, _refresh_cnv, once=True)

        def _on_plot_bin(ev):
            st = launcher._cnv_state.setdefault(str(sample_dir), {})
            v = getattr(ev, "args", None) if hasattr(ev, "args") else getattr(ev, "value", None)
            if v is None and hasattr(ev, "value"):
                v = ev.value
            # NiceGUI may pass a dict {'value': index, 'label': '...'} instead of the key
            if isinstance(v, dict):
                keys_ordered = list(_PLOT_BIN_KEY_TO_BP.keys())
                idx = v.get("value", 0)
                key = keys_ordered[idx] if 0 <= idx < len(keys_ordered) else _PLOT_BIN_KEY_DEFAULT
            else:
                key = v
            st["plot_bin_width"] = _PLOT_BIN_KEY_TO_BP.get(
                key, _PLOT_BIN_KEY_TO_BP[_PLOT_BIN_KEY_DEFAULT]
            )
            st["_force_chrom_refresh"] = True  # force re-render with new bin width
            ui.timer(0.1, _refresh_cnv, once=True)

        def _on_assumed_ploidy_display() -> None:
            st = launcher._cnv_state.setdefault(str(sample_dir), {})
            try:
                v = float(cnv_assumed_ploidy.value)
                st["plot_assumed_ploidy"] = v if v > 0 else 2.0
            except (TypeError, ValueError):
                st["plot_assumed_ploidy"] = 2.0
            st["_force_color_refresh"] = True
            ui.timer(0.1, _refresh_cnv, once=True)

        def _on_run_ichor_py(_ev=None) -> None:
            """Run ichor_py on CNV3 from GUI state or from CNV3.npy on disk (no full reanalysis)."""
            st = launcher._cnv_state.setdefault(str(sample_dir), {})
            cnv3_map, bw = _cnv3_map_and_binwidth_for_ichor(st, sample_dir)
            if not cnv3_map or bw <= 0:
                try:
                    ui.notify(
                        "Need CNV3 data and bin width (CNV3.npy and CNV_dict.npy in the sample folder, "
                        "or load the sample in this panel first).",
                        type="negative",
                    )
                except Exception:
                    pass
                return
            cnv_ichor_py.set_text("Tumour fraction/ploidy (ichor_py): computing…")
            log = logging.getLogger("robin.gui.cnv")

            def _after_run(ok: bool) -> None:
                try:
                    cnv_run_ichor_py.enable()
                except Exception:
                    pass
                st["_last_refresh"] = 0
                if ok:
                    st["ichor_params_m"] = 0
                    st["ichor_seg_m"] = 0
                _refresh_cnv()
                if ok:
                    try:
                        ui.notify("ichor_py results saved to sample folder.", type="positive")
                    except Exception:
                        pass
                else:
                    try:
                        ui.notify("ichor_py did not complete (see logs).", type="warning")
                    except Exception:
                        pass

            try:
                cnv_run_ichor_py.disable()
            except Exception:
                pass

            def _run_blocking() -> bool:
                return persist_ichor_py_outputs(str(sample_dir), cnv3_map, bw, log)

            try:
                asyncio.get_running_loop()
            except RuntimeError:
                _after_run(_run_blocking())
                return

            async def _async_run() -> None:
                ok = await asyncio.to_thread(_run_blocking)
                _after_run(ok)

            asyncio.create_task(_async_run())

        # Bind both native change and model-value updates for robustness
        cnv_chrom_select.on("change", _on_chrom)
        cnv_chrom_select.on("update:model-value", _on_chrom)

        # Gene selection zoom
        def _on_gene(ev):
            st = launcher._cnv_state.setdefault(str(sample_dir), {})
            selected_gene = _val(ev, "All") or "All"
            st["selected_gene"] = selected_gene
            st["_force_gene_refresh"] = True  # Force refresh for gene selection
            # Trigger immediate refresh to update all UI elements
            ui.timer(0.1, _refresh_cnv, once=True)

        cnv_gene_select.on("change", _on_gene)
        cnv_gene_select.on("update:model-value", _on_gene)
        cnv_scale.on("change", _on_scale)
        cnv_scale.on("update:model-value", _on_scale)
        cnv_plot_bin.on("change", _on_plot_bin)
        cnv_plot_bin.on("update:model-value", _on_plot_bin)
        cnv_segments.on("change", _on_segments)
        cnv_segments.on("update:model-value", _on_segments)
        cnv_segment_source.on("change", _on_segment_source)
        cnv_segment_source.on("update:model-value", _on_segment_source)
        cnv_bp.on("change", _on_bp)
        cnv_bp.on("update:model-value", _on_bp)
        cnv_color.on("change", _on_color)
        cnv_color.on("update:model-value", _on_color)
        cnv_assumed_ploidy.on("change", _on_assumed_ploidy_display)
        cnv_assumed_ploidy.on("update:model-value", _on_assumed_ploidy_display)
        cnv_run_ichor_py.on("click", _on_run_ichor_py)
    except Exception:
        pass

    def _sync_cnv_echarts_theme_if_needed() -> None:
        """Re-apply axis/tooltip/title colours when the user toggles light/dark mode."""
        try:
            dark = _is_dark_mode()
        except Exception:
            dark = False
        key = str(sample_dir)
        st = launcher._cnv_state.get(key, {})
        if st.get("cnv_plot_theme_dark") == dark:
            return
        _apply_cnv_echart_chrome(cnv_abs, dark)
        _apply_cnv_echart_chrome(cnv_diff, dark)
        try:
            cnv_abs.update()
            cnv_diff.update()
        except Exception:
            pass
        launcher._cnv_state.setdefault(key, {})["cnv_plot_theme_dark"] = dark

    # Start the refresh timer (every 30 seconds)
    refresh_timer = ui.timer(30.0, _refresh_cnv, active=True, immediate=False)
    ui.timer(0.5, _refresh_cnv, once=True)
    unregister_cnv_theme_sync = register_theme_sync_callback(
        _sync_cnv_echarts_theme_if_needed,
        element=cnv_abs,
        interval_s=0.5,
        immediate=True,
    )
    try:
        ui.context.client.on_disconnect(
            lambda: (refresh_timer.deactivate(), unregister_cnv_theme_sync())
        )
    except Exception:
        pass
