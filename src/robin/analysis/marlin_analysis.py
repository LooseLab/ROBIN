"""
MARLIN methylation classification for ROBIN.

MARLIN (Methylation- and AI-guided Rapid Leukemia Subtype Inference) classifies
acute leukemia methylation subtypes from sparse CpG beta profiles. This module
mirrors the Random Forest / Sturgeon workflow job pattern:

  bed_conversion parquet → probe-mapped betas → Keras model → marlin_scores.csv

TensorFlow is an optional dependency (``pip install robin[marlin]``). The trained
model is downloaded from Zenodo on first use via ``robin.utils.marlin_manager``.
"""

from __future__ import annotations

import gzip
import logging
import os
import time
from dataclasses import dataclass, field
from typing import Any, Dict, List, Mapping, Optional

import pandas as pd

from robin.logging_config import get_job_logger
from robin.utils.marlin_manager import (
    DEFAULT_GENOME_BUILD,
    ensure_marlin_assets,
)

try:
    from robin.analysis.utilities.merge_bedmethyl import load_minimal_modkit_data
except ImportError:  # pragma: no cover - optional in minimal envs
    load_minimal_modkit_data = None


def _is_fail_only_expected(job) -> bool:
    """True when errors are expected for fail-only BAM submissions."""
    try:
        if not bool(job.context.metadata.get("fail_only_bam_submission", False)):
            return False
        fp = getattr(job.context, "filepath", "") or ""
        base = os.path.basename(fp).lower()
        return ("fail" in base) and ("pass" not in base)
    except Exception:
        return False


@dataclass
class MarlinMetadata:
    """Container for MARLIN analysis metadata and results."""

    sample_id: str
    parquet_path: str
    analysis_timestamp: float
    batch_number: int
    scores_file_path: Optional[str] = None
    bed_file_path: Optional[str] = None
    covered_cpgs: int = 0
    top_class: Optional[str] = None
    top_score: Optional[float] = None
    processing_steps: List[str] = field(default_factory=list)
    error_message: Optional[str] = None
    results: Dict[str, Any] = field(default_factory=dict)


def _normalize_chrom(value: object) -> str:
    text = str(value).strip()
    if text.lower().startswith("chr"):
        return text
    if text.upper() in {"M", "MT"}:
        return "chrM"
    return f"chr{text}"


def load_probe_position_map(
    probes_bed_path: str | os.PathLike,
) -> Dict[tuple[str, int], str]:
    """
    Build (chrom, ref_position) → probe_id from a MARLIN probe BED.

    Matches MARLIN's R realtime pipeline, which joins modkit pileup on chrom + V2
    (0-based start) against the probe BED.
    """
    mapping: Dict[tuple[str, int], str] = {}
    open_fn = gzip.open if str(probes_bed_path).endswith(".gz") else open
    with open_fn(probes_bed_path, "rt", encoding="utf-8") as handle:
        for line in handle:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split("\t")
            if len(parts) < 4:
                continue
            chrom = _normalize_chrom(parts[0])
            try:
                start = int(parts[1])
            except ValueError:
                continue
            probe_id = parts[3].strip()
            if probe_id:
                mapping[(chrom, start)] = probe_id
    if not mapping:
        raise ValueError(f"No probe positions loaded from {probes_bed_path}")
    return mapping


def bedmethyl_to_probe_values(
    modkit_df: pd.DataFrame,
    position_map: Mapping[tuple[str, int], str],
) -> Dict[str, float]:
    """
    Map modkit/bedmethyl rows to MARLIN probe beta values.

    ``percent_modified`` may be 0–100 or 0–1; values > 1 are treated as percent.
    Multiple strand/position hits for the same probe are averaged.
    """
    if modkit_df is None or modkit_df.empty:
        return {}

    df = modkit_df.copy()
    chrom_col = "chrom" if "chrom" in df.columns else "chr"
    start_col = (
        "chromStart"
        if "chromStart" in df.columns
        else "start_pos" if "start_pos" in df.columns else "start"
    )
    value_col = (
        "percent_modified"
        if "percent_modified" in df.columns
        else (
            "fraction"
            if "fraction" in df.columns
            else "score" if "score" in df.columns else None
        )
    )
    if value_col is None:
        raise ValueError(
            "Methylation dataframe missing percent_modified/fraction/score column"
        )

    sums: Dict[str, float] = {}
    counts: Dict[str, int] = {}
    for chrom, start, raw in zip(
        df[chrom_col].tolist(), df[start_col].tolist(), df[value_col].tolist()
    ):
        try:
            start_i = int(start)
            beta = float(raw)
        except (TypeError, ValueError):
            continue
        if beta > 1.0:
            beta = beta / 100.0
        probe_id = position_map.get((_normalize_chrom(chrom), start_i))
        if probe_id is None:
            continue
        sums[probe_id] = sums.get(probe_id, 0.0) + beta
        counts[probe_id] = counts.get(probe_id, 0) + 1

    return {probe: sums[probe] / counts[probe] for probe in sums}


def write_marlin_probe_bed(
    probe_values: Mapping[str, float],
    position_map: Mapping[tuple[str, int], str],
    output_path: str,
) -> str:
    """Write chrom/start/end/beta/probe_id BED for debugging / MARLIN predict_bed."""
    # Invert map to one representative interval per probe (min start).
    probe_coords: Dict[str, tuple[str, int]] = {}
    for (chrom, start), probe_id in position_map.items():
        existing = probe_coords.get(probe_id)
        if existing is None or start < existing[1]:
            probe_coords[probe_id] = (chrom, start)

    rows = []
    for probe_id, beta in probe_values.items():
        coords = probe_coords.get(probe_id)
        if coords is None:
            continue
        chrom, start = coords
        rows.append((chrom, start, start + 1, beta, probe_id))
    rows.sort(key=lambda r: (r[0], r[1]))

    os.makedirs(os.path.dirname(output_path) or ".", exist_ok=True)
    with open(output_path, "w", encoding="utf-8") as handle:
        for chrom, start, end, beta, probe_id in rows:
            handle.write(f"{chrom}\t{start}\t{end}\t{beta}\t{probe_id}\n")
    return output_path


def append_marlin_scores(
    scores_path: str,
    prediction,
    *,
    timestamp_ms: float,
) -> None:
    """Append one prediction row to ``marlin_scores.csv`` (fractional softmax scores)."""
    row: Dict[str, Any] = {
        "timestamp": timestamp_ms,
        "covered_cpgs": prediction.covered_cpgs,
        "number_probes": prediction.covered_cpgs,
    }
    row.update({name: float(score) for name, score in prediction.scores.items()})
    new_df = pd.DataFrame([row])

    if os.path.exists(scores_path):
        try:
            existing = pd.read_csv(scores_path)
            combined = pd.concat([existing, new_df], ignore_index=True)
        except Exception:
            combined = new_df
    else:
        combined = new_df

    # Keep timestamp first, then covered_cpgs, then class columns.
    cols = list(combined.columns)
    ordered = [c for c in ("timestamp", "covered_cpgs", "number_probes") if c in cols]
    ordered.extend(sorted(c for c in cols if c not in ordered))
    combined[ordered].to_csv(scores_path, index=False)


class MarlinAnalysis:
    """MARLIN analysis worker."""

    def __init__(
        self,
        work_dir: Optional[str] = None,
        *,
        genome_build: str = DEFAULT_GENOME_BUILD,
        model_path: Optional[str] = None,
        download_if_missing: bool = True,
    ):
        self.work_dir = work_dir or os.getcwd()
        self.genome_build = genome_build or DEFAULT_GENOME_BUILD
        self.model_path = model_path
        self.download_if_missing = download_if_missing
        self.bambatch: Dict[str, int] = {}
        self._predictor = None
        self._position_map: Optional[Dict[tuple[str, int], str]] = None
        self._assets: Optional[Dict[str, Any]] = None

        logger = logging.getLogger("robin.marlin")
        logger.info("Marlin Analysis initialized (genome=%s)", self.genome_build)

    def _ensure_predictor(self):
        if self._predictor is not None:
            return self._predictor

        from robin.analysis.marlin_python import MARLINPredictor

        assets = ensure_marlin_assets(
            genome_build=self.genome_build,
            model_path=self.model_path,
            download_if_missing=self.download_if_missing,
        )
        self._assets = assets
        self._position_map = load_probe_position_map(assets["probes"])
        self._predictor = MARLINPredictor.from_paths(
            model_path=assets["model"],
            feature_path=assets["features"],
            annotation_path=assets["annotations"],
        )
        return self._predictor

    def process_parquet_file(self, parquet_path: str, sample_id: str) -> MarlinMetadata:
        logger = logging.getLogger("robin.marlin")
        start_time = time.time()

        if sample_id not in self.bambatch:
            self.bambatch[sample_id] = 1

        result = MarlinMetadata(
            sample_id=sample_id,
            parquet_path=parquet_path,
            analysis_timestamp=start_time,
            batch_number=self.bambatch[sample_id],
        )

        try:
            if load_minimal_modkit_data is None:
                raise ImportError("merge_bedmethyl helpers are unavailable")
            if not os.path.exists(parquet_path):
                raise FileNotFoundError(f"Parquet file not found: {parquet_path}")
            result.processing_steps.append("file_validation")

            predictor = self._ensure_predictor()
            assert self._position_map is not None
            result.processing_steps.append("model_loaded")

            modkit_df = load_minimal_modkit_data(parquet_path)
            result.processing_steps.append("modkit_data_loaded")

            probe_values = bedmethyl_to_probe_values(modkit_df, self._position_map)
            if not probe_values:
                raise ValueError(
                    "No MARLIN probes overlapped methylation calls in parquet "
                    f"(genome_build={self.genome_build})"
                )
            result.processing_steps.append("probes_mapped")

            sample_dir = os.path.join(self.work_dir, sample_id)
            os.makedirs(sample_dir, exist_ok=True)
            bed_path = os.path.join(sample_dir, "MarlinBed.bed")
            write_marlin_probe_bed(probe_values, self._position_map, bed_path)
            result.bed_file_path = bed_path
            result.processing_steps.append("bed_file_created")

            prediction = predictor.predict_from_probe_values(probe_values)
            result.covered_cpgs = prediction.covered_cpgs
            result.top_class = prediction.top_class
            result.top_score = prediction.top_score
            result.processing_steps.append("prediction_complete")

            scores_path = os.path.join(sample_dir, "marlin_scores.csv")
            append_marlin_scores(
                scores_path, prediction, timestamp_ms=start_time * 1000
            )
            result.scores_file_path = scores_path
            result.processing_steps.append("results_saved")

            result.results = {
                "status": "success",
                "sample_id": sample_id,
                "analysis_time": time.time() - start_time,
                "batch_number": result.batch_number,
                "covered_cpgs": prediction.covered_cpgs,
                "top_class": prediction.top_class,
                "top_score": prediction.top_score,
                "scores_file": scores_path,
                "bed_file": bed_path,
                "genome_build": self.genome_build,
                "processing_steps": result.processing_steps.copy(),
            }
            self.bambatch[sample_id] += 1
            logger.info(
                "MARLIN %s → %s (%.3f) covered_cpgs=%s",
                sample_id,
                prediction.top_class,
                prediction.top_score,
                prediction.covered_cpgs,
            )
            return result

        except Exception as exc:
            logger.error(
                "MARLIN analysis failed for %s: %s", sample_id, exc, exc_info=True
            )
            result.error_message = str(exc)
            result.processing_steps.append("analysis_failed")
            return result


def process_multiple_files(
    parquet_paths,
    metadata_list,
    work_dir,
    logger,
    *,
    genome_build: str = DEFAULT_GENOME_BUILD,
    model_path: Optional[str] = None,
    job_id=None,
):
    """Process multiple parquet files for one sample into one marlin_scores.csv."""
    if not parquet_paths or not metadata_list:
        raise ValueError("parquet_paths and metadata_list must not be empty")
    if len(parquet_paths) != len(metadata_list):
        raise ValueError("parquet_paths and metadata_list must have the same length")

    sample_id = metadata_list[0].get("sample_id", "unknown")
    analysis_result = {
        "sample_id": sample_id,
        "parquet_paths": parquet_paths,
        "analysis_timestamp": time.time(),
        "processing_steps": [],
        "error_message": None,
        "files_processed": 0,
        "total_files": len(parquet_paths),
        "total_batches": 0,
        "scores_file_path": None,
        "bed_file_path": None,
    }

    try:
        sample_dir = os.path.join(work_dir, sample_id)
        os.makedirs(sample_dir, exist_ok=True)
        analyzer = MarlinAnalysis(
            work_dir=work_dir,
            genome_build=genome_build,
            model_path=model_path,
        )
        analysis_result["processing_steps"].append("analyzer_created")

        for i, parquet_path in enumerate(parquet_paths):
            logger.info(
                "MARLIN file %s/%s: %s",
                i + 1,
                len(parquet_paths),
                os.path.basename(parquet_path),
            )
            if not os.path.exists(parquet_path):
                logger.warning("Parquet not found: %s", parquet_path)
                continue
            marlin_result = analyzer.process_parquet_file(parquet_path, sample_id)
            if marlin_result.error_message:
                logger.warning(
                    "MARLIN skipped %s: %s",
                    os.path.basename(parquet_path),
                    marlin_result.error_message,
                )
                continue
            analysis_result["files_processed"] += 1
            analysis_result["total_batches"] += 1
            analysis_result["scores_file_path"] = marlin_result.scores_file_path
            analysis_result["bed_file_path"] = marlin_result.bed_file_path
            try:
                from robin.workflow_ray import notify_coordinator_files_completed

                notify_coordinator_files_completed("marlin", 1, job_id=job_id)
            except Exception:
                pass

        if analysis_result["files_processed"] == 0:
            analysis_result["error_message"] = (
                "No files could be processed successfully"
            )
            analysis_result["processing_steps"].append("no_files_processed")
            return analysis_result

        analysis_result["processing_steps"].append("analysis_complete")
        return analysis_result
    except Exception as exc:
        logger.error("Multi-file MARLIN analysis failed for %s: %s", sample_id, exc)
        analysis_result["error_message"] = str(exc)
        analysis_result["processing_steps"].append("analysis_failed")
        return analysis_result


def _genome_build_from_job(job) -> str:
    try:
        meta = job.context.metadata or {}
        for key in ("marlin_genome_build", "reference_genome", "genome_build"):
            value = meta.get(key)
            if value:
                return str(value)
        bam_meta = meta.get("bam_metadata") or {}
        for key in ("reference_genome", "genome_build"):
            value = bam_meta.get(key)
            if value:
                return str(value)
    except Exception:
        pass
    return DEFAULT_GENOME_BUILD


def marlin_handler(job, work_dir=None):
    """Workflow handler for MARLIN classification jobs."""
    logger = get_job_logger(str(job.job_id), job.job_type, job.context.filepath)
    suppress_expected = _is_fail_only_expected(job)
    genome_build = _genome_build_from_job(job)

    try:
        batched_job = job.context.metadata.get("_batched_job")
        if batched_job:
            batch_size = batched_job.get_file_count()
            sample_id = batched_job.get_sample_id()
            batch_id = batched_job.batch_id
            logger.info(
                "Processing MARLIN batch: %s files for sample '%s' (batch_id: %s)",
                batch_size,
                sample_id,
                batch_id,
            )

            metadata_list = []
            parquet_paths = []
            for i, bam_path in enumerate(batched_job.get_filepaths()):
                file_context = batched_job.contexts[i]
                file_metadata = dict(file_context.metadata.get("bam_metadata", {}))
                file_sample_id = file_context.get_sample_id()
                file_metadata["sample_id"] = (
                    file_sample_id if file_sample_id != "unknown" else sample_id
                )
                parquet_path = file_context.results.get("bed_conversion", {}).get(
                    "parquet_path"
                )
                if parquet_path:
                    parquet_paths.append(parquet_path)
                    metadata_list.append(file_metadata)
                else:
                    logger.warning(
                        "No parquet path for batch file %s: %s",
                        i + 1,
                        os.path.basename(bam_path),
                    )

            if not parquet_paths:
                error_msg = (
                    "No parquet paths found from bed conversion results in batch"
                )
                if suppress_expected:
                    logger.warning(
                        "%s (expected for fail-only BAM submission)", error_msg
                    )
                    job.context.add_result(
                        "marlin_analysis",
                        {"status": "expected_failure", "reason": error_msg},
                    )
                else:
                    logger.error(error_msg)
                    job.context.add_error("marlin_analysis", error_msg)
                return

            batch_work_dir = (
                os.path.dirname(parquet_paths[0]) if work_dir is None else work_dir
            )
            if work_dir is not None:
                os.makedirs(work_dir, exist_ok=True)

            batch_result = process_multiple_files(
                parquet_paths=parquet_paths,
                metadata_list=metadata_list,
                work_dir=batch_work_dir,
                logger=logger,
                genome_build=genome_build,
                job_id=job.job_id,
            )
            job.context.add_metadata(
                "marlin_analysis",
                {
                    "batch_result": batch_result,
                    "batch_size": batch_size,
                    "sample_id": sample_id,
                    "batch_id": batch_id,
                    "files_processed": batch_result.get("files_processed", 0),
                    "total_files": batch_result.get("total_files", 0),
                },
            )
            if batch_result.get("error_message"):
                if suppress_expected:
                    job.context.add_result(
                        "marlin_analysis",
                        {
                            "status": "expected_failure",
                            "error_message": batch_result["error_message"],
                            "sample_id": sample_id,
                        },
                    )
                else:
                    job.context.add_error(
                        "marlin_analysis", batch_result["error_message"]
                    )
            else:
                job.context.add_result(
                    "marlin_analysis",
                    {
                        "status": "success",
                        "sample_id": sample_id,
                        "analysis_time": batch_result.get("analysis_timestamp", 0),
                        "processing_steps": batch_result.get("processing_steps", []),
                        "scores_file": batch_result.get("scores_file_path", ""),
                        "bed_file": batch_result.get("bed_file_path", ""),
                        "files_processed": batch_result.get("files_processed", 0),
                        "total_files": batch_result.get("total_files", 0),
                        "genome_build": genome_build,
                    },
                )
            return

        bed_conversion_results = job.context.results.get("bed_conversion", {})
        parquet_path = bed_conversion_results.get("parquet_path")
        if not parquet_path:
            raise ValueError("No parquet file path found in bed_conversion results")

        bam_metadata = job.context.metadata.get("bam_metadata", {})
        sample_id = bam_metadata.get("sample_id", "unknown")
        if work_dir is None:
            work_dir = os.path.dirname(parquet_path)
        else:
            os.makedirs(work_dir, exist_ok=True)

        analyzer = MarlinAnalysis(work_dir=work_dir, genome_build=genome_build)
        result = analyzer.process_parquet_file(parquet_path, sample_id)
        job.context.add_metadata("marlin_analysis", result.results)
        job.context.add_metadata("marlin_processing_steps", result.processing_steps)

        if result.error_message:
            if suppress_expected:
                job.context.add_result(
                    "marlin_analysis",
                    {
                        "status": "expected_failure",
                        "error_message": result.error_message,
                    },
                )
            else:
                job.context.add_error("marlin_analysis", result.error_message)
                logger.error("MARLIN analysis failed: %s", result.error_message)
        else:
            job.context.add_result(
                "marlin_analysis",
                {
                    "status": "success",
                    "sample_id": result.sample_id,
                    "analysis_time": result.results.get("analysis_time", 0),
                    "batch_number": result.batch_number,
                    "processing_steps": result.processing_steps,
                    "scores_file": result.scores_file_path,
                    "bed_file": result.bed_file_path,
                    "covered_cpgs": result.covered_cpgs,
                    "top_class": result.top_class,
                    "top_score": result.top_score,
                    "genome_build": genome_build,
                },
            )
            try:
                from robin.workflow_ray import notify_coordinator_files_completed

                notify_coordinator_files_completed("marlin", 1, job_id=job.job_id)
            except Exception:
                pass
            logger.info(
                "MARLIN complete for %s: %s (%.3f)",
                sample_id,
                result.top_class,
                result.top_score or 0.0,
            )

    except Exception as exc:
        if suppress_expected:
            logger.warning(
                "Expected MARLIN failure for fail-only BAM submission: %s", exc
            )
            job.context.add_result(
                "marlin_analysis",
                {"status": "expected_failure", "error_message": str(exc)},
            )
            return
        job.context.add_error("marlin_analysis", str(exc))
        logger.error("Error in MARLIN analysis for %s: %s", job.context.filepath, exc)
