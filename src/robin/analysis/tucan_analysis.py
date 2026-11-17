"""
Tucan methylation classification for ROBIN.

Tucan (solid Tumor Classification using Nanopore) classifies pediatric solid
tumors and lymphomas from sparse CpG methylation calls. This module mirrors the
MARLIN / Sturgeon workflow job pattern:

  bed_conversion parquet → exact probe-mapped BED → tucan.predict → tucan_scores.csv

Probe mapping uses exact chrom/position match by default (margin=0). Sturgeon's
historic ±25 bp window is opt-in via ``ROBIN_TUCAN_PROBE_MARGIN``.

The ``tucan`` package and Hugging Face model are optional
(``pip install robin[tucan]``). See https://github.com/UMCUGenetics/tucan
"""

from __future__ import annotations

import logging
import os
import tempfile
import time
from dataclasses import dataclass, field
from typing import Any, Dict, List, Optional

import pandas as pd

from robin.logging_config import get_job_logger
from robin.utils.tucan_manager import (
    DEFAULT_NUM_CPGS,
    ensure_tucan_assets,
    get_default_num_cpgs,
    get_default_num_samplings,
    get_probe_margin,
    require_tucan_package,
)

try:
    from robin.analysis.utilities.merge_bedmethyl import (
        load_modkit_data,
        modkit_pileup_file_to_bed,
    )
except ImportError:  # pragma: no cover - optional in minimal envs
    load_modkit_data = None
    modkit_pileup_file_to_bed = None

logger = logging.getLogger("robin.tucan")

SCORE_META_COLUMNS = frozenset(
    {
        "timestamp",
        "number_probes",
        "covered_cpgs",
        "probes",
    }
)


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
class TucanMetadata:
    """Container for Tucan analysis metadata and results."""

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


def binarize_methylation_calls(bed_df: pd.DataFrame) -> pd.DataFrame:
    """Ensure ``methylation_call`` is 0/1 as expected by Tucan."""
    if bed_df is None or bed_df.empty:
        return bed_df
    out = bed_df.copy()
    if "methylation_call" not in out.columns:
        raise ValueError("Probe BED missing methylation_call column")
    out["methylation_call"] = (out["methylation_call"].astype(float) > 0).astype(int)
    return out


def append_tucan_scores(
    scores_path: str,
    prediction_df: pd.DataFrame,
    *,
    timestamp_ms: float,
) -> None:
    """Append one (or more) Tucan prediction rows to ``tucan_scores.csv``."""
    if prediction_df is None or prediction_df.empty:
        raise ValueError("Tucan prediction DataFrame is empty")

    rows = prediction_df.copy()
    if "probes" in rows.columns and "number_probes" not in rows.columns:
        rows = rows.rename(columns={"probes": "number_probes"})
    if "number_probes" in rows.columns and "covered_cpgs" not in rows.columns:
        rows["covered_cpgs"] = rows["number_probes"]
    rows["timestamp"] = timestamp_ms

    # Prefer timestamp / coverage meta first, then class columns.
    cols = list(rows.columns)
    ordered = [
        c
        for c in ("timestamp", "covered_cpgs", "number_probes", "probes")
        if c in cols
    ]
    ordered.extend(c for c in cols if c not in ordered)
    new_df = rows[ordered]

    if os.path.exists(scores_path):
        try:
            existing = pd.read_csv(scores_path)
            combined = pd.concat([existing, new_df], ignore_index=True)
        except Exception:
            combined = new_df
    else:
        combined = new_df

    os.makedirs(os.path.dirname(scores_path) or ".", exist_ok=True)
    combined.to_csv(scores_path, index=False)


def _top_prediction(prediction_df: pd.DataFrame) -> tuple[Optional[str], Optional[float], int]:
    """Return (top_class, top_score, covered_cpgs) from a Tucan output frame."""
    if prediction_df is None or prediction_df.empty:
        return None, None, 0
    row = prediction_df.iloc[-1]
    covered = 0
    for key in ("number_probes", "covered_cpgs", "probes"):
        if key in prediction_df.columns:
            try:
                covered = int(float(row[key]))
                break
            except Exception:
                pass
    best_class = None
    best_score = None
    for col, value in row.items():
        if str(col).strip().lower() in SCORE_META_COLUMNS:
            continue
        try:
            score = float(value)
        except Exception:
            continue
        if best_score is None or score > best_score:
            best_score = score
            best_class = str(col)
    return best_class, best_score, covered


def run_tucan_predict(
    bed_path: str,
    model_zip: str,
    output_path: str,
    *,
    num_cpgs: Optional[int] = None,
    num_samplings: Optional[int] = None,
) -> pd.DataFrame:
    """Call ``tucan.cli.predict`` and return the written scores DataFrame."""
    require_tucan_package()
    from tucan.cli import predict

    n = get_default_num_cpgs() if num_cpgs is None else int(num_cpgs)
    s = get_default_num_samplings() if num_samplings is None else int(num_samplings)

    predict(
        input_file=bed_path,
        model_zip_path=model_zip,
        n=n,
        output_file=output_path,
        num_samples=s,
        file_type="bed",
    )
    return pd.read_csv(output_path)


class TucanAnalysis:
    """Tucan analysis worker."""

    def __init__(
        self,
        work_dir: Optional[str] = None,
        *,
        model_path: Optional[str] = None,
        download_if_missing: bool = True,
        num_cpgs: Optional[int] = None,
        num_samplings: Optional[int] = None,
        probe_margin: Optional[int] = None,
    ):
        self.work_dir = work_dir or os.getcwd()
        self.model_path = model_path
        self.download_if_missing = download_if_missing
        self.num_cpgs = num_cpgs
        self.num_samplings = num_samplings
        self.probe_margin = (
            get_probe_margin() if probe_margin is None else int(probe_margin)
        )
        self.bambatch: Dict[str, int] = {}
        self._assets: Optional[Dict[str, Any]] = None

        logger.info(
            "Tucan Analysis initialized (probe_margin=%s)", self.probe_margin
        )

    def _ensure_assets(self) -> Dict[str, Any]:
        if self._assets is not None:
            return self._assets
        require_tucan_package()
        self._assets = ensure_tucan_assets(
            model_path=self.model_path,
            download_if_missing=self.download_if_missing,
        )
        return self._assets

    def process_parquet_file(
        self, parquet_path: str, sample_id: str
    ) -> TucanMetadata:
        start_time = time.time()

        if sample_id not in self.bambatch:
            self.bambatch[sample_id] = 1

        result = TucanMetadata(
            sample_id=sample_id,
            parquet_path=parquet_path,
            analysis_timestamp=start_time,
            batch_number=self.bambatch[sample_id],
        )

        try:
            if load_modkit_data is None or modkit_pileup_file_to_bed is None:
                raise ImportError("merge_bedmethyl helpers are unavailable")
            if not os.path.exists(parquet_path):
                raise FileNotFoundError(f"Parquet file not found: {parquet_path}")
            result.processing_steps.append("file_validation")

            assets = self._ensure_assets()
            result.processing_steps.append("model_loaded")

            modkit_df = load_modkit_data(parquet_path)
            result.processing_steps.append("modkit_data_loaded")

            sample_dir = os.path.join(self.work_dir, sample_id)
            os.makedirs(sample_dir, exist_ok=True)
            bed_path = os.path.join(sample_dir, "TucanBed.bed")

            with tempfile.NamedTemporaryFile(
                dir=sample_dir, suffix=".bed.tmp", delete=False
            ) as tmp:
                tmp_path = tmp.name
            try:
                bed_df = modkit_pileup_file_to_bed(
                    modkit_df,
                    tmp_path,
                    str(assets["mapping_probes"]),
                    margin=self.probe_margin,
                )
            finally:
                try:
                    if os.path.exists(tmp_path):
                        os.unlink(tmp_path)
                except OSError:
                    pass

            if bed_df is None or bed_df.empty:
                raise ValueError(
                    "No Tucan probes overlapped methylation calls in parquet"
                )
            bed_df = binarize_methylation_calls(bed_df)
            bed_df.to_csv(bed_path, sep="\t", index=False, header=True)
            result.bed_file_path = bed_path
            result.processing_steps.append("bed_file_created")

            raw_out = os.path.join(sample_dir, "tucan_raw_prediction.csv")
            prediction_df = run_tucan_predict(
                bed_path,
                str(assets["model_zip"]),
                raw_out,
                num_cpgs=self.num_cpgs,
                num_samplings=self.num_samplings,
            )
            result.processing_steps.append("prediction_complete")

            top_class, top_score, covered = _top_prediction(prediction_df)
            result.top_class = top_class
            result.top_score = top_score
            result.covered_cpgs = covered

            scores_path = os.path.join(sample_dir, "tucan_scores.csv")
            append_tucan_scores(
                scores_path, prediction_df, timestamp_ms=start_time * 1000
            )
            result.scores_file_path = scores_path
            result.processing_steps.append("results_saved")

            result.results = {
                "status": "success",
                "sample_id": sample_id,
                "analysis_time": time.time() - start_time,
                "batch_number": result.batch_number,
                "covered_cpgs": covered,
                "top_class": top_class,
                "top_score": top_score,
                "scores_file": scores_path,
                "bed_file": bed_path,
                "num_cpgs": self.num_cpgs
                if self.num_cpgs is not None
                else DEFAULT_NUM_CPGS,
                "probe_margin": self.probe_margin,
                "processing_steps": result.processing_steps.copy(),
            }
            self.bambatch[sample_id] += 1
            logger.info(
                "Tucan %s → %s (%.3f) covered_cpgs=%s",
                sample_id,
                top_class,
                top_score or 0.0,
                covered,
            )
            return result

        except Exception as exc:
            logger.error(
                "Tucan analysis failed for %s: %s", sample_id, exc, exc_info=True
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
    model_path: Optional[str] = None,
    num_cpgs: Optional[int] = None,
    num_samplings: Optional[int] = None,
    probe_margin: Optional[int] = None,
    job_id=None,
):
    """Process multiple parquet files for one sample into one tucan_scores.csv."""
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
        analyzer = TucanAnalysis(
            work_dir=work_dir,
            model_path=model_path,
            num_cpgs=num_cpgs,
            num_samplings=num_samplings,
            probe_margin=probe_margin,
        )
        analysis_result["processing_steps"].append("analyzer_created")

        for i, parquet_path in enumerate(parquet_paths):
            logger.info(
                "Tucan file %s/%s: %s",
                i + 1,
                len(parquet_paths),
                os.path.basename(parquet_path),
            )
            if not os.path.exists(parquet_path):
                logger.warning("Parquet not found: %s", parquet_path)
                continue
            tucan_result = analyzer.process_parquet_file(parquet_path, sample_id)
            if tucan_result.error_message:
                logger.warning(
                    "Tucan skipped %s: %s",
                    os.path.basename(parquet_path),
                    tucan_result.error_message,
                )
                continue
            analysis_result["files_processed"] += 1
            analysis_result["total_batches"] += 1
            analysis_result["scores_file_path"] = tucan_result.scores_file_path
            analysis_result["bed_file_path"] = tucan_result.bed_file_path
            try:
                from robin.workflow_ray import notify_coordinator_files_completed

                notify_coordinator_files_completed("tucan", 1, job_id=job_id)
            except Exception:
                pass

        if analysis_result["files_processed"] == 0:
            analysis_result["error_message"] = "No files could be processed successfully"
            analysis_result["processing_steps"].append("no_files_processed")
            return analysis_result

        analysis_result["processing_steps"].append("analysis_complete")
        return analysis_result
    except Exception as exc:
        logger.error("Multi-file Tucan analysis failed for %s: %s", sample_id, exc)
        analysis_result["error_message"] = str(exc)
        analysis_result["processing_steps"].append("analysis_failed")
        return analysis_result


def tucan_handler(job, work_dir=None):
    """Workflow handler for Tucan classification jobs."""
    logger = get_job_logger(str(job.job_id), job.job_type, job.context.filepath)
    suppress_expected = _is_fail_only_expected(job)

    try:
        batched_job = job.context.metadata.get("_batched_job")
        if batched_job:
            batch_size = batched_job.get_file_count()
            sample_id = batched_job.get_sample_id()
            batch_id = batched_job.batch_id
            logger.info(
                "Processing Tucan batch: %s files for sample '%s' (batch_id: %s)",
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
                error_msg = "No parquet paths found from bed conversion results in batch"
                if suppress_expected:
                    logger.warning("%s (expected for fail-only BAM submission)", error_msg)
                    job.context.add_result(
                        "tucan_analysis",
                        {"status": "expected_failure", "reason": error_msg},
                    )
                else:
                    logger.error(error_msg)
                    job.context.add_error("tucan_analysis", error_msg)
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
                job_id=job.job_id,
            )
            job.context.add_metadata(
                "tucan_analysis",
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
                        "tucan_analysis",
                        {
                            "status": "expected_failure",
                            "error_message": batch_result["error_message"],
                            "sample_id": sample_id,
                        },
                    )
                else:
                    job.context.add_error(
                        "tucan_analysis", batch_result["error_message"]
                    )
            else:
                job.context.add_result(
                    "tucan_analysis",
                    {
                        "status": "success",
                        "sample_id": sample_id,
                        "analysis_time": batch_result.get("analysis_timestamp", 0),
                        "processing_steps": batch_result.get("processing_steps", []),
                        "scores_file": batch_result.get("scores_file_path", ""),
                        "bed_file": batch_result.get("bed_file_path", ""),
                        "files_processed": batch_result.get("files_processed", 0),
                        "total_files": batch_result.get("total_files", 0),
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

        analyzer = TucanAnalysis(work_dir=work_dir)
        result = analyzer.process_parquet_file(parquet_path, sample_id)
        job.context.add_metadata("tucan_analysis", result.results)
        job.context.add_metadata("tucan_processing_steps", result.processing_steps)

        if result.error_message:
            if suppress_expected:
                job.context.add_result(
                    "tucan_analysis",
                    {"status": "expected_failure", "error_message": result.error_message},
                )
            else:
                job.context.add_error("tucan_analysis", result.error_message)
                logger.error("Tucan analysis failed: %s", result.error_message)
        else:
            job.context.add_result(
                "tucan_analysis",
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
                },
            )
            try:
                from robin.workflow_ray import notify_coordinator_files_completed

                notify_coordinator_files_completed("tucan", 1, job_id=job.job_id)
            except Exception:
                pass
            logger.info(
                "Tucan complete for %s: %s (%.3f)",
                sample_id,
                result.top_class,
                result.top_score or 0.0,
            )

    except Exception as exc:
        if suppress_expected:
            logger.warning("Expected Tucan failure for fail-only BAM submission: %s", exc)
            job.context.add_result(
                "tucan_analysis",
                {"status": "expected_failure", "error_message": str(exc)},
            )
            return
        job.context.add_error("tucan_analysis", str(exc))
        logger.error("Error in Tucan analysis for %s: %s", job.context.filepath, exc)
