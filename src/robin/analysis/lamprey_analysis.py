"""
Lamprey methylation classification for ROBIN (research / evaluation only).

Lamprey classifies hematological malignancy subtypes from sparse CpG calls via
an ONNX model. ROBIN does **not** vendor Lamprey; install the upstream package
separately and set ``ROBIN_LAMPREY_RESEARCH_ACK=1`` before first model download.

Pipeline (mirrors MARLIN / Sturgeon)::

  bed_conversion parquet → hg38 probe map → binarized feature vector → ONNX
  → lamprey_scores.csv
"""

from __future__ import annotations

import logging
import os
import time
from dataclasses import dataclass, field
from typing import Any, Dict, List, Mapping, Optional, Sequence

import numpy as np
import pandas as pd

from robin.logging_config import get_job_logger
from robin.utils.lamprey_manager import (
    DEFAULT_GENOME_BUILD,
    RESEARCH_DISCLAIMER,
    ensure_lamprey_assets,
    require_research_ack,
)

try:
    from robin.analysis.utilities.merge_bedmethyl import load_minimal_modkit_data
except ImportError:  # pragma: no cover
    load_minimal_modkit_data = None


def _is_fail_only_expected(job) -> bool:
    try:
        if not bool(job.context.metadata.get("fail_only_bam_submission", False)):
            return False
        fp = getattr(job.context, "filepath", "") or ""
        base = os.path.basename(fp).lower()
        return ("fail" in base) and ("pass" not in base)
    except Exception:
        return False


@dataclass(frozen=True)
class LampreyPrediction:
    scores: Dict[str, float]
    covered_cpgs: int
    temperature: float
    top_class: str
    top_score: float
    diagnostic: bool  # Lamprey's native >=0.90 flag (informational)


@dataclass
class LampreyMetadata:
    sample_id: str
    parquet_path: str
    analysis_timestamp: float
    batch_number: int
    scores_file_path: Optional[str] = None
    covered_cpgs: int = 0
    top_class: Optional[str] = None
    top_score: Optional[float] = None
    diagnostic: Optional[bool] = None
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


def _softmax(x: np.ndarray, axis: int = 1) -> np.ndarray:
    x = x - np.max(x, axis=axis, keepdims=True)
    exp_x = np.exp(x)
    return exp_x / np.sum(exp_x, axis=axis, keepdims=True)


def load_probe_names(bed_path: str | os.PathLike) -> List[str]:
    """Ordered probe IDs from a Lamprey probe BED (skips header)."""
    names: List[str] = []
    with open(bed_path, "r", encoding="utf-8") as handle:
        header = handle.readline()
        if not header:
            raise ValueError(f"Empty probe BED: {bed_path}")
        for line in handle:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split("\t")
            if len(parts) < 4:
                continue
            # Skip accidental re-read of header-like rows
            if parts[0] == "chrom":
                continue
            names.append(parts[3].strip())
    if not names:
        raise ValueError(f"No probe names in {bed_path}")
    return names


def load_probe_position_map(
    probes_bed_path: str | os.PathLike,
) -> Dict[tuple[str, int], str]:
    """
    Build (chrom, position) → probe_id from Lamprey ``probe_hg38.bed``.

    Lamprey intervals are 2 bp (CpG). Map both bases so single-base bedmethyl
    rows still match.
    """
    mapping: Dict[tuple[str, int], str] = {}
    with open(probes_bed_path, "r", encoding="utf-8") as handle:
        header = handle.readline()
        for line in handle:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split("\t")
            if len(parts) < 4 or parts[0] == "chrom":
                continue
            chrom = _normalize_chrom(parts[0])
            try:
                start = int(parts[1])
                end = int(parts[2]) if len(parts) > 2 else start + 2
            except ValueError:
                continue
            probe_id = parts[3].strip()
            if not probe_id:
                continue
            for pos in range(start, max(end, start + 1)):
                mapping[(chrom, pos)] = probe_id
    if not mapping:
        raise ValueError(f"No probe positions loaded from {probes_bed_path}")
    return mapping


def bedmethyl_to_probe_calls(
    modkit_df: pd.DataFrame,
    position_map: Mapping[tuple[str, int], str],
) -> Dict[str, int]:
    """
    Map bedmethyl rows to Lamprey methylation calls ``{-1, 1}``.

    Matches upstream ``pileup_to_lamprey``: percent > 50 → 1, percent < 50 → -1,
    percent == 50 dropped. Multiple hits for one probe: majority / last-write
    after averaging percent then thresholding.
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
            "Methylation dataframe missing percent_modified/fraction/score"
        )

    sums: Dict[str, float] = {}
    counts: Dict[str, int] = {}
    for chrom, start, raw in zip(
        df[chrom_col].tolist(), df[start_col].tolist(), df[value_col].tolist()
    ):
        try:
            start_i = int(start)
            value = float(raw)
        except (TypeError, ValueError):
            continue
        # Normalize to percent 0–100
        percent = value * 100.0 if value <= 1.0 else value
        probe_id = position_map.get((_normalize_chrom(chrom), start_i))
        if probe_id is None:
            continue
        sums[probe_id] = sums.get(probe_id, 0.0) + percent
        counts[probe_id] = counts.get(probe_id, 0) + 1

    calls: Dict[str, int] = {}
    for probe_id, total in sums.items():
        avg = total / counts[probe_id]
        if avg == 50.0:
            continue
        calls[probe_id] = 1 if avg > 50.0 else -1
    return calls


def build_feature_vector(
    probe_names: Sequence[str],
    probe_calls: Mapping[str, int],
) -> tuple[np.ndarray, int]:
    """Return float32 feature vector in model order and covered-probe count."""
    vector = np.zeros(len(probe_names), dtype=np.float32)
    for idx, name in enumerate(probe_names):
        call = probe_calls.get(name)
        if call is None:
            continue
        # Upstream maps 0 → -1; we already emit ±1 only.
        vector[idx] = float(-1 if call == 0 else call)
    n_used = int(np.count_nonzero(vector))
    return vector, n_used


def append_lamprey_scores(
    scores_path: str,
    prediction: LampreyPrediction,
    *,
    timestamp_ms: float,
) -> None:
    row: Dict[str, Any] = {
        "timestamp": timestamp_ms,
        "covered_cpgs": prediction.covered_cpgs,
        "number_probes": prediction.covered_cpgs,
        "temperature": prediction.temperature,
        "diagnostic": int(prediction.diagnostic),
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

    cols = list(combined.columns)
    ordered = [
        c
        for c in (
            "timestamp",
            "covered_cpgs",
            "number_probes",
            "temperature",
            "diagnostic",
        )
        if c in cols
    ]
    ordered.extend(c for c in cols if c not in ordered)
    combined[ordered].to_csv(scores_path, index=False)


class LampreyPredictor:
    """In-memory ONNX predictor mirroring Lamprey nanopore inference."""

    def __init__(
        self,
        session,
        probe_names: Sequence[str],
        class_names: Sequence[str],
        bin_centers: np.ndarray,
        temps: np.ndarray,
        *,
        display_names: Optional[Mapping[str, str]] = None,
        input_name: str,
        output_name: str,
    ):
        self.session = session
        self.probe_names = list(probe_names)
        self.class_names = list(class_names)
        self.bin_centers = np.asarray(bin_centers)
        self.temps = np.asarray(temps)
        self.display_names = dict(display_names or {})
        self.input_name = input_name
        self.output_name = output_name

    @classmethod
    def from_assets(cls, assets: Mapping[str, os.PathLike | str]) -> "LampreyPredictor":
        import yaml

        try:
            import onnxruntime as ort
        except ModuleNotFoundError as exc:
            raise ModuleNotFoundError(
                "onnxruntime is required for Lamprey. "
                "Install with: pip install 'robin[lamprey]'"
            ) from exc

        model_dir = assets["model_dir"]
        onnx_path = os.path.join(str(model_dir), "model.onnx")
        if not os.path.isfile(onnx_path):
            raise FileNotFoundError(f"model.onnx not found in {model_dir}")

        with open(assets["classification_yaml"], "r", encoding="utf-8") as handle:
            classification = yaml.safe_load(handle)
        class_names = list(classification["encoder"]["type"].keys())
        display_names = classification.get("names", {}).get("type", {})

        probe_names = load_probe_names(assets["feature_order_bed"])
        bin_centers = np.load(os.path.join(str(model_dir), "bin_centers.npy"))
        temps = np.load(os.path.join(str(model_dir), "temps.npy"))

        so = ort.SessionOptions()
        so.intra_op_num_threads = 1
        so.inter_op_num_threads = 1
        sess = ort.InferenceSession(
            onnx_path, sess_options=so, providers=["CPUExecutionProvider"]
        )
        return cls(
            session=sess,
            probe_names=probe_names,
            class_names=class_names,
            bin_centers=bin_centers,
            temps=temps,
            display_names=display_names,
            input_name=sess.get_inputs()[0].name,
            output_name=sess.get_outputs()[0].name,
        )

    def predict_from_probe_calls(
        self, probe_calls: Mapping[str, int]
    ) -> LampreyPrediction:
        vector, n_used = build_feature_vector(self.probe_names, probe_calls)
        temperature = float(self.temps[np.argmin(np.abs(self.bin_centers - n_used))])
        batch = vector.reshape(1, -1)
        outputs = self.session.run([self.output_name], {self.input_name: batch})[0]
        outputs = outputs / np.exp(temperature)
        outputs = _softmax(outputs, axis=1)
        probs = outputs[0]
        scores = {name: float(score) for name, score in zip(self.class_names, probs)}
        top_idx = int(np.argmax(probs))
        top_score = float(probs[top_idx])
        return LampreyPrediction(
            scores=scores,
            covered_cpgs=n_used,
            temperature=temperature,
            top_class=self.class_names[top_idx],
            top_score=top_score,
            diagnostic=top_score >= 0.9,
        )


class LampreyAnalysis:
    """Lamprey analysis worker (research / evaluation only)."""

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
        self._predictor: Optional[LampreyPredictor] = None
        self._position_map: Optional[Dict[tuple[str, int], str]] = None
        logging.getLogger("robin.lamprey").info(
            "Lamprey Analysis initialized (genome=%s). %s",
            self.genome_build,
            RESEARCH_DISCLAIMER,
        )

    def _ensure_predictor(self) -> LampreyPredictor:
        if self._predictor is not None:
            return self._predictor
        require_research_ack()
        assets = ensure_lamprey_assets(
            genome_build=self.genome_build,
            model_path=self.model_path,
            download_if_missing=self.download_if_missing,
        )
        self._position_map = load_probe_position_map(assets["probes_bed"])
        self._predictor = LampreyPredictor.from_assets(assets)
        return self._predictor

    def process_parquet_file(
        self, parquet_path: str, sample_id: str
    ) -> LampreyMetadata:
        logger = logging.getLogger("robin.lamprey")
        start_time = time.time()
        if sample_id not in self.bambatch:
            self.bambatch[sample_id] = 1

        result = LampreyMetadata(
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

            probe_calls = bedmethyl_to_probe_calls(modkit_df, self._position_map)
            if not probe_calls:
                raise ValueError(
                    "No Lamprey probes overlapped methylation calls in parquet "
                    f"(genome_build={self.genome_build})"
                )
            result.processing_steps.append("probes_mapped")

            prediction = predictor.predict_from_probe_calls(probe_calls)
            result.covered_cpgs = prediction.covered_cpgs
            result.top_class = prediction.top_class
            result.top_score = prediction.top_score
            result.diagnostic = prediction.diagnostic
            result.processing_steps.append("prediction_complete")

            sample_dir = os.path.join(self.work_dir, sample_id)
            os.makedirs(sample_dir, exist_ok=True)
            scores_path = os.path.join(sample_dir, "lamprey_scores.csv")
            append_lamprey_scores(
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
                "diagnostic": prediction.diagnostic,
                "scores_file": scores_path,
                "genome_build": self.genome_build,
                "research_use_only": True,
                "processing_steps": result.processing_steps.copy(),
            }
            self.bambatch[sample_id] += 1
            logger.info(
                "Lamprey %s → %s (%.3f) covered=%s diagnostic=%s [research/evaluation only]",
                sample_id,
                prediction.top_class,
                prediction.top_score,
                prediction.covered_cpgs,
                prediction.diagnostic,
            )
            return result
        except Exception as exc:
            logger.error(
                "Lamprey analysis failed for %s: %s", sample_id, exc, exc_info=True
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
    }
    try:
        os.makedirs(os.path.join(work_dir, sample_id), exist_ok=True)
        analyzer = LampreyAnalysis(
            work_dir=work_dir,
            genome_build=genome_build,
            model_path=model_path,
        )
        analysis_result["processing_steps"].append("analyzer_created")
        for i, parquet_path in enumerate(parquet_paths):
            logger.info(
                "Lamprey file %s/%s: %s",
                i + 1,
                len(parquet_paths),
                os.path.basename(parquet_path),
            )
            if not os.path.exists(parquet_path):
                logger.warning("Parquet not found: %s", parquet_path)
                continue
            lamprey_result = analyzer.process_parquet_file(parquet_path, sample_id)
            if lamprey_result.error_message:
                logger.warning(
                    "Lamprey skipped %s: %s",
                    os.path.basename(parquet_path),
                    lamprey_result.error_message,
                )
                continue
            analysis_result["files_processed"] += 1
            analysis_result["total_batches"] += 1
            analysis_result["scores_file_path"] = lamprey_result.scores_file_path
            try:
                from robin.workflow_ray import notify_coordinator_files_completed

                notify_coordinator_files_completed("lamprey", 1, job_id=job_id)
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
        logger.error("Multi-file Lamprey analysis failed for %s: %s", sample_id, exc)
        analysis_result["error_message"] = str(exc)
        analysis_result["processing_steps"].append("analysis_failed")
        return analysis_result


def lamprey_handler(job, work_dir=None):
    """Workflow handler for Lamprey classification (research / evaluation only)."""
    logger = get_job_logger(str(job.job_id), job.job_type, job.context.filepath)
    suppress_expected = _is_fail_only_expected(job)
    logger.warning(RESEARCH_DISCLAIMER)

    try:
        batched_job = job.context.metadata.get("_batched_job")
        if batched_job:
            batch_size = batched_job.get_file_count()
            sample_id = batched_job.get_sample_id()
            batch_id = batched_job.batch_id
            logger.info(
                "Processing Lamprey batch: %s files for sample '%s' (batch_id: %s)",
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
                    job.context.add_result(
                        "lamprey_analysis",
                        {"status": "expected_failure", "reason": error_msg},
                    )
                else:
                    job.context.add_error("lamprey_analysis", error_msg)
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
                genome_build=DEFAULT_GENOME_BUILD,
                job_id=job.job_id,
            )
            job.context.add_metadata(
                "lamprey_analysis",
                {
                    "batch_result": batch_result,
                    "batch_size": batch_size,
                    "sample_id": sample_id,
                    "batch_id": batch_id,
                    "files_processed": batch_result.get("files_processed", 0),
                    "total_files": batch_result.get("total_files", 0),
                    "research_use_only": True,
                },
            )
            if batch_result.get("error_message"):
                if suppress_expected:
                    job.context.add_result(
                        "lamprey_analysis",
                        {
                            "status": "expected_failure",
                            "error_message": batch_result["error_message"],
                            "sample_id": sample_id,
                        },
                    )
                else:
                    job.context.add_error(
                        "lamprey_analysis", batch_result["error_message"]
                    )
            else:
                job.context.add_result(
                    "lamprey_analysis",
                    {
                        "status": "success",
                        "sample_id": sample_id,
                        "scores_file": batch_result.get("scores_file_path", ""),
                        "files_processed": batch_result.get("files_processed", 0),
                        "total_files": batch_result.get("total_files", 0),
                        "genome_build": DEFAULT_GENOME_BUILD,
                        "research_use_only": True,
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

        analyzer = LampreyAnalysis(work_dir=work_dir, genome_build=DEFAULT_GENOME_BUILD)
        result = analyzer.process_parquet_file(parquet_path, sample_id)
        job.context.add_metadata("lamprey_analysis", result.results)
        job.context.add_metadata("lamprey_processing_steps", result.processing_steps)

        if result.error_message:
            if suppress_expected:
                job.context.add_result(
                    "lamprey_analysis",
                    {
                        "status": "expected_failure",
                        "error_message": result.error_message,
                    },
                )
            else:
                job.context.add_error("lamprey_analysis", result.error_message)
        else:
            job.context.add_result(
                "lamprey_analysis",
                {
                    "status": "success",
                    "sample_id": result.sample_id,
                    "scores_file": result.scores_file_path,
                    "covered_cpgs": result.covered_cpgs,
                    "top_class": result.top_class,
                    "top_score": result.top_score,
                    "diagnostic": result.diagnostic,
                    "genome_build": DEFAULT_GENOME_BUILD,
                    "research_use_only": True,
                },
            )
            try:
                from robin.workflow_ray import notify_coordinator_files_completed

                notify_coordinator_files_completed("lamprey", 1, job_id=job.job_id)
            except Exception:
                pass
    except Exception as exc:
        if suppress_expected:
            job.context.add_result(
                "lamprey_analysis",
                {"status": "expected_failure", "error_message": str(exc)},
            )
            return
        job.context.add_error("lamprey_analysis", str(exc))
        logger.error("Error in Lamprey analysis for %s: %s", job.context.filepath, exc)
