"""Run MNP-Flex locally via Docker and adapt outputs for ROBIN."""

from __future__ import annotations

import csv
import json
import logging
import re
import shutil
import subprocess
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

from robin.analysis.mnpflex_config import MNPFlexConfig

logger = logging.getLogger(__name__)

_QC_SUMMARY_SUFFIX = "_qc_summary.csv"
_PNG_MAPPINGS = (
    ("_coverage.png", "qc_coverage_plot.png"),
    ("_methylation_density.png", "qc_methylation_density_plot.png"),
    ("_mgmt_region.png", "mgmt_region_plot.png"),
)


def is_docker_daemon_available(binary: str = "docker") -> Tuple[bool, str]:
    try:
        result = subprocess.run(
            [binary, "info"],
            capture_output=True,
            text=True,
            timeout=30,
            check=False,
        )
    except FileNotFoundError:
        return False, f"{binary} executable not found on PATH."
    except subprocess.TimeoutExpired:
        return False, f"{binary} info timed out."
    if result.returncode != 0:
        detail = (result.stderr or result.stdout or "").strip()
        return False, f"{binary} daemon is not available: {detail}"
    return True, ""


def is_docker_image_available(image: str, binary: str = "docker") -> Tuple[bool, str]:
    try:
        result = subprocess.run(
            [binary, "image", "inspect", image],
            capture_output=True,
            text=True,
            timeout=30,
            check=False,
        )
    except FileNotFoundError:
        return False, f"{binary} executable not found on PATH."
    except subprocess.TimeoutExpired:
        return False, f"{binary} image inspect timed out."
    if result.returncode != 0:
        detail = (result.stderr or result.stdout or "").strip()
        return (
            False,
            f"Docker image {image!r} is not available locally. "
            f"Load or pull the image first. {detail}",
        )
    return True, ""


def validate_docker_runtime(config: MNPFlexConfig) -> None:
    err = config.validation_error()
    if err:
        raise RuntimeError(err)
    ok, message = is_docker_daemon_available(config.docker_binary)
    if not ok:
        raise RuntimeError(message)
    ok, message = is_docker_image_available(
        config.docker_image or "", config.docker_binary
    )
    if not ok:
        raise RuntimeError(message)


def _read_metric_value_csv(path: Path) -> Dict[str, str]:
    out: Dict[str, str] = {}
    with path.open("r", newline="", encoding="utf-8") as fh:
        reader = csv.DictReader(fh)
        for row in reader:
            metric = (row.get("metric") or "").strip()
            if metric:
                out[metric] = (row.get("value") or "").strip()
    return out


def _read_single_row_csv(path: Path) -> Dict[str, str]:
    with path.open("r", newline="", encoding="utf-8") as fh:
        reader = csv.DictReader(fh)
        row = next(reader, None)
    return dict(row or {})


def _parse_classifier_label(label: str) -> Dict[str, str]:
    text = (label or "").strip()
    if not text:
        return {"name": "Unknown", "version": "Unknown", "classifier_type": "docker"}
    match = re.match(r"^(?P<name>.+?)\s+v(?P<version>\S+)$", text)
    if match:
        return {
            "name": match.group("name").strip(),
            "version": match.group("version").strip(),
            "classifier_type": "docker",
        }
    return {"name": text, "version": "Unknown", "classifier_type": "docker"}


def _safe_float(value: Any) -> Optional[float]:
    if value is None or value == "":
        return None
    try:
        return float(value)
    except (TypeError, ValueError):
        return None


def _build_hierarchy_from_lims(lims_row: Dict[str, str], description: str) -> List[Dict[str, Any]]:
    levels = [
        (lims_row.get("Super Family"), lims_row.get("Score")),
        (lims_row.get("Family"), lims_row.get("Score.1")),
        (lims_row.get("Class"), lims_row.get("Score.2")),
        (lims_row.get("Subclass"), lims_row.get("Score.3")),
    ]

    node: Optional[Dict[str, Any]] = None
    for idx, (group, score) in enumerate(reversed(levels)):
        group_text = (group or "").strip()
        if not group_text:
            continue
        current: Dict[str, Any] = {
            "group": group_text,
            "score": _safe_float(score),
            "description": description if idx == 0 else "",
            "members": [node] if node else [],
        }
        node = current
    return [node] if node else []


def _build_scores_from_cal(path: Path) -> List[Dict[str, Any]]:
    scores: List[Dict[str, Any]] = []
    with path.open("r", newline="", encoding="utf-8") as fh:
        reader = csv.reader(fh)
        header = next(reader, None)
        for row in reader:
            if len(row) < 2:
                continue
            subclass = (row[0] or "").strip()
            score = _safe_float(row[1])
            if not subclass or score is None:
                continue
            scores.append(
                {
                    "score": score,
                    "reference_group": {
                        "molecular_subclass": subclass,
                        "name": subclass,
                    },
                }
            )
    scores.sort(key=lambda item: float(item.get("score") or 0), reverse=True)
    return scores


def _discover_file_prefix(docker_dir: Path) -> str:
    for path in sorted(docker_dir.glob(f"*{_QC_SUMMARY_SUFFIX}")):
        return path.name[: -len(_QC_SUMMARY_SUFFIX)]
    raise FileNotFoundError(
        f"No Docker QC summary file (*{_QC_SUMMARY_SUFFIX}) found in {docker_dir}"
    )


def build_bundle_summary_from_docker_dir(
    docker_dir: Path,
    *,
    docker_image: Optional[str] = None,
) -> Dict[str, Any]:
    prefix = _discover_file_prefix(docker_dir)
    qc_path = docker_dir / f"{prefix}{_QC_SUMMARY_SUFFIX}"
    mgmt_path = docker_dir / f"{prefix}_mgmt_status.csv"
    lims_path = docker_dir / f"{prefix}_lims.csv"
    annotation_path = docker_dir / f"{prefix}_annotation.csv"
    scores_path = docker_dir / f"{prefix}_scores_cal.csv"

    if not qc_path.exists():
        raise FileNotFoundError(f"Missing Docker QC summary: {qc_path}")

    qc_metrics = _read_metric_value_csv(qc_path)
    mgmt_row = _read_single_row_csv(mgmt_path) if mgmt_path.exists() else {}
    lims_row = _read_single_row_csv(lims_path) if lims_path.exists() else {}
    annotation_row = (
        _read_single_row_csv(annotation_path) if annotation_path.exists() else {}
    )
    description = (annotation_row.get("Description") or "").strip()
    classifier = _parse_classifier_label(lims_row.get("Classifier", ""))

    hierarchy = _build_hierarchy_from_lims(lims_row, description)
    scores = _build_scores_from_cal(scores_path) if scores_path.exists() else []

    return {
        "qc": {
            "status": qc_metrics.get("status", "Unknown"),
            "avg_coverage": _safe_float(qc_metrics.get("average_coverage")),
            "missing_site_count": _safe_float(qc_metrics.get("missing_sites_count")),
        },
        "mgmt": {
            "status": mgmt_row.get("status", "Unknown"),
            "average": _safe_float(mgmt_row.get("average")),
            "site_count": _safe_float(mgmt_row.get("sites")),
        },
        "classifier_summary": {
            "classifier": classifier,
            "scores": scores,
            "summary_hierarchical": hierarchy,
        },
        "source": "docker",
        "docker_image": docker_image,
        "docker_file_prefix": prefix,
    }


def find_docker_output_dir(search_root: Path, bed_stem: str) -> Path:
    """Locate the directory containing Docker CSV outputs."""
    candidates = [
        search_root / bed_stem,
        search_root,
    ]
    for candidate in candidates:
        if candidate.is_dir() and list(candidate.glob(f"*{_QC_SUMMARY_SUFFIX}")):
            return candidate
    matches = sorted(search_root.rglob(f"*{_QC_SUMMARY_SUFFIX}"))
    if matches:
        return matches[0].parent
    raise FileNotFoundError(
        f"Could not find Docker MNP-Flex outputs under {search_root} "
        f"for input stem {bed_stem!r}."
    )


def adapt_docker_outputs(
    docker_dir: Path,
    output_dir: Path,
    *,
    docker_image: Optional[str] = None,
) -> Path:
    """Write bundle_summary.json and standard plot names into output_dir."""
    output_dir.mkdir(parents=True, exist_ok=True)
    raw_dir = output_dir / "docker_raw"
    if docker_dir.resolve() != raw_dir.resolve():
        raw_dir.mkdir(parents=True, exist_ok=True)
        for path in docker_dir.iterdir():
            dest = raw_dir / path.name
            if path.is_dir():
                if dest.exists():
                    shutil.rmtree(dest)
                shutil.copytree(path, dest)
            else:
                shutil.copy2(path, dest)
        docker_dir = raw_dir

    summary = build_bundle_summary_from_docker_dir(
        docker_dir, docker_image=docker_image
    )
    prefix = summary.get("docker_file_prefix", "")
    summary_path = output_dir / "bundle_summary.json"
    with summary_path.open("w", encoding="utf-8") as fh:
        json.dump(summary, fh, indent=2)

    for src_suffix, dest_name in _PNG_MAPPINGS:
        src = docker_dir / f"{prefix}{src_suffix}"
        if src.exists():
            shutil.copy2(src, output_dir / dest_name)

    return summary_path


def run_docker_mnpflex(
    *,
    config: MNPFlexConfig,
    bed_path: Path,
    sample_id: str,
    output_dir: Path,
) -> Path:
    """Execute the configured Docker image and adapt outputs for ROBIN."""
    validate_docker_runtime(config)
    if not bed_path.exists():
        raise FileNotFoundError(f"MNP-Flex input BED not found: {bed_path}")

    docker_workspace = output_dir / "docker_workspace"
    docker_workspace.mkdir(parents=True, exist_ok=True)
    input_mount = bed_path.parent.resolve()
    output_mount = docker_workspace.resolve()
    container_input = f"/input/{bed_path.name}"

    cmd = [
        config.docker_binary,
        "run",
        "--rm",
        "-v",
        f"{input_mount}:/input:ro",
        "-v",
        f"{output_mount}:/output",
        *config.docker_extra_args,
        config.docker_image or "",
        "--input",
        container_input,
        "--sample",
        sample_id,
    ]
    logger.info("[MNPFlex] Running Docker: %s", " ".join(cmd))
    try:
        result = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            timeout=config.docker_timeout_s,
            check=False,
        )
    except subprocess.TimeoutExpired as exc:
        raise TimeoutError(
            f"MNP-Flex Docker run timed out after {config.docker_timeout_s}s."
        ) from exc

    if result.returncode != 0:
        stderr = (result.stderr or "").strip()
        stdout = (result.stdout or "").strip()
        detail = stderr or stdout or f"exit code {result.returncode}"
        raise RuntimeError(f"MNP-Flex Docker run failed: {detail}")

    docker_output_dir = find_docker_output_dir(output_mount, bed_path.stem)
    return adapt_docker_outputs(
        docker_output_dir,
        output_dir,
        docker_image=config.docker_image,
    )
