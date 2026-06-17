"""Unified entry point for MNP-Flex analysis (Docker or Epignostix API)."""

from __future__ import annotations

import logging
from pathlib import Path
from typing import Any, Dict, Optional

from robin.analysis.mnpflex_bed import select_input_bed_for_config
from robin.analysis.mnpflex_config import MNPFlexConfig, load_mnpflex_config
from robin.analysis.mnpflex_docker import run_docker_mnpflex
from robin.utils.mnpflex_client_standalone import MNPFlexClient

logger = logging.getLogger(__name__)


def run_mnpflex_analysis(
    *,
    sample_dir: Path,
    sample_id: str,
    output_dir: Path,
    config: Optional[MNPFlexConfig] = None,
) -> Dict[str, Any]:
    """Run MNP-Flex for a sample and write results into output_dir."""
    cfg = config or load_mnpflex_config()
    err = cfg.validation_error()
    if err:
        raise RuntimeError(err)
    if cfg.backend == "disabled":
        raise RuntimeError(
            "MNP-Flex is disabled. Set MNPFLEX_BACKEND to docker or api."
        )

    output_dir.mkdir(parents=True, exist_ok=True)
    bed_path = select_input_bed_for_config(sample_dir, sample_id, cfg)
    logger.info(
        "[MNPFlex] Running sample=%s backend=%s bed=%s output_dir=%s",
        sample_id,
        cfg.backend,
        bed_path,
        output_dir,
    )

    if cfg.backend == "docker":
        summary_path = run_docker_mnpflex(
            config=cfg,
            bed_path=bed_path,
            sample_id=sample_id,
            output_dir=output_dir,
        )
        return {
            "backend": "docker",
            "bed_path": str(bed_path),
            "bundle_summary_path": str(summary_path),
            "docker_image": cfg.docker_image,
        }

    if cfg.backend == "api":
        client = MNPFlexClient(
            base_url=cfg.base_url,
            username=cfg.username or "",
            password=cfg.password or "",
            verify_ssl=False,
            client_id=cfg.client_id,
            client_secret=cfg.client_secret,
            scope=cfg.scope,
        )
        result = client.upload_retrieve_cleanup(
            bed_file_path=str(bed_path),
            sample_identifier=sample_id,
            workflow_id=cfg.workflow_id,
            output_dir=str(output_dir),
        )
        return {
            "backend": "api",
            "bed_path": str(bed_path),
            **result,
        }

    raise RuntimeError(f"Unsupported MNP-Flex backend: {cfg.backend}")
