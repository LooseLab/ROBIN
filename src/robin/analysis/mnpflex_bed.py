"""Shared helpers to prepare MNP-Flex BED inputs from sample parquet."""

from __future__ import annotations

import os
from pathlib import Path
from typing import Optional

from robin import resources
from robin.analysis.mnpflex_config import MNPFlexConfig, MNPFlexDockerInput
from robin.analysis.utilities.matkit import reconstruct_full_bedmethyl_for_mnpflex
from robin.analysis.utilities.mnp_flex import APIClient as MnpFlexApiClient


def find_parquet_file(sample_dir: Path, sample_id: str) -> Optional[Path]:
    preferred = sample_dir / f"{sample_id}.parquet"
    if preferred.exists():
        return preferred
    matches = list(sample_dir.glob("*.parquet"))
    return matches[0] if matches else None


def build_bed_file_from_parquet(sample_dir: Path, sample_id: str) -> Path:
    parquet_path = find_parquet_file(sample_dir, sample_id)
    if not parquet_path or not parquet_path.exists():
        raise RuntimeError("No parquet data found for this sample.")
    bed_path = sample_dir / f"{sample_id}.mnpflex.bed"
    bed_df = reconstruct_full_bedmethyl_for_mnpflex(str(parquet_path))
    bed_df.to_csv(bed_path, sep="\t", index=False, header=False)
    return bed_path


def build_subset_bed_from_parquet(sample_dir: Path, sample_id: str) -> Path:
    bed_path = build_bed_file_from_parquet(sample_dir, sample_id)
    subset_path = sample_dir / f"{sample_id}.MNPFlex.subset.bed"
    reference_bed = os.path.join(
        os.path.dirname(os.path.abspath(resources.__file__)),
        "mnp_flex_sample_clean.bed",
    )
    api_client = MnpFlexApiClient(base_url="https://mnp-flex.org", verify_ssl=False)
    api_client.process_streaming(reference_bed, str(bed_path), str(subset_path))
    return subset_path


def select_input_bed(
    sample_dir: Path,
    sample_id: str,
    *,
    docker_input: MNPFlexDockerInput = "full",
    for_api: bool = False,
) -> Path:
    """Return the BED file to submit for MNP-Flex analysis."""
    if for_api:
        return build_subset_bed_from_parquet(sample_dir, sample_id)
    if docker_input == "subset":
        return build_subset_bed_from_parquet(sample_dir, sample_id)
    return build_bed_file_from_parquet(sample_dir, sample_id)


def select_input_bed_for_config(
    sample_dir: Path,
    sample_id: str,
    config: MNPFlexConfig,
) -> Path:
    if config.backend == "api":
        return select_input_bed(sample_dir, sample_id, for_api=True)
    return select_input_bed(sample_dir, sample_id, docker_input=config.docker_input)
