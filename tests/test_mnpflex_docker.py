"""Tests for MNP-Flex Docker output adaptation."""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from robin.analysis.mnpflex_docker import (
    adapt_docker_outputs,
    build_bundle_summary_from_docker_dir,
    find_docker_output_dir,
    format_mnpflex_runtime_error,
    run_docker_mnpflex,
)
from robin.analysis.mnpflex_config import MNPFlexConfig


FIXTURE_DIR = (
    Path(__file__).resolve().parent / "fixtures" / "mnpflex_docker" / "26D22147.MNPFlex"
)


def test_format_mnpflex_runtime_error_docker_daemon() -> None:
    raw = (
        "docker daemon is not available: Cannot connect to the Docker daemon "
        "at unix:///Users/me/.docker/run/docker.sock"
    )
    msg = format_mnpflex_runtime_error(raw)
    assert "Docker is not running" in msg
    assert "Docker Desktop" in msg


def test_build_bundle_summary_from_fixture() -> None:
    summary = build_bundle_summary_from_docker_dir(
        FIXTURE_DIR, docker_image="mnpflex-synnovis:1.0.0"
    )
    assert summary["source"] == "docker"
    assert summary["qc"]["missing_site_count"] == 1438.0
    assert summary["qc"]["avg_coverage"] == 3.25
    assert "WARNING" in str(summary["qc"]["status"])
    assert summary["mgmt"]["status"] == "Unmethylated"
    assert summary["mgmt"]["site_count"] == 136.0
    assert summary["classifier_summary"]["classifier"]["name"] == "Brain"
    assert summary["classifier_summary"]["classifier"]["version"] == "12.8"
    hierarchy = summary["classifier_summary"]["summary_hierarchical"]
    assert hierarchy
    assert hierarchy[0]["group"] == "Adult-Type Diffuse Gliomas"
    scores = summary["classifier_summary"]["scores"]
    assert len(scores) >= 180
    assert scores[0]["reference_group"]["molecular_subclass"]


def test_adapt_docker_outputs_writes_bundle_summary(tmp_path: Path) -> None:
    output_dir = tmp_path / "mnpflex_results_sample"
    summary_path = adapt_docker_outputs(
        FIXTURE_DIR,
        output_dir,
        docker_image="mnpflex-synnovis:1.0.0",
    )
    assert summary_path == output_dir / "bundle_summary.json"
    assert summary_path.exists()
    payload = json.loads(summary_path.read_text(encoding="utf-8"))
    assert payload["source"] == "docker"
    assert (output_dir / "docker_raw").exists()


def test_find_docker_output_dir_supports_nested_layout(tmp_path: Path) -> None:
    nested = tmp_path / "docker_workspace" / "sample.bedstem"
    nested.mkdir(parents=True)
    for path in FIXTURE_DIR.glob("*.csv"):
        (nested / path.name).write_text(path.read_text(encoding="utf-8"), encoding="utf-8")
    found = find_docker_output_dir(tmp_path / "docker_workspace", "sample.bedstem")
    assert found == nested


def test_run_docker_mnpflex_invokes_container(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    bed_path = tmp_path / "sample.mnpflex.bed"
    bed_path.write_text("chr1\t1\t2\n", encoding="utf-8")
    output_dir = tmp_path / "results"
    docker_dir = output_dir / "docker_workspace" / bed_path.stem
    docker_dir.mkdir(parents=True)

    for path in FIXTURE_DIR.glob("*.csv"):
        (docker_dir / path.name).write_text(path.read_text(encoding="utf-8"), encoding="utf-8")

    captured: dict = {}

    def fake_run(cmd, **kwargs):
        captured["cmd"] = cmd
        class Result:
            returncode = 0
            stdout = ""
            stderr = ""

        return Result()

    monkeypatch.setattr("robin.analysis.mnpflex_docker.subprocess.run", fake_run)
    monkeypatch.setattr(
        "robin.analysis.mnpflex_docker.validate_docker_runtime", lambda config: None
    )

    config = MNPFlexConfig(
        backend="docker",
        docker_image="registry.example.com/mnpflex:2.0.0",
        docker_timeout_s=60,
        docker_binary="docker",
        docker_extra_args=["--memory=8g"],
        docker_input="full",
        username=None,
        password=None,
        base_url="https://app.epignostix.com",
        workflow_id=18,
        client_id="ROBIN",
        client_secret="SECRET",
        scope="",
    )
    summary_path = run_docker_mnpflex(
        config=config,
        bed_path=bed_path,
        sample_id="sample",
        output_dir=output_dir,
    )
    assert summary_path.exists()
    cmd = captured["cmd"]
    assert "registry.example.com/mnpflex:2.0.0" in cmd
    assert "--memory=8g" in cmd
    assert "/input/sample.mnpflex.bed" in cmd
    assert "--sample" in cmd
