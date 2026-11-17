"""Tests for workflow TOML configuration."""

from __future__ import annotations

from pathlib import Path

import click
import pytest
from click.testing import CliRunner

from robin.cli import main
from robin.workflow_config import (
    WORKFLOW_CONFIG_EXAMPLE,
    load_workflow_toml,
    merge_workflow_params,
)


def test_load_workflow_toml_parses_example(tmp_path: Path) -> None:
    config_path = tmp_path / "workflow.toml"
    config_path.write_text(
        """
path = "incoming_bams"
workflow = ["cnv", "mgmt", "sturgeon"]
center = "NUH"
target_panel = "rCNS2"
work_dir = "~/runs/demo"
reference = "/refs/genome.fa"
log_level = "INFO"
deduplicate_jobs = ["sturgeon", "mgmt"]
with_gui = false
""",
        encoding="utf-8",
    )

    config = load_workflow_toml(config_path)

    assert config["path"] == Path("incoming_bams")
    assert config["workflow"] == "cnv,mgmt,sturgeon"
    assert config["center"] == "NUH"
    assert config["target_panel"] == "rCNS2"
    assert config["work_dir"] == Path.home() / "runs/demo"
    assert config["reference"] == Path("/refs/genome.fa")
    assert config["log_level"] == "INFO"
    assert config["deduplicate_jobs"] == ("sturgeon", "mgmt")
    assert config["with_gui"] is False


def test_merge_workflow_params_uses_toml_defaults() -> None:
    @click.command()
    @click.pass_context
    @click.option("--toml", "toml_path", type=click.Path(path_type=Path), default=None)
    @click.option("--workflow", "-w", default=None)
    @click.option("--center", default=None)
    @click.option("--target-panel", "target_panel", default=None)
    @click.argument("path", required=False, type=click.Path(path_type=Path))
    def _cmd(
        ctx: click.Context,
        toml_path: Path | None,
        path: Path | None,
        workflow: str | None,
        center: str | None,
        target_panel: str | None,
    ) -> None:
        merged = merge_workflow_params(
            ctx,
            toml_path,
            {
                "path": path,
                "workflow": workflow,
                "center": center,
                "target_panel": target_panel,
            },
        )
        click.echo(
            f"{merged['path']}|{merged['workflow']}|{merged['center']}|{merged['target_panel']}"
        )

    runner = CliRunner()

    def _fake_load(path: Path) -> dict[str, object]:
        return {
            "path": Path("/data/bams"),
            "workflow": "mgmt,sturgeon",
            "center": "NUH",
            "target_panel": "rCNS2",
        }

    with pytest.MonkeyPatch.context() as monkeypatch:
        monkeypatch.setattr("robin.workflow_config.load_workflow_toml", _fake_load)
        result = runner.invoke(_cmd, ["--toml", "settings.toml"])

    assert result.exit_code == 0
    assert result.output.strip() == "/data/bams|mgmt,sturgeon|NUH|rCNS2"


def test_merge_workflow_params_cli_overrides_toml() -> None:
    @click.command()
    @click.pass_context
    @click.option("--toml", "toml_path", type=click.Path(path_type=Path), default=None)
    @click.option("--workflow", "-w", default=None)
    @click.option("--center", default=None)
    @click.option("--target-panel", "target_panel", default=None)
    @click.argument("path", required=False, type=click.Path(path_type=Path))
    def _cmd(
        ctx: click.Context,
        toml_path: Path | None,
        path: Path | None,
        workflow: str | None,
        center: str | None,
        target_panel: str | None,
    ) -> None:
        merged = merge_workflow_params(
            ctx,
            toml_path,
            {
                "path": path,
                "workflow": workflow,
                "center": center,
                "target_panel": target_panel,
            },
        )
        click.echo(f"{merged['path']}|{merged['workflow']}|{merged['center']}")

    runner = CliRunner()

    def _fake_load(path: Path) -> dict[str, object]:
        return {
            "path": Path("/data/bams"),
            "workflow": "mgmt,sturgeon",
            "center": "NUH",
            "target_panel": "rCNS2",
        }

    with pytest.MonkeyPatch.context() as monkeypatch:
        monkeypatch.setattr("robin.workflow_config.load_workflow_toml", _fake_load)
        result = runner.invoke(
            _cmd,
            [
                "--toml",
                "settings.toml",
                "/override/path",
                "-w",
                "cnv",
                "--center",
                "Auckland",
            ],
        )

    assert result.exit_code == 0
    assert result.output.strip() == "/override/path|cnv|Auckland"


def test_merge_workflow_params_requires_required_fields() -> None:
    @click.command()
    @click.pass_context
    @click.option("--workflow", "-w", default=None)
    @click.option("--center", default=None)
    @click.option("--target-panel", "target_panel", default=None)
    @click.argument("path", required=False, type=click.Path(path_type=Path))
    def _cmd(
        ctx: click.Context,
        path: Path | None,
        workflow: str | None,
        center: str | None,
        target_panel: str | None,
    ) -> None:
        merge_workflow_params(
            ctx,
            None,
            {
                "path": path,
                "workflow": workflow,
                "center": center,
                "target_panel": target_panel,
            },
        )

    runner = CliRunner()
    result = runner.invoke(_cmd, [])
    assert result.exit_code != 0
    assert "Missing required workflow setting" in result.output


def test_workflow_example_toml_is_valid() -> None:
    example = Path(__file__).resolve().parents[1] / "examples" / "workflow.example.toml"
    config = load_workflow_toml(example)
    assert config["center"] == "NUH"
    assert "cnv" in config["workflow"]


def test_workflow_config_example_documents_keys() -> None:
    assert "path =" in WORKFLOW_CONFIG_EXAMPLE
    assert "target_panel" in WORKFLOW_CONFIG_EXAMPLE
