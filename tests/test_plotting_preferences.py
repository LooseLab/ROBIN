from __future__ import annotations

from pathlib import Path

from robin.gui.plotting_preferences import (
    CNV_REPORT_SCALE_NORMALIZED_DIFFERENCE,
    CNV_REPORT_SCALE_PLOIDY,
    PLOTTING_PREFERENCES_KEY,
    PlottingPreferencesConfig,
    cnv_report_plot_caption,
    resolve_cnv_report_scale,
    resolve_cnv_summary_normalized,
)
from robin.reference_contigs import DEFAULT_REFERENCE_CONTIG_SCOPE
from robin.security import SecurityStore


def test_plotting_preferences_defaults() -> None:
    config = PlottingPreferencesConfig()
    assert config.cnv_report_scale == CNV_REPORT_SCALE_PLOIDY
    assert config.reference_contig_scope == DEFAULT_REFERENCE_CONTIG_SCOPE


def test_plotting_preferences_legacy_log2_alias() -> None:
    config = PlottingPreferencesConfig.from_dict({"cnv_report_scale": "log2_ratio"})
    assert config.cnv_report_scale == CNV_REPORT_SCALE_NORMALIZED_DIFFERENCE


def test_cnv_report_plot_caption_normalized() -> None:
    caption = cnv_report_plot_caption(CNV_REPORT_SCALE_NORMALIZED_DIFFERENCE)
    assert "log" in caption.lower()


def test_security_store_plotting_preferences_roundtrip(tmp_path: Path) -> None:
    store = SecurityStore(db_path=tmp_path / "security.db")
    payload = PlottingPreferencesConfig(
        cnv_report_scale=CNV_REPORT_SCALE_NORMALIZED_DIFFERENCE,
    ).to_dict()
    store.set_gui_setting(PLOTTING_PREFERENCES_KEY, payload, updated_by_user_id=None)
    loaded = store.get_gui_setting(PLOTTING_PREFERENCES_KEY)
    assert loaded is not None
    restored = PlottingPreferencesConfig.from_dict(loaded)
    assert resolve_cnv_report_scale(restored.cnv_report_scale) == (
        CNV_REPORT_SCALE_NORMALIZED_DIFFERENCE
    )


def test_resolve_cnv_summary_normalized_uses_admin_preference(tmp_path: Path) -> None:
    store = SecurityStore(db_path=tmp_path / "security.db")
    store.set_gui_setting(
        PLOTTING_PREFERENCES_KEY,
        PlottingPreferencesConfig(
            cnv_report_scale=CNV_REPORT_SCALE_NORMALIZED_DIFFERENCE,
        ).to_dict(),
        updated_by_user_id=None,
    )
    prefs = PlottingPreferencesConfig.from_dict(
        store.get_gui_setting(PLOTTING_PREFERENCES_KEY)
    )
    assert resolve_cnv_summary_normalized(None, plotting_preferences=prefs) is True
    assert resolve_cnv_summary_normalized(False, plotting_preferences=prefs) is False
    assert resolve_cnv_summary_normalized(True, plotting_preferences=prefs) is True


def test_resolve_cnv_summary_normalized_defaults_to_ploidy() -> None:
    assert (
        resolve_cnv_summary_normalized(
            None, plotting_preferences=PlottingPreferencesConfig()
        )
        is False
    )


def test_robin_report_loads_admin_plotting_preferences(
    tmp_path: Path, monkeypatch
) -> None:
    """When plotting_preferences is omitted, RobinReport must load the store default."""
    from robin.reporting.report import RobinReport

    prefs = PlottingPreferencesConfig(
        cnv_report_scale=CNV_REPORT_SCALE_NORMALIZED_DIFFERENCE,
    )
    monkeypatch.setattr(
        "robin.gui.plotting_preferences.load_plotting_preferences",
        lambda store=None: prefs,
    )

    sample_dir = tmp_path / "sample"
    sample_dir.mkdir()
    (sample_dir / "master.csv").write_text(
        "read_id,channel,mux,start_time,duration,passed_filter\n"
    )

    report = RobinReport(
        filename=str(tmp_path / "test_report.pdf"),
        output=str(sample_dir),
        center="test",
        plotting_preferences=None,
    )
    assert report.cnv_summary_normalized is True
    assert (
        report.plotting_preferences.cnv_report_scale
        == CNV_REPORT_SCALE_NORMALIZED_DIFFERENCE
    )


def test_robin_report_respects_explicit_ploidy_override(
    tmp_path: Path, monkeypatch
) -> None:
    from robin.reporting.report import RobinReport

    prefs = PlottingPreferencesConfig(
        cnv_report_scale=CNV_REPORT_SCALE_NORMALIZED_DIFFERENCE,
    )
    monkeypatch.setattr(
        "robin.gui.plotting_preferences.load_plotting_preferences",
        lambda store=None: prefs,
    )

    sample_dir = tmp_path / "sample"
    sample_dir.mkdir()
    (sample_dir / "master.csv").write_text(
        "read_id,channel,mux,start_time,duration,passed_filter\n"
    )

    report = RobinReport(
        filename=str(tmp_path / "test_report.pdf"),
        output=str(sample_dir),
        center="test",
        cnv_summary_normalized=False,
        plotting_preferences=None,
    )
    assert report.cnv_summary_normalized is False


def test_robin_report_empty_config_does_not_load_store(
    tmp_path: Path, monkeypatch
) -> None:
    """Explicit empty PlottingPreferencesConfig must not be replaced by store load."""
    from robin.reporting.report import RobinReport

    called = {"count": 0}

    def _boom(store=None):
        called["count"] += 1
        raise AssertionError("load_plotting_preferences should not be called")

    monkeypatch.setattr(
        "robin.gui.plotting_preferences.load_plotting_preferences",
        _boom,
    )

    sample_dir = tmp_path / "sample"
    sample_dir.mkdir()
    (sample_dir / "master.csv").write_text(
        "read_id,channel,mux,start_time,duration,passed_filter\n"
    )

    report = RobinReport(
        filename=str(tmp_path / "test_report.pdf"),
        output=str(sample_dir),
        center="test",
        plotting_preferences=PlottingPreferencesConfig(),
    )
    assert report.cnv_summary_normalized is False
    assert called["count"] == 0
