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
from robin.security import SecurityStore


def test_plotting_preferences_defaults() -> None:
    config = PlottingPreferencesConfig()
    assert config.cnv_report_scale == CNV_REPORT_SCALE_PLOIDY


def test_plotting_preferences_legacy_log2_alias() -> None:
    config = PlottingPreferencesConfig.from_dict(
        {"cnv_report_scale": "log2_ratio"}
    )
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
    assert resolve_cnv_summary_normalized(None, plotting_preferences=PlottingPreferencesConfig()) is False
