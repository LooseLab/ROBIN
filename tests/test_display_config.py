from __future__ import annotations

from pathlib import Path

from robin.gui.display_config import (
    SAMPLE_DISPLAY_KEY,
    SampleDisplayConfig,
    any_sample_details_visible,
    config_from_workflow_steps,
    get_visible_classification_steps,
    is_section_visible,
)
from robin.security import SecurityStore


def test_display_config_hides_admin_override() -> None:
    config = SampleDisplayConfig(sections={"sturgeon": False})
    assert is_section_visible(
        "sturgeon",
        workflow_steps=["sturgeon", "cnv"],
        display_config=config,
    ) is False
    assert is_section_visible(
        "cnv",
        workflow_steps=["sturgeon", "cnv"],
        display_config=config,
    ) is True


def test_display_config_respects_workflow_steps() -> None:
    assert is_section_visible("cnv", workflow_steps=["mgmt"]) is False
    assert is_section_visible("cnv", workflow_steps=["cnv", "mgmt"]) is True


def test_fusion_child_hidden_when_parent_hidden() -> None:
    config = SampleDisplayConfig(sections={"fusion": False})
    assert is_section_visible(
        "fusion_target",
        workflow_steps=["fusion"],
        display_config=config,
    ) is False


def test_classification_visible_set() -> None:
    config = SampleDisplayConfig(
        sections={"sturgeon": True, "nanodx": False, "pannanodx": False}
    )
    visible = get_visible_classification_steps(
        ["sturgeon", "nanodx", "pannanodx", "random_forest"],
        config,
    )
    assert visible == {"sturgeon", "random_forest"}


def test_config_from_workflow_steps() -> None:
    config = config_from_workflow_steps(["cnv", "sturgeon"], role="user")
    assert config.role_sections["user"]["cnv"] is True
    assert config.role_sections["user"]["mgmt"] is False
    assert config.role_sections["user"]["sturgeon"] is True
    assert config.role_sections["user"]["nanodx"] is False


def test_role_sections_use_viewer_role() -> None:
    config = SampleDisplayConfig(
        role_sections={
            "user": {"cnv": False, "sturgeon": True},
            "admin": {"cnv": True},
        }
    )
    assert (
        is_section_visible(
            "cnv",
            workflow_steps=["cnv"],
            display_config=config,
            viewer_role="user",
        )
        is False
    )
    assert (
        is_section_visible(
            "cnv",
            workflow_steps=["cnv"],
            display_config=config,
            viewer_role="admin",
        )
        is True
    )


def test_legacy_sections_migrate_to_user_role() -> None:
    config = SampleDisplayConfig.from_dict(
        {"schema_version": 1, "sections": {"nanodx": False}}
    )
    assert config.role_sections["user"]["nanodx"] is False
    assert config.role_sections["admin"] == {}


def test_any_sample_details_visible() -> None:
    config = SampleDisplayConfig(sections={"target": False, "snp": False, "fusion": False})
    assert any_sample_details_visible(["target", "fusion", "snp_analysis"], config) is False

    config2 = SampleDisplayConfig(sections={"fusion": True, "fusion_target": True})
    assert any_sample_details_visible(["fusion"], config2) is True

    assert is_section_visible(
        "snp",
        workflow_steps=["snp_analysis"],
        surface="sample_details",
    ) is True
    assert is_section_visible(
        "snp",
        workflow_steps=["cnv"],
        surface="sample_details",
    ) is False


def test_security_store_gui_settings_roundtrip(tmp_path: Path) -> None:
    store = SecurityStore(db_path=tmp_path / "security.db")
    payload = SampleDisplayConfig(sections={"cnv": False}).to_dict()
    store.set_gui_setting(SAMPLE_DISPLAY_KEY, payload, updated_by_user_id=None)
    loaded = store.get_gui_setting(SAMPLE_DISPLAY_KEY)
    assert loaded is not None
    assert loaded["sections"]["cnv"] is False
    restored = SampleDisplayConfig.from_dict(loaded)
    assert is_section_visible(
        "cnv",
        workflow_steps=["cnv"],
        display_config=restored,
    ) is False
