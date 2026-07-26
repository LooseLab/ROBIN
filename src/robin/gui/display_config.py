"""Sample page display configuration: registry, persistence model, and visibility."""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any, Dict, FrozenSet, List, Optional, Set

CLASSIFICATION_STEP_IDS = frozenset(
    {"sturgeon", "nanodx", "pannanodx", "random_forest", "marlin", "lamprey", "tucan"}
)

SAMPLE_DISPLAY_KEY = "sample_display"
SAMPLE_DETAILS_SURFACE = "sample_details"
DISPLAY_CONFIG_SCHEMA_VERSION = 3
DISPLAY_ROLES = ("user", "admin")
DISPLAY_ROLE_LABELS = {
    "user": "Standard users",
    "admin": "Administrators",
}
DEFAULT_VIEWER_ROLE = "user"


@dataclass(frozen=True)
class DisplaySection:
    """A block that can be shown or hidden on sample pages (and optionally reports)."""

    id: str
    label: str
    group: str
    workflow_step: Optional[str] = None
    parent_id: Optional[str] = None
    default_visible: bool = True
    surfaces: FrozenSet[str] = frozenset({"sample_page", "report"})


DISPLAY_SECTIONS: Dict[str, DisplaySection] = {
    "sturgeon": DisplaySection(
        "sturgeon", "Sturgeon", "classification", workflow_step="sturgeon"
    ),
    "nanodx": DisplaySection(
        "nanodx", "NanoDX", "classification", workflow_step="nanodx"
    ),
    "pannanodx": DisplaySection(
        "pannanodx", "PanNanoDX", "classification", workflow_step="pannanodx"
    ),
    "random_forest": DisplaySection(
        "random_forest",
        "Random Forest",
        "classification",
        workflow_step="random_forest",
    ),
    "marlin": DisplaySection(
        "marlin",
        "MARLIN",
        "classification",
        workflow_step="marlin",
    ),
    "lamprey": DisplaySection(
        "lamprey",
        "Lamprey (research)",
        "classification",
        workflow_step="lamprey",
    ),
    "tucan": DisplaySection(
        "tucan",
        "Tucan",
        "classification",
        workflow_step="tucan",
    ),
    "target": DisplaySection(
        "target",
        "Coverage / targets (+ IGV & genes on More details)",
        "analysis",
        workflow_step="target",
        surfaces=frozenset({"sample_page", "sample_details", "report"}),
    ),
    "cnv": DisplaySection("cnv", "CNV", "analysis", workflow_step="cnv"),
    "mgmt": DisplaySection("mgmt", "MGMT", "analysis", workflow_step="mgmt"),
    "fusion": DisplaySection(
        "fusion",
        "Fusion analysis (+ fusion pairs on More details)",
        "analysis",
        workflow_step="fusion",
        surfaces=frozenset({"sample_page", "sample_details", "report"}),
    ),
    "itd": DisplaySection(
        "itd",
        "ITDs / insertions",
        "analysis",
        workflow_step="itd",
    ),
    "fusion_target": DisplaySection(
        "fusion_target",
        "Fusion — target panel",
        "analysis",
        workflow_step="fusion",
        parent_id="fusion",
        surfaces=frozenset({"sample_page", "sample_details", "report"}),
    ),
    "fusion_genome": DisplaySection(
        "fusion_genome",
        "Fusion — genome-wide",
        "analysis",
        workflow_step="fusion",
        parent_id="fusion",
        surfaces=frozenset({"sample_page", "sample_details", "report"}),
    ),
    "bed_coverage": DisplaySection(
        "bed_coverage", "BED coverage", "analysis", workflow_step="fusion"
    ),
    "mnpflex": DisplaySection("mnpflex", "MNP-Flex", "v12_classifier"),
    "snp": DisplaySection(
        "snp",
        "SNP analysis",
        "sample_details",
        surfaces=frozenset({"sample_details"}),
    ),
    "output_files": DisplaySection(
        "output_files", "Output files", "other", surfaces=frozenset({"sample_page"})
    ),
}

DISPLAY_GROUP_ORDER = ("classification", "v12_classifier", "analysis", "sample_details", "other")
DISPLAY_GROUP_LABELS = {
    "classification": "Classification",
    "v12_classifier": "V12 Classifier",
    "analysis": "Analysis",
    "sample_details": "More details page",
    "other": "Other",
}


@dataclass
class SampleDisplayConfig:
    """Admin-controlled visibility for sample page sections, per role."""

    schema_version: int = DISPLAY_CONFIG_SCHEMA_VERSION
    sections: Dict[str, bool] = field(default_factory=dict)
    role_sections: Dict[str, Dict[str, bool]] = field(default_factory=dict)
    report_sections: Optional[Dict[str, bool]] = None
    role_report_sections: Optional[Dict[str, Dict[str, bool]]] = None
    updated_at: Optional[str] = None
    updated_by: Optional[str] = None

    def to_dict(self) -> Dict[str, Any]:
        role_sections = {
            role: dict(self.role_sections.get(role) or {})
            for role in DISPLAY_ROLES
        }
        out: Dict[str, Any] = {
            "schema_version": self.schema_version,
            "role_sections": role_sections,
            # Legacy single-map field: standard-user settings for older readers.
            "sections": dict(role_sections.get("user") or self.sections),
        }
        if self.report_sections is not None:
            out["report_sections"] = dict(self.report_sections)
        if self.role_report_sections is not None:
            out["role_report_sections"] = {
                role: dict(overrides)
                for role, overrides in self.role_report_sections.items()
            }
        if self.updated_at:
            out["updated_at"] = self.updated_at
        if self.updated_by:
            out["updated_by"] = self.updated_by
        return out

    @classmethod
    def _migrate_role_sections(
        cls,
        role_sections: Dict[str, Dict[str, bool]],
        schema_version: int,
    ) -> Dict[str, Dict[str, bool]]:
        if schema_version >= DISPLAY_CONFIG_SCHEMA_VERSION:
            return role_sections
        # v3: SNP was incorrectly forced off by admin UI / workflow bugs; restore default.
        migrated = {role: dict(mapping) for role, mapping in role_sections.items()}
        for role in DISPLAY_ROLES:
            mapping = migrated.setdefault(role, {})
            if mapping.get("snp") is False:
                del mapping["snp"]
        return migrated

    @classmethod
    def from_dict(cls, data: Optional[Dict[str, Any]]) -> "SampleDisplayConfig":
        if not data:
            return cls()
        report_sections = data.get("report_sections")
        role_report_raw = data.get("role_report_sections")
        role_sections: Dict[str, Dict[str, bool]] = {}
        raw_role_sections = data.get("role_sections")
        if isinstance(raw_role_sections, dict):
            for role, mapping in raw_role_sections.items():
                if isinstance(mapping, dict):
                    role_sections[str(role)] = {
                        str(k): bool(v) for k, v in mapping.items()
                    }
        legacy_sections = {str(k): bool(v) for k, v in (data.get("sections") or {}).items()}
        if legacy_sections and not role_sections.get("user"):
            role_sections["user"] = dict(legacy_sections)
        if "admin" not in role_sections:
            role_sections["admin"] = {}
        schema_version = int(data.get("schema_version") or 1)
        role_sections = cls._migrate_role_sections(role_sections, schema_version)
        if schema_version < DISPLAY_CONFIG_SCHEMA_VERSION:
            if legacy_sections.get("snp") is False:
                legacy_sections.pop("snp", None)
            user_map = role_sections.get("user")
            if user_map is not None and user_map.get("snp") is False:
                user_map.pop("snp", None)
        sections = dict(role_sections.get("user") or legacy_sections)
        return cls(
            schema_version=DISPLAY_CONFIG_SCHEMA_VERSION,
            sections=sections,
            role_sections=role_sections,
            report_sections=(
                {str(k): bool(v) for k, v in report_sections.items()}
                if isinstance(report_sections, dict)
                else None
            ),
            role_report_sections=(
                {
                    str(role): {str(k): bool(v) for k, v in overrides.items()}
                    for role, overrides in role_report_raw.items()
                    if isinstance(overrides, dict)
                }
                if isinstance(role_report_raw, dict)
                else None
            ),
            updated_at=data.get("updated_at"),
            updated_by=data.get("updated_by"),
        )

    def section_flag(
        self,
        section_id: str,
        *,
        surface: str = "sample_page",
        role: str = DEFAULT_VIEWER_ROLE,
    ) -> Optional[bool]:
        """Explicit admin override for a section and role, if set."""
        role_key = role if role in DISPLAY_ROLES else DEFAULT_VIEWER_ROLE
        role_map = self.role_sections.get(role_key) or {}
        if surface == "report":
            if self.role_report_sections and role_key in self.role_report_sections:
                role_report = self.role_report_sections[role_key]
                if section_id in role_report:
                    return bool(role_report[section_id])
            if self.report_sections is not None and section_id in self.report_sections:
                return bool(self.report_sections[section_id])
        if section_id in role_map:
            return bool(role_map[section_id])
        if role_key == DEFAULT_VIEWER_ROLE and section_id in self.sections:
            return bool(self.sections[section_id])
        return None

    def with_role_updates(
        self,
        role: str,
        sections: Dict[str, bool],
        *,
        updated_at: Optional[str] = None,
        updated_by: Optional[str] = None,
    ) -> "SampleDisplayConfig":
        role_key = role if role in DISPLAY_ROLES else DEFAULT_VIEWER_ROLE
        merged_roles = {
            r: dict(self.role_sections.get(r) or {}) for r in DISPLAY_ROLES
        }
        role_map = dict(merged_roles.get(role_key) or {})
        role_map.update(sections)
        merged_roles[role_key] = role_map
        legacy_sections = (
            dict(merged_roles["user"])
            if role_key == "user"
            else dict(self.sections)
        )
        return SampleDisplayConfig(
            schema_version=self.schema_version,
            sections=legacy_sections,
            role_sections=merged_roles,
            report_sections=self.report_sections,
            role_report_sections=self.role_report_sections,
            updated_at=updated_at or self.updated_at,
            updated_by=updated_by or self.updated_by,
        )

    def with_updates(
        self,
        sections: Dict[str, bool],
        *,
        updated_at: Optional[str] = None,
        updated_by: Optional[str] = None,
    ) -> "SampleDisplayConfig":
        return self.with_role_updates(
            DEFAULT_VIEWER_ROLE,
            sections,
            updated_at=updated_at,
            updated_by=updated_by,
        )


def _normalize_workflow_steps(workflow_steps: Optional[List[str]]) -> Set[str]:
    if not workflow_steps:
        return set()
    out: Set[str] = set()
    for step in workflow_steps:
        name = step.split(":")[-1] if ":" in step else step
        out.add(name)
    return out


def _workflow_allows_section(
    section: DisplaySection, workflow_steps: Optional[List[str]]
) -> bool:
    if not workflow_steps:
        return True
    if section.workflow_step is None:
        return True
    enabled = _normalize_workflow_steps(workflow_steps)
    return section.workflow_step in enabled


def is_section_visible(
    section_id: str,
    *,
    workflow_steps: Optional[List[str]] = None,
    display_config: Optional[SampleDisplayConfig] = None,
    surface: str = "sample_page",
    viewer_role: str = DEFAULT_VIEWER_ROLE,
) -> bool:
    """Return True when a section should appear on the given surface."""
    section = DISPLAY_SECTIONS.get(section_id)
    if section is None:
        return True
    if surface not in section.surfaces:
        return False

    if section.parent_id:
        if not is_section_visible(
            section.parent_id,
            workflow_steps=workflow_steps,
            display_config=display_config,
            surface=surface,
            viewer_role=viewer_role,
        ):
            return False

    if not _workflow_allows_section(section, workflow_steps):
        return False

    if display_config is not None:
        override = display_config.section_flag(
            section_id, surface=surface, role=viewer_role
        )
        if override is not None:
            return override

    return section.default_visible


def get_visible_classification_steps(
    workflow_steps: Optional[List[str]] = None,
    display_config: Optional[SampleDisplayConfig] = None,
    *,
    surface: str = "sample_page",
    viewer_role: str = DEFAULT_VIEWER_ROLE,
) -> Set[str]:
    return {
        step
        for step in CLASSIFICATION_STEP_IDS
        if is_section_visible(
            step,
            workflow_steps=workflow_steps,
            display_config=display_config,
            surface=surface,
            viewer_role=viewer_role,
        )
    }


def any_classification_visible(
    workflow_steps: Optional[List[str]] = None,
    display_config: Optional[SampleDisplayConfig] = None,
    *,
    surface: str = "sample_page",
    viewer_role: str = DEFAULT_VIEWER_ROLE,
) -> bool:
    return bool(
        get_visible_classification_steps(
            workflow_steps,
            display_config,
            surface=surface,
            viewer_role=viewer_role,
        )
    )


def any_sample_details_visible(
    workflow_steps: Optional[List[str]] = None,
    display_config: Optional[SampleDisplayConfig] = None,
    *,
    viewer_role: str = DEFAULT_VIEWER_ROLE,
) -> bool:
    """True when the More details page would show at least one analysis block."""
    if is_section_visible(
        "target",
        workflow_steps=workflow_steps,
        display_config=display_config,
        surface=SAMPLE_DETAILS_SURFACE,
        viewer_role=viewer_role,
    ):
        return True
    if is_section_visible(
        "snp",
        workflow_steps=workflow_steps,
        display_config=display_config,
        surface=SAMPLE_DETAILS_SURFACE,
        viewer_role=viewer_role,
    ):
        return True
    if is_section_visible(
        "fusion",
        workflow_steps=workflow_steps,
        display_config=display_config,
        surface=SAMPLE_DETAILS_SURFACE,
        viewer_role=viewer_role,
    ) and (
        is_section_visible(
            "fusion_target",
            workflow_steps=workflow_steps,
            display_config=display_config,
            surface=SAMPLE_DETAILS_SURFACE,
            viewer_role=viewer_role,
        )
        or is_section_visible(
            "fusion_genome",
            workflow_steps=workflow_steps,
            display_config=display_config,
            surface=SAMPLE_DETAILS_SURFACE,
            viewer_role=viewer_role,
        )
    ):
        return True
    return False


def config_from_workflow_steps(
    workflow_steps: Optional[List[str]],
    *,
    role: str = DEFAULT_VIEWER_ROLE,
) -> SampleDisplayConfig:
    """Build a display config that mirrors the active workflow steps for one role."""
    sections: Dict[str, bool] = {}
    for section_id, section in DISPLAY_SECTIONS.items():
        if section.workflow_step is None:
            sections[section_id] = section.default_visible
        else:
            sections[section_id] = _workflow_allows_section(section, workflow_steps)
    base = SampleDisplayConfig()
    return base.with_role_updates(role, sections)


def section_admin_surface(section: DisplaySection) -> str:
    """Pick the surface used to resolve admin checkbox state for a section."""
    if "sample_page" in section.surfaces:
        return "sample_page"
    if SAMPLE_DETAILS_SURFACE in section.surfaces:
        return SAMPLE_DETAILS_SURFACE
    if "report" in section.surfaces:
        return "report"
    return "sample_page"


def effective_section_map(
    workflow_steps: Optional[List[str]] = None,
    display_config: Optional[SampleDisplayConfig] = None,
    *,
    surface: str = "sample_page",
    viewer_role: str = DEFAULT_VIEWER_ROLE,
) -> Dict[str, bool]:
    """Resolve visibility for every registered section (for admin UI)."""
    return {
        section_id: is_section_visible(
            section_id,
            workflow_steps=workflow_steps,
            display_config=display_config,
            surface=section_admin_surface(DISPLAY_SECTIONS[section_id]),
            viewer_role=viewer_role,
        )
        for section_id in DISPLAY_SECTIONS
    }


def sections_for_group(group: str) -> List[DisplaySection]:
    return [
        s
        for s in DISPLAY_SECTIONS.values()
        if s.group == group and s.parent_id is None
    ] + [
        s
        for s in DISPLAY_SECTIONS.values()
        if s.group == group and s.parent_id is not None
    ]

