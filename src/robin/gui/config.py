"""
Configuration and helper functions for GUI section visibility.
"""

from typing import TYPE_CHECKING, List, Optional, Set

if TYPE_CHECKING:
    from robin.gui.display_config import SampleDisplayConfig

# Re-export get_confidence_level from classification_config for backward compatibility
try:
    from robin.classification_config import get_confidence_level
except ImportError:
    # Fallback if classification_config is not available
    def get_confidence_level(classifier: str, confidence: float) -> str:
        """Fallback confidence level function."""
        if confidence >= 80:
            return "High confidence"
        elif confidence >= 50:
            return "Medium confidence"
        elif confidence >= 20:
            return "Low confidence"
        else:
            return "Very low confidence"


# Map workflow step names to GUI section names
WORKFLOW_STEP_TO_SECTION = {
    "target": "target",
    "mgmt": "mgmt",
    "fusion": "fusion",
    "itd": "itd",
    "cnv": "cnv",
    "sturgeon": "sturgeon",
    "nanodx": "nanodx",
    "random_forest": "random_forest",
    "pannanodx": "pannanodx",
    "marlin": "marlin",
    "lamprey": "lamprey",
    "tucan": "tucan",
}

# Map workflow steps to their display names in classification section
CLASSIFICATION_STEPS = {
    "sturgeon": "Sturgeon",
    "nanodx": "NanoDX",
    "random_forest": "Random Forest",
    "pannanodx": "PanNanoDX",
    "marlin": "MARLIN",
    "lamprey": "Lamprey (research)",
    "tucan": "Tucan",
}


def get_enabled_sections(workflow_steps: Optional[List[str]]) -> Set[str]:
    """
    Determine which sections should be enabled based on workflow steps.

    Args:
        workflow_steps: List of workflow step names (e.g., ['target', 'mgmt', 'sturgeon'])

    Returns:
        Set of enabled section names
    """
    if not workflow_steps:
        # If no workflow steps specified, show all sections (backward compatibility)
        return set(WORKFLOW_STEP_TO_SECTION.values())

    enabled = set()
    for step in workflow_steps:
        # Handle workflow steps that might have queue prefixes (e.g., "classification:sturgeon")
        step_name = step.split(":")[-1] if ":" in step else step
        if step_name in WORKFLOW_STEP_TO_SECTION:
            enabled.add(WORKFLOW_STEP_TO_SECTION[step_name])

    return enabled


def is_section_enabled(section_name: str, workflow_steps: Optional[List[str]]) -> bool:
    """
    Check if a specific section should be enabled.

    Args:
        section_name: Name of the section to check
        workflow_steps: List of workflow step names

    Returns:
        True if section should be enabled, False otherwise
    """
    enabled_sections = get_enabled_sections(workflow_steps)

    # If no workflow steps specified, show all sections (backward compatibility)
    if not workflow_steps:
        return True

    return section_name in enabled_sections


def get_enabled_classification_steps(workflow_steps: Optional[List[str]]) -> Set[str]:
    """
    Get the set of enabled classification steps.

    Args:
        workflow_steps: List of workflow step names

    Returns:
        Set of enabled classification step names (e.g., {'sturgeon', 'nanodx'})
    """
    enabled_sections = get_enabled_sections(workflow_steps)
    return enabled_sections.intersection(set(CLASSIFICATION_STEPS.keys()))


def is_section_visible(
    section_id: str,
    *,
    workflow_steps: Optional[List[str]] = None,
    display_config: Optional["SampleDisplayConfig"] = None,
    surface: str = "sample_page",
    viewer_role: Optional[str] = None,
) -> bool:
    """Check workflow and admin display config for section visibility."""
    from robin.gui.display_config import (
        DEFAULT_VIEWER_ROLE,
    )
    from robin.gui.display_config import is_section_visible as _resolve

    return _resolve(
        section_id,
        workflow_steps=workflow_steps,
        display_config=display_config,
        surface=surface,
        viewer_role=viewer_role or DEFAULT_VIEWER_ROLE,
    )


def get_visible_classification_steps(
    workflow_steps: Optional[List[str]] = None,
    display_config: Optional["SampleDisplayConfig"] = None,
    *,
    surface: str = "sample_page",
    viewer_role: Optional[str] = None,
) -> Set[str]:
    """Classification steps visible on the sample page or in reports."""
    from robin.gui.display_config import (
        DEFAULT_VIEWER_ROLE,
    )
    from robin.gui.display_config import get_visible_classification_steps as _resolve

    return _resolve(
        workflow_steps,
        display_config,
        surface=surface,
        viewer_role=viewer_role or DEFAULT_VIEWER_ROLE,
    )


def any_classification_visible(
    workflow_steps: Optional[List[str]] = None,
    display_config: Optional["SampleDisplayConfig"] = None,
    *,
    surface: str = "sample_page",
    viewer_role: Optional[str] = None,
) -> bool:
    from robin.gui.display_config import (
        DEFAULT_VIEWER_ROLE,
    )
    from robin.gui.display_config import any_classification_visible as _resolve

    return _resolve(
        workflow_steps,
        display_config,
        surface=surface,
        viewer_role=viewer_role or DEFAULT_VIEWER_ROLE,
    )


def resolve_viewer_role(launcher: object) -> str:
    """Map the signed-in GUI user to a display-config role key."""
    from robin.gui.display_config import DEFAULT_VIEWER_ROLE

    if launcher is None:
        return DEFAULT_VIEWER_ROLE
    get_user_id = getattr(launcher, "_get_current_user_id", None)
    store = getattr(launcher, "security_store", None)
    if not callable(get_user_id) or store is None:
        return DEFAULT_VIEWER_ROLE
    user_id = get_user_id()
    if user_id is not None and store.user_has_role(int(user_id), "admin"):
        return "admin"
    return "user"


def launcher_visibility_context(launcher: object, *, surface: str = "sample_page"):
    """Return (workflow_steps, display_config, viewer_role) from a launcher."""
    workflow_steps = (
        launcher.workflow_steps
        if launcher is not None and hasattr(launcher, "workflow_steps")
        else None
    )
    display_config = (
        launcher.display_config
        if launcher is not None and hasattr(launcher, "display_config")
        else None
    )
    return workflow_steps, display_config, resolve_viewer_role(launcher)


def any_sample_details_visible_for_launcher(launcher: object) -> bool:
    from robin.gui.display_config import any_sample_details_visible

    workflow_steps, display_config, viewer_role = launcher_visibility_context(launcher)
    return any_sample_details_visible(
        workflow_steps, display_config, viewer_role=viewer_role
    )
