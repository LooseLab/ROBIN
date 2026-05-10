"""
Module: theme

This module defines the theme and layout for the entire application, following the **Editorial Bioinformatics**
design system (see ``design.md``): clear hierarchy, slate neutrals, emerald accents, and semantic visualization colors.

It includes:

- A context manager `frame` to create a custom page frame with navigation, header, and footer.
- Utility functions to handle dark mode and remote access toggling.
- Material Design 3 color tokens, typography scale, and spacing system.
- Modern component styling with proper elevation and surface treatments.
- Responsive design patterns and accessibility improvements.

Functions:

- frame(navtitle: str): Context manager for creating a consistent page layout with header, footer, and navigation.
- cleanup_and_exit(): Handles cleanup operations before shutting down the application.
- dark_mode(event: events.ValueChangeEventArguments): Toggles dark mode based on the event argument.
- use_on_air(args: events.ValueChangeEventArguments): Toggles remote access based on the event argument.

Constants:

- IMAGEFILE: Path to the image file used in the header and footer.
- HEADER_HTML: HTML content for the header.
- STYLE_CSS, M3_COMPONENTS_CSS, MOSAIC_COMPONENTS_CSS: CSS bundles for the app shell.
- EDITORIAL_FONTS_HTML: Google Fonts link block (Manrope, Inter, JetBrains Mono).

External Dependencies:

- contextlib.contextmanager
- nicegui (ui, app, events, core, air)
- pathlib.Path
- robin.images
- os
- psutil
- platform
"""

from contextlib import contextmanager
from signal import siginterrupt
from packaging import version
import requests
import asyncio
import logging
import subprocess
import importlib.metadata
import weakref
import time
from typing import Callable, Optional, Any, Dict, List


from nicegui import ui, app, events, core, run, background_tasks
import nicegui.air

from pathlib import Path

# These will be set by the get_imagefile() and get_version() functions
IMAGEFILE = None
__about__ = None


import os
import psutil
import platform

# Check if we're in development mode
is_development_mode = os.environ.get("ROBIN_DEV_MODE", "").lower() in ("1", "true", "yes", "on")

# Process large BAMs individually (do not use alongside live runs)
_process_large_bams_enabled = os.environ.get("ROBIN_PROCESS_LARGE_BAMS", "0").strip().lower() in ("1", "true", "yes", "on")

# Single global theme sync timer state
_theme_sync_callbacks: Dict[int, Dict[str, Any]] = {}
_theme_sync_next_id = 0
_theme_sync_timer_started = False
_THEME_SYNC_TICK_SECONDS = 0.5


def ui_element_exists(element: Any) -> bool:
    """Best-effort check that a NiceGUI element still exists."""
    if element is None:
        return True
    try:
        if getattr(element, "_deleted", False):
            return False
        if getattr(element, "client", None) is None:
            return False
        # Access id to force resolution for stale/disconnected elements.
        _ = element.id
        return True
    except Exception:
        return False


def _run_theme_sync_callbacks() -> None:
    """Run registered theme sync callbacks from a single global timer."""
    now = time.monotonic()
    dead_ids: List[int] = []
    for callback_id, item in list(_theme_sync_callbacks.items()):
        element_ref = item.get("element_ref")
        interval_s = float(item.get("interval_s", 1.0))
        last_run = float(item.get("last_run", 0.0))

        if now - last_run < interval_s:
            continue

        target_element = element_ref() if element_ref else None
        if element_ref and not ui_element_exists(target_element):
            dead_ids.append(callback_id)
            continue

        cb = item.get("callback")
        if cb is None:
            dead_ids.append(callback_id)
            continue

        try:
            cb()
            item["last_run"] = now
        except Exception:
            # Keep UI responsive; callback owners can re-register if needed.
            pass

    for callback_id in dead_ids:
        _theme_sync_callbacks.pop(callback_id, None)


def ensure_theme_sync_timer() -> None:
    """Ensure the single global theme sync timer is started."""
    global _theme_sync_timer_started
    if _theme_sync_timer_started:
        return
    ui.timer(_THEME_SYNC_TICK_SECONDS, _run_theme_sync_callbacks, active=True)
    _theme_sync_timer_started = True


def register_theme_sync_callback(
    callback: Callable[[], None],
    *,
    element: Optional[Any] = None,
    interval_s: float = 1.0,
    immediate: bool = False,
) -> Callable[[], None]:
    """
    Register a callback to run from the single global theme timer.

    If `element` is provided, callback execution stops automatically once that
    element no longer exists.
    """
    global _theme_sync_next_id
    ensure_theme_sync_timer()
    _theme_sync_next_id += 1
    callback_id = _theme_sync_next_id
    _theme_sync_callbacks[callback_id] = {
        "callback": callback,
        "interval_s": max(0.1, float(interval_s)),
        "last_run": 0.0,
        "element_ref": weakref.ref(element) if element is not None else None,
    }

    if immediate:
        try:
            callback()
        except Exception:
            pass

    def _unregister() -> None:
        _theme_sync_callbacks.pop(callback_id, None)

    return _unregister


def get_imagefile():
    """Get the path to the ROBIN logo image file."""
    global IMAGEFILE
    if IMAGEFILE is None:
        try:
            from robin.gui import images
            IMAGEFILE = os.path.join(
                os.path.dirname(os.path.abspath(images.__file__)), "ROBIN_logo_small.png"
            )
        except (ImportError, AttributeError):
            # Fallback path when running standalone
            IMAGEFILE = os.path.join(
                os.path.dirname(os.path.abspath(__file__)), "images", "ROBIN_logo_small.png"
            )
    return IMAGEFILE


def get_about():
    """Get the __about__ module with version information."""
    global __about__
    if __about__ is None:
        try:
            from robin import __about__
        except ImportError:
            # Fallback when running standalone - create a minimal __about__ object
            class MockAbout:
                __version__ = "standalone-test"
            __about__ = MockAbout()
    return __about__


def get_version_from_github():
    response = requests.get(
        "https://raw.githubusercontent.com/LooseLab/ROBIN/main/src/robin/__about__.py"
    )
    response.raise_for_status()
    remote_version_str = None
    for line in response.text.split("\n"):
        if line.startswith("__version__"):
            remote_version_str = line.split("=")[1].strip().strip('"').strip("'")
            break
    return remote_version_str


def styled_table(*, columns, rows=None, pagination=20, class_size="table-xs", **kwargs):
    """Create a NiceGUI table with Material Design 3 styling and compact layout.

    Args:
        columns: columns definition passed to ui.table
        rows: initial rows
        pagination: rows per page (0 or None disables pagination)
        class_size: table size class (e.g., "table-xs", "table-sm")
        **kwargs: forwarded to ui.table

    Returns:
        Tuple of (container, table) where container is the overflow wrapper column and table is the ui.table instance.
    """
    # Add CSS to hide pagination for tables with no-pagination class (only once)
    if not hasattr(styled_table, '_pagination_css_added'):
        ui.add_head_html("""
            <style>
                .no-pagination .q-table__bottom {
                    display: none !important;
                }
            </style>
        """)
        styled_table._pagination_css_added = True
    
    # Outer container with M3 styling and mobile touch scrolling support
    container_classes = "w-full overflow-x-auto compact-table elevation-1 rounded-lg"
    if pagination == 0 or pagination is None:
        container_classes += " no-pagination"
    
    container = ui.column().classes(container_classes)
    with container:
        # If pagination is 0 or None, disable pagination completely
        if pagination == 0 or pagination is None:
            table = ui.table(
                columns=columns, rows=rows or [], pagination=None, **kwargs
            )
        else:
            table = ui.table(
                columns=columns, rows=rows or [], pagination=pagination, **kwargs
            )
        try:
            table.classes(replace=f"table w-full {class_size} text-xs")
        except Exception:
            table.classes(f"table w-full {class_size} text-xs")
        # Use Quasar's dense mode with M3 styling for maximum compactness
        try:
            if pagination == 0 or pagination is None:
                table.props("dense flat wrap-cells hide-pagination")
            else:
                table.props("dense flat wrap-cells")
        except Exception:
            pass
    return container, table


def unpack_qtable_request_pagination(event_args: Any) -> Optional[Dict[str, Any]]:
    """Extract the pagination dict from a NiceGUI Quasar QTable ``request`` event payload."""
    if event_args is None:
        return None
    if isinstance(event_args, list):
        if not event_args:
            return None
        event_args = event_args[0]
    if not isinstance(event_args, dict):
        return None
    pag = event_args.get("pagination")
    if isinstance(pag, dict):
        return pag
    return event_args


def clamp_qtable_server_pagination(
    pagination: Dict[str, Any],
    *,
    rows_number: int,
    rows_per_page_default: int = 100,
) -> Dict[str, Any]:
    """Normalize ``page`` / ``rowsPerPage`` / ``rowsNumber`` for server-side QTable paging.

    Treats ``rowsPerPage <= 0`` (Quasar \"All\") as ``max(1, rows_number)`` so we never
    materialize zero pages or rely on client-side \"all rows\" behavior.
    """
    try:
        rpp = int(pagination.get("rowsPerPage"))
    except (TypeError, ValueError):
        rpp = rows_per_page_default
    if rpp <= 0:
        rpp = max(1, int(rows_number))

    total = max(0, int(rows_number))
    max_page = max(1, (total + rpp - 1) // rpp) if total else 1
    try:
        page = int(pagination.get("page") or 1)
    except (TypeError, ValueError):
        page = 1
    page = max(1, min(page, max_page))

    out = dict(pagination)
    out["page"] = page
    out["rowsPerPage"] = rpp
    out["rowsNumber"] = total
    out.setdefault("sortBy", None)
    out.setdefault("descending", False)
    return out


def wire_qtable_server_pagination_handlers(
    table: Any,
    refill: Callable[[Dict[str, Any]], None],
) -> None:
    """Drive server-side row slicing from Quasar footer controls.

    Changing **page** or **records per page** updates the ``pagination`` model via
    NiceGUI's ``update:pagination`` event; some QTable setups also emit ``request``.
    Listening to **both** keeps rows-per-page and page navigation in sync with Python.
    """
    def _on_request(e: Any) -> None:
        pag = unpack_qtable_request_pagination(getattr(e, "args", None))
        if pag is not None:
            refill(pag)

    def _on_pagination_change(e: Any) -> None:
        val = getattr(e, "value", None)
        if isinstance(val, dict):
            refill(val)

    table.on("request", _on_request)
    table.on_pagination_change(_on_pagination_change)


def styled_server_paged_table(
    *,
    columns: List[Dict[str, Any]],
    rows=None,
    pagination: Dict[str, Any],
    class_size: str = "table-xs",
    row_key: str = "__row_idx",
    rows_per_page_options: Optional[List[int]] = None,
    **kwargs: Any,
):
    """Styled QTable with footer pagination for server-driven row slicing.

    After creating the table, call :func:`wire_qtable_server_pagination_handlers`
    so **page** and **rows-per-page** changes refetch rows (``request`` alone is not
    always emitted). Use ``row_key`` that exists on every rendered row (often ``__row_idx``).
    """
    import json

    container_classes = "w-full overflow-x-auto compact-table elevation-1 rounded-lg"
    container = ui.column().classes(container_classes)
    with container:
        table = ui.table(
            columns=columns,
            rows=rows or [],
            row_key=row_key,
            pagination=pagination,
            **kwargs,
        )
        try:
            table.classes(replace=f"table w-full {class_size} text-xs")
        except Exception:
            table.classes(f"table w-full {class_size} text-xs")
        opts = rows_per_page_options or [25, 50, 100, 250]
        opt_str = json.dumps(opts, separators=(",", ":"))
        table.props(f'dense flat wrap-cells rows-per-page-options="{opt_str}"')
    return container, table


async def check_version():
    """
    Check the current version against the remote version on GitHub.
    Shows a notification or dialog to the user about their version status.
    """
    # Skip version check in development mode
    if is_development_mode:
        return
        
    # Check if version has already been checked in this app session
    try:
        if app.storage.general.get("version_checked", False):
            return
    except RuntimeError:
        # Storage not available in this context, continue with version check
        pass

    try:
        remote_version_str = await run.io_bound(get_version_from_github)

        if not remote_version_str:
            with ui.dialog() as dialog, ui.card().classes(
                "robin-dialog-surface p-4 md:p-5 min-w-[18rem] max-w-md"
            ):
                ui.label("Version check").classes(
                    "classification-insight-heading text-headline-small q-mb-sm"
                )
                ui.label(
                    "Could not determine remote version. Please check manually."
                ).classes("classification-insight-foot q-mb-md")
                ui.button("OK", on_click=dialog.close).props("color=primary no-caps")
            dialog.open()
            return

        local_version = version.parse(get_about().__version__)
        remote_version = version.parse(remote_version_str)

        if local_version == remote_version:
            ui.notify("Your ROBIN installation is up to date!", type="positive")
        elif local_version < remote_version:
            with ui.dialog() as dialog, ui.card().classes(
                "robin-dialog-surface p-4 md:p-5 min-w-[18rem] max-w-md"
            ):
                ui.label("Update available").classes(
                    "classification-insight-heading text-headline-small q-mb-sm"
                )
                ui.label(f"Your version: {local_version}").classes(
                    "classification-insight-foot"
                )
                ui.label(f"Latest version: {remote_version}").classes(
                    "classification-insight-foot"
                )
                ui.label(
                    "Would you like to visit the GitHub repository to update?"
                ).classes("classification-insight-foot q-mb-md")
                with ui.row().classes("w-full justify-center gap-2 flex-wrap"):
                    ui.button(
                        "Continue with current version",
                        on_click=dialog.close,
                    ).props("flat no-caps outline")
                    ui.button(
                        "Visit GitHub",
                        on_click=lambda: ui.open("https://github.com/LooseLab/ROBIN"),
                        icon="open_in_new",
                    ).props("color=primary no-caps")
            dialog.open()
        else:
            with ui.dialog() as dialog, ui.card().classes(
                "robin-dialog-surface p-4 md:p-5 min-w-[18rem] max-w-md"
            ):
                ui.label("Development version").classes(
                    "classification-insight-heading text-headline-small q-mb-sm"
                )
                ui.label(
                    f"You are running a development version ({local_version})."
                ).classes("classification-insight-foot")
                ui.label(f"Latest release: {remote_version}").classes(
                    "classification-insight-foot"
                )
                ui.label(
                    "This version may be unstable and is only for testing. "
                    "It is not recommended for production use."
                ).classes("classification-insight-foot q-mb-md")
                ui.label(
                    "Please consider using the latest release instead."
                ).classes("classification-insight-foot q-mb-md")
                ui.button("OK", on_click=dialog.close).props("color=primary no-caps")
            dialog.open()

    except requests.RequestException:
        with ui.dialog() as dialog, ui.card().classes(
            "robin-dialog-surface p-4 md:p-5 min-w-[18rem] max-w-md"
        ):
            ui.label("Connection error").classes(
                "classification-insight-heading text-headline-small q-mb-sm"
            )
            ui.label("Could not check for updates.").classes("classification-insight-foot")
            ui.label(
                "Either you are not connected to the internet or you cannot access "
                "https://www.github.com/looselab/robin."
            ).classes("classification-insight-foot q-mb-md")
            ui.label("Please manually check for updates.").classes(
                "classification-insight-foot q-mb-md"
            )
            ui.button("OK", on_click=dialog.close).props("color=primary no-caps")
        dialog.open()
    except Exception as e:
        with ui.dialog() as dialog, ui.card().classes(
            "robin-dialog-surface p-4 md:p-5 min-w-[18rem] max-w-md"
        ):
            ui.label("Error").classes(
                "classification-insight-heading text-headline-small q-mb-sm"
            )
            ui.label(f"Error checking version: {str(e)}").classes(
                "classification-insight-foot"
            )
            ui.button("OK", on_click=dialog.close).props("color=primary no-caps")
        dialog.open()

    # Mark version as checked for this app session
    try:
        app.storage.general["version_checked"] = True
    except RuntimeError:
        # Storage not available in this context, skip setting the flag
        pass


# IMAGEFILE is now defined above with fallback handling

# Module-level variables
quitdialog = None
_remote_access_initialized = False
_remote_access_target = False
_remote_access_task = None
_remote_access_request_id = 0

MENU_BREAKPOINT = 1200


def _is_local_client() -> bool:
    """Return True if the current request is from localhost or 127.0.0.1 (same device as the server)."""
    try:
        # Host is set by auth middleware on each request
        host = (app.storage.user.get("_request_host") or "").strip().lower()
        if not host:
            return False
        # Strip port if present (e.g. "localhost:8080" -> "localhost")
        host_part = host.split(":")[0].lower()
        return host_part in ("localhost", "127.0.0.1")
    except Exception:
        return False


class GlobalSystemMetrics:
    """Global system metrics singleton that provides CPU and RAM usage data."""
    
    _instance = None
    _timer_active = False
    
    def __new__(cls):
        if cls._instance is None:
            cls._instance = super().__new__(cls)
            cls._instance.cpu = 0
            cls._instance.ram = 0
        return cls._instance
    
    def start_timer(self):
        """Start the global metrics timer if not already running."""
        if not self._timer_active:
            self._timer_active = True
            ui.timer(1.0, self.update_metrics)
    
    def update_metrics(self):
        """Update CPU and RAM metrics."""
        try:
            self.cpu = round(psutil.getloadavg()[1] / os.cpu_count() * 100, 1)
            self.ram = round(psutil.virtual_memory()[2], 1)
        except (OSError, AttributeError):
            # Handle cases where psutil might fail or os.cpu_count() returns None
            self.cpu = 0
            self.ram = 0


# Global metrics instance
global_metrics = GlobalSystemMetrics()

# Read the HTML content for the header
HEADER_HTML = (Path(__file__).parent / "static" / "header.html").read_text()

# Read the CSS styles for the application
STYLE_CSS = (Path(__file__).parent / "static" / "styles.css").read_text()
M3_COMPONENTS_CSS = (Path(__file__).parent / "static" / "m3-components.css").read_text()
MOSAIC_COMPONENTS_CSS = (Path(__file__).parent / "static" / "mosaic-components.css").read_text()

# Google Fonts: Manrope (headlines), Inter (UI/data labels), JetBrains Mono (technical strings)
EDITORIAL_FONTS_HTML = """
<link rel="preconnect" href="https://fonts.googleapis.com">
<link rel="preconnect" href="https://fonts.gstatic.com" crossorigin>
<link href="https://fonts.googleapis.com/css2?family=Inter:wght@400;500;600;700&family=JetBrains+Mono:wght@400;500;600&family=Manrope:wght@600;700;800&display=swap" rel="stylesheet">
"""


@contextmanager
def frame(navtitle: str, batphone=False, smalltitle=None, center: str = None, setup_notifications=None):
    """
    Context manager to create a custom page frame with Material Design 3 styling and consistent behavior across all pages.

    Args:
        navtitle (str): The title to display in the navigation header.
        batphone (bool): Whether to show the BATMAN mode title.
        smalltitle (str): The title to display on small screens.
        center (str): Center ID running the analysis.
        setup_notifications (callable): Optional callback to set up notification system with the notification container.

    Yields:
        None
    """
    global quitdialog
    if batphone:
        navtitle = f"BATMAN & {navtitle}"

    # Initialize remote access policy once at app launch (global runtime state).
    initialize_remote_access_once()

    # Store center in app storage if provided
    if center:
        app.storage.general["center"] = center

    # Add custom HTML and CSS to the head of the page
    ui.add_head_html(
        '<script src="https://cdn.jsdelivr.net/npm/igv@3.7.0/dist/igv.min.js"></script>'
    )
    ui.add_head_html(EDITORIAL_FONTS_HTML)
    ui.add_head_html(
        HEADER_HTML + f"<style>{STYLE_CSS}</style><style>{M3_COMPONENTS_CSS}</style><style>{MOSAIC_COMPONENTS_CSS}</style>"
    )
    # Add mobile-specific responsive CSS
    ui.add_head_html("""
        <style>
        /* Mobile-first responsive design for classification cards */
        @media (max-width: 768px) {
            .classification-cards {
                flex-direction: column !important;
            }
            .classification-card {
                flex: none !important;
                width: 100% !important;
                margin-bottom: 1rem;
            }
        }
        
        /* Ensure proper spacing on mobile */
        @media (max-width: 480px) {
            .mobile-padding {
                padding: 0.375rem !important;
            }
            .mobile-text {
                font-size: 0.75rem !important;
            }
        }
        
        /* Better touch targets for mobile */
        @media (max-width: 768px) {
            .mobile-button {
                min-height: 44px !important;
                min-width: 44px !important;
            }
        }
        
        /* Responsive breakpoints for better mobile experience */
        @media (max-width: 640px) {
            .text-headline-large {
                font-size: 1.5rem !important;
                line-height: 1.2 !important;
            }
            .text-headline-medium {
                font-size: 1.125rem !important;
            }
        }
        
        /* Better spacing for mobile footer */
        @media (max-width: 640px) {
            .footer-mobile {
                flex-direction: column !important;
                gap: 0.5rem !important;
                text-align: center !important;
            }
        }
        
        /* Ultra-compact footer for mobile */
        @media (max-width: 768px) {
            .footer-compact {
                padding: 0.0625rem 0.125rem !important;
                min-height: auto !important;
            }
            .footer-compact .mobile-button {
                min-height: 24px !important;
                min-width: 24px !important;
                padding: 0.0625rem 0.125rem !important;
                margin: 0 0.0625rem !important;
                font-size: 0.6875rem !important;
            }
        }
        
        /* Desktop footer with more space */
        @media (min-width: 769px) {
            .footer-compact {
                padding: 0.5rem 1rem !important;
            }
        }
        
        /* Ultra-compact footer buttons for landscape phones */
        @media (max-width: 896px) and (orientation: landscape) {
            .footer-compact .mobile-button {
                min-height: 20px !important;
                min-width: 20px !important;
                padding: 0.03125rem 0.0625rem !important;
                margin: 0 0.03125rem !important;
                font-size: 0.625rem !important;
            }
        }
        
        @media (max-width: 667px) and (orientation: landscape) {
            .footer-compact .mobile-button {
                min-height: 18px !important;
                min-width: 18px !important;
                padding: 0.015625rem 0.03125rem !important;
                margin: 0 0.015625rem !important;
                font-size: 0.5625rem !important;
            }
        }
        
        /* Smaller labels for activity monitors on mobile */
        @media (max-width: 768px) {
            .text-body-small {
                font-size: 0.75rem !important;
            }
        }
        
        /* Mobile dashboard cards - stack vertically on small screens */
        @media (max-width: 768px) {
            .mobile-run-details {
                flex-direction: column !important;
            }
            .mobile-classification-details {
                flex-direction: column !important;
            }
            .mobile-analysis-details {
                flex-direction: column !important;
            }
            .mobile-dashboard-card {
                flex: none !important;
                width: 100% !important;
                margin-bottom: 0.5rem;
            }
        }
        
        /* Better spacing for dashboard cards on mobile */
        @media (max-width: 480px) {
            .mobile-dashboard-card {
                padding: 0.75rem !important;
            }
            .mobile-dashboard-card .mosaic-card__content {
                font-size: 0.875rem !important;
                line-height: 1.2 !important;
            }
        }
        
        /* Ultra-compact header and footer on mobile */
        @media (max-width: 768px) {
            .q-header {
                min-height: 32px !important;
                padding: 0.125rem 0.25rem !important;
            }
            .q-footer {
                min-height: 32px !important;
                padding: 0.125rem 0.25rem !important;
            }
        }
        
        /* Even more compact on very small screens */
        @media (max-width: 480px) {
            .q-header {
                min-height: 28px !important;
                padding: 0.0625rem 0.125rem !important;
            }
            .q-footer {
                min-height: 28px !important;
                padding: 0.0625rem 0.125rem !important;
            }
        }
        
        /* Ultra-compact for landscape phones */
        @media (max-width: 896px) and (orientation: landscape) {
            .q-header {
                min-height: 24px !important;
                padding: 0.03125rem 0.0625rem !important;
            }
            .q-footer {
                min-height: 24px !important;
                padding: 0.03125rem 0.0625rem !important;
            }
        }
        
        /* Even more compact for small landscape phones */
        @media (max-width: 667px) and (orientation: landscape) {
            .q-header {
                min-height: 20px !important;
                padding: 0.015625rem 0.03125rem !important;
            }
            .q-footer {
                min-height: 20px !important;
                padding: 0.015625rem 0.03125rem !important;
            }
        }
        
        /* Ensure header and footer stay at page extremities */
        @media (max-width: 768px) {
            .q-header {
                position: fixed !important;
                top: 0 !important;
                left: 0 !important;
                right: 0 !important;
                z-index: 1000 !important;
            }
            .q-footer {
                position: fixed !important;
                bottom: 0 !important;
                left: 0 !important;
                right: 0 !important;
                z-index: 1000 !important;
            }
            /* Add padding to main content to account for fixed header/footer */
            .q-page {
                padding-top: 24px !important;
                padding-bottom: 24px !important;
            }
        }
        
        @media (max-width: 480px) {
            .q-page {
                padding-top: 20px !important;
                padding-bottom: 20px !important;
            }
        }
        
        /* Landscape phone content padding */
        @media (max-width: 896px) and (orientation: landscape) {
            .q-page {
                padding-top: 18px !important;
                padding-bottom: 18px !important;
            }
        }
        
        @media (max-width: 667px) and (orientation: landscape) {
            .q-page {
                padding-top: 16px !important;
                padding-bottom: 16px !important;
            }
        }
        
        /* Ultra-compact header text on mobile */
        @media (max-width: 768px) {
            .q-header .text-headline-medium {
                font-size: 1.25rem !important;
                line-height: 1.2 !important;
            }
        }
        
        @media (max-width: 480px) {
            .q-header .text-headline-medium {
                font-size: 1.125rem !important;
                line-height: 1.1 !important;
            }
        }
        
        /* Ultra-compact header text for landscape phones */
        @media (max-width: 896px) and (orientation: landscape) {
            .q-header .text-headline-medium {
                font-size: 1rem !important;
                line-height: 1 !important;
            }
        }
        
        @media (max-width: 667px) and (orientation: landscape) {
            .q-header .text-headline-medium {
                font-size: 0.875rem !important;
                line-height: 1 !important;
            }
        }
        
        /* Ultra-compact logos for landscape phones */
        @media (max-width: 896px) and (orientation: landscape) {
            .q-header img {
                width: 24px !important;
            }
            .q-footer img {
                width: 16px !important;
            }
        }
        
        @media (max-width: 667px) and (orientation: landscape) {
            .q-header img {
                width: 20px !important;
            }
            .q-footer img {
                width: 14px !important;
            }
        }
        
        /* Ensure proper centering of main content */
        @media (max-width: 768px) {
            .main-content-card {
                margin-left: 0.25rem !important;
                margin-right: 0.25rem !important;
                max-width: calc(100% - 0.5rem) !important;
            }
        }

        /* Better centering on very small screens */
        @media (max-width: 480px) {
            .main-content-card {
                margin-left: 0.125rem !important;
                margin-right: 0.125rem !important;
                max-width: calc(100% - 0.25rem) !important;
            }
        }
        </style>
    """)
    ui.add_head_html(
        """
        <script>
        function emitSize() {
            emitEvent('resize', {
                width: document.body.offsetWidth,
                height: document.body.offsetHeight,
            });
        }
        window.onload = emitSize;
        window.onresize = emitSize;
        </script>
    """
    )

    # Create disclaimer dialog that appears on first visit with M3 styling
    async def show_disclaimer():
        # Skip disclaimer in development mode
        if is_development_mode:
            return
            
        try:
            disclaimer_acknowledged = app.storage.general.get("disclaimer_acknowledged", False)
        except RuntimeError:
            # Storage not available in this context, show disclaimer
            disclaimer_acknowledged = False
        
        if not disclaimer_acknowledged:
            with ui.dialog().props(
                "persistent"
            ) as disclaimer_dialog, ui.card().classes(
                "robin-dialog-surface p-4 md:p-6 w-full max-w-2xl max-h-[90vh] overflow-y-auto"
            ):
                ui.label("Disclaimer").classes(
                    "classification-insight-heading text-headline-small q-mb-sm"
                )
                ui.label(
                    "This tool and the data generated by it are intended for research use only and should not be used for "
                    "direct diagnostic purposes. The methylation-based classifications and other analyses provided here may "
                    "be considered by neuropathologists as supplementary information in the context of comprehensive "
                    "diagnostic assessment, which should include clinical history, radiological findings, and complete "
                    "histopathological and molecular evaluation. The final interpretation and diagnosis should always be "
                    "made by qualified healthcare professionals based on all available information."
                ).classes("classification-insight-foot q-mb-md")

                if batphone:
                    ui.label("BATMAN mode").classes(
                        "target-coverage-panel__meta-label q-mb-sm"
                    )
                    ui.label(
                        "You are running this tool in BATMAN mode. "
                        "This is a beta version of the tool and may not be fully functional. "
                        "BATMAN means: Breakpoint Adaptive Targeting alongside Methylation Analysis on Nanopore. "
                        "Target regions update in real time based on detected breakpoints. "
                        "This code only works with ReadFish at this time."
                    ).classes("classification-insight-foot q-mb-md")

                if _process_large_bams_enabled:
                    ui.label("Process large BAMs").classes(
                        "target-coverage-panel__meta-label q-mb-sm"
                    )
                    ui.label(
                        "ROBIN_PROCESS_LARGE_BAMS is enabled. Do not use this option alongside live runs."
                    ).classes("text-amber-700 q-mb-md")

                def acknowledge():
                    try:
                        app.storage.general["disclaimer_acknowledged"] = True
                    except RuntimeError:
                        # Storage not available in this context, skip setting the flag
                        pass
                    disclaimer_dialog.close()

                ui.button("I agree", on_click=acknowledge, icon="check_circle").props(
                    "color=primary no-caps"
                )
            disclaimer_dialog.open()

    ui.timer(0.5, show_disclaimer, once=True)

    # Warn at launch when process-large-BAMs mode is on (do not use with live runs)
    def show_large_bam_launch_warning():
        if _process_large_bams_enabled:
            ui.notify(
                "ROBIN_PROCESS_LARGE_BAMS is enabled. Do not use this option alongside live runs.",
                type="warning",
                timeout=10000,
            )

    ui.timer(1.5, show_large_bam_launch_warning, once=True)

    # Add version check timer
    ui.timer(1.0, check_version, once=True)

    # Create a persistent dialog for quitting the app with M3 styling
    quitdialog = ui.dialog().props("persistent")

    async def quit_app():
        quitdialog.close()
        await cleanup_and_exit()

    def logout_user():
        """Clear auth session and return to login page."""
        try:
            app.storage.user.clear()
        except RuntimeError:
            # Storage may be unavailable in some contexts; still redirect.
            pass
        ui.navigate.to("/login")

    with quitdialog, ui.card().classes(
        "robin-dialog-surface p-4 md:p-5 min-w-[18rem] max-w-md"
    ):
        ui.label("Quit R.O.B.I.N?").classes(
            "classification-insight-heading text-headline-small q-mb-sm"
        )
        ui.label(
            "Quitting the app will stop running methylation analysis."
        ).classes("classification-insight-foot")
        ui.label("If you want to keep analysis running, click Cancel.").classes(
            "classification-insight-foot"
        )
        ui.label(
            "You can safely close this window and analysis will keep running in the background."
        ).classes("classification-insight-foot q-mb-md")
        with ui.row().classes("w-full justify-center gap-2 flex-wrap"):
            ui.button("Cancel", on_click=quitdialog.close).props("flat no-caps outline")
            ui.button("Really quit", icon="logout", on_click=quit_app).props(
                "color=negative no-caps"
            )

    # Create a header with navigation title and menu using M3 styling
    header_classes = "items-center duration-200 p-0 px-2 no-wrap elevation-1"
    if batphone:
        header_classes += " batphone"

    with ui.header(elevated=True).classes(header_classes):
        # Use flexbox layout instead of grid to prevent overlap
        with ui.row().classes("w-full items-center justify-between px-1 py-0.5 sm:px-3 sm:py-1.5"):
            # Left: Hamburger, then title (responsive)
            with ui.row().classes("items-center gap-1 sm:gap-2 min-w-0 flex-1"):
                with ui.button(icon="menu").classes("rounded-md flex-shrink-0"):
                    with ui.menu() as menu:
                        ui.menu_item("Home", lambda: ui.navigate.to("/")).classes(
                            "text-body-medium"
                        )
                        ui.menu_item(
                            "View Samples", lambda: ui.navigate.to("/live_data")
                        ).classes("text-body-medium")
                        ui.menu_item(
                            "Generate Sample ID",
                            lambda: ui.navigate.to("/sample_id_generator"),
                        ).classes("text-body-medium")
                        ui.menu_item(
                            "Watched Folders",
                            lambda: ui.navigate.to("/watched_folders"),
                        ).classes("text-body-medium")
                        ui.menu_item(
                            "Activity Monitor",
                            lambda: ui.navigate.to("/robin"),
                        ).classes("text-body-medium")
                        ui.menu_item(
                            "Workflow",
                            lambda: ui.navigate.to("/workflow"),
                        ).classes("text-body-medium")
                        ui.menu_item(
                            "Documentation",
                            lambda: ui.navigate.to(
                                "https://looselab.github.io/ROBIN/"
                            ),
                        ).classes("text-body-medium")
                        ui.separator()
                        ui.switch(
                            "Allow Remote Access", on_change=use_on_air
                        ).classes("ml-4 bg-transparent").props(
                            'color="primary"'
                        ).bind_value(app.storage.general, "use_on_air")
                        ui.separator()
                        ui.switch("Dark Mode").classes("ml-4 bg-transparent").props(
                            'color="primary"'
                        ).bind_value(app.storage.user, "dark_mode")
                        ui.dark_mode().bind_value(app.storage.user, "dark_mode")
                        ui.separator()
                        ui.menu_item("Close", menu.close).classes(
                            "text-body-medium"
                        )
                        ui.button(
                            "LOG OUT", icon="logout", on_click=logout_user
                        ).classes("bg-error text-white rounded-md")
                        if _is_local_client():
                            ui.button(
                                "Quit", icon="logout", on_click=quitdialog.open
                            ).classes("bg-error text-white rounded-md")

                with ui.row().classes("items-center min-w-0 shrink"):
                    ui.html(navtitle, sanitize=False).classes(
                        f"max-[{MENU_BREAKPOINT}px]:hidden text-headline-medium drop-shadow font-bold truncate"
                    ).style("font-weight: 700; font-family: var(--font-display)")
                    ui.html(smalltitle, sanitize=False).classes(
                        f"min-[{MENU_BREAKPOINT+1}px]:hidden text-headline-medium drop-shadow font-bold truncate"
                    ).style("font-weight: 700; font-family: var(--font-display)")

            # Right: Metrics and logo
            with ui.row().classes("items-center gap-2 flex-shrink-0"):
                ui.label(f"Viewing: {platform.node()}").classes(
                    f"max-[{MENU_BREAKPOINT}px]:hidden text-body-medium"
                )
                ui.label("CPU").classes(
                    f"max-[{MENU_BREAKPOINT}px]:hidden text-body-small"
                )
                cpu_activity = ui.circular_progress(max=100).classes(
                    f"max-[{MENU_BREAKPOINT}px]:hidden"
                )

                ui.label("RAM").classes(
                    f"max-[{MENU_BREAKPOINT}px]:hidden text-body-small"
                )
                ram_utilisation = ui.circular_progress(max=100).classes(
                    f"max-[{MENU_BREAKPOINT}px]:hidden"
                )

                global_metrics.start_timer()
                cpu_activity.bind_value_from(global_metrics, "cpu")
                ram_utilisation.bind_value_from(global_metrics, "ram")

                ui.image(get_imagefile()).classes(
                    "flex-shrink-0 w-10 sm:w-14 object-contain"
                )

    with ui.column().classes(
        "w-full h-full max-w-full overflow-hidden flex flex-col items-center "
        "px-1"
    ) as main_content:
        pass

    # Create a global notification container for progress updates
    with ui.column() as notification_container:
        pass  # This will hold our notifications
    
    # Set up notification system if callback provided
    if setup_notifications:
        setup_notifications(notification_container)

    # Create a footer with useful information and quit button using M3 styling
    footer_classes = "items-center duration-200 p-0 px-2 no-wrap elevation-1"
    if batphone:
        footer_classes += " batphone"
    with ui.footer().classes(footer_classes):
        with ui.dialog() as dialog, ui.card().classes(
            "robin-dialog-surface p-4 md:p-5 min-w-[16rem] max-w-sm"
        ):
            ui.label("Links").classes(
                "classification-insight-heading text-headline-small q-mb-sm"
            )
            ui.separator().classes("mgmt-detail-separator")
            with ui.column().classes("gap-2"):
                ui.link("Code on GitHub", "https://github.com/looselab/robin").classes(
                    "classification-insight-foot"
                )
                ui.link(
                    "Rapid CNS2 Paper",
                    "https://link.springer.com/article/10.1007/s00401-022-02415-6",
                ).classes("classification-insight-foot")
                ui.link(
                    "Sturgeon Classifier",
                    "https://www.nature.com/articles/s41586-023-06615-2",
                ).classes("classification-insight-foot")
                ui.link(
                    "Protocol",
                    "https://www.protocols.io/view/intra-operative-nanopore-sequencing-to-classify-br-c65qzg5w",
                ).classes("classification-insight-foot")
                ui.link("Oxford Nanopore", "https://nanoporetech.com/").classes(
                    "classification-insight-foot"
                )
                ui.link("epi2me labs", "https://labs.epi2me.io/").classes(
                    "classification-insight-foot"
                )
                ui.link("Looselab", "https://looselab.github.io/").classes(
                    "classification-insight-foot"
                )
            ui.button("Close", on_click=dialog.close).props("color=primary no-caps")
        
        # Footer content - ultra-compact on mobile
        with ui.row().classes(
            "w-full items-center justify-between px-1 py-0.5 sm:px-3 sm:py-1.5 gap-0.5 no-wrap footer-compact"
        ):
            # Left side: Logo (ultra-small on mobile)
            ui.image(get_imagefile()).classes(
                "flex-shrink-0 w-7 sm:w-10 object-contain"
            )
            
            # Center: Buttons with proper spacing
            with ui.row().classes("items-center gap-2 flex-shrink-0"):
                ui.colors(primary="#16A34A")  # Emerald 600 (Editorial Bioinformatics primary)
                ui.button("Links", on_click=dialog.open).classes("rounded-md mobile-button text-xs px-2 py-1")

                with ui.button(icon="info").classes("rounded-md mobile-button px-2 py-1"):
                    with ui.menu() as menu:
                        ui.label().bind_text_from(
                            app, "urls", backward=lambda n: f"Available urls: {n}"
                        ).classes("text-body-medium")
                        ui.label("Version: " + get_about().__version__).classes(
                            "text-body-medium"
                        )
            
            # Right side: Compact copyright (mobile only)
            ui.label("©Looselab").classes(
                "text-xs text-weight-italic flex-shrink-0"
            )
            
            # Desktop-only additional info
            ui.label("Not for diagnostic use.").classes(
                f"min-[{MENU_BREAKPOINT+1}px]:block hidden text-xs text-weight-italic flex-shrink-0"
            )

    with main_content:
        yield


async def cleanup_and_exit():
    """
    Handle any necessary cleanup operations before exiting the application and then shut down the application.

    Returns:
        None

    Example:
        >>> cleanup_and_exit()
        None
    """
    logging.info("User initiated shutdown via UI")

    # Create and show shutdown modal with M3 styling
    with ui.dialog().props("persistent") as shutdown_dialog, ui.card().classes(
        "robin-dialog-surface p-4 md:p-5 w-full max-w-sm"
    ):
        ui.label("Shutting down").classes(
            "classification-insight-heading text-headline-small q-mb-sm"
        )
        ui.label(
            "R.O.B.I.N is shutting down. Please wait while we clean up…"
        ).classes("classification-insight-foot q-mb-md")
        with ui.row().classes("w-full justify-center"):
            ui.spinner(size="lg", color="primary")

    # Close the quit dialog if it's open
    if quitdialog:
        quitdialog.close()

    # Show the shutdown dialog
    shutdown_dialog.open()
    logging.info("Shutdown dialog opened")

    # Wait a moment to ensure the dialog is visible
    await asyncio.sleep(1.0)
    logging.info("Initial wait complete")

    # Perform cleanup
    logging.info("Performing cleanup operations...")
    logging.info("Shutting down ROBIN... from theme.py")
    logging.info(
        "Here we need to do some very graceful shutdown to make sure we don't leave any threads running and we don't leave any files open."
    )

    # Close the shutdown dialog
    shutdown_dialog.close()
    logging.info("Shutdown dialog closed")

    # Shutdown the application
    logging.info("Application shutdown initiated")
    app.shutdown()


def set_remote_access(enabled: bool) -> None:
    print(f"[ROBIN] Setting remote access to {enabled}")
    """Queue a global NiceGUI On Air state change."""
    global _remote_access_target, _remote_access_task, _remote_access_request_id
    _remote_access_target = bool(enabled)
    _remote_access_request_id += 1
    request_id = _remote_access_request_id
    if _remote_access_task is not None and not _remote_access_task.done():
        _remote_access_task.cancel()
    _remote_access_task = background_tasks.create(
        _apply_remote_access_state(request_id=request_id, target=_remote_access_target)
    )


async def _apply_remote_access_state(request_id: int, target: bool) -> None:
    """Apply one requested On Air state; stale requests are ignored."""
    global _remote_access_task
    try:
        if request_id != _remote_access_request_id:
            print("[ROBIN] Ignoring stale On Air request")
            return

        if target:
            if core.air is None:
                print("[ROBIN] Activating On Air: creating Air instance")
                core.air = nicegui.air.Air("")
            # Air has an internal reconnect timer; re-enable connect path when toggled on.
            core.air.connecting = False
            print("[ROBIN] Activating On Air: connecting")
            await core.air.connect()
            return

        air_instance = core.air
        if air_instance is not None:
            # Air starts an internal timer(5, self.connect); keep connecting=True to block auto-reconnect.
            air_instance.connecting = True
            print("[ROBIN] Removing On Air: disconnecting")
            await air_instance.disconnect()
            print("[ROBIN] Removing On Air: auto-reconnect blocked")
        else:
            print("[ROBIN] On Air already inactive")
    except asyncio.CancelledError:
        print("[ROBIN] Cancelled stale On Air transition")
        raise
    finally:
        _remote_access_task = None


def initialize_remote_access_once() -> None:
    """Set initial remote access state once per process launch."""
    global _remote_access_initialized
    if _remote_access_initialized:
        return
    _remote_access_initialized = True
    set_remote_access(False)
    try:
        app.storage.general["use_on_air"] = False
    except RuntimeError:
        # Storage may not be available in some contexts.
        pass


def use_on_air(args: events.ValueChangeEventArguments):
    """
    Enable or disable remote access based on the value of the event argument.

    Args:
        args (events.ValueChangeEventArguments): The event argument containing the value for remote access toggle.

    Returns:
        None

    Example:
        >>> args = events.ValueChangeEventArguments(value=True)
        >>> use_on_air(args)
        None
    """
    value = args.value
    if isinstance(value, str):
        enabled = value.strip().lower() in {"1", "true", "yes", "on"}
    else:
        enabled = bool(value)
    print(f"[ROBIN] Allow Remote Access toggle is {'ON' if enabled else 'OFF'}")
    set_remote_access(enabled)


def create_home_page():
    """Create the home page content."""
    with frame(
        "<strong>R</strong>apid nanop<strong>O</strong>re <strong>B</strong>rain intraoperat<strong>I</strong>ve classificatio<strong>N</strong>",
        smalltitle="<strong>R.O.B.I.N</strong>",
    ):
        
        ui.label("Welcome to the Application").classes(
            "text-headline-large text-center px-3"
        )
        with ui.row().classes('items-center m-auto'):
            ui.circular_progress(value=0.1, show_value=False, size="xs")
            ui.circular_progress(value=0.1, show_value=False, size="xl") 
        with ui.row().classes('items-center m-auto'):
            with ui.circular_progress(value=0.1, show_value=False, size="sm") as progress:
                ui.button(
                    icon='star',
                    on_click=lambda: progress.set_value(progress.value + 0.1)
                ).props('flat round')
            ui.label('click to increase progress')
        with ui.card().classes("w-full max-w-4xl mx-auto mobile-padding main-content-card").style("border: 2px solid var(--md-primary)"):
            with ui.row().classes("w-full flex justify-between items-center flex-wrap gap-2"):
                ui.label('Sample Name').classes("text-headline-medium flex-shrink-0")
                ui.button("Button").classes("bg-primary text-white rounded-md mobile-button")
            ui.separator().classes().style("border: 1px solid var(--md-primary)")
            with ui.card().classes("w-full bg-gradient-to-r from-blue-50 to-indigo-50 mobile-padding"):
                ui.label("Run Information").classes("text-lg font-semibold mb-3 text-blue-800")
                with ui.row().classes("w-full gap-2 sm:gap-6 items-center flex-wrap"):
                    ui.label("Run").classes("text-body-medium text-xs sm:text-sm")
                    ui.label("Model").classes("text-body-medium text-xs sm:text-sm")
                    ui.label("Device").classes("text-body-medium text-xs sm:text-sm")
                    ui.label("Flow Cell").classes("text-body-medium text-xs sm:text-sm")
                    ui.label("Sample").classes("text-body-medium text-xs sm:text-sm")


            with ui.card().classes("w-full mobile-padding"):
                ui.label("Classification Results").classes("text-lg font-semibold mb-3 text-blue-800")
                # Use responsive grid: 2 columns on desktop, 1 column on mobile
                with ui.row().classes("w-full gap-2 sm:gap-3 flex-wrap classification-cards"):
                    # Sturgeon Classification
                    with ui.card().classes("flex-1 min-w-0 elevation-4 rounded-xl bg-gradient-to-br from-blue-50 to-indigo-50 border-l-4 border-blue-500 classification-card"):
                        ui.label("Sturgeon Classification").classes("font-bold text-blue-800 mb-2 text-sm sm:text-base")
                        ui.label("Class: --").classes("font-bold text-medium text-blue-600 text-xs sm:text-sm")
                        ui.label("Confidence: --%").classes("text-xs sm:text-sm text-blue-600")
                        ui.label("Probes: --").classes("text-xs sm:text-sm text-blue-600")
                        ui.label("Model: --").classes("text-xs sm:text-sm text-blue-600")
                        ui.label("Features: --").classes("text-xs sm:text-sm text-blue-600")

                    # NanoDX Classification
                    with ui.card().classes("flex-1 min-w-0 elevation-4 rounded-xl bg-gradient-to-br from-green-50 to-green-100 border-l-4 border-green-500 classification-card"):
                        ui.label("NanoDX Classification").classes("font-bold text-green-800 mb-2 text-sm sm:text-base")
                        ui.label("Class: --").classes("font-bold text-medium text-green-600 text-xs sm:text-sm")
                        ui.label("Confidence: --%").classes("text-xs sm:text-sm text-green-600")
                        ui.label("Probes: --").classes("text-xs sm:text-sm text-green-600")
                        ui.label("Model: --").classes("text-xs sm:text-sm text-green-600")
                        ui.label("Features: --").classes("text-xs sm:text-sm text-green-600")

                    # PanNanoDX Classification
                    with ui.card().classes("flex-1 min-w-0 elevation-4 rounded-xl bg-gradient-to-br from-purple-50 to-purple-100 border-l-4 border-purple-500 classification-card"):
                        ui.label("PanNanoDX Classification").classes("font-bold text-purple-800 mb-2 text-sm sm:text-base")
                        ui.label("Class: --").classes("font-bold text-medium text-purple-600 text-xs sm:text-sm")
                        ui.label("Confidence: --%").classes("text-xs sm:text-sm text-purple-600")
                        ui.label("Probes: --").classes("text-xs sm:text-sm text-purple-600")
                        ui.label("Model: --").classes("text-xs sm:text-sm text-purple-600")
                        ui.label("Features: --").classes("text-xs sm:text-sm text-purple-600")

                    # Random Forest Classification
                    with ui.card().classes("flex-1 min-w-0 elevation-4 rounded-xl bg-gradient-to-br from-orange-50 to-orange-100 border-l-4 border-orange-500 classification-card"):
                        ui.label("Random Forest Classification").classes("font-bold text-orange-800 mb-2 text-sm sm:text-base")
                        ui.label("Class: --").classes("font-bold text-medium text-orange-600 text-xs sm:text-sm")
                        ui.label("Confidence: --%").classes("text-xs sm:text-sm text-orange-600")
                        ui.label("Probes: --").classes("text-xs sm:text-sm text-orange-600")
                        ui.label("Model: --").classes("text-xs sm:text-sm text-orange-600")
                        ui.label("Features: --").classes("text-xs sm:text-sm text-orange-600")


            ui.label('text below').classes("text-body-medium px-3")
            with ui.card_section().classes("px-3"):
                ui.image('https://picsum.photos/id/684/640/360').classes("w-full h-auto rounded-lg")
                ui.label('Lorem ipsum dolor sit amet, consectetur adipiscing elit, ...').classes("text-body-medium mobile-text")


def create_standalone_page():
    """Create a simplified standalone page without storage dependencies."""
    # Add custom HTML and CSS to the head of the page
    ui.add_head_html(
        '<script src="https://cdn.jsdelivr.net/npm/igv@3.7.0/dist/igv.min.js"></script>'
    )
    ui.add_head_html(EDITORIAL_FONTS_HTML)
    ui.add_head_html(
        HEADER_HTML + f"<style>{STYLE_CSS}</style><style>{M3_COMPONENTS_CSS}</style><style>{MOSAIC_COMPONENTS_CSS}</style>"
    )
    
    # Create a simple header (same shell padding pattern as frame())
    with ui.header(elevated=True).classes("items-center duration-200 p-0 px-2 no-wrap elevation-1"):
        with ui.row().classes("w-full items-center justify-between px-1 py-0.5 sm:px-3 sm:py-1.5"):
            ui.html("<strong>R.O.B.I.N</strong>", sanitize=False).classes("text-headline-medium drop-shadow font-bold").style(
                "font-weight: 700; font-family: var(--font-display)"
            )
            ui.image(get_imagefile()).style("width: 50px").classes("ml-auto")
    
    # Create main content
    with ui.column().classes("w-full h-full max-w-full overflow-hidden p-6"):
        ui.label("Welcome to ROBIN Theme Test").classes("text-headline-large text-center")
        ui.label("This is a standalone test of the ROBIN theme system.").classes("text-body-large text-center mt-4")
        ui.label(f"Version: {get_about().__version__}").classes("text-body-medium text-center mt-2")
    
    # Create a simple footer (same shell padding pattern as frame())
    with ui.footer().classes("items-center duration-200 p-0 px-2 no-wrap elevation-1"):
        with ui.row().classes("w-full items-center justify-between px-1 py-0.5 sm:px-3 sm:py-1.5"):
            ui.image(get_imagefile()).style("width: 40px")
            ui.label("ROBIN Theme Test - Standalone Mode").classes("text-body-small")


def create_workflow_page():
    """Display the ROBIN workflow diagram with editorial insight styling."""
    with frame(
        "ROBIN Workflow",
        smalltitle="Workflow",
    ):
        # Get versions for the diagram
        sturgeon_version = get_sturgeon_version()
        crossnn_version = get_crossnn_version()
        cnv_from_bam_version = get_cnv_from_bam_version()

        with ui.element("div").classes("w-full min-w-0").props(
            "id=workflow-diagram-page"
        ):
            with ui.column().classes(
                "w-full max-w-6xl mx-auto gap-3 p-2 md:p-3"
            ):
                with ui.element("div").classes(
                    "classification-insight-shell w-full min-w-0"
                ):
                    ui.label("Pipeline overview").classes(
                        "classification-insight-heading text-headline-small"
                    )
                    ui.label(
                        "From MinKNOW through R.O.B.I.N: sequencing, alignment, "
                        "preprocessing, classifiers, analyses, and report generation."
                    ).classes("classification-insight-foot")

                with ui.element("div").classes(
                    "classification-insight-card w-full min-w-0"
                ):
                    with ui.column().classes(
                        "w-full min-w-0 gap-3 p-2 md:p-3"
                    ):
                        with ui.row().classes(
                            "items-center gap-2 min-w-0"
                        ):
                            ui.icon("account_tree").classes(
                                "classification-insight-icon"
                            )
                            ui.label("Workflow diagram").classes(
                                "classification-insight-model flex-1 min-w-0"
                            )
                        ui.label(
                            "Interactive chart: pan and zoom if your browser supports it."
                        ).classes("classification-insight-foot")

                        with ui.element("div").classes(
                            "w-full min-w-0 overflow-x-auto workflow-diagram-scroll"
                        ):
                            ui.mermaid(
            f"""
flowchart TD
    %% Style definitions with M3 color palette
    classDef minKNOW fill:#E8DEF8,stroke:#6750A4,stroke-width:2px,color:#1D192B,font-size:14px,font-weight:500
    classDef robin fill:#EADDFF,stroke:#6750A4,stroke-width:2px,color:#21005D,font-size:14px,font-weight:500
    classDef classifier fill:#FFD8E4,stroke:#7D5260,stroke-width:2px,color:#31111D,font-size:14px,font-weight:500
    classDef analysis fill:#FFD8E4,stroke:#7D5260,stroke-width:2px,color:#31111D,font-size:14px,font-weight:500
    classDef output fill:#EADDFF,stroke:#6750A4,stroke-width:2px,color:#21005D,font-size:14px,font-weight:500
    classDef headerLabel fill:#ffffff00,stroke:#ffffff00,color:#1C1B1F,font-size:18px,font-weight:700

    %% Custom header nodes
    robinLabel["<b>R.O.B.I.N v{get_about().__version__}</b>"]:::headerLabel
    minKNOWLabel["<b>MinKNOW Pipeline</b>"]:::headerLabel

    %% MinKNOW pipeline
    subgraph MinKNOW[" "]
        sequencing["Sequencing"]
        adaptive["Adaptive Sampling"]
        alignment["Alignment"]
        bam["BAM Files"]
    end
    minKNOWLabel -.-> MinKNOW

    %% R.O.B.I.N. pipeline
    subgraph ROBIN[" "]
        runInfo["Extract Run Information"]
        groupRuns["Group by Run"]
        mergeBam["Merge BAM Files"]
        methylation["Extract Methylation<br>(modkit retired)"]
        individualBam["Process BAMs Individually"]
        cnv["CNV Analysis<br>(CNV from BAM v{cnv_from_bam_version})"]
        fusion["Fusion Analysis"]
        mgmt["MGMT Analysis"]
        coverage["Coverage Analysis"]
        snv["Optional SNP/V Analysis<br>(ClairS-To)"]
        crossnnCNS["CrossNN - CNS<br>({crossnn_version})"]
        crossnnPan["CrossNN - PanCancer<br>({crossnn_version})"]
        sturgeon["Sturgeon<br>({sturgeon_version})"]
        randomForest["Random Forest"]
        integrate["Integrate Results"]
        report["Report Generation"]
    end
    robinLabel -.-> ROBIN

    %% Connections
    sequencing --> adaptive
    adaptive --> alignment
    alignment --> bam
    bam --> runInfo
    runInfo --> groupRuns

    %% Split after groupRuns
    groupRuns -->|"For methylation"| mergeBam
    groupRuns -->|"For other analyses"| individualBam

    mergeBam --> methylation
    methylation --> crossnnCNS & crossnnPan & sturgeon & randomForest
    crossnnCNS --> integrate
    crossnnPan --> integrate
    sturgeon --> integrate
    randomForest --> integrate

    individualBam --> cnv & fusion & mgmt & coverage
    coverage --> snv

    integrate --> report
    cnv --> report
    fusion --> report
    mgmt --> report
    snv --> report

    %% Styling
    class sequencing,adaptive,alignment,bam minKNOW
    class runInfo,groupRuns,mergeBam,individualBam,integrate robin
    class crossnnCNS,crossnnPan,sturgeon,randomForest classifier
    class cnv,fusion,mgmt,coverage,snv analysis
    class report output
    %% Subgraph background coloring with M3 colors
    style MinKNOW fill:#E8DEF8,stroke:#6750A4,stroke-width:2px
    style ROBIN fill:#EADDFF,stroke:#6750A4,stroke-width:2px
""",
            config={
                "theme": "redux",
                "look": "neo",
                "flowchart": {"curve": "basis", "defaultRenderer": "elk"},
            },
        ).classes("w-full min-w-0 workflow-diagram-mermaid")


def register_theme_pages():
    """Register the theme pages. This function should be called when the module is imported."""
    @ui.page("/")
    def home_page():
        create_home_page()

    @ui.page("/workflow")
    def workflow_page():
        create_workflow_page()


def get_modkit_version():
    """Get the version of modkit installed."""
    try:
        result = subprocess.run(["modkit", "--version"], capture_output=True, text=True)
        if result.returncode == 0:
            return result.stdout.strip()
        return "unknown"
    except (subprocess.SubprocessError, FileNotFoundError, OSError):
        return "unknown"


def get_sturgeon_version():
    """Get the version of sturgeon installed."""
    try:
        result = subprocess.run(
            ["sturgeon", "--version"], capture_output=True, text=True
        )
        if result.returncode == 0:
            return result.stdout.strip()
        return "unknown"
    except (subprocess.SubprocessError, FileNotFoundError, OSError):
        try:
            return importlib.metadata.version("sturgeon")
        except importlib.metadata.PackageNotFoundError:
            return "unknown"


def get_crossnn_version():
    """Get the version of CrossNN installed."""
    try:
        return importlib.metadata.version("nanoDX")
    except importlib.metadata.PackageNotFoundError:
        try:
            result = subprocess.run(
                ["git", "describe", "--tags"],
                cwd=os.path.join(
                    os.path.dirname(os.path.abspath(__file__)),
                    "submodules/nanoDX",
                ),
                capture_output=True,
                text=True,
            )
            if result.returncode == 0:
                return result.stdout.strip()
            return "unknown"
        except (subprocess.SubprocessError, FileNotFoundError, OSError):
            return "unknown"


def get_cnv_from_bam_version():
    """Get the version of CNV from BAM installed."""
    try:
        return importlib.metadata.version("cnv_from_bam")
    except importlib.metadata.PackageNotFoundError:
        try:
            result = subprocess.run(
                ["git", "describe", "--tags"],
                cwd=os.path.join(
                    os.path.dirname(os.path.abspath(__file__)),
                    "submodules/cnv_from_bam",
                ),
                capture_output=True,
                text=True,
            )
            if result.returncode == 0:
                return result.stdout.strip()
            return "unknown"
        except (subprocess.SubprocessError, FileNotFoundError, OSError):
            return "unknown"


def main():
    """
    Main function to test the theme by creating a simple page using the frame context manager.

    Example:
        >>> main()
        None
    """
    ui.add_head_html(EDITORIAL_FONTS_HTML)

    # Register some fonts that we might need later on (guarded)
    fonts_dir = Path(__file__).parent / "fonts"
    if fonts_dir.exists():
        app.add_static_files("/fonts", str(fonts_dir))

    # For standalone execution, create UI directly without @ui.page decorators
    # This completely avoids the global scope issue with NiceGUI
    create_standalone_page()

    ui.run(storage_secret="robin")


# Register theme pages when module is imported (but not when run as main)
if __name__ not in {"__main__", "__mp_main__"}:
    register_theme_pages()

if __name__ in {"__main__", "__mp_main__"}:
    main()


def get_process_ram_usage():
    """Get RAM usage for Python and R processes."""
    python_ram = 0
    r_ram = 0

    for proc in psutil.process_iter(["pid", "name", "memory_info"]):
        try:
            # Get process name and memory info
            name = proc.info["name"].lower()
            mem_info = proc.info["memory_info"]

            # Skip if memory_info is None
            if mem_info is None:
                continue

            # Calculate RAM usage in GB
            ram_gb = mem_info.rss / (1024 * 1024 * 1024)  # Convert bytes to GB

            # Only count processes that are actually using memory
            if ram_gb > 0:
                if "python" in name:
                    python_ram += ram_gb
                elif "r" in name or "Rscript" in name:
                    r_ram += ram_gb

        except (psutil.NoSuchProcess, psutil.AccessDenied, psutil.ZombieProcess):
            continue
        except Exception as e:
            logging.debug(
                f"Error getting memory info for process {proc.info.get('name', 'unknown')}: {str(e)}"
            )
            continue

    return round(python_ram, 2), round(
        r_ram, 2
    )  # Round to 2 decimal places for cleaner display


def debug_process_tree():
    """
    Debug function to understand why we're picking up incorrect child processes.
    Shows the full process tree and explains why each process is being counted.
    """
    try:
        main_process = psutil.Process(os.getpid())

        logging.debug("=== Process Tree Debug ===")
        logging.debug(
            f"Main ROBIN process: {main_process.name()} (PID: {main_process.pid})"
        )
        logging.debug(f"Command line: {' '.join(main_process.cmdline())}")
        logging.debug(f"Parent PID: {main_process.ppid()}")

        # Get all processes in the system
        all_processes = []
        for proc in psutil.process_iter(["pid", "name", "ppid", "cmdline"]):
            try:
                all_processes.append(proc.info)
            except (psutil.NoSuchProcess, psutil.AccessDenied):
                continue

        # Find all processes that could be considered "children"
        potential_children = []
        for proc_info in all_processes:
            pid = proc_info["pid"]
            ppid = proc_info["ppid"]

            # Check if this process is in our process tree
            if _is_in_process_tree(pid, main_process.pid, all_processes):
                potential_children.append(proc_info)

        logging.debug(
            f"\nFound {len(potential_children)} processes in ROBIN's process tree:"
        )

        for i, proc_info in enumerate(potential_children):
            pid = proc_info["pid"]
            name = proc_info["name"]
            ppid = proc_info["ppid"]
            cmdline = " ".join(proc_info["cmdline"]) if proc_info["cmdline"] else "N/A"

            # Determine why this process is being counted
            reason = _explain_process_inclusion(proc_info, main_process.pid)

            logging.debug(f"\n{i+1}. {name} (PID: {pid}, PPID: {ppid})")
            logging.debug(f"   Command: {cmdline}")
            logging.debug(f"   Reason: {reason}")

            # Show memory usage if available
            try:
                proc = psutil.Process(pid)
                rss = proc.memory_info().rss / (1024 * 1024 * 1024)
                logging.debug(f"   Memory: {rss:.2f}GB RSS")
            except Exception as e:
                logging.debug(f"   Memory: <access denied>: {e}")
                
        logging.debug("\n=== End Process Tree Debug ===")

    except Exception as e:
        logging.debug(f"Process tree debug failed: {e}")


def _is_in_process_tree(target_pid, root_pid, all_processes):
    """
    Check if a process is in the process tree starting from root_pid.
    """
    if target_pid == root_pid:
        return True

    # Find the process
    target_proc = None
    for proc_info in all_processes:
        if proc_info["pid"] == target_pid:
            target_proc = proc_info
            break

    if not target_proc:
        return False

    # Recursively check if parent is in the tree
    return _is_in_process_tree(target_proc["ppid"], root_pid, all_processes)


def _explain_process_inclusion(proc_info, main_pid):
    """
    Explain why a process is being included in our child process count.
    """
    ppid = proc_info["ppid"]
    name = proc_info["name"]
    cmdline = proc_info["cmdline"]

    reasons = []

    # Direct child
    if ppid == main_pid:
        reasons.append("Direct child of main ROBIN process")

    # Python process with ROBIN in command line
    if name in ["python", "python3", "python3.9", "robin"]:
        if cmdline and any("robin" in arg.lower() for arg in cmdline):
            reasons.append("Python process with ROBIN in command line")

    # Analysis tool
    if name in ["modkit", "matkit", "sturgeon", "R", "Rscript", "bedtools"]:
        reasons.append("Analysis tool spawned by ROBIN")

    # Process in tree but not direct child
    if ppid != main_pid and ppid != 1:
        reasons.append("Process in ROBIN's process tree (indirect child)")

    # Reparented process
    if ppid == 1:
        reasons.append("Reparented to init/systemd (orphaned process)")

    # Shared library or dependency
    if name in ["libc", "libpython", "libssl", "libcrypto"]:
        reasons.append("Shared library process")

    if not reasons:
        reasons.append("Unknown reason - process in tree but unclear relationship")

    return "; ".join(reasons)
