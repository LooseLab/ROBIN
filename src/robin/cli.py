"""
Command-line interface for robin.

This module provides a comprehensive CLI for running bioinformatics workflows on BAM files.
It supports both simplified and legacy workflow formats, with options for distributed
computing using Ray and traditional threading.

Key Features:
- Simplified workflow syntax: 'mgmt,sturgeon' (auto-queue assignment)
- Legacy workflow syntax: 'preprocessing:bed_conversion,mgmt:mgmt,classification:sturgeon'
- Automatic dependency management (e.g., bed_conversion for sturgeon analysis)
- Ray-based distributed computing support
- Configurable logging levels per job type
- Job deduplication by sample ID
- Progress tracking and verbose output options

Code Quality Improvements:
- Modular function design with single responsibilities
- Comprehensive input validation and error handling
- Constants for configuration values
- Better type hints and documentation
- Graceful error handling with informative messages
- Consistent error reporting to stderr
- Proper cleanup on interruption
"""

# Suppress pkg_resources deprecation warnings from sorted_nearest
import warnings
warnings.filterwarnings(
    "ignore", message="pkg_resources is deprecated", category=UserWarning
)
# Suppress matplotlib tight_layout warnings
warnings.filterwarnings(
    "ignore", message="The figure layout has changed to tight", category=UserWarning
)

import os
import sys
import csv
import shutil
import tempfile
from pathlib import Path
from typing import Optional, List, Dict, Tuple, Any, Iterable

import click
import logging

# Check if we're in development mode
is_development_mode = os.environ.get("ROBIN_DEV_MODE", "").lower() in ("1", "true", "yes", "on")

from robin.workflow_simple import default_file_classifier, Job

# Many analysis handlers have optional third-party dependencies. Import them lazily
# so lightweight commands (e.g. `robin utils update-models`) still work.
_analysis_import_error: Optional[BaseException] = None
try:
    from robin.analysis.bam_preprocessor import bam_preprocessing_handler
    from robin.analysis.mgmt_analysis import mgmt_handler
    from robin.analysis.cnv_analysis import cnv_handler
    from robin.analysis.bed_conversion import bed_conversion_handler
    from robin.analysis.sturgeon_analysis import sturgeon_handler
    from robin.analysis.nanodx_analysis import nanodx_handler, pannanodx_handler
    from robin.analysis.random_forest_analysis import random_forest_handler
    from robin.analysis.target_analysis import target_handler
    from robin.analysis.fusion_analysis import fusion_handler
    from robin.analysis.utilities.matkit import run_matkit
    from robin.analysis.mgmt_analysis import extract_mgmt_site_rows_from_bed
except Exception as e:
    _analysis_import_error = e
    bam_preprocessing_handler = None  # type: ignore[assignment]
    mgmt_handler = None  # type: ignore[assignment]
    cnv_handler = None  # type: ignore[assignment]
    bed_conversion_handler = None  # type: ignore[assignment]
    sturgeon_handler = None  # type: ignore[assignment]
    nanodx_handler = None  # type: ignore[assignment]
    pannanodx_handler = None  # type: ignore[assignment]
    random_forest_handler = None  # type: ignore[assignment]
    target_handler = None  # type: ignore[assignment]
    fusion_handler = None  # type: ignore[assignment]
    run_matkit = None  # type: ignore[assignment]
    extract_mgmt_site_rows_from_bed = None  # type: ignore[assignment]
from robin.logging_config import (
    configure_logging,
)
from robin.utils.sequencing_files import (
    DEFAULT_GRCH38_REFERENCE_URL,
    copy_panel_bed_to,
    copy_panel_source_bed_if_present,
    describe_reference_action,
    materialize_reference,
    panel_bed_filename,
    panel_source_available,
    panel_source_filename,
)


def _download_missing_models(missing_files, models_dir):
    """Download missing model files (same asset manifest logic as ``robin utils update-models``)."""
    import json
    import hashlib
    import urllib.request
    import urllib.error
    import os
    
    print("\n🔄 Attempting to download missing models...")
    
    # Find project root from models_dir location
    # models_dir is src/robin/models, so project root is 3 levels up
    project_root = models_dir.parent.parent.parent
    
    # Load assets manifest
    try:
        assets_file = project_root / "assets.json"
        if not assets_file.exists():
            print(f"❌ assets.json not found at {assets_file}. Cannot download models automatically.")
            return False
        
        with open(assets_file, 'r') as f:
            manifest = json.load(f)
    except Exception as e:
        print(f"❌ Failed to load assets manifest: {e}")
        return False
    
    # Asset name mapping
    asset_mapping = {
        "general.zip": "general_model",
        "Capper_et_al_NN_v2.pkl": "capper_model", 
        "pancan_devel_v5i_NN_v2.pkl": "pancan_model"
    }
    
    github_token = os.getenv('GITHUB_TOKEN')
    if not github_token:
        print("ℹ️  No GITHUB_TOKEN found. Trying public download...")
    
    success_count = 0
    for filename in missing_files:
        if filename not in asset_mapping:
            print(f"⚠️  Unknown model file: {filename}")
            continue
            
        asset_name = asset_mapping[filename]
        if asset_name not in manifest["assets"]:
            print(f"❌ Asset '{asset_name}' not found in manifest")
            continue
            
        asset_info = manifest["assets"][asset_name]
        asset_url = asset_info["url"]
        expected_sha256 = asset_info["sha256"]
        
        target_path = models_dir / filename
        
        try:
            print(f"\n📥 Downloading {filename}...")
            
            # Download the file
            headers = {}
            if github_token:
                headers["Authorization"] = f"Bearer {github_token}"
            
            request = urllib.request.Request(asset_url, headers=headers)
            
            with urllib.request.urlopen(request) as response:
                with open(target_path, 'wb') as f:
                    f.write(response.read())
            
            # Verify checksum
            print("🔍 Verifying checksum...")
            sha256_hash = hashlib.sha256()
            with open(target_path, "rb") as f:
                for chunk in iter(lambda: f.read(4096), b""):
                    sha256_hash.update(chunk)
            calculated_sha256 = sha256_hash.hexdigest()
            
            if calculated_sha256 != expected_sha256:
                print(f"❌ Checksum mismatch for {filename}")
                print(f"Expected: {expected_sha256}")
                print(f"Got:      {calculated_sha256}")
                target_path.unlink()
                continue
            
            print(f"✅ Successfully downloaded {filename}")
            success_count += 1
            
        except urllib.error.HTTPError as e:
            if e.code == 401:
                print(f"❌ Authentication failed for {filename}. Need GitHub token.")
            elif e.code == 404:
                print(f"❌ Asset not found: {filename}")
            else:
                print(f"❌ HTTP error {e.code} downloading {filename}: {e.reason}")
        except Exception as e:
            print(f"❌ Failed to download {filename}: {e}")
    
    return success_count == len(missing_files)


def _check_models_or_exit():
    """Ensure required runtime assets exist; auto-download missing ones."""
    try:
        from robin.utils.model_checker import get_models_directory, check_model_files
        from robin.utils.model_updater import update_models as _update_models
        from robin.utils.clinvar_manager import ensure_clinvar_files
    except Exception as e:
        click.echo(f"❌ Could not load asset bootstrap helpers: {e}", err=True)
        sys.exit(1)

    # 1) Ensure required models are available (auto-download missing on first run).
    try:
        models_dir = get_models_directory()
    except Exception as e:
        click.echo(f"❌ Could not locate ROBIN models directory: {e}", err=True)
        sys.exit(1)

    all_present, missing_files, _present_files = check_model_files()
    if not all_present:
        click.echo("Missing required model files detected. Downloading now...")
        ok, msgs = _update_models(models_dir=models_dir, overwrite=False)
        for m in msgs:
            click.echo(m)
        if not ok:
            click.echo(
                "❌ Failed to download required model files. "
                "You can retry with: robin utils update-models",
                err=True,
            )
            sys.exit(1)

        all_present_after, missing_after, _present_after = check_model_files()
        if not all_present_after:
            click.echo(
                f"❌ Required model files are still missing after download: {', '.join(missing_after)}",
                err=True,
            )
            sys.exit(1)

    # 2) Ensure ClinVar VCF resources are available (download or generate as needed).
    try:
        ensure_clinvar_files(download_if_missing=True)
    except Exception as e:
        click.echo(
            "❌ Failed to set up ClinVar VCF files automatically. "
            f"You can retry with: robin utils update-clinvar\nReason: {e}",
            err=True,
        )
        sys.exit(1)


# Constants
VALID_LOG_LEVELS = {"DEBUG", "INFO", "WARNING", "ERROR"}
VALID_JOB_TYPES = {
    "preprocessing",
    "bed_conversion",
    "mgmt",
    "cnv",
    "target",
    "fusion",
    "sturgeon",
    "nanodx",
    "pannanodx",
    "random_forest",
}
DEFAULT_LOG_LEVEL = "ERROR"
DEFAULT_ANALYSIS_WORKERS = 1
DEFAULT_TIMEOUT = 5.0

# Queue mapping for job types
QUEUE_MAPPING = {
    "preprocessing": "preprocessing",
    "bed_conversion": "bed_conversion",
    "mgmt": "mgmt",
    "cnv": "cnv",
    "target": "target",
    "fusion": "fusion",
    "sturgeon": "classification",
    "nanodx": "classification",
    "pannanodx": "classification",
    "random_forest": "slow",
}

# Jobs that require bed_conversion as a dependency
JOBS_REQUIRING_BED_CONVERSION = {"sturgeon", "nanodx", "pannanodx", "random_forest"}

from robin.reporting.sections.disclaimer_text import EXTENDED_DISCLAIMER_TEXT

try:
    from rich.console import Console
    from rich.panel import Panel

    _RICH_AVAILABLE = True
except Exception:
    _RICH_AVAILABLE = False

_RICH_CONSOLE = Console() if _RICH_AVAILABLE else None


def _echo_styled(message: str, level: str = "info") -> None:
    """Print with optional rich styling (no emojis)."""
    if not _RICH_AVAILABLE or _RICH_CONSOLE is None:
        click.echo(message)
        return

    styles = {
        "info": "cyan",
        "success": "green",
        "warn": "yellow",
        "error": "red",
        "header": "bold magenta",
        "muted": "dim",
    }
    style = styles.get(level, "white")
    _RICH_CONSOLE.print(message, style=style)

# Disclaimer text for user acknowledgment
DISCLAIMER_TEXT = EXTENDED_DISCLAIMER_TEXT

# Handler configurations
HANDLER_CONFIGS = [
    # (queue_type, job_type, handler_func, legacy_queue_type, needs_work_dir)
    ("preprocessing", "preprocessing", bam_preprocessing_handler, None, False),
    ("bed_conversion", "bed_conversion", bed_conversion_handler, None, True),
    ("mgmt", "mgmt", mgmt_handler, "analysis", True),
    ("cnv", "cnv", cnv_handler, "analysis", True),
    ("target", "target", target_handler, "analysis", True),
    ("fusion", "fusion", fusion_handler, "analysis", True),
    ("classification", "sturgeon", sturgeon_handler, None, True),
    ("classification", "nanodx", nanodx_handler, None, True),
    ("classification", "pannanodx", pannanodx_handler, None, True),
    ("slow", "random_forest", random_forest_handler, None, True),
]


def _get_user_acknowledgment() -> bool:
    """Display a disclaimer and require explicit user acknowledgment.

    Returns:
        bool: True if the user types 'I agree' exactly, otherwise False.
    """
    # Skip disclaimer in development mode
    if is_development_mode:
        return True
        
    if _RICH_AVAILABLE and _RICH_CONSOLE is not None:
        _RICH_CONSOLE.print(
            Panel(
                DISCLAIMER_TEXT,
                title="DISCLAIMER",
                title_align="left",
                border_style="yellow",
                style="yellow",
                expand=False,
            )
        )
    else:
        click.echo("\n" + "=" * 70)
        click.echo("DISCLAIMER:")
        click.echo("=" * 70)
        click.echo(DISCLAIMER_TEXT)
        click.echo("=" * 70)
    _echo_styled("\nTo proceed, please type 'I agree' (exactly as shown):", level="warn")
    try:
        response = input().strip()
    except (KeyboardInterrupt, EOFError):
        click.echo("\nAcknowledgment interrupted. Exiting.", err=True)
        return False

    if response == "I agree":
        return True

    click.echo(
        "Incorrect acknowledgment. Please run the command again and type 'I agree' to acknowledge.",
        err=True,
    )
    return False


def _warn_if_process_large_bams() -> None:
    """If ROBIN_PROCESS_LARGE_BAMS is set, print a warning not to use with live runs."""
    if os.environ.get("ROBIN_PROCESS_LARGE_BAMS", "0").strip().lower() in ("1", "true", "yes", "on"):
        _echo_styled(
            "Warning: ROBIN_PROCESS_LARGE_BAMS is enabled. Do not use this option alongside live runs.",
            level="warn",
        )


@click.group()
@click.version_option()
def main() -> None:
    """Robin now uses Little John - his second in command who kept the merry men in line."""
    pass


def _iter_mgmt_bams(root: Path, recursive: bool) -> Iterable[Path]:
    """Yield mgmt_sorted.bam paths under a root directory."""
    if recursive:
        yield from root.rglob("mgmt_sorted.bam")
        return

    direct = root / "mgmt_sorted.bam"
    if direct.exists():
        yield direct

    for entry in root.iterdir():
        if not entry.is_dir():
            continue
        candidate = entry / "mgmt_sorted.bam"
        if candidate.exists():
            yield candidate


@main.group()
def password() -> None:
    """Manage the GUI login password (stored as a hash)."""
    pass


@password.command("set")
def password_set() -> None:
    """Set or replace the GUI password. Prompts twice for confirmation; input is never echoed."""
    try:
        from robin.gui_launcher import set_gui_password_interactive
    except ImportError as e:
        click.echo(f"GUI password module not available: {e}", err=True)
        sys.exit(1)
    if not set_gui_password_interactive():
        sys.exit(1)
    click.echo("GUI password set successfully.")


@main.group()
def utils() -> None:
    """Utility helpers for inspecting robin outputs."""
    pass


@utils.command()
@click.argument("output_dir", type=click.Path(exists=True, path_type=Path))
@click.option(
    "--recursive",
    "-r",
    is_flag=True,
    help="Search for mgmt_sorted.bam recursively under output_dir.",
)
@click.option(
    "--out",
    "-o",
    "out_path",
    type=click.Path(path_type=Path),
    default="-",
    help="Output TSV path (default: stdout).",
)
def mgmt(output_dir: Path, recursive: bool, out_path: Path) -> None:
    """Summarize MGMT CpG site methylation counts from mgmt_sorted.bam files."""
    if not output_dir.exists():
        click.echo(f"Error: output_dir not found: {output_dir}", err=True)
        sys.exit(1)

    bam_paths = list(_iter_mgmt_bams(output_dir, recursive))
    if not bam_paths:
        click.echo(
            f"No mgmt_sorted.bam files found under {output_dir} (recursive={recursive}).",
            err=True,
        )
        sys.exit(1)

    output_stream = sys.stdout if str(out_path) == "-" else open(out_path, "w", newline="")
    try:
        writer = csv.writer(output_stream, delimiter="\t")
        writer.writerow(
            [
                "sample_id",
                "run_path",
                "site",
                "chr",
                "pos",
                "cov_fwd",
                "cov_rev",
                "cov_total",
                "meth_fwd",
                "meth_rev",
                "meth_total",
                "meth_pct",
            ]
        )

        for bam_path in sorted(bam_paths):
            run_dir = bam_path.parent
            sample_id = run_dir.name
            rel_run_path = os.path.relpath(run_dir, output_dir)

            with tempfile.NamedTemporaryFile(
                suffix=".bed", delete=False
            ) as temp_bed:
                temp_bed_path = temp_bed.name

            try:
                run_matkit(str(bam_path), temp_bed_path)
                site_rows = extract_mgmt_site_rows_from_bed(temp_bed_path)
            except Exception as exc:
                click.echo(
                    f"Failed to process {bam_path}: {exc}",
                    err=True,
                )
                continue
            finally:
                try:
                    os.remove(temp_bed_path)
                except OSError:
                    pass

            if not site_rows:
                click.echo(
                    f"No MGMT site data found in {bam_path}",
                    err=True,
                )
                continue

            for row in site_rows:
                meth_fwd = int(row.get("meth_fwd", 0))
                meth_rev = int(row.get("meth_rev", 0))
                meth_total = meth_fwd + meth_rev
                cov_total = int(row.get("cov_total", 0))
                meth_pct = round((meth_total / cov_total) * 100.0, 2) if cov_total else 0.0
                site_label = str(row.get("site", "")).split(" ")[0] if row.get("site") else ""
                writer.writerow(
                    [
                        sample_id,
                        rel_run_path,
                        site_label,
                        row.get("chr", ""),
                        row.get("pos", ""),
                        row.get("cov_fwd", 0),
                        row.get("cov_rev", 0),
                        cov_total,
                        meth_fwd,
                        meth_rev,
                        meth_total,
                        meth_pct,
                    ]
                )
    finally:
        if output_stream is not sys.stdout:
            output_stream.close()


@utils.command("update-clinvar")
def update_clinvar() -> None:
    """Update ClinVar to the newest available NCBI version (best-effort)."""

    try:
        from robin.utils.clinvar_manager import update_clinvar_if_newer

        updated = update_clinvar_if_newer(download_if_missing=True)
        if updated:
            click.echo("ClinVar updated successfully.")
        else:
            click.echo("ClinVar is already up to date.")
    except Exception as e:
        click.echo(f"Failed to update ClinVar: {e}", err=True)
        sys.exit(1)


@utils.command("update-models")
@click.option(
    "--models-dir",
    type=click.Path(path_type=Path),
    default=None,
    help="Directory to place model files (defaults to ROBIN's models directory).",
)
@click.option(
    "--manifest",
    "manifest_path",
    type=click.Path(path_type=Path),
    default=None,
    help="Path to assets.json manifest (defaults to repo-root assets.json or ROBIN_ASSETS_MANIFEST).",
)
@click.option(
    "--overwrite",
    is_flag=True,
    help="Overwrite existing model files (default: skip existing).",
)
def update_models(models_dir: Optional[Path], manifest_path: Optional[Path], overwrite: bool) -> None:
    """Download/update ROBIN model files using the assets manifest."""
    try:
        from robin.utils.model_checker import get_models_directory
        from robin.utils.model_updater import update_models as _update_models

        target_dir = models_dir or get_models_directory()
        ok, msgs = _update_models(
            models_dir=target_dir,
            manifest_path=str(manifest_path) if manifest_path else None,
            overwrite=overwrite,
        )
        for m in msgs:
            click.echo(m)
        if not ok:
            sys.exit(1)
        click.echo("Model update complete.")
    except Exception as e:
        click.echo(f"Failed to update models: {e}", err=True)
        sys.exit(1)


def _remove_panel_from_system(panel_name: str) -> bool:
    """Remove a custom panel from ROBIN system."""
    try:
        # Check if it's a built-in panel
        built_in_panels = {"rCNS2", "AML"}
        if panel_name in built_in_panels:
            click.echo(f"Error: Cannot remove built-in panel '{panel_name}'", err=True)
            click.echo("Built-in panels (rCNS2, AML) cannot be removed.", err=True)
            return False
        
        # Get resources directory
        try:
            from robin import resources
            resources_dir = Path(resources.__file__).parent
        except ImportError:
            click.echo("Error: Could not locate ROBIN resources directory", err=True)
            return False
        
        # Check if panel exists
        panel_filename = f"{panel_name}_panel_name_uniq.bed"
        panel_path = resources_dir / panel_filename
        
        if not panel_path.exists():
            click.echo(f"Error: Panel '{panel_name}' not found", err=True)
            click.echo(f"Expected file: {panel_path}", err=True)
            return False
        
        # Confirm removal
        source_preview = resources_dir / panel_source_filename(panel_name)
        click.echo(f"Panel '{panel_name}' will be removed:")
        click.echo(f"  Processed BED: {panel_path}")
        click.echo(f"  Size: {panel_path.stat().st_size} bytes")
        if source_preview.exists():
            click.echo(f"  Original upload: {source_preview}")
            click.echo(f"  Size: {source_preview.stat().st_size} bytes")
        
        # Ask for confirmation
        try:
            confirm = input(f"\nAre you sure you want to remove panel '{panel_name}'? Type 'yes' to confirm: ").strip()
            if confirm.lower() != 'yes':
                click.echo("Panel removal cancelled.")
                return False
        except (KeyboardInterrupt, EOFError):
            click.echo("\nPanel removal cancelled.")
            return False
        
        # Remove processed BED and optional stored original upload
        panel_path.unlink()
        source_path = resources_dir / panel_source_filename(panel_name)
        if source_path.exists():
            source_path.unlink()
            click.echo(f"Removed file: {source_path}")

        click.echo(f"Panel '{panel_name}' removed successfully")
        click.echo(f"Removed file: {panel_path}")

        return True
        
    except Exception as e:
        click.echo(f"Error removing panel: {e}", err=True)
        return False


@main.command()
@click.argument("panel_name", type=str)
@click.option(
    "--force",
    "-f",
    is_flag=True,
    help="Skip confirmation prompt (use with caution)"
)
def remove_panel(panel_name: str, force: bool) -> None:
    """Remove a custom panel from ROBIN.
    
    PANEL_NAME: Name of the panel to remove
    
    Built-in panels (rCNS2, AML) cannot be removed.
    """
    if not _get_user_acknowledgment():
        sys.exit(1)
    
    # Validate panel name
    if not panel_name or not panel_name.strip():
        click.echo("Error: Panel name cannot be empty", err=True)
        sys.exit(1)
    
    panel_name = panel_name.strip()
    
    # Check if it's a built-in panel
    built_in_panels = {"rCNS2", "AML"}
    if panel_name in built_in_panels:
        click.echo(f"Error: Cannot remove built-in panel '{panel_name}'", err=True)
        click.echo("Built-in panels (rCNS2, AML) cannot be removed.", err=True)
        sys.exit(1)
    
    # Get resources directory
    try:
        from robin import resources
        resources_dir = Path(resources.__file__).parent
    except ImportError:
        click.echo("Error: Could not locate ROBIN resources directory", err=True)
        sys.exit(1)
    
    # Check if panel exists
    panel_filename = f"{panel_name}_panel_name_uniq.bed"
    panel_path = resources_dir / panel_filename
    source_path = resources_dir / panel_source_filename(panel_name)

    if not panel_path.exists():
        click.echo(f"Error: Panel '{panel_name}' not found", err=True)
        click.echo(f"Expected file: {panel_path}", err=True)
        click.echo("\nUse 'robin list-panels' to see available panels.", err=True)
        sys.exit(1)

    # Show panel info
    click.echo(f"Panel '{panel_name}' found:")
    click.echo(f"  Processed BED: {panel_path}")
    click.echo(f"  Size: {panel_path.stat().st_size} bytes")
    if source_path.exists():
        click.echo(f"  Original upload: {source_path}")
        click.echo(f"  Size: {source_path.stat().st_size} bytes")
    
    # Confirm removal unless --force is used
    if not force:
        try:
            confirm = input(f"\nAre you sure you want to remove panel '{panel_name}'? Type 'yes' to confirm: ").strip()
            if confirm.lower() != 'yes':
                click.echo("Panel removal cancelled.")
                return
        except (KeyboardInterrupt, EOFError):
            click.echo("\nPanel removal cancelled.")
            return
    
    # Remove processed BED and optional stored original
    try:
        panel_path.unlink()
        click.echo(f"\nRemoved file: {panel_path}")
        if source_path.exists():
            source_path.unlink()
            click.echo(f"Removed file: {source_path}")
        click.echo(f"\nPanel '{panel_name}' removed successfully!")
        click.echo("Use 'robin list-panels' to see remaining panels")
    except Exception as e:
        click.echo(f"Error removing panel file: {e}", err=True)
        sys.exit(1)


@main.command()
def list_panels() -> None:
    """List all available panels in ROBIN."""
    if not _get_user_acknowledgment():
        sys.exit(1)
    
    panels = _get_available_panels()
    
    click.echo("Available panels in ROBIN:\n")
    
    # Built-in panels
    built_in_panels = ["rCNS2", "AML"]
    custom_panels = [p for p in panels if p not in built_in_panels]
    
    click.echo("BUILT-IN PANELS:")
    for panel in built_in_panels:
        click.echo(f"  • {panel}")
    
    if custom_panels:
        click.echo("\nCUSTOM PANELS:")
        for panel in custom_panels:
            click.echo(f"  • {panel}")
    else:
        click.echo("\nCUSTOM PANELS:")
        click.echo("  (none)")
    
    click.echo(f"\nTotal panels: {len(panels)}")
    click.echo("\nUsage: Use --target-panel <panel_name> in workflow commands")
    click.echo("Example: robin workflow /path/to/bams --workflow mgmt,target --target-panel rCNS2")
    click.echo("\nPanel management:")
    click.echo("  • Add panel: robin add-panel <bed_file> <panel_name>")
    click.echo("  • Remove panel: robin remove-panel <panel_name>")


@main.command()
def list_job_types() -> None:
    """List all available job types organized by queue category."""
    if not _get_user_acknowledgment():
        sys.exit(1)

    job_types = {
        "preprocessing": ["preprocessing - Extract metadata from BAM files"],
        "bed_conversion": ["bed_conversion - Convert BAM files to BED format"],
        "mgmt": ["mgmt - MGMT methylation analysis"],
        "cnv": ["cnv - Copy number variation analysis"],
        "target": ["target - Target analysis"],
        "fusion": ["fusion - Fusion detection analysis"],
        "classification": [
            "sturgeon - Sturgeon classification analysis",
            "nanodx - NanoDX analysis",
            "pannanodx - PanNanoDX analysis",
        ],
        "slow": ["random_forest - Random Forest analysis"],
    }

    click.echo("Available job types in robin:\n")

    for queue, jobs in job_types.items():
        click.echo(f"{queue.upper()} QUEUE:")
        for job in jobs:
            click.echo(f"  • {job}")
        click.echo()

    click.echo("Usage examples:")
    click.echo(
        "  • Simplified format (recommended): 'mgmt,sturgeon' (bed_conversion auto-added)"
    )
    click.echo(
        "  • Full pipeline (simplified): 'mgmt,cnv,target,fusion,sturgeon,nanodx,pannanodx,random_forest' (bed_conversion auto-added)"
    )
    click.echo(
        "  • Legacy format with queue prefixes: 'preprocessing:bed_conversion,mgmt:mgmt,classification:sturgeon'"
    )
    click.echo(
        "\nNote: 'preprocessing' is automatically added as the first step if not specified."
    )
    click.echo(
        "Note: 'bed_conversion' is automatically added when needed for sturgeon, nanodx, pannanodx, or random_forest jobs."
    )
    click.echo(
        "Note: Each analysis type (mgmt, cnv, target, fusion) has its own queue for parallel processing."
    )
    click.echo(
        "Note: The system automatically determines the appropriate queue for each job type in simplified format."
    )


def _validate_bed_file(bed_path: Path) -> Tuple[bool, List[str]]:
    """Validate BED file format and return (is_valid, error_messages)."""
    errors = []
    
    if not bed_path.exists():
        errors.append(f"BED file does not exist: {bed_path}")
        return False, errors
    
    if not bed_path.is_file():
        errors.append(f"Path is not a file: {bed_path}")
        return False, errors
    
    try:
        with open(bed_path, 'r') as f:
            line_count = 0
            for line_num, line in enumerate(f, 1):
                line = line.strip()
                if not line or line.startswith('#'):
                    continue
                
                line_count += 1
                parts = line.split('\t')
                
                if len(parts) < 3:
                    errors.append(f"Line {line_num}: Invalid BED format - must have at least 3 columns (chromosome, start, end)")
                    continue
                
                # Validate chromosome
                chrom = parts[0]
                if not chrom.startswith('chr'):
                    errors.append(f"Line {line_num}: Chromosome must start with 'chr': {chrom}")
                
                # Validate start and end positions
                try:
                    start = int(parts[1])
                    end = int(parts[2])
                    if start < 0 or end < 0:
                        errors.append(f"Line {line_num}: Start and end positions must be non-negative")
                    if start >= end:
                        errors.append(f"Line {line_num}: Start position must be less than end position")
                except ValueError:
                    errors.append(f"Line {line_num}: Start and end positions must be integers")
                
                # Check for gene name in 4th column (required for ROBIN panels)
                if len(parts) < 4 or not parts[3].strip():
                    errors.append(f"Line {line_num}: Missing gene name in 4th column")
                
                # Optional: validate 6-column BED format if present
                if len(parts) >= 6:
                    # Validate score (5th column) - should be numeric or "."
                    score = parts[4].strip()
                    if score != ".":
                        try:
                            score_int = int(score)
                            if score_int < 0:
                                errors.append(f"Line {line_num}: Score must be non-negative")
                        except ValueError:
                            errors.append(f"Line {line_num}: Score must be an integer or '.'")
                    
                    # Validate strand (6th column) - should be + or -
                    strand = parts[5].strip()
                    if strand not in ['+', '-']:
                        errors.append(f"Line {line_num}: Strand must be '+' or '-', got: {strand}")
                
                # Limit error reporting to first 10 errors
                if len(errors) >= 10:
                    errors.append("... (additional errors truncated)")
                    break
            
            if line_count == 0:
                errors.append("BED file contains no valid data lines")
    
    except Exception as e:
        errors.append(f"Error reading BED file: {e}")
    
    return len(errors) == 0, errors


def _generate_unique_gene_bed(input_bed_path: Path, output_bed_path: Path) -> bool:
    """Generate unique gene BED file from input BED file."""
    try:
        import pandas as pd
        
        # Read the input BED file - handle both 4-column and 6-column BED formats
        try:
            # Try 6-column format first (chrom, start, end, gene, score, strand)
            df = pd.read_csv(
                input_bed_path,
                sep='\t',
                header=None,
                names=['chrom', 'start', 'end', 'gene', 'score', 'strand'],
                comment='#'
            )
        except ValueError:
            # Fallback to 4-column format (chrom, start, end, gene)
            df = pd.read_csv(
                input_bed_path,
                sep='\t',
                header=None,
                names=['chrom', 'start', 'end', 'gene'],
                comment='#'
            )
        
        # Process gene names - handle comma-separated genes
        processed_regions = []
        for _, row in df.iterrows():
            genes = [g.strip() for g in str(row['gene']).split(',')]
            for gene in genes:
                if gene and gene != 'nan':  # Skip empty gene names and NaN values
                    processed_regions.append({
                        'chrom': row['chrom'],
                        'start': int(row['start']),
                        'end': int(row['end']),
                        'gene': gene
                    })
        
        # Convert to DataFrame and remove duplicates
        processed_df = pd.DataFrame(processed_regions)
        
        # Remove duplicates based on gene name (keep first occurrence)
        processed_df = processed_df.drop_duplicates(subset=['gene'], keep='first')
        
        # Sort by chromosome and position
        processed_df = processed_df.sort_values(['chrom', 'start', 'end'])
        
        # Write to output file in standard 4-column BED format
        processed_df[['chrom', 'start', 'end', 'gene']].to_csv(
            output_bed_path,
            sep='\t',
            header=False,
            index=False
        )
        
        return True
        
    except Exception as e:
        click.echo(f"Error generating unique gene BED file: {e}", err=True)
        return False


def _get_available_panels() -> List[str]:
    """Get list of available panels from resources directory."""
    panels = ["rCNS2", "AML"]  # Built-in panels
    
    try:
        # Try to find the resources directory without importing robin module
        # Look for the resources directory relative to this file
        current_file = Path(__file__)
        resources_dir = current_file.parent.parent / "robin" / "resources"
        
        if resources_dir.exists():
            # Look for custom panels (files ending with _panel_name_uniq.bed)
            for bed_file in resources_dir.glob("*_panel_name_uniq.bed"):
                panel_name = bed_file.stem.replace("_panel_name_uniq", "")
                if panel_name not in panels:
                    panels.append(panel_name)
            
            panels.sort()
        
    except Exception:
        # Fallback to built-in panels only - don't fail on any import or other errors
        pass
    
    return panels


@utils.command("sequencing-files")
@click.option(
    "--panel",
    "-p",
    "panel",
    type=click.Choice(_get_available_panels()),
    required=True,
    help="Target gene panel (same as ROBIN workflow --target-panel).",
)
@click.option(
    "--reference",
    "-r",
    "reference",
    type=str,
    default=DEFAULT_GRCH38_REFERENCE_URL,
    show_default=True,
    help=(
        "Reference genome: HTTPS URL to download, or path to a local FASTA "
        "(e.g. .fa or .fa.gz)."
    ),
)
@click.option(
    "--output-dir",
    "-o",
    "output_dir",
    type=click.Path(path_type=Path),
    default=None,
    help="Directory for outputs (default: ./reference_files in the current working directory).",
)
@click.option(
    "--yes",
    "-y",
    is_flag=True,
    help="Proceed without interactive confirmation (for scripts).",
)
def sequencing_files(
    panel: str,
    reference: str,
    output_dir: Optional[Path],
    yes: bool,
) -> None:
    """Copy panel BED(s) and reference genome for sequencing.

    Copies the processed ``*_panel_name_uniq.bed`` and, when packaged, ``*_panel_source.bed``
    (original / unprocessed design) into one folder for adaptive sampling and aligners.
    The default reference is NCBI GRCh38 (no-alt analysis set, UCSC-style names), a large download.
    """
    out = (output_dir or (Path.cwd() / "reference_files")).resolve()
    bed_name = panel_bed_filename(panel)
    source_name = panel_source_filename(panel)
    has_source = panel_source_available(panel)

    click.echo("")
    click.echo("Planned actions:")
    click.echo(f"  Output directory: {out}")
    step = 1
    click.echo(
        f"  {step}. Processed panel BED: copy '{bed_name}' from ROBIN "
        "(unique genes; used for ROBIN analyses)."
    )
    step += 1
    if has_source:
        click.echo(
            f"  {step}. Source panel BED: copy '{source_name}' "
            "(unprocessed design; shipped with ROBIN or from `robin add-panel`)."
        )
        step += 1
    else:
        click.echo(
            f"  — No '{source_name}' in ROBIN for this panel "
            "(optional file; add it under `robin/resources` to include it here)."
        )
    click.echo(f"  {step}. Reference genome:")
    for line in describe_reference_action(reference, out).split("\n"):
        click.echo(f"     {line}")

    if not yes:
        click.echo("")
        if not click.confirm(
            "Proceed with copying the panel BED(s) and fetching the reference as described above?",
            default=False,
        ):
            raise click.Abort()

    try:
        bed_dest = copy_panel_bed_to(panel, out)
        click.echo(f"Copied processed panel BED to: {bed_dest}")
    except FileNotFoundError as e:
        click.echo(f"Error: {e}", err=True)
        sys.exit(1)

    src_dest = copy_panel_source_bed_if_present(panel, out)
    if src_dest:
        click.echo(f"Copied source panel BED to: {src_dest}")

    try:
        ref_dest = materialize_reference(reference, out)
        click.echo(f"Reference genome at: {ref_dest}")
    except (FileNotFoundError, RuntimeError) as e:
        click.echo(f"Error: {e}", err=True)
        sys.exit(1)

    click.echo("")
    click.echo(
        "Done. Use the source or processed BED as needed for adaptive sampling; "
        "use the reference FASTA for alignment."
    )


def _register_panel_in_system(panel_name: str, bed_path: Path) -> bool:
    """Register the new panel in ROBIN system by updating relevant files."""
    try:
        # Create a simple registration by copying the BED file to resources
        # and updating any necessary configuration files
        
        # For now, we'll just ensure the file is in the right location
        # The actual registration happens when the panel is referenced by name
        # in analysis modules like target_analysis.py and fusion_work.py
        
        click.echo(f"Panel '{panel_name}' registered successfully")
        click.echo(f"BED file location: {bed_path}")
        click.echo(f"Panel can now be used with --target-panel {panel_name}")
        
        return True
        
    except Exception as e:
        click.echo(f"Error registering panel: {e}", err=True)
        return False


@main.command()
@click.argument("bed_file", type=click.Path(exists=True, path_type=Path))
@click.argument("panel_name", type=str)
@click.option(
    "--validate-only",
    is_flag=True,
    help="Only validate the BED file format without adding the panel"
)
def add_panel(bed_file: Path, panel_name: str, validate_only: bool) -> None:
    """Add a custom panel to ROBIN.
    
    BED_FILE: Path to the BED file containing panel regions
    PANEL_NAME: Name for the panel (e.g., 'CustomPanel', 'MyPanel')
    
    The BED file should be in standard format with at least 4 columns:
    chromosome, start, end, gene_name(s) [, score, strand]
    
    Supported formats:
    - 4-column: chr1, 1000000, 2000000, GENE1
    - 6-column: chr1, 1000000, 2000000, GENE1, 0, +
    
    Gene names can be comma-separated for regions covering multiple genes.
    The output will be a unique gene list with one entry per gene.

    ROBIN also stores an unmodified copy of your upload as
    '{panel_name}_panel_source.bed' alongside the processed file.
    """
    if not _get_user_acknowledgment():
        sys.exit(1)
    
    # Validate panel name
    if not panel_name or not panel_name.strip():
        click.echo("Error: Panel name cannot be empty", err=True)
        sys.exit(1)
    
    panel_name = panel_name.strip()
    
    # Check for reserved panel names
    reserved_names = {"rCNS2", "AML"}
    if panel_name in reserved_names:
        click.echo(f"Error: Panel name '{panel_name}' is reserved. Please choose a different name.", err=True)
        sys.exit(1)
    
    click.echo(f"Validating BED file: {bed_file}")
    
    # Validate BED file format
    is_valid, errors = _validate_bed_file(bed_file)
    
    if not is_valid:
        click.echo("BED file validation failed:", err=True)
        for error in errors:
            click.echo(f"  • {error}", err=True)
        sys.exit(1)
    
    click.echo("BED file format validation passed")
    
    if validate_only:
        click.echo("Validation complete. Use without --validate-only to add the panel.")
        return
    
    # Generate output paths
    try:
        from robin import resources
        resources_dir = Path(resources.__file__).parent
    except ImportError:
        click.echo("Error: Could not locate ROBIN resources directory", err=True)
        sys.exit(1)
    
    output_filename = f"{panel_name}_panel_name_uniq.bed"
    output_path = resources_dir / output_filename
    source_path = resources_dir / panel_source_filename(panel_name)

    # Check if panel already exists
    if output_path.exists():
        click.echo(f"Error: Panel '{panel_name}' already exists at {output_path}", err=True)
        click.echo("Please choose a different panel name or remove the existing panel first.", err=True)
        sys.exit(1)
    if source_path.exists():
        click.echo(f"Error: File already exists: {source_path}", err=True)
        click.echo("Remove it or choose a different panel name.", err=True)
        sys.exit(1)

    click.echo(f"Storing original BED (unprocessed): {source_path}")
    try:
        shutil.copy2(bed_file, source_path)
    except OSError as e:
        click.echo(f"Error copying original BED file: {e}", err=True)
        sys.exit(1)

    click.echo(f"Generating unique gene BED file: {output_path}")

    # Generate unique gene BED file
    if not _generate_unique_gene_bed(bed_file, output_path):
        click.echo("Failed to generate unique gene BED file", err=True)
        try:
            source_path.unlink(missing_ok=True)
        except OSError:
            pass
        sys.exit(1)

    click.echo("Unique gene BED file generated successfully")

    # Register panel in system
    if not _register_panel_in_system(panel_name, output_path):
        click.echo("Failed to register panel in system", err=True)
        try:
            output_path.unlink(missing_ok=True)
            source_path.unlink(missing_ok=True)
        except OSError:
            pass
        sys.exit(1)

    click.echo(f"\nPanel '{panel_name}' added successfully!")
    click.echo(f"  Processed BED (unique genes): {output_path}")
    click.echo(f"  Original upload: {source_path}")
    click.echo(f"Usage: Use --target-panel {panel_name} in workflow commands")
    click.echo(f"Example: robin workflow /path/to/bams --workflow mgmt,target --target-panel {panel_name}")
    click.echo(f"Remove: robin remove-panel {panel_name}")


def _parse_job_log_levels(job_log_level: Tuple[str, ...]) -> Dict[str, str]:
    """Parse job-specific log level specifications."""
    job_levels = {}

    for job_level_spec in job_log_level:
        if ":" in job_level_spec:
            job_type, level = job_level_spec.split(":", 1)
            job_type = job_type.strip()
            level = level.strip().upper()

            if level not in VALID_LOG_LEVELS:
                click.echo(
                    f"Warning: Invalid log level '{level}' for job '{job_type}'. Valid levels: {', '.join(VALID_LOG_LEVELS)}",
                    err=True,
                )
                continue

            job_levels[job_type] = level
        else:
            click.echo(
                f"Warning: Invalid job log level format '{job_level_spec}'. Use format 'job_type:level'",
                err=True,
            )
    return job_levels


def _parse_command_mappings(commands: Tuple[str, ...]) -> Dict[str, str]:
    """Parse command mapping specifications."""
    command_map = {}
    for cmd_mapping in commands:
        if ":" in cmd_mapping:
            job_type, command = cmd_mapping.split(":", 1)
            job_type = job_type.strip()
            command = command.strip()

            if not job_type or not command:
                click.echo(
                    f"Warning: Empty job type or command in '{cmd_mapping}'", err=True
                )
                continue

            command_map[job_type] = command
        else:
            click.echo(
                f"Warning: Invalid command mapping format '{cmd_mapping}'. Use format 'job_type:command'",
                err=True,
            )
    return command_map


def _validate_workflow_steps(workflow_steps: List[str]) -> List[str]:
    """Validate and clean workflow steps."""
    cleaned_steps = []
    for step in workflow_steps:
        step = step.strip()
        if not step:
            continue

        if ":" in step:
            queue_type, job_type = step.split(":", 1)
            queue_type = queue_type.strip()
            job_type = job_type.strip()

            if job_type not in VALID_JOB_TYPES:
                click.echo(
                    f"Warning: Unknown job type '{job_type}' in step '{step}'", err=True
                )
                continue

            cleaned_steps.append(f"{queue_type}:{job_type}")
        else:
            # Simplified format - job type only
            if step not in VALID_JOB_TYPES:
                click.echo(f"Warning: Unknown job type '{step}'", err=True)
                continue
            cleaned_steps.append(step)

    return cleaned_steps


def _ensure_preprocessing_step(workflow_steps: List[str]) -> List[str]:
    """Ensure preprocessing is the first step in the workflow."""
    if not workflow_steps or not workflow_steps[0].endswith(":preprocessing"):
        workflow_steps.insert(0, "preprocessing:preprocessing")
    return workflow_steps


def _convert_simplified_workflow(job_types: List[str]) -> List[str]:
    """Convert simplified job types to standard workflow format with queue prefixes."""
    # Check if any jobs require bed_conversion as a dependency
    needs_bed_conversion = any(
        job in JOBS_REQUIRING_BED_CONVERSION for job in job_types
    )

    # Add bed_conversion if needed and not already present
    if needs_bed_conversion and "bed_conversion" not in job_types:
        job_types = ["bed_conversion"] + job_types

    converted_steps = []
    for job_type in job_types:
        queue_type = QUEUE_MAPPING.get(job_type)
        if queue_type:
            converted_steps.append(f"{queue_type}:{job_type}")
        else:
            # Unknown job type, default to slow queue
            converted_steps.append(f"slow:{job_type}")

    return converted_steps


def _create_ray_workflow_runner(
    verbose: bool,
    analysis_workers: int,
    legacy_analysis_queue: bool,
    log_level: str,
    target_panel: str,
    preprocessing_workers: int = 1,
    bed_workers: int = 1,
    reference: Optional[Path] = None,
) -> Any:
    """Create and configure Ray-based workflow runner (Ray Core)."""
    try:
        # Prefer Ray Core implementation
        import asyncio
        from robin import workflow_ray as wrn

        class _RayCoreWrapper:
            def __init__(self, reference: Optional[Path] = None, target_panel: str = None):
                self.manager = type(
                    "_DummyManager", (), {"get_priority_info": lambda _self: {}}
                )()
                self.coordinator = None  # Will be set when workflow starts
                self.reference = reference  # Store reference genome for GUI access
                self.target_panel = target_panel  # Store target panel for GUI access

                # Debug logging for reference genome (only in verbose mode)
                if self.reference:
                    pass  # Reference genome stored successfully
                else:
                    pass  # No reference genome stored

            def register_handler(self, *args, **kwargs):
                # Handlers are wired inside workflow_ray; keep API parity
                return None

            def register_command_handler(self, *args, **kwargs):
                return None

            def submit_sample_job(
                self,
                sample_dir: str,
                job_type: str,
                sample_id: str = None,
                force_regenerate: bool = False,
            ) -> bool:
                """Submit a job for an existing sample directory using Ray workflow."""
                try:
                    # Try to get coordinator with retries
                    import time
                    from robin import workflow_ray as wrn

                    max_retries = 5
                    retry_delay = 1.0

                    for attempt in range(max_retries):
                        # Try to get coordinator
                        if self.coordinator is None:
                            try:
                                self.coordinator = wrn.get_coordinator_sync()
                            except Exception:
                                pass

                        if self.coordinator is not None:
                            break

                        if attempt < max_retries - 1:
                            time.sleep(retry_delay)
                            retry_delay *= 1.5  # Exponential backoff

                    if self.coordinator is None:
                        return False

                    # Submit the job using the coordinator
                    # Ray remote calls are synchronous from the caller's perspective
                    result = self.coordinator.submit_sample_job.remote(
                        sample_dir, job_type, sample_id, force_regenerate
                    )

                    # Wait for the result (this is the actual async part)
                    import ray

                    try:
                        final_result = ray.get(result, timeout=30.0)
                        return final_result
                    except ray.exceptions.GetTimeoutError:
                        # Submission may still have reached the coordinator; avoid false-negative UI status.
                        return True
                    except Exception:
                        return False

                except Exception:
                    return False

            def submit_snp_analysis_job(
                self,
                sample_dir: str,
                sample_id: str = None,
                reference: str = None,
                threads: int = 4,
                force_regenerate: bool = False,
            ) -> bool:
                """Submit a SNP analysis job for an existing sample directory using Ray workflow."""
                try:
                    sid_for_log = sample_id or Path(sample_dir).name
                    click.echo(
                        f"[SNP] Submitting SNP analysis for sample '{sid_for_log}'..."
                    )
                    # Try to get coordinator with retries
                    import time
                    from robin import workflow_ray as wrn

                    max_retries = 5
                    retry_delay = 1.0

                    for attempt in range(max_retries):
                        # Try to get coordinator
                        if self.coordinator is None:
                            try:
                                self.coordinator = wrn.get_coordinator_sync()
                            except Exception:
                                pass

                        if self.coordinator is not None:
                            break

                        if attempt < max_retries - 1:
                            time.sleep(retry_delay)
                            retry_delay *= 1.5  # Exponential backoff

                    if self.coordinator is None:
                        click.echo(
                            f"[SNP] Failed to get coordinator for sample '{sid_for_log}'.",
                            err=True,
                        )
                        return False

                    # Submit the SNP analysis job using the coordinator
                    # Ray remote calls are synchronous from the caller's perspective
                    result = self.coordinator.submit_snp_analysis_job.remote(
                        sample_dir, sample_id, reference, threads, force_regenerate
                    )

                    # Wait for the result (this is the actual async part)
                    import ray

                    try:
                        final_result = ray.get(result, timeout=30.0)
                        if final_result:
                            click.echo(
                                f"[SNP] Queued SNP analysis for sample '{sid_for_log}'."
                            )
                        else:
                            click.echo(
                                f"[SNP] Submission returned False for sample '{sid_for_log}'.",
                                err=True,
                            )
                        return final_result
                    except ray.exceptions.GetTimeoutError:
                        click.echo(
                            f"[SNP] Submit confirmation timed out for sample '{sid_for_log}', treating as queued.",
                            err=True,
                        )
                        return True
                    except Exception as e:
                        click.echo(
                            f"[SNP] Failed waiting for submit confirmation for sample '{sid_for_log}': {e}",
                            err=True,
                        )
                        return False

                except Exception as e:
                    click.echo(
                        f"[SNP] Unexpected submission error for sample '{sample_id or Path(sample_dir).name}': {e}",
                        err=True,
                    )
                    return False

            def submit_target_bam_finalize_job(
                self,
                sample_dir: str,
                sample_id: str = None,
                target_panel: str = None,
                reference: str = None,
            ) -> bool:
                """Submit a target BAM finalization job for an existing sample directory using Ray workflow."""
                try:
                    sid_for_log = sample_id or Path(sample_dir).name
                    click.echo(
                        f"[Finalize] Submitting target BAM finalization for sample '{sid_for_log}'..."
                    )
                    import time
                    from robin import workflow_ray as wrn

                    max_retries = 5
                    retry_delay = 1.0

                    for attempt in range(max_retries):
                        if self.coordinator is None:
                            try:
                                self.coordinator = wrn.get_coordinator_sync()
                            except Exception:
                                pass

                        if self.coordinator is not None:
                            break

                        if attempt < max_retries - 1:
                            time.sleep(retry_delay)
                            retry_delay *= 1.5

                    if self.coordinator is None:
                        click.echo(
                            f"[Finalize] Failed to get coordinator for sample '{sid_for_log}'.",
                            err=True,
                        )
                        return False

                    result = self.coordinator.submit_target_bam_finalize_job.remote(
                        sample_dir, sample_id, target_panel, reference
                    )

                    import ray

                    try:
                        final_result = ray.get(result, timeout=30.0)
                        if final_result:
                            click.echo(
                                f"[Finalize] Queued target BAM finalization for sample '{sid_for_log}'."
                            )
                        else:
                            click.echo(
                                f"[Finalize] Submission returned False for sample '{sid_for_log}'.",
                                err=True,
                            )
                        return final_result
                    except ray.exceptions.GetTimeoutError:
                        click.echo(
                            f"[Finalize] Submit confirmation timed out for sample '{sid_for_log}', treating as queued.",
                            err=True,
                        )
                        return True
                    except Exception as e:
                        click.echo(
                            f"[Finalize] Failed waiting for submit confirmation for sample '{sid_for_log}': {e}",
                            err=True,
                        )
                        return False

                except Exception as e:
                    click.echo(
                        f"[Finalize] Unexpected submission error for sample '{sample_id or Path(sample_dir).name}': {e}",
                        err=True,
                    )
                    return False

            def is_sample_ready_for_snp_analysis(
                self, sample_dir: str
            ) -> tuple[bool, list[str]]:
                """Check if a sample directory is ready for SNP analysis."""
                import os

                required_files = ["target.bam", "targets_exceeding_threshold.bed"]

                missing_files = []
                for filename in required_files:
                    file_path = os.path.join(sample_dir, filename)
                    if not os.path.exists(file_path):
                        missing_files.append(filename)

                is_ready = len(missing_files) == 0
                return is_ready, missing_files

            def run_workflow(
                self,
                path: Path,
                workflow_steps: List[str],
                no_process_existing: bool,
                no_progress: bool,
                no_watch: bool,
                log_level_local: str,
                analysis_workers_local: int,
                preset: Optional[str] = None,
            ) -> None:
                # Store reference to self for coordinator access
                self_ref = self

                async def _run_with_coordinator():

                    # Run the workflow first to create the coordinator
                    await wrn.run(
                        plan=workflow_steps,
                        paths=[str(path)],
                        target_panel=target_panel,
                        analysis_workers=analysis_workers_local,
                        process_existing=not no_process_existing,
                        monitor=not no_progress,
                        watch=(not no_watch),
                        patterns=["*.bam"],
                        ignore_patterns=None,
                        recursive=True,
                        work_dir=None,
                        log_level=log_level_local,
                        preset=preset,
                        workflow_runner=self_ref,
                        reference=str(reference) if reference else None,
                    )

                    # After the workflow starts, try to get the coordinator reference
                    try:
                        coord = wrn.get_coordinator_sync()
                        if coord:
                            self_ref.coordinator = coord
                    except Exception:
                        pass

                asyncio.run(_run_with_coordinator())

        return _RayCoreWrapper(reference=reference, target_panel=target_panel)
    except ImportError as e:
        click.echo(
            f"Warning: Ray Core not available ({e}). Falling back to threading-based workflow.",
            err=True,
        )
        return None


def _initialize_ray(num_cpus: Optional[int], include_dashboard: bool = True) -> None:
    """Initialize Ray with specified CPU count and dashboard option."""
    if num_cpus is not None:
        if num_cpus <= 0:
            click.echo(
                f"Warning: Invalid CPU count {num_cpus}. Must be positive.", err=True
            )
            return

        import ray

        if not ray.is_initialized():
            # Suppress Ray logging to prevent interference with progress bars
            import logging
            import json

            ray_logger = logging.getLogger("ray")
            ray_logger.setLevel(logging.ERROR)
            logging.getLogger("ray.worker").setLevel(logging.ERROR)
            logging.getLogger("ray.remote").setLevel(logging.ERROR)
            logging.getLogger("ray.actor").setLevel(logging.ERROR)
            logging.getLogger("ray.util").setLevel(logging.ERROR)

            try:
                # Set Ray environment variables to reduce verbose output
                os.environ["RAY_DISABLE_IMPORT_WARNING"] = "1"
                os.environ["RAY_DISABLE_DEPRECATION_WARNING"] = "1"
                os.environ["RAY_OBJECT_STORE_ALLOW_SLOW_STORAGE"] = "1"
                # Raise dashboard/State API job list limit (default 10k); set only if not already set
                if "RAY_MAX_LIMIT_FROM_DATA_SOURCE" not in os.environ:
                    os.environ["RAY_MAX_LIMIT_FROM_DATA_SOURCE"] = "100000"

                # Bind dashboard to 0.0.0.0 when supported so it's reachable off-host
                init_kwargs = {
                    "num_cpus": num_cpus,
                    "ignore_reinit_error": True,
                    "object_store_memory": 1000000000,  # 1GB object store
                    "temp_dir": "/tmp/ray",  # Use /tmp for temporary files
                    "_system_config": {
                        "object_spilling_config": json.dumps(
                            {
                                "type": "filesystem",
                                "params": {"directory_path": "/tmp/ray/spill"},
                            }
                        ),
                        "max_direct_call_object_size": 1000000,  # 1MB
                        "object_store_full_delay_ms": 100,
                        "object_store_full_max_retries": 0,
                    },
                }

                if include_dashboard:
                    init_kwargs.update(
                        {
                            "include_dashboard": True,
                            "dashboard_host": os.environ.get(
                                "RAY_DASHBOARD_HOST", "0.0.0.0"
                            ),
                        }
                    )

                try:
                    ray.init(**init_kwargs)
                except TypeError:
                    # Older Ray versions may not support dashboard args
                    # Remove dashboard-specific args and retry
                    init_kwargs.pop("include_dashboard", None)
                    init_kwargs.pop("dashboard_host", None)
                    ray.init(**init_kwargs)
                pass  # Ray initialized successfully
            except Exception as e:
                click.echo(f"Warning: Failed to initialize Ray: {e}", err=True)
        else:
            pass  # Ray already initialized


def _configure_queue_priorities(runner: Any, queue_priority: Tuple[str, ...]) -> None:
    """Configure queue priorities for Ray workflow."""
    if not queue_priority:
        return

    click.echo("Configuring queue priorities:")
    for priority_spec in queue_priority:
        if ":" in priority_spec:
            queue_type, priority_str = priority_spec.split(":", 1)
            queue_type = queue_type.strip()
            priority_str = priority_str.strip()

            try:
                priority = int(priority_str)
                if priority < 0:
                    click.echo(
                        f"Warning: Priority {priority} for queue '{queue_type}' is negative. Using 0 instead.",
                        err=True,
                    )
                    priority = 0

                runner.manager.set_queue_priority(queue_type, priority)
                click.echo(f"  {queue_type}: {priority}")
            except ValueError:
                click.echo(
                    f"Warning: Invalid priority value '{priority_str}' for queue '{queue_type}'. Must be an integer.",
                    err=True,
                )
        else:
            click.echo(
                f"Warning: Invalid priority specification '{priority_spec}'. Use format 'queue:priority'",
                err=True,
            )


def _create_workflow_runner(
    use_ray: bool,
    verbose: bool,
    analysis_workers: int,
    legacy_analysis_queue: bool,
    log_level: str,
    target_panel: str,
    preprocessing_workers: int = 1,
    bed_workers: int = 1,
    reference: Optional[Path] = None,
    center: str = None,
) -> Any:
    """Create the appropriate workflow runner based on configuration."""
    if analysis_workers < 1:
        click.echo(
            f"Warning: Invalid analysis_workers value {analysis_workers}. Using {DEFAULT_ANALYSIS_WORKERS} instead.",
            err=True,
        )
        analysis_workers = DEFAULT_ANALYSIS_WORKERS

    if use_ray:
        runner = _create_ray_workflow_runner(
            verbose,
            analysis_workers,
            legacy_analysis_queue,
            log_level,
            target_panel,  # Pass target_panel parameter to Ray workflow runner
            preprocessing_workers=preprocessing_workers,
            bed_workers=bed_workers,
            reference=reference,  # Pass reference parameter to Ray workflow runner
        )
        if runner is None:
            # Fallback to threading-based workflow
            from robin.workflow_simple import WorkflowRunner

            runner = WorkflowRunner(
                target_panel=target_panel,
                verbose=verbose,
                analysis_workers=analysis_workers,
                use_separate_analysis_queues=not legacy_analysis_queue,
                preprocessing_workers=preprocessing_workers,
                bed_workers=bed_workers,
                reference=reference,
                center=center,
            )
        return runner
    else:
        from robin.workflow_simple import WorkflowRunner

        return WorkflowRunner(
            target_panel=target_panel,
            verbose=verbose,
            analysis_workers=analysis_workers,
            use_separate_analysis_queues=not legacy_analysis_queue,
            preprocessing_workers=preprocessing_workers,
            bed_workers=bed_workers,
            reference=reference,
            center=center,
        )


def _register_handlers(
    runner: Any,
    legacy_analysis_queue: bool,
    work_dir: Optional[Path],
    target_panel: str,
    reference: Optional[Path] = None,
    center: str = None,
) -> None:
    """Register all workflow handlers with the runner."""
    
    # Validate target panel
    available_panels = _get_available_panels()
    if target_panel not in available_panels:
        click.echo(
            f"Warning: Target panel '{target_panel}' is not in detected panels: {available_panels}. "
            f"This may be a custom panel - proceeding with specified panel.",
            err=True,
        )
        # Don't reset to rCNS2 - allow custom panels to be used
    else:
        click.echo(f"Using target panel: {target_panel}")
    
    # Track handlers that should accept target_panel
    handlers_requiring_panel = {"target", "fusion", "cnv"}
    registered_panel_handlers = set()
    
    for (
        queue_type,
        job_type,
        handler_func,
        legacy_queue_type,
        needs_work_dir,
    ) in HANDLER_CONFIGS:
        # Determine the actual queue type based on legacy mode
        actual_queue_type = (
            legacy_queue_type
            if legacy_analysis_queue and legacy_queue_type
            else queue_type
        )

        # Create handler with work directory and/or reference genome if specified and needed
        if work_dir and needs_work_dir:
            if reference and job_type == "target":
                # Special handling for target analysis with reference genome and target panel
                def create_handler_with_work_dir_and_ref(
                    handler, work_dir_path, ref_path, center_param, panel_param
                ):
                    return lambda job: handler(
                        job, work_dir=str(work_dir_path), reference=str(ref_path), target_panel=job.context.metadata.get("target_panel", panel_param)
                    )

                final_handler = create_handler_with_work_dir_and_ref(
                    handler_func, work_dir, reference, center, target_panel
                )
            elif reference and job_type == "mgmt":
                # Special handling for MGMT analysis with reference genome
                def create_mgmt_handler_with_work_dir_and_ref(
                    handler, work_dir_path, ref_path
                ):
                    return lambda job: handler(
                        job, work_dir=str(work_dir_path), reference=str(ref_path)
                    )

                final_handler = create_mgmt_handler_with_work_dir_and_ref(
                    handler_func, work_dir, reference
                )
            elif job_type == "fusion":
                # Special handling for fusion analysis with target panel
                def create_fusion_handler_with_work_dir(
                    handler, work_dir_path, panel_param
                ):
                    return lambda job: handler(
                        job, work_dir=str(work_dir_path), target_panel=job.context.metadata.get("target_panel", panel_param)
                    )

                final_handler = create_fusion_handler_with_work_dir(
                    handler_func, work_dir, target_panel
                )
            else:
                # Standard work directory handling
                if job_type in ["target", "cnv"]:
                    # Analysis with work_dir and target_panel
                    def create_analysis_handler_with_work_dir(handler, work_dir_path, panel_param):
                        return lambda job: handler(job, work_dir=str(work_dir_path), target_panel=job.context.metadata.get("target_panel", panel_param))
                    
                    final_handler = create_analysis_handler_with_work_dir(handler_func, work_dir, target_panel)
                else:
                    # Standard work directory handling for other job types
                    def create_handler_with_work_dir(handler, work_dir_path, center_param):
                        return lambda job: handler(job, work_dir=str(work_dir_path))

                    final_handler = create_handler_with_work_dir(handler_func, work_dir, center)
        elif reference and job_type == "target":
            # Reference genome only (no work_dir needed) with target panel
            def create_handler_with_ref(handler, ref_path, center_param, panel_param):
                return lambda job: handler(job, reference=str(ref_path), target_panel=job.context.metadata.get("target_panel", panel_param))

            final_handler = create_handler_with_ref(handler_func, reference, center, target_panel)
        elif job_type in ["fusion", "cnv"]:
            # Analysis with target panel only (no work_dir needed)
            def create_analysis_handler_with_panel(handler, panel_param):
                return lambda job: handler(job, target_panel=job.context.metadata.get("target_panel", panel_param))

            final_handler = create_analysis_handler_with_panel(handler_func, target_panel)
        elif job_type == "preprocessing":
            # Special handling for preprocessing to pass center
            def create_preprocessing_handler(handler, center_param):
                return lambda job: handler(job, center=center_param)

            final_handler = create_preprocessing_handler(handler_func, center)
        else:
            final_handler = handler_func

        # Track handlers that require target_panel
        if job_type in handlers_requiring_panel:
            import inspect
            sig = inspect.signature(handler_func)
            if "target_panel" in sig.parameters:
                registered_panel_handlers.add(job_type)
            else:
                click.echo(
                    f"Warning: Handler for {job_type} does not accept target_panel parameter. "
                    f"Panel information will not be passed to this handler.",
                    err=True,
                )

        try:
            runner.register_handler(actual_queue_type, job_type, final_handler)
        except Exception as e:
            click.echo(
                f"Warning: Failed to register handler for {queue_type}:{job_type}: {e}",
                err=True,
            )
    
    # Report on panel handler registration
    missing_panel_handlers = handlers_requiring_panel - registered_panel_handlers
    if missing_panel_handlers:
        click.echo(
            f"Warning: The following handlers do not support target_panel parameter: {missing_panel_handlers}",
            err=True,
        )
    else:
        click.echo(f"Successfully registered panel-aware handlers: {registered_panel_handlers}")


def _register_command_handlers(
    runner: Any, command_map: Dict[str, str], legacy_analysis_queue: bool
) -> None:
    """Register command handlers with the runner."""
    for job_type, command in command_map.items():
        if ":" in job_type:
            queue_type, actual_job_type = job_type.split(":", 1)

            # Handle legacy analysis queue mapping
            if legacy_analysis_queue and queue_type in [
                "mgmt",
                "cnv",
                "target",
                "fusion",
            ]:
                actual_queue_type = "analysis"
            else:
                actual_queue_type = queue_type

            try:
                runner.register_command_handler(
                    actual_queue_type, actual_job_type, command
                )
            except Exception as e:
                click.echo(
                    f"Warning: Failed to register command handler for {queue_type}:{actual_job_type}: {e}",
                    err=True,
                )
        else:
            # Default to preprocessing queue
            try:
                runner.register_command_handler("preprocessing", job_type, command)
            except Exception as e:
                click.echo(
                    f"Warning: Failed to register command handler for preprocessing:{job_type}: {e}",
                    err=True,
                )


def _create_classifier_with_work_dir(work_dir: Path, workflow_steps: List[str], target_panel: str):
    """Create a classifier function that includes work directory in job context."""

    def classifier_with_work_dir(filepath: str) -> List[Job]:
        jobs = default_file_classifier(filepath, workflow_steps, target_panel)
        for job in jobs:
            job.context.add_metadata("work_dir", str(work_dir))
        return jobs

    return classifier_with_work_dir


def _validate_inputs(
    path: Path, workflow: str, analysis_workers: int, ray_num_cpus: Optional[int]
) -> None:
    """Validate input parameters and provide helpful error messages."""
    if not path.exists():
        raise click.BadParameter(f"Path '{path}' does not exist")

    if not path.is_dir():
        raise click.BadParameter(f"Path '{path}' is not a directory")

    if not workflow.strip():
        raise click.BadParameter("Workflow cannot be empty")

    if analysis_workers < 1:
        raise click.BadParameter("analysis_workers must be at least 1")

    if ray_num_cpus is not None and ray_num_cpus < 1:
        raise click.BadParameter("ray_num_cpus must be at least 1")


def _display_workflow_config(
    path: Path,
    center: str,
    work_dir: Optional[Path],
    workflow_steps: List[str],
    command_map: dict,
    log_level: str,
    job_levels: dict,
    deduplicate_jobs: tuple,
    legacy_analysis_queue: bool,
    analysis_workers: int,
    no_process_existing: bool,
    uses_simplified_format: bool = False,
    original_workflow: str = "",
    use_ray: bool = False,
    ray_num_cpus: Optional[int] = None,
    queue_priority: tuple = (),
) -> None:
    """Display workflow configuration information."""
    _echo_styled(f"Center: {center}", level="info")
    
    if no_process_existing:
            _echo_styled(
            f"Starting workflow on {path} for BAM files (skipping existing files)..."
            )
    else:
            _echo_styled(
            f"Starting workflow on {path} for BAM files (will process existing files first)..."
            )

    if work_dir:
        _echo_styled(f"Output directory: {work_dir}", level="info")
    else:
        _echo_styled(
            "Output directory: Not specified (using input directory)", level="warn"
        )

    if uses_simplified_format:
        # Extract job types from the final workflow steps
        final_job_types = [step.split(":")[1] for step in workflow_steps if ":" in step]
        _echo_styled(f"Workflow plan (simplified): {final_job_types}", level="info")
        _echo_styled(f"Auto-assigned queues: {workflow_steps}", level="info")

        # Check if bed_conversion was auto-added
        if original_workflow and isinstance(original_workflow, str):
            original_jobs = [step.strip() for step in original_workflow.split(",")]
            if (
                "bed_conversion" in final_job_types
                and "bed_conversion" not in original_jobs
            ):
                _echo_styled(
                    "Note: bed_conversion was automatically added as it's required for other jobs",
                    level="warn",
                )
    else:
        _echo_styled(f"Workflow plan: {workflow_steps}", level="info")

    _echo_styled(f"Commands: {command_map}", level="info")
    _echo_styled(f"Log_level: {log_level}", level="info")

    if job_levels:
        click.echo(f"Job log levels: {job_levels}")
    if deduplicate_jobs:
        click.echo(f"Job deduplication: {list(deduplicate_jobs)}")

    # Display Ray configuration if enabled
    if use_ray:
        _echo_styled("Distributed computing: Ray (experimental)", level="info")
        if ray_num_cpus:
            click.echo(f"  - Ray CPUs: {ray_num_cpus}")
        else:
            click.echo("  - Ray CPUs: auto-detect")
        click.echo(
            f"  - Analysis workers per type: {analysis_workers} (distributed across Ray cluster)"
        )
        click.echo(f"  - Log level: {log_level} (applied to all Ray actors)")

        # Display priority configuration
        if queue_priority:
            click.echo(f"  - Queue priorities: {list(queue_priority)}")
    else:
        _echo_styled(
            "Distributed computing: Disabled (using threading)", level="warn"
        )
        click.echo("Worker configuration:")
        if legacy_analysis_queue:
            click.echo(
                "  - Analysis queue mode: Legacy (single queue for all analysis types)"
            )
            click.echo(f"  - Analysis workers: {analysis_workers}")
        else:
            click.echo("  - Analysis queue mode: Separate queues per analysis type")
            click.echo(
                f"  - Analysis workers per type: {analysis_workers} (MGMT, CNV, Target, Fusion each get {analysis_workers} workers)"
            )

        click.echo("  - Other queues: 1 worker each (fixed)")

    _echo_styled("Press Ctrl+C to stop", level="warn")
    _echo_styled("Running the workflow...", level="success")


@main.command()
@click.argument("path", type=click.Path(exists=True, path_type=Path))
@click.option(
    "--workflow",
    "-w",
    required=True,
    help="Workflow plan. Can be specified in two formats:\n1. With queue prefixes: 'preprocessing:bed_conversion,mgmt:mgmt,classification:sturgeon'\n2. Simplified (auto-queue): 'mgmt,sturgeon' - system automatically determines appropriate queue for each job type and adds bed_conversion when needed",
)
@click.option(
    "--center",
    required=True,
    help="Center ID running the analysis (e.g., 'Sherwood', 'Auckland', 'New York')",
)
@click.option(
    "--commands",
    "-c",
    multiple=True,
    help="Command mappings (e.g., 'index:samtools index {file}')",
)
@click.option(
    "--verbose",
    "-v",
    is_flag=True,
    help="Enable verbose output and detailed error traces",
)
@click.option(
    "--no-process-existing",
    is_flag=True,
    help="Skip processing existing files, only watch for new changes",
)
@click.option(
    "--work-dir",
    "-d",
    type=click.Path(path_type=Path),
    help="Base output directory for analysis results",
)
@click.option(
    "--log-level",
    default=DEFAULT_LOG_LEVEL,
    type=click.Choice(list(VALID_LOG_LEVELS)),
    help=f"Global log level (default: {DEFAULT_LOG_LEVEL})",
)
@click.option(
    "--job-log-level",
    multiple=True,
    help=f"Set log level for specific job (e.g., 'preprocessing:DEBUG', 'mgmt:WARNING'). Valid levels: {', '.join(VALID_LOG_LEVELS)}",
)
@click.option(
    "--deduplicate-jobs",
    multiple=True,
    help="Job types to deduplicate by sample ID (e.g., 'sturgeon', 'mgmt'). Jobs of these types will only run once per sample, even if multiple upstream jobs complete simultaneously.",
)
@click.option(
    "--no-progress", is_flag=True, help="Disable progress bars for file processing"
)
@click.option(
    "--analysis-workers",
    type=int,
    default=DEFAULT_ANALYSIS_WORKERS,
    help=f"Number of analysis workers per analysis queue (default: {DEFAULT_ANALYSIS_WORKERS})",
)
@click.option(
    "--preprocessing-workers",
    type=int,
    default=1,
    help="Number of preprocessing workers (default: 1)",
)
@click.option(
    "--bed-workers",
    type=int,
    default=1,
    help="Number of bed_conversion workers (default: 1)",
)
@click.option(
    "--legacy-analysis-queue",
    is_flag=True,
    help="Use legacy single analysis queue instead of separate queues per analysis type",
)
@click.option(
    "--use-ray/--no-use-ray",
    default=True,
    help="Enable Ray distributed computing (default: on). Disable with --no-use-ray.",
)
@click.option(
    "--use-ray-core",
    is_flag=True,
    help="Deprecated: Ray Core is now the default when --use-ray is provided.",
)
@click.option(
    "--ray-num-cpus",
    type=int,
    default=None,
    help="Number of CPUs to use for Ray (default: auto-detect). Only used when --use-ray is specified.",
)
@click.option(
    "--queue-priority",
    multiple=True,
    help="Set priority for a queue (e.g., 'preprocessing:10', 'bed_conversion:9', 'mgmt:5'). Only used when --use-ray is specified.",
)
@click.option(
    "--show-priorities",
    is_flag=True,
    help="Show current queue priorities and exit. Only used when --use-ray is specified.",
)
@click.option(
    "--reference",
    "-r",
    type=click.Path(exists=True, path_type=Path),
    help="Path to reference genome (FASTA format). Required for SNP calling and some other analyses.",
)
@click.option(
    "--no-watch",
    is_flag=True,
    help="Do not watch directories for new files (default: watch enabled).",
)
@click.option(
    "--with-gui/--no-gui",
    default=True,
    help="Launch NiceGUI workflow monitor (default: on). Disable with --no-gui.",
)
@click.option(
    "--gui-host",
    default="0.0.0.0",
    show_default=True,
    help="Host interface for the GUI server (e.g., 0.0.0.0 to listen on all interfaces).",
)
@click.option(
    "--gui-port",
    type=int,
    default=8081,
    show_default=True,
    help="Port for the GUI server.",
)
@click.option(
    "--preset",
    type=click.Choice(["p2i", "standard", "high"]),
    default="standard",
    help="Execution preset for Ray Core: 'p2i' (2 CPUs cap, 4 grouped actors, concurrency 1), 'standard' (default; 6 CPUs cap, grouped actors, analysis uses analysis_workers), 'high' (per-job-type actors).",
)
@click.option(
    "--ray-dashboard/--no-ray-dashboard",
    default=True,
    help="Enable Ray dashboard (default: on). Disable with --no-ray-dashboard. Only used when --use-ray is specified.",
)
@click.option(
    "--target-panel",
    type=click.Choice(_get_available_panels()),
    required=True,
    help="Target gene panel for fusion analysis. Use 'robin add-panel' to add custom panels.",
)
def workflow(
    path: Path,
    workflow: str,
    center: str,
    commands: tuple[str, ...],
    verbose: bool,
    no_process_existing: bool,
    work_dir: Optional[Path],
    log_level: str,
    job_log_level: tuple[str, ...],
    deduplicate_jobs: tuple[str, ...],
    no_progress: bool,
    analysis_workers: int,
    preprocessing_workers: int,
    bed_workers: int,
    legacy_analysis_queue: bool,
    use_ray: bool,
    use_ray_core: bool,
    ray_num_cpus: Optional[int],
    queue_priority: tuple[str, ...],
    show_priorities: bool,
    reference: Optional[Path],
    with_gui: bool,
    gui_host: str,
    gui_port: int,
    no_watch: bool,
    preset: Optional[str],
    ray_dashboard: bool,
    target_panel: str,
) -> None:
    """Run various operations on BAM files in a directory. Preprocessing is automatically included as the first step."""
    try:
        # Check for required model files first
        _check_models_or_exit()
        
        # Validate reference genome if provided
        if reference:
            try:
                # Import the validation function from matkit
                from robin.analysis.utilities.matkit import _ensure_fasta_index
                
                # Convert Path to string for the validation function
                ref_path = str(reference) if isinstance(reference, Path) else reference
                
                # Use click.echo for visibility even when log level is ERROR
                _echo_styled(f"Validating reference genome: {reference}", level="info")
                
                # Validate and ensure index exists
                _ensure_fasta_index(ref_path)
                
                _echo_styled(
                    f"Reference genome validated and indexed: {reference}",
                    level="success",
                )
            except Exception as e:
                error_msg = (
                    f"Failed to validate reference genome: {e}\n"
                    f"Please ensure the file exists and is a valid FASTA file, "
                    f"or remove the --reference parameter if not needed for this workflow."
                )
                click.echo(f"❌ {error_msg}", err=True)
                sys.exit(1)
        
        # Require user acknowledgment before proceeding
        if not _get_user_acknowledgment():
            sys.exit(1)

        _warn_if_process_large_bams()

        # Validate input parameters
        _validate_inputs(path, workflow, analysis_workers, ray_num_cpus)

        # Parse and validate inputs
        job_levels = _parse_job_log_levels(job_log_level)
        command_map = _parse_command_mappings(commands)

        # Configure logging
        configure_logging(global_level=log_level, job_levels=job_levels)

        # Parse workflow plan and convert to standard format if needed
        original_workflow_string = workflow  # Store the original workflow string
        workflow_steps = [step.strip() for step in workflow.split(",")]

        # Validate and clean workflow steps
        workflow_steps = _validate_workflow_steps(workflow_steps)

        # Check if workflow uses simplified format (no queue prefixes)
        uses_simplified_format = all(":" not in step for step in workflow_steps)

        if uses_simplified_format:
            # Convert simplified format to standard format with automatic queue assignment
            workflow_steps = _convert_simplified_workflow(workflow_steps)

        # Ensure preprocessing is the first step to extract metadata
        workflow_steps = _ensure_preprocessing_step(workflow_steps)

        # Validate that we have at least some workflow steps
        if not workflow_steps:
            raise click.BadParameter("No valid workflow steps found after validation")

        # If using Ray and CPU count specified, initialize Ray BEFORE creating runner
        # so the runner respects the requested CPU resources
        if use_ray and ray_num_cpus is not None:
            _initialize_ray(ray_num_cpus, ray_dashboard)

        # Ray Core engine path (default when --use-ray is used)
        if use_ray:
            # Debug: Log reference genome status
            # Ensure Ray is initialized if CPUs not specified
            try:
                import ray

                # Raise dashboard/State API job list limit (default 10k) if not set
                if "RAY_MAX_LIMIT_FROM_DATA_SOURCE" not in os.environ:
                    os.environ["RAY_MAX_LIMIT_FROM_DATA_SOURCE"] = "100000"

                if not ray.is_initialized():
                    # Apply preset CPU caps for Ray Core if provided
                    init_kwargs = {
                        "ignore_reinit_error": True,
                    }

                    if ray_dashboard:
                        init_kwargs.update(
                            {
                                "include_dashboard": True,
                                "dashboard_host": os.environ.get(
                                    "RAY_DASHBOARD_HOST", "0.0.0.0"
                                ),
                            }
                        )

                    if preset in {"p2i", "standard"}:
                        init_kwargs["num_cpus"] = 2 if preset == "p2i" else 6  # Increased from 4 to 6 for standard
                    try:
                        ray.init(**init_kwargs)
                    except TypeError:
                        # Older Ray versions may not support dashboard args
                        init_kwargs.pop("include_dashboard", None)
                        init_kwargs.pop("dashboard_host", None)
                        ray.init(**init_kwargs)
            except Exception:
                pass

            # Display configuration
            _display_workflow_config(
                path=path,
                center=center,
                work_dir=work_dir,
                workflow_steps=workflow_steps,
                command_map=command_map,
                log_level=log_level,
                job_levels=job_levels,
                deduplicate_jobs=deduplicate_jobs,
                legacy_analysis_queue=legacy_analysis_queue,
                analysis_workers=analysis_workers,
                no_process_existing=no_process_existing,
                uses_simplified_format=uses_simplified_format,
                original_workflow=workflow,
                use_ray=True,
                ray_num_cpus=ray_num_cpus,
                queue_priority=queue_priority,
            )

            # Create workflow runner for Ray workflow (needed for GUI integration)
            runner = _create_workflow_runner(
                True,  # use_ray=True
                verbose,
                analysis_workers,
                legacy_analysis_queue,
                log_level,
                preprocessing_workers=preprocessing_workers,
                bed_workers=bed_workers,
                reference=reference,  # Add reference parameter for Ray workflow too
                center=center,
                target_panel=target_panel,
            )

            # Run Ray Core implementation
            try:
                import asyncio
                from robin import workflow_ray as wrn

                asyncio.run(
                    wrn.run(
                        plan=workflow_steps,
                        paths=[str(path)],
                        target_panel=target_panel,
                        analysis_workers=analysis_workers,
                        process_existing=not no_process_existing,
                        monitor=not no_progress,
                        watch=(not no_watch),
                        patterns=["*.bam"],
                        ignore_patterns=None,
                        recursive=True,
                        work_dir=str(work_dir) if work_dir else None,
                        log_level=log_level,
                        preset=preset,
                        workflow_runner=runner,
                        reference=str(reference) if reference else None,
                        gui_host=gui_host,
                        gui_port=gui_port,
                        center=center,
                    )
                )
            except KeyboardInterrupt:
                print("\n[SHUTDOWN] Interrupted by user (Ctrl-C)")
                print("[SHUTDOWN] Initiating graceful shutdown...")
                # Attempt to shutdown Ray gracefully
                try:
                    import ray
                    if ray.is_initialized():
                        print("[SHUTDOWN] Shutting down Ray coordinator...")
                        # Get the coordinator and shutdown gracefully
                        try:
                            coord = ray.get_actor("robin_coordinator")
                            ray.kill(coord)
                            print("[SHUTDOWN] Ray coordinator stopped")
                        except Exception:
                            pass
                        # Shutdown Ray
                        print("[SHUTDOWN] Shutting down Ray...")
                        ray.shutdown()
                        print("[SHUTDOWN] Ray shutdown complete")
                except Exception as e:
                    print(f"[SHUTDOWN] Warning: Error shutting down Ray: {e}")
                # Attempt to shutdown GUI server if running
                try:
                    from robin.gui.app import get_gui_launcher as _get  # type: ignore

                    gl = _get()
                    if gl is not None:
                        print("[SHUTDOWN] Shutting down GUI server...")
                        try:
                            from nicegui import app as ng_app  # type: ignore

                            ng_app.shutdown()
                            print("[SHUTDOWN] GUI server shutdown complete")
                        except Exception:
                            pass
                except Exception:
                    pass
                print("[SHUTDOWN] Workflow stopped gracefully")
            return

        # Create workflow runner
        runner = _create_workflow_runner(
            use_ray,
            verbose,
            analysis_workers,
            legacy_analysis_queue,
            log_level,
            preprocessing_workers=preprocessing_workers,
            bed_workers=bed_workers,
            reference=reference,
            center=center,
            target_panel=target_panel,
        )

        # Handle Ray-specific configuration
        if use_ray and hasattr(runner, "manager"):
            # Handle show-priorities option
            if show_priorities:
                click.echo("Current Queue Priorities:")
                priority_info = runner.manager.get_priority_info()
                for queue_type, priority in sorted(
                    priority_info["queue_priorities"].items(),
                    key=lambda x: x[1],
                    reverse=True,
                ):
                    click.echo(f"  {queue_type}: {priority}")
                return

            # Configure queue priorities
            _configure_queue_priorities(runner, queue_priority)

            # Ray already initialized above when ray_num_cpus is provided

        # Configure job deduplication
        if deduplicate_jobs:
            for job_type in deduplicate_jobs:
                job_type = job_type.strip()
                if job_type:  # Skip empty job types
                    try:
                        runner.manager.add_deduplication_job_type(job_type)
                    except Exception as e:
                        click.echo(
                            f"Warning: Failed to enable deduplication for '{job_type}': {e}",
                            err=True,
                        )
            valid_dedup_jobs = [j.strip() for j in deduplicate_jobs if j.strip()]
            if valid_dedup_jobs:
                click.echo(f"Job deduplication enabled for: {valid_dedup_jobs}")

        # Register handlers and command handlers
        _register_handlers(runner, legacy_analysis_queue, work_dir, target_panel, reference, center)
        _register_command_handlers(runner, command_map, legacy_analysis_queue)

        # For Ray workflow, reinitialize processors after handlers are registered
        if use_ray and hasattr(runner, "manager"):
            try:
                runner.manager._reinitialize_processors()
            except Exception as e:
                click.echo(
                    f"Warning: Failed to reinitialize Ray processors: {e}", err=True
                )

        # Create classifier function
        classifier_func = None
        if work_dir:
            classifier_func = _create_classifier_with_work_dir(work_dir, workflow_steps, target_panel)

        # Display configuration
        _display_workflow_config(
            path=path,
            center=center,
            work_dir=work_dir,
            workflow_steps=workflow_steps,
            command_map=command_map,
            log_level=log_level,
            job_levels=job_levels,
            deduplicate_jobs=deduplicate_jobs,
            legacy_analysis_queue=legacy_analysis_queue,
            analysis_workers=analysis_workers,
            no_process_existing=no_process_existing,
            uses_simplified_format=uses_simplified_format,
            original_workflow=original_workflow_string,
            use_ray=use_ray,
            ray_num_cpus=ray_num_cpus,
            queue_priority=queue_priority,
        )

        # Launch GUI if requested
        gui_launcher = None
        if with_gui:
            try:
                print("Launching NiceGUI workflow monitor with vertical tabs...")

                # Import the new GUI launcher
                try:
                    from robin.gui.app import launch_gui

                    # Launch GUI using the new launcher FIRST so the global sender is ready
                    gui_launcher = launch_gui(
                        host=gui_host,
                        port=gui_port,
                        show=False,
                        workflow_runner=runner,
                        workflow_steps=workflow_steps,
                        monitored_directory=str(work_dir) if work_dir else str(path),
                        center=center,
                    )

                    # Now install workflow hooks for real-time monitoring
                    try:
                        from robin.workflow_hooks import install_workflow_hooks

                        install_workflow_hooks(runner, workflow_steps, str(path))
                        click.echo(
                            "Workflow state hooks installed for real-time monitoring"
                        )
                    except Exception as e:
                        click.echo(
                            f"Failed to install workflow hooks: {e}. GUI will show static information only."
                        )

                    base_url = f"http://{gui_host}:{gui_port}"
                    print(f"GUI launched successfully on {base_url}")
                    print(f"   Welcome page: {base_url}/")
                    print(f"   Workflow monitor: {base_url}/robin")
                    print(f"   Sample tracking: {base_url}/live_data")
                    print(f"   Individual samples: {base_url}/live_data/sampleID/")
                    print(
                        "   Open your browser to monitor the workflow with the new navigation structure"
                    )
                    print(
                        "   The GUI will run in the background while the workflow executes"
                    )

                except ImportError as e:
                    click.echo(f"Failed to import GUI launcher: {e}")
                    click.echo("   Please ensure the gui_launcher module is available")

            except Exception as e:
                click.echo(
                    f"Failed to launch GUI: {e}. Continuing with workflow only.",
                    err=True,
                )

        # Run the workflow
        runner.run_workflow(
            watch_dir=str(path),
            workflow_plan=workflow_steps,
            recursive=True,
            patterns=["*.bam"],
            ignore_patterns=None,
            classifier_func=classifier_func,
            process_existing=not no_process_existing,
            show_progress=not no_progress,
        )

    except KeyboardInterrupt:
        print("\n[SHUTDOWN] Interrupted by user (Ctrl-C)")
        print("[SHUTDOWN] Initiating graceful shutdown...")
        if "runner" in locals() and hasattr(runner, "manager"):
            print("[SHUTDOWN] Stopping workflow manager...")
            try:
                runner.manager.stop(timeout=DEFAULT_TIMEOUT)
                print("[SHUTDOWN] Workflow manager stopped")
            except Exception as e:
                print(f"[SHUTDOWN] Warning: Error during shutdown: {e}")
                click.echo(f"Warning: Error during shutdown: {e}", err=True)
        # Attempt to shutdown GUI server if running
        try:
            from robin.gui.app import get_gui_launcher as _get  # type: ignore

            gl = _get()
            if gl is not None:
                print("[SHUTDOWN] Shutting down GUI server...")
                try:
                    from nicegui import app as ng_app  # type: ignore

                    ng_app.shutdown()
                    print("[SHUTDOWN] GUI server shutdown complete")
                except Exception:
                    pass
        except Exception:
            pass
        print("[SHUTDOWN] Workflow stopped gracefully")
        click.echo("\nWorkflow stopped by user")
    except click.BadParameter as e:
        click.echo(f"Parameter error: {e}", err=True)
        sys.exit(1)
    except Exception as e:
        click.echo(f"Unexpected error: {e}", err=True)
        if verbose:
            import traceback

            click.echo(f"Traceback: {traceback.format_exc()}", err=True)
        sys.exit(1)


if __name__ == "__main__":
    main()
