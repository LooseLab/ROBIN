# methylartist_locus_capture.py
import os
import runpy
import shutil
import sys
from contextlib import contextmanager, redirect_stderr
from typing import Any, Dict, List, Optional

from matplotlib.figure import Figure


@contextmanager
def _patch_fig_savefig(capture_dict):
    orig = Figure.savefig

    def _capture(self, *args, **kwargs):
        capture_dict["fig"] = self
        # swallow file I/O
        return None

    Figure.savefig = _capture
    try:
        yield
    finally:
        Figure.savefig = orig


@contextmanager
def _patch_argv(argv):
    old = sys.argv
    sys.argv = argv
    try:
        yield
    finally:
        sys.argv = old


@contextmanager
def _suppress_methylartist_warnings():
    """Temporarily suppress WARNING-level logs and stderr noise from methylartist."""
    import logging

    previous_disable = logging.root.manager.disable
    # Disable WARNING and below (keeps ERROR/CRITICAL visible)
    logging.disable(logging.WARNING)
    devnull = open(os.devnull, "w")
    try:
        with redirect_stderr(devnull):
            yield
    finally:
        logging.disable(previous_disable)
        try:
            devnull.close()
        except Exception:
            pass


@contextmanager
def _prepend_conda_env_lib_to_ld_library_path():
    """
    Prepend this Python env's lib/ to LD_LIBRARY_PATH so native deps (SciPy, ICU
    libicui18n, etc.) resolve libstdc++.so.6 from conda (libstdcxx-ng) instead of
    /lib/x86_64-linux-gnu/ (fixes CXXABI_1.3.15 on older hosts when methylartist
    loads via runpy).
    """
    prefix = (
        os.environ.get("CONDA_PREFIX")
        or getattr(sys, "base_prefix", None)
        or sys.prefix
    )
    lib_dir = os.path.join(prefix, "lib")
    if not os.path.isdir(lib_dir):
        yield
        return
    previous = os.environ.get("LD_LIBRARY_PATH")
    os.environ["LD_LIBRARY_PATH"] = f"{lib_dir}:{previous}" if previous else lib_dir
    try:
        yield
    finally:
        if previous is None:
            os.environ.pop("LD_LIBRARY_PATH", None)
        else:
            os.environ["LD_LIBRARY_PATH"] = previous


def _find_cli_or_raise() -> str:
    cli = shutil.which("methylartist")
    if not cli:
        raise RuntimeError(
            "Couldn't find 'methylartist' on PATH. "
            "Activate the environment where you installed it, "
            "or add the env's bin/Scripts directory to PATH."
        )
    return cli


def has_bam_index(bam_path: str) -> bool:
    """
    Check if a BAM file has a valid index (.bai or .csi) that supports fetch.
    Returns False if the file has no index or the index is invalid.
    """
    import os

    bai_path = f"{bam_path}.bai"
    csi_path = f"{bam_path}.csi"
    if not os.path.exists(bai_path) and not os.path.exists(csi_path):
        return False
    try:
        import pysam

        with pysam.AlignmentFile(bam_path, "rb") as bam:
            if hasattr(bam, "has_index") and not bam.has_index:
                return False
            if not bam.references:
                return False
            # Actually try a fetch to ensure index works
            first_ref = bam.references[0]
            ref_len = bam.get_reference_length(first_ref) or 1000
            list(bam.fetch(first_ref, 0, min(1000, ref_len)))
            return True
    except Exception:
        return False


def _choose_main_axis(fig: Figure):
    axes = [ax for ax in getattr(fig, "axes", []) if hasattr(ax, "get_xlim")]
    if not axes:
        return fig.gca()
    spans = []
    for ax in axes:
        try:
            x0, x1 = ax.get_xlim()
            spans.append((ax, abs(x1 - x0)))
        except Exception:
            continue
    if not spans:
        return fig.gca()
    spans.sort(key=lambda t: t[1], reverse=True)
    return spans[0][0]


def locus_figure(
    interval: str,
    bam_path: str,
    motif: str = "CG",
    mods: str = "m",
    extra_cli: Optional[List[str]] = None,
    site_rows: Optional[List[Dict[str, Any]]] = None,
):
    """
    Run the installed 'methylartist locus' command but return the Matplotlib Figure
    object instead of writing an image to disk.

    Parameters
    ----------
    interval : str          e.g. 'chr10:129466536-129467536'
    bam_path : str          path to BAM with MM/ML tags
    motif : str             e.g. 'CG' (default)
    mods : str              e.g. 'm' for 5mC (default)
    extra_cli : list[str]   any extra CLI flags, e.g. ["--minqual","20","--reads","2000"]

    Returns
    -------
    matplotlib.figure.Figure
    """
    cli = _find_cli_or_raise()

    argv = [
        "methylartist",
        "locus",
        "-i",
        str(interval),
        "-b",
        str(bam_path),
        "-o",
        "ignored.png",  # required by CLI; we'll intercept savefig anyway
        "--motif",
        str(motif),
        "--mods",
        str(mods),
    ]
    if extra_cli:
        argv.extend([str(a) for a in extra_cli])

    # Final safeguard: ensure all argv entries are plain strings
    argv = [str(a) for a in argv]

    # Check if BAM file is indexed and the index is valid before calling methylartist
    import logging
    import os

    import pysam

    if not has_bam_index(bam_path):
        raise RuntimeError(
            f"BAM file is not indexed or index is invalid: {bam_path}. "
            f"Index file (.bai or .csi) not found or cannot be used for fetch."
        )
    bai_path = f"{bam_path}.bai"

    # Check if index file has content (not empty)
    try:
        bai_size = os.path.getsize(bai_path)
        if bai_size == 0:
            raise RuntimeError(
                f"BAM index file is empty: {bai_path}. "
                f"The index file exists but has no content. "
                f"This may indicate the index creation was interrupted. "
                f"Try regenerating the index."
            )
    except OSError as e:
        logging.warning(f"[MGMT] Could not check index file size: {bai_path} - {e}")
        # Continue anyway - might be a permission issue

    # Verify the index is actually readable and valid
    try:
        with pysam.AlignmentFile(bam_path, "rb") as test_bam:
            # Try to access the header first
            _ = test_bam.header

            # Verify file has references
            if not test_bam.references:
                raise RuntimeError(
                    f"BAM file has no references: {bam_path}. "
                    f"The file may be empty or corrupted."
                )

            # Verify the index is actually readable
            try:
                test_bam.check_index()
            except AttributeError:
                # Older pysam versions might not have check_index
                # Try to actually use the index by fetching a small region
                first_chrom = test_bam.references[0]
                # Try fetching a small region - this will fail if index is invalid
                list(
                    test_bam.fetch(
                        first_chrom,
                        0,
                        min(1000, test_bam.get_reference_length(first_chrom) or 1000),
                    )
                )
            except Exception as idx_error:
                logging.warning(
                    f"[MGMT] BAM index exists but is not readable: {bai_path} - {idx_error}"
                )
                raise RuntimeError(
                    f"BAM file index is invalid or corrupted: {bam_path}. "
                    f"The index file exists but cannot be read. "
                    f"Error: {idx_error}. "
                    f"This may indicate the index is incomplete or corrupted. "
                    f"Try regenerating the index with: pysam.index('{bam_path}')"
                )
    except RuntimeError:
        # Re-raise our custom RuntimeError
        raise
    except Exception as e:
        # Other errors (file not readable, etc.)
        logging.warning(f"[MGMT] Failed to validate BAM index: {bam_path} - {e}")
        raise RuntimeError(
            f"BAM file index validation failed: {bam_path}. "
            f"Error: {e}. "
            f"The index file exists but may be invalid or the BAM file may be corrupted."
        )

    captured = {"fig": None}

    with (
        _prepend_conda_env_lib_to_ld_library_path(),
        _patch_fig_savefig(captured),
        _patch_argv(argv),
        _suppress_methylartist_warnings(),
    ):
        # Execute the installed CLI script as __main__
        # Catch SystemExit exceptions from methylartist (e.g., when BAM is not indexed)
        try:
            runpy.run_path(cli, run_name="__main__")
        except SystemExit as e:
            # SystemExit(0) is normal for successful completion
            # Only raise error if exit code is non-zero
            exit_code = e.code if hasattr(e, "code") else 0
            if exit_code != 0:
                error_msg = (
                    str(e) if str(e) else f"methylartist exited with code {exit_code}"
                )
                raise RuntimeError(f"methylartist failed: {error_msg}")

    # Try to get figure from capture first
    fig = captured["fig"]

    # Fallback: if savefig wasn't called, try to get the figure from matplotlib's figure manager
    if fig is None:
        import logging

        import matplotlib.pyplot as plt

        # Try to get the current figure
        try:
            # Get all figure numbers
            figure_numbers = plt.get_fignums()
            if figure_numbers:
                # Get the most recent figure (highest number)
                latest_fig_num = max(figure_numbers)
                fig = plt.figure(latest_fig_num)
                logging.debug(
                    f"[MGMT] Captured figure from matplotlib figure manager (figure {latest_fig_num})"
                )
            else:
                # Try to get current figure (might be None)
                fig = plt.gcf()
                if fig and len(fig.get_axes()) > 0:
                    logging.debug(f"[MGMT] Captured current figure from matplotlib")
                else:
                    fig = None
        except Exception as e:
            logging.debug(
                f"[MGMT] Failed to get figure from matplotlib figure manager: {e}"
            )
            fig = None

    if fig is None:
        # Provide more helpful error message
        import logging

        logging.error(
            f"[MGMT] methylartist locus did not produce a figure. "
            f"BAM: {bam_path}, Interval: {interval}"
        )
        raise RuntimeError(
            f"methylartist locus did not produce a figure. "
            f"This may indicate that the BAM file has no reads in the specified interval "
            f"({interval}) or methylartist encountered an error. "
            f"Check that the BAM file contains methylation data (MM/ML tags) for this region."
        )

    # Verify the figure has content (axes with data)
    try:
        axes = fig.get_axes()
        if not axes:
            import logging

            logging.warning(
                f"[MGMT] Figure has no axes - methylartist may have produced an empty figure"
            )
            # Still return the figure - let the caller handle empty figures
        else:
            # Check if at least one axis has data
            has_data = False
            for ax in axes:
                if hasattr(ax, "has_data") and ax.has_data():
                    has_data = True
                    break
                # Check if axis has any artists (lines, patches, etc.)
                if len(ax.get_children()) > 0:
                    has_data = True
                    break
            if not has_data:
                import logging

                logging.debug(
                    f"[MGMT] Figure axes exist but contain no data - this may indicate no reads in interval"
                )
    except Exception as e:
        import logging

        logging.debug(f"[MGMT] Could not verify figure content: {e}")
        # Continue anyway - figure might still be valid

    # Always add MGMT CpG site annotations if we're in the MGMT interval
    # Check if interval matches MGMT region (chr10:129466536-129467536)
    try:
        import logging

        # Determine which sites to annotate
        mgmt_sites = []
        if "chr10" in interval and "129466" in interval:
            # Use site_rows if available, otherwise use default positions
            if site_rows:
                for row in site_rows:
                    pos_str = str(row.get("pos", ""))
                    site_label = str(row.get("site", ""))
                    # Extract just the number (e.g., "1" from "Site 1" or "1" from "Site 1 (CpG ...)")
                    site_num = (
                        site_label.split("(")[0].strip().replace("Site ", "").strip()
                    )
                    if not site_num:
                        # Fallback: try to extract number from the label
                        import re

                        match = re.search(r"\d+", site_label)
                        site_num = match.group(0) if match else site_label.strip()

                    parts = [
                        p.strip() for p in pos_str.split("/") if p.strip().isdigit()
                    ]
                    if len(parts) >= 2:
                        try:
                            p1 = float(parts[0])
                            p2 = float(parts[1])
                            mgmt_sites.append((p1, p2, site_num))
                        except (ValueError, IndexError):
                            continue

            # If no sites from site_rows, use default MGMT positions
            if not mgmt_sites:
                mgmt_sites = [
                    (129467255, 129467256, "1"),
                    (129467258, 129467259, "2"),
                    (129467262, 129467263, "3"),
                    (129467272, 129467273, "4"),
                ]

            # Get all axes in the figure
            axes = [ax for ax in fig.get_axes() if hasattr(ax, "get_xlim")]

            if not axes:
                axes = [fig.gca()]

            # Parse the interval to get the start position (methylartist uses relative coordinates)
            interval_start = None
            try:
                # Parse interval like "chr10:129466536-129467536"
                if ":" in interval:
                    chrom, coords = interval.split(":", 1)
                    if "-" in coords:
                        start_str, end_str = coords.split("-", 1)
                        interval_start = int(start_str)
                        interval_end = int(end_str)
            except Exception:
                pass

            # Convert absolute genomic positions to relative positions (relative to interval start)
            if interval_start:
                relative_mgmt_sites = []
                for p1, p2, label in mgmt_sites:
                    rel_p1 = p1 - interval_start
                    rel_p2 = p2 - interval_start
                    relative_mgmt_sites.append((rel_p1, rel_p2, label))
                mgmt_sites = relative_mgmt_sites

            # Find the main axis for text annotations (only calculate once)
            main_ax = _choose_main_axis(fig)

            # Pre-calculate all label positions with proper staggering
            # First, collect all sites that need annotations
            label_positions_dict = {}
            if main_ax:
                ylim = main_ax.get_ylim()
                y_range = ylim[1] - ylim[0]
                xlim = main_ax.get_xlim()
                x_range = xlim[1] - xlim[0]

                # Estimate label width in data coordinates (approximate)
                # Font size 9 with padding means roughly 0.02-0.03 of x_range per character
                max_label_len = (
                    max([len(l) for _, _, l in mgmt_sites]) if mgmt_sites else 1
                )
                label_width_data = (max_label_len + 2) * 0.015 * x_range

                # Calculate center positions for all sites
                site_centers = []
                for p1, p2, label in mgmt_sites:
                    center_x = (p1 + p2) / 2.0
                    if xlim[0] <= center_x <= xlim[1]:
                        site_centers.append((center_x, label))

                # Sort by x position to process left to right
                site_centers.sort(key=lambda x: x[0])

                # Calculate staggered y positions to avoid overlap
                base_y_offset = 0.02
                vertical_spacing = 0.025  # Fraction of y_range for spacing

                for idx, (center_x, label) in enumerate(site_centers):
                    # Start with base position
                    annotation_y = ylim[1] - (base_y_offset * y_range)

                    # Check for horizontal overlap with previously placed labels
                    overlap_threshold = label_width_data * 0.8  # 80% of label width
                    stagger_level = 0

                    # Find all overlapping labels and determine best stagger position
                    overlapping_labels = []
                    for existing_label, (
                        existing_x,
                        existing_y,
                        existing_stagger,
                    ) in label_positions_dict.items():
                        if abs(center_x - existing_x) < overlap_threshold:
                            overlapping_labels.append(
                                (existing_x, existing_y, existing_stagger)
                            )

                    if overlapping_labels:
                        # Find the maximum stagger level among overlapping labels
                        max_stagger = max([s for _, _, s in overlapping_labels])
                        stagger_level = max_stagger + 1

                        # Alternate direction: even stagger levels go up, odd levels go down
                        # This creates a zigzag pattern
                        if stagger_level % 2 == 0:
                            # Go up from the highest overlapping label
                            max_y = max([y for _, y, _ in overlapping_labels])
                            annotation_y = max_y + (vertical_spacing * y_range)
                        else:
                            # Go down from the lowest overlapping label
                            min_y = min([y for _, y, _ in overlapping_labels])
                            annotation_y = min_y - (vertical_spacing * y_range)

                        # Ensure we don't go outside the plot area
                        annotation_y = max(
                            ylim[0] + 0.05 * y_range,
                            min(ylim[1] - 0.02 * y_range, annotation_y),
                        )

                    # Store the position with stagger level
                    label_positions_dict[label] = (
                        center_x,
                        annotation_y,
                        stagger_level,
                    )

            # Add vertical lines to ALL axes (all subplots)
            for ax_idx, ax in enumerate(axes):
                try:
                    xlim = ax.get_xlim()

                    # For each MGMT site, draw vertical lines at both CpG positions
                    for p1, p2, label in mgmt_sites:
                        # Only draw if positions are within the visible x-axis range
                        if xlim[0] <= p1 <= xlim[1] or xlim[0] <= p2 <= xlim[1]:
                            # Draw dashed vertical lines - more subtle, lower zorder so data points appear on top
                            ax.axvline(
                                p1,
                                color="crimson",
                                linestyle="--",
                                linewidth=1.0,
                                alpha=0.5,
                                zorder=5,
                            )
                            ax.axvline(
                                p2,
                                color="crimson",
                                linestyle="--",
                                linewidth=1.0,
                                alpha=0.5,
                                zorder=5,
                            )

                            # Add text annotation on the main axis only
                            if ax == main_ax and label in label_positions_dict:
                                center_x, annotation_y, _ = label_positions_dict[label]

                                ax.text(
                                    center_x,
                                    annotation_y,
                                    label,
                                    ha="center",
                                    va="bottom",
                                    fontsize=9,
                                    fontweight="bold",
                                    color="crimson",
                                    bbox=dict(
                                        boxstyle="round,pad=0.3",
                                        facecolor="white",
                                        edgecolor="crimson",
                                        alpha=0.9,
                                        linewidth=1.5,
                                    ),
                                    zorder=10,
                                )
                except Exception as e:
                    logging.debug(f"[MGMT] Failed to add annotations to axis: {e}")
                    continue

            # Methylartist builds axes via a standalone GridSpec (no figure=), so
            # enabling constrained_layout here triggers:
            # "There are no gridspecs with layoutgrids..."
            # Annotations (axvline/text) do not need a layout engine pass.
            try:
                if hasattr(fig, "set_layout_engine"):
                    fig.set_layout_engine(None)
                elif hasattr(fig, "set_constrained_layout"):
                    fig.set_constrained_layout(False)
            except Exception:
                pass

            logging.debug(
                f"[MGMT] Successfully added annotations for {len(mgmt_sites)} sites"
            )
        else:
            logging.debug(
                f"[MGMT] Interval {interval} is not MGMT region, skipping annotations"
            )
    except Exception as e:
        import logging

        logging.exception(f"[MGMT] Failed to add annotations: {e}")

    return fig


# --- Optional: persistence helpers (same-Python/Matplotlib only) ---
def figure_is_renderable(fig: Figure) -> bool:
    """
    Return True if *fig* can be drawn/saved with the current matplotlib.

    Pickled figures often unpickle successfully but fail at render time after a
    matplotlib upgrade (e.g. missing ``_axis_map`` on Axes).
    """
    if fig is None:
        return False
    try:
        axes = fig.get_axes()
        if not axes:
            return False
        for ax in axes:
            if not hasattr(ax, "_axis_map"):
                return False
        import io
        import warnings

        # Older pickled MGMT figures may still have constrained_layout enabled
        # against a non-figure GridSpec; clear it before the probe save.
        try:
            if hasattr(fig, "get_constrained_layout") and fig.get_constrained_layout():
                if hasattr(fig, "set_layout_engine"):
                    fig.set_layout_engine(None)
                else:
                    fig.set_constrained_layout(False)
        except Exception:
            pass

        buf = io.BytesIO()
        with warnings.catch_warnings():
            warnings.filterwarnings(
                "ignore",
                message=r".*no gridspecs with layoutgrids.*",
                category=UserWarning,
            )
            fig.savefig(buf, format="svg")
        return True
    except Exception:
        return False


def save_figure_pickle(fig: Figure, path: str) -> None:
    """Save a Matplotlib Figure object for later re-use (not portable across versions)."""
    import pickle

    with open(path, "wb") as f:
        pickle.dump(fig, f)


def load_figure_pickle(path: str) -> Figure:
    """Reload a pickled Matplotlib Figure object."""
    import pickle

    with open(path, "rb") as f:
        return pickle.load(f)


def try_load_figure_pickle(
    path: str, *, remove_if_invalid: bool = True
) -> Optional[Figure]:
    """
    Load a pickled figure when it is renderable in the current matplotlib.

    Returns None if unpickling fails, axes are incompatible, or savefig would fail.
    """
    import logging

    import matplotlib.pyplot as plt

    try:
        fig = load_figure_pickle(path)
    except Exception as exc:
        logging.debug(f"[MGMT] Could not unpickle figure from {path}: {exc}")
        return None

    if figure_is_renderable(fig):
        return fig

    logging.warning(
        f"[MGMT] Pickled figure at {path} is not renderable in the current "
        "matplotlib environment (often a version mismatch after upgrade); "
        "will regenerate from BAM"
    )
    try:
        plt.close(fig)
    except Exception:
        pass
    if remove_if_invalid:
        try:
            os.remove(path)
        except OSError:
            pass
    return None


# --- Optional: portable-ish JSON via mpld3 (fidelity may vary) ---
def save_figure_mpld3(fig: Figure, path: str) -> None:
    import mpld3

    with open(path, "w") as f:
        f.write(mpld3.fig_to_json(fig))


def load_figure_mpld3(path: str) -> Figure:
    import mpld3

    with open(path, "r") as f:
        return mpld3.json_to_fig(f.read())
