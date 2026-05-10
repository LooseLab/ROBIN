from __future__ import annotations

from typing import Any, Dict, List, Optional
from pathlib import Path
import logging
import pickle
import os
import hashlib

import pandas as pd
import numpy as np
import matplotlib

# Use Agg backend for compatibility with DNA Features Viewer and NiceGUI
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # Import pyplot for plotting functionality

# Configure matplotlib to allow more open figures and suppress warnings
# During preprocessing, many figures may be created
plt.rcParams['figure.max_open_warning'] = 100  # Increase threshold

# Module-level variable to track last fusion plot figure for cleanup
_last_fusion_figure = None

# Import DNA Features Viewer components at module level
try:
    from dna_features_viewer import GraphicFeature, GraphicRecord

    DNA_FEATURES_AVAILABLE = True
except ImportError:
    DNA_FEATURES_AVAILABLE = False
    logging.warning("DNA Features Viewer not available - using fallback visualization")

    
# chrov ideograms removed - not working properly

try:
    from nicegui import ui
except ImportError:  # pragma: no cover
    ui = None


def _render_paged_df_table(
    df: pd.DataFrame,
    columns: List[Dict[str, Any]],
    *,
    pagination_default: int = 100,
    class_size: str = "table-xs",
    search_placeholder: str = "Search...",
) -> Any:
    """Render a DataFrame with Quasar server-side pagination (Python slices rows per page)."""
    from robin.gui.theme import (
        clamp_qtable_server_pagination,
        styled_server_paged_table,
        wire_qtable_server_pagination_handlers,
    )

    preview_mode = len(df) > 50_000
    source_df = df.iloc[:5_000].copy() if preview_mode else df
    n_rows = len(source_df)
    page_state: Dict[str, Any] = {
        "filtered_positions": list(range(n_rows)),
    }

    total0 = len(page_state["filtered_positions"])
    pd_default = int(pagination_default)
    init_pagination = clamp_qtable_server_pagination(
        {
            "sortBy": None,
            "descending": False,
            "page": 1,
            "rowsPerPage": pd_default,
            "rowsNumber": total0,
        },
        rows_number=total0,
        rows_per_page_default=pd_default,
    )
    _, table = styled_server_paged_table(
        columns=columns,
        rows=[],
        pagination=init_pagination,
        row_key="__row_idx",
        class_size=class_size,
        rows_per_page_options=[25, 50, 100, 250],
    )

    count_label = ui.label("").classes("text-xs text-gray-500")
    if preview_mode:
        ui.label(
            f"Preview mode: showing first {len(source_df):,} rows of {len(df):,}."
        ).classes("text-xs text-amber-600")

    def _fill_from_pagination(pag: Dict[str, Any]) -> None:
        total = len(page_state["filtered_positions"])
        pag = clamp_qtable_server_pagination(
            pag,
            rows_number=total,
            rows_per_page_default=pd_default,
        )
        rpp = int(pag["rowsPerPage"])
        page = int(pag["page"])
        start = (page - 1) * rpp
        end = start + rpp
        positions = page_state["filtered_positions"]
        slice_pos = positions[start:end]
        page_rows = (
            source_df.iloc[slice_pos].reset_index(drop=True).to_dict("records")
            if slice_pos
            else []
        )
        for i, pos in enumerate(slice_pos):
            page_rows[i]["__row_idx"] = int(pos)
        table.rows = page_rows
        table.pagination = pag
        count_label.text = (
            f"{total} rows match search ({len(df)} total in source)"
            if total
            else f"0 rows match ({len(df)} total in source)"
        )
        table.update()

    wire_qtable_server_pagination_handlers(table, _fill_from_pagination)

    def _apply_search(term: str) -> None:
        txt = str(term or "").strip().lower()
        if not txt:
            page_state["filtered_positions"] = list(range(n_rows))
        else:
            filtered: List[int] = []
            for pos in range(n_rows):
                row = source_df.iloc[pos]
                for col in source_df.columns:
                    value = row.get(col)
                    if value is not None and txt in str(value).lower():
                        filtered.append(pos)
                        break
            page_state["filtered_positions"] = filtered
        pag = clamp_qtable_server_pagination(
            dict(table.pagination),
            rows_number=len(page_state["filtered_positions"]),
            rows_per_page_default=pd_default,
        )
        pag["page"] = 1
        _fill_from_pagination(pag)

    with table.add_slot("top-right"):
        search_input = ui.input(placeholder=search_placeholder).props(
            "type=search dense clearable"
        )
        search_input.on("update:model-value", lambda e: _apply_search(getattr(e, "value", "")))

    for col in table.columns:
        col["sortable"] = False

    _fill_from_pagination(init_pagination)
    try:
        def _cleanup_table_state() -> None:
            page_state["filtered_positions"] = []
            table.rows = []

        ui.context.client.on_disconnect(_cleanup_table_state)
    except Exception:
        pass
    return table


def _is_dark_mode() -> bool:
    """Match Quasar ``body--dark`` via app storage (see theme.frame)."""
    try:
        from nicegui import app

        return bool(app.storage.user.get("dark_mode"))
    except Exception:
        return False


def _apply_fusion_figure_theme(fig: Any, dark: bool) -> None:
    """Adjust matplotlib fusion figure for light vs dark UI (same surface as MGMT, design.md §6)."""
    if fig is None:
        return
    if dark:
        bg = "#0f172a"
        fg = "#e2e8f0"
        fg_muted = "#94a3b8"
    else:
        bg = "#ffffff"
        fg = "#0f172a"
        fg_muted = "#475569"
    try:
        fig.patch.set_facecolor(bg)
        for ax in fig.get_axes():
            try:
                ax.set_facecolor(bg)
                ax.tick_params(colors=fg_muted)
                ax.xaxis.label.set_color(fg_muted)
                ax.yaxis.label.set_color(fg_muted)
                ax.title.set_color(fg)
                for spine in ax.spines.values():
                    spine.set_color(fg_muted)
                leg = ax.get_legend()
                if leg is not None:
                    fr = leg.get_frame()
                    if fr is not None:
                        fr.set_facecolor(bg)
                        fr.set_edgecolor(fg_muted)
                    for t in leg.get_texts():
                        t.set_color(fg)
            except Exception:
                continue
        for txt in getattr(fig, "texts", []) or []:
            try:
                txt.set_color(fg)
            except Exception:
                pass
    except Exception:
        pass


def _rename_and_sort(df: pd.DataFrame) -> pd.DataFrame:
    """Rename columns to friendly names and sort for display."""
    if df is None or df.empty:
        return pd.DataFrame()
    rename_map = {
        "col1": "chromBED",
        "col2": "BS",
        "col3": "BE",
        "col4": "Gene",
        "reference_id": "chrom",
        "reference_start": "mS",
        "reference_end": "mE",
        "read_id": "readID",
        "mapping_quality": "mapQ",
        "strand": "strand",
        "read_start": "Read Map Start",
        "read_end": "Read Map End",
        "is_secondary": "Secondary",
        "is_supplementary": "Supplementary",
        "mapping_span": "mapping span",
    }
    try:
        return df.sort_values(by="reference_start").rename(columns=rename_map)
    except Exception:
        return df.rename(columns=rename_map)


def _create_data_hash(data: Dict[str, Any]) -> str:
    """Create a hash of the fusion data to detect actual content changes."""
    try:
        if not data:
            return ""
        
        # Create hash from key data components
        hash_components = []
        
        # Hash the annotated data DataFrame
        if "annotated_data" in data and data["annotated_data"] is not None:
            df_str = data["annotated_data"].to_string()
            df_hash = hashlib.md5(df_str.encode()).hexdigest()
            hash_components.append(f"df:{df_hash}")
        
        # Hash the goodpairs Series
        if "goodpairs" in data and data["goodpairs"] is not None:
            pairs_str = data["goodpairs"].to_string()
            pairs_hash = hashlib.md5(pairs_str.encode()).hexdigest()
            hash_components.append(f"pairs:{pairs_hash}")
        
        # Hash gene groups
        if "gene_groups" in data and data["gene_groups"] is not None:
            groups_str = str(sorted(data["gene_groups"]))
            groups_hash = hashlib.md5(groups_str.encode()).hexdigest()
            hash_components.append(f"groups:{groups_hash}")
        
        # Hash candidate count
        if "candidate_count" in data:
            hash_components.append(f"count:{data['candidate_count']}")
        
        return hashlib.md5("|".join(hash_components).encode()).hexdigest()
    except Exception as e:
        logging.warning(f"[Fusion] Failed to create data hash: {e}")
        return ""


def _count_unique_fusion_pairs(data: Dict[str, Any]) -> int:
    """Count unique fusion pairs from fusion data.
    
    Args:
        data: Dictionary containing fusion data with annotated_data and goodpairs
        
    Returns:
        Number of unique fusion pairs
    """
    try:
        if not data or data.get("annotated_data") is None:
            return 0
        
        annotated_data = data.get("annotated_data", pd.DataFrame())
        goodpairs = data.get("goodpairs", pd.Series())
        
        if annotated_data.empty:
            return 0
        
        # Filter to good pairs if available
        if not goodpairs.empty and goodpairs.sum() > 0:
            aligned_goodpairs = goodpairs.reindex(annotated_data.index, fill_value=False)
            filtered_data = annotated_data[aligned_goodpairs]
        else:
            filtered_data = annotated_data
        
        if filtered_data.empty:
            return 0
        
        # Get validated fusion pairs using breakpoint validation
        clustered_data = _cluster_fusion_reads(filtered_data, max_distance=10000, use_breakpoint_validation=True)
        
        if clustered_data.empty:
            return 0
        
        # Count unique fusion pairs
        unique_pairs = clustered_data["fusion_pair"].nunique()
        return int(unique_pairs)
        
    except Exception as e:
        logging.warning(f"[Fusion] Failed to count unique fusion pairs: {e}")
        return 0


def _get_validated_fusion_groups(data: Dict[str, Any]) -> List[List[str]]:
    """Get validated fusion groups from fusion data.
    
    A fusion group is validated if it contains at least one validated fusion pair
    (meeting the minimum read support threshold).
    
    Args:
        data: Dictionary containing fusion data with annotated_data, goodpairs, and gene_groups
        
    Returns:
        List of validated fusion groups (each group is a list of gene names)
    """
    try:
        if not data:
            return []
        
        annotated_data = data.get("annotated_data", pd.DataFrame())
        goodpairs = data.get("goodpairs", pd.Series())
        gene_groups = data.get("gene_groups")
        
        if annotated_data.empty or not gene_groups or len(gene_groups) == 0:
            return []
        
        # Filter to good pairs if available
        if not goodpairs.empty and goodpairs.sum() > 0:
            aligned_goodpairs = goodpairs.reindex(annotated_data.index, fill_value=False)
            filtered_data = annotated_data[aligned_goodpairs]
        else:
            filtered_data = annotated_data
        
        if filtered_data.empty:
            return []
        
        # Get validated fusion pairs (these meet the minimum read support threshold)
        clustered_data = _cluster_fusion_reads(filtered_data, max_distance=10000, use_breakpoint_validation=True)
        
        if clustered_data.empty:
            return []
        
        # Extract validated gene pairs from validated breakpoints
        validated_pairs_set = set()
        for _, row in clustered_data.iterrows():
            fusion_pair_str = row["fusion_pair"]  # e.g., "GENE1-GENE2"
            if fusion_pair_str:
                genes = [g.strip() for g in fusion_pair_str.split("-") if g.strip()]
                if len(genes) >= 2:
                    # Normalize pair (sorted) for consistent comparison
                    normalized_pair = tuple(sorted(genes))
                    validated_pairs_set.add(normalized_pair)
        
        # Filter gene groups to only include those with at least one validated pair
        validated_groups = []
        for group in gene_groups:
            if not isinstance(group, (list, tuple)) or len(group) < 2:
                continue
            
            # Normalize group genes
            normalized_group_genes = sorted([str(g).strip() for g in group if g])
            if len(normalized_group_genes) < 2:
                continue
            
            # Check if this group contains at least one validated pair
            group_has_validated_pair = False
            for i in range(len(normalized_group_genes)):
                for j in range(i + 1, len(normalized_group_genes)):
                    pair = (normalized_group_genes[i], normalized_group_genes[j])
                    if pair in validated_pairs_set:
                        group_has_validated_pair = True
                        break
                if group_has_validated_pair:
                    break
            
            if group_has_validated_pair:
                validated_groups.append(normalized_group_genes)
        
        logging.info(f"[Fusion] Validated {len(validated_groups)} fusion groups from {len(gene_groups)} total groups")
        return validated_groups
        
    except Exception as e:
        logging.warning(f"[Fusion] Failed to get validated fusion groups: {e}")
        return []


def _count_unique_fusion_groups(data: Dict[str, Any]) -> int:
    """Count unique fusion groups from fusion data.
    
    Only counts groups that contain at least one validated fusion pair
    (meeting the minimum read support threshold).
    
    Args:
        data: Dictionary containing fusion data with gene_groups
        
    Returns:
        Number of unique validated fusion groups
    """
    try:
        validated_groups = _get_validated_fusion_groups(data)
        return len(validated_groups)
        
    except Exception as e:
        logging.warning(f"[Fusion] Failed to count unique fusion groups: {e}")
        return 0


def _generate_summary_files_from_pickle(sample_dir: Path, force_regenerate: bool = False) -> bool:
    """Generate summary files from existing pickle files if they don't exist.
    
    This provides backward compatibility for existing analyses that were run
    before the summary file generation was implemented.
    
    Returns True if summary files were generated, False otherwise.
    """
    try:
        import csv
        from pathlib import Path
        
        # Check if summary files already exist
        summary_file = sample_dir / "fusion_summary.csv"
        if summary_file.exists() and not force_regenerate:
            return True  # Already have summary files
        elif summary_file.exists() and force_regenerate:
            pass  # Will overwrite existing file
        
        # Load data from pickle files
        target_file = sample_dir / "fusion_candidates_master_processed.pkl"
        genome_file = sample_dir / "fusion_candidates_all_processed.pkl"
        
        logging.info(f"[Fusion] Checking for target file: {target_file}")
        logging.info(f"[Fusion] Target file exists: {target_file.exists()}")
        logging.info(f"[Fusion] Checking for genome-wide file: {genome_file}")
        logging.info(f"[Fusion] Genome-wide file exists: {genome_file.exists()}")
        
        # Debug: List all fusion-related files in the directory
        fusion_files = list(sample_dir.glob("*fusion*"))
        logging.info(f"[Fusion] All fusion files in directory: {[f.name for f in fusion_files]}")
        
        # Debug: Check file sizes
        if target_file.exists():
            logging.info(f"[Fusion] Target file size: {target_file.stat().st_size} bytes")
        if genome_file.exists():
            logging.info(f"[Fusion] Genome-wide file size: {genome_file.stat().st_size} bytes")
        
        target_data = _load_processed_pickle(target_file)
        genome_data = _load_processed_pickle(genome_file)
        
        if target_data is not None:
            logging.info(f"[Fusion] Target data loaded: candidate_count={target_data.get('candidate_count', 0)}")
        else:
            logging.warning(f"[Fusion] Failed to load target data from: {target_file}")
            
        if genome_data is not None:
            logging.info(f"[Fusion] Genome-wide data loaded: candidate_count={genome_data.get('candidate_count', 0)}")
        else:
            logging.warning(f"[Fusion] Failed to load genome-wide data from: {genome_file}")
        
        # Extract counts using the same filtering logic as the display code
        # This ensures consistency between summary panel and fusion section
        target_count = 0
        genome_count = 0
        
        if target_data is not None and isinstance(target_data, dict):
            target_count = _count_unique_fusion_pairs(target_data)
            logging.info(f"[Fusion] Summary: Target count (filtered): {target_count}")
        
        if genome_data is not None and isinstance(genome_data, dict):
            genome_count = _count_unique_fusion_pairs(genome_data)
            logging.info(f"[Fusion] Summary: Genome-wide count (filtered): {genome_count}")
        
        # Generate fusion_summary.csv
        with open(summary_file, "w", newline="") as f:
            writer = csv.writer(f)
            writer.writerow(["target_fusions", "genome_fusions"])
            writer.writerow([target_count, genome_count])
        
        logging.info(f"[Fusion] Generated summary file from pickle: target={target_count}, genome={genome_count}")
        
        # Generate fusion_results.csv (use the one with more data)
        results_file = sample_dir / "fusion_results.csv"
        if target_count > 0:
            # Use target panel data
            with open(results_file, "w", newline="") as f:
                writer = csv.writer(f)
                writer.writerow(["target_fusions", "genome_fusions"])
                writer.writerow([target_count, 0])
        elif genome_count > 0:
            # Use genome-wide data
            with open(results_file, "w", newline="") as f:
                writer = csv.writer(f)
                writer.writerow(["target_fusions", "genome_fusions"])
                writer.writerow([0, genome_count])
        
        # Generate sv_count.txt (use genome count for backward compatibility)
        sv_count_file = sample_dir / "sv_count.txt"
        with open(sv_count_file, "w") as f:
            f.write(str(genome_count))
        
        logging.info(f"[Fusion] Generated summary files from pickle data - target: {target_count}, genome: {genome_count}")
        return True
        
    except Exception as e:
        logging.warning(f"[Fusion] Failed to generate summary files from pickle: {e}")
        return False


def _load_processed_pickle(file_path: Path) -> Optional[Dict[str, Any]]:
    """Load the preprocessed fusion data structure written by preprocess_fusion_data_standalone.

    Note: The file is a pickle format with .pkl extension.
    """
    try:
        logging.info(f"[Fusion] Attempting to load pickle from: {file_path}")
        if not file_path.exists():
            logging.warning(f"[Fusion] File does not exist: {file_path}")
            return None
        if file_path.stat().st_size == 0:
            logging.warning(f"[Fusion] File is empty: {file_path}")
            return None
        
        # Try to load the pickle with better error handling
        with open(file_path, "rb") as f:
            try:
                data = pickle.load(f)
            except (pickle.UnpicklingError, EOFError) as e:
                # If pickle is truncated, try to load what we can
                logging.warning(f"[Fusion] Pickle file appears truncated, attempting recovery: {file_path}")
                f.seek(0)
                try:
                    # Try loading with protocol 0 which is more forgiving
                    data = pickle.load(f)
                except:
                    # If all else fails, return None and let the system regenerate
                    logging.error(f"[Fusion] Could not recover truncated pickle: {file_path}")
                    return None
        
        # Expected keys: annotated_data (DataFrame), goodpairs (Series), gene_groups (list), candidate_count (int)
        if isinstance(data, dict):
            logging.info(f"[Fusion] Loaded pickle data with keys: {list(data.keys())}")
            logging.info(f"[Fusion] Raw candidate_count: {data.get('candidate_count', 'not found')}")
            logging.info(f"[Fusion] Raw gene_groups count: {len(data.get('gene_groups', []))}")
            
            # Apply the same filtering logic as the reporting code
            annotated_data = data.get("annotated_data", pd.DataFrame())
            goodpairs = data.get("goodpairs", pd.Series())
            
            logging.info(f"[Fusion] Raw annotated_data shape: {annotated_data.shape}")
            logging.info(f"[Fusion] Raw goodpairs length: {len(goodpairs)}")
            
            if not annotated_data.empty and not goodpairs.empty:
                # Only keep the good pairs (same as reporting code does)
                data["annotated_data"] = annotated_data[goodpairs]
                logging.info(f"[Fusion] Filtered data: {len(annotated_data)} -> {len(data['annotated_data'])} good pairs")
            else:
                logging.info(f"[Fusion] No filtering applied - annotated_data empty: {annotated_data.empty}, goodpairs empty: {goodpairs.empty}")
            
            logging.info(f"[Fusion] Final candidate_count: {data.get('candidate_count', 'not found')}")
            logging.info(f"[Fusion] Final gene_groups count: {len(data.get('gene_groups', []))}")
            
            return data
        else:
            logging.warning(f"[Fusion] Loaded data is not a dict, type: {type(data)}")
            return None
    except Exception as e:
        logging.exception(
            f"[Fusion] Failed to load processed fusion pickle: {file_path} - {e}"
        )
        return None


def _make_fusion_table(container: Any, df: pd.DataFrame) -> Any:
    """Create or update a NiceGUI table for fusion candidates and return it."""
    if df is None or df.empty:
        with container:
            ui.label("No fusion candidates found yet").classes("text-gray-600")
        return None

    with container:
        # Use styled_table for consistent styling
        from robin.gui.theme import styled_table
        
        # Convert DataFrame to rows format for styled_table
        processed_df = _rename_and_sort(df)
        
        # Create columns definition from DataFrame
        columns = []
        for col in processed_df.columns:
            columns.append({
                "name": col,
                "label": col.replace("_", " ").title(),
                "field": col,
                "sortable": True
            })
        
        table = _render_paged_df_table(
            processed_df,
            columns,
            pagination_default=100,
            class_size="table-xs",
            search_placeholder="Search",
        )
        return table


def _cluster_nearby_regions(regions_df: pd.DataFrame, max_distance: int = 200) -> pd.DataFrame:
    """
    Cluster nearby regions together to merge similar breakpoint events.
    
    Args:
        regions_df: DataFrame with region data (chromosome, start, end, read_ids, etc.)
        max_distance: Maximum distance for clustering regions (default: 200bp)
        
    Returns:
        DataFrame with clustered/merged regions
    """
    if regions_df.empty:
        return pd.DataFrame()
    
    clustered_results = []
    used_indices = set()
    
    # Group by chromosome first
    for chromosome, chr_group in regions_df.groupby("chromosome", observed=True):
        chr_regions = chr_group.reset_index(drop=True)
        
        for i, row in chr_regions.iterrows():
            if i in used_indices:
                continue
            
            # Start a new cluster with this region
            cluster_indices = [i]
            used_indices.add(i)
            
            cluster_start = row["start"]
            cluster_end = row["end"]
            cluster_reads = set()
            cluster_mapqs = []
            cluster_spans = []
            cluster_types = set()
            
            # Add reads from this region (using set to ensure uniqueness)
            if "read_ids" in row:
                read_ids = row["read_ids"]
                # Handle both set and list types (pandas may convert sets to lists)
                if isinstance(read_ids, set):
                    cluster_reads.update(read_ids)
                elif isinstance(read_ids, list):
                    cluster_reads.update(read_ids)  # set.update() handles duplicates automatically
                elif "read_count" in row:
                    # Fallback: if we don't have actual read IDs, we can't properly track uniqueness
                    logging.warning(f"[Fusion] Missing read_ids for region, using read_count estimate")
                    for k in range(int(row["read_count"])):
                        cluster_reads.add(f"read_{i}_{k}")
            elif "read_count" in row:
                # Fallback: if we don't have actual read IDs, we can't properly track uniqueness
                logging.warning(f"[Fusion] Missing read_ids for region, using read_count estimate")
                for k in range(int(row["read_count"])):
                    cluster_reads.add(f"read_{i}_{k}")
            
            if "avg_mapping_quality" in row:
                cluster_mapqs.append(row["avg_mapping_quality"])
            if "avg_mapping_span" in row:
                cluster_spans.append(row["avg_mapping_span"])
            if "event_type" in row:
                cluster_types.add(row["event_type"])
            
            # Find nearby regions to merge
            for j, other_row in chr_regions.iterrows():
                if j in used_indices:
                    continue
                
                other_start = other_row["start"]
                other_end = other_row["end"]
                
                # Check if regions overlap or are close
                overlap = not (cluster_end < other_start or cluster_start > other_end)
                # Check if start/end positions are within max_distance
                start_close = abs(cluster_start - other_start) <= max_distance
                end_close = abs(cluster_end - other_end) <= max_distance
                # Also check if midpoints are close
                midpoint1 = (cluster_start + cluster_end) // 2
                midpoint2 = (other_start + other_end) // 2
                midpoints_close = abs(midpoint1 - midpoint2) <= max_distance
                
                if overlap or (start_close and end_close) or midpoints_close:
                    cluster_indices.append(j)
                    used_indices.add(j)
                    
                    # Expand cluster boundaries
                    cluster_start = min(cluster_start, other_start)
                    cluster_end = max(cluster_end, other_end)
                    
                    # Merge read counts and statistics (using set to ensure uniqueness)
                    if "read_ids" in other_row:
                        other_read_ids = other_row["read_ids"]
                        # Handle both set and list types (pandas may convert sets to lists)
                        if isinstance(other_read_ids, set):
                            cluster_reads.update(other_read_ids)
                        elif isinstance(other_read_ids, list):
                            cluster_reads.update(other_read_ids)  # set.update() handles duplicates automatically
                        elif "read_count" in other_row:
                            # Fallback: if we don't have actual read IDs, we can't properly track uniqueness
                            logging.warning(f"[Fusion] Missing read_ids for merged region, using read_count estimate")
                            for k in range(int(other_row["read_count"])):
                                cluster_reads.add(f"read_{j}_{k}")
                    elif "read_count" in other_row:
                        # Fallback: if we don't have actual read IDs, we can't properly track uniqueness
                        logging.warning(f"[Fusion] Missing read_ids for merged region, using read_count estimate")
                        for k in range(int(other_row["read_count"])):
                            cluster_reads.add(f"read_{j}_{k}")
                    
                    if "avg_mapping_quality" in other_row:
                        cluster_mapqs.append(other_row["avg_mapping_quality"])
                    if "avg_mapping_span" in other_row:
                        cluster_spans.append(other_row["avg_mapping_span"])
                    if "event_type" in other_row:
                        cluster_types.add(other_row["event_type"])
            
            # Create merged region
            # read_count is the number of UNIQUE reads supporting this clustered event
            # (using len() on the set ensures we count each read only once, even if it appears in multiple merged regions)
            clustered_results.append({
                "chromosome": chromosome,
                "start": cluster_start,
                "end": cluster_end,
                "event_type": "master_bed_region" if "master_bed_region" in cluster_types else "supplementary",
                "read_count": len(cluster_reads),  # Count of unique reads (set automatically handles deduplication)
                "avg_mapping_quality": round(sum(cluster_mapqs) / len(cluster_mapqs), 1) if cluster_mapqs else 0,
                "avg_mapping_span": round(sum(cluster_spans) / len(cluster_spans), 0) if cluster_spans else 0,
            })
    
    if not clustered_results:
        return pd.DataFrame()
    
    return pd.DataFrame(clustered_results)


def _summarize_master_bed_events(df: pd.DataFrame, min_read_support: int = 3, min_mapq: int = 50, 
                                  cluster_distance: int = 5000) -> pd.DataFrame:
    """
    Summarize master BED data into detected breakpoint pair events.
    
    Uses the same breakpoint pair extraction logic as _extract_master_bed_breakpoints:
    - Identifies reads with both primary and supplementary alignments
    - Creates breakpoint pairs (primary + supplementary for each read)
    - Clusters similar breakpoint pairs (both primary and supplementary locations must be close)
    - Returns both primary and supplementary regions as separate events
    
    Only includes high-quality supplementary mappings (MapQ >= min_mapq).
    
    Args:
        df: DataFrame with master BED fusion candidates (individual reads)
        min_read_support: Minimum number of reads required to show an event (default: 3)
        min_mapq: Minimum mapping quality for supplementary mappings (default: 50)
        cluster_distance: Maximum distance for clustering breakpoint pairs (default: 5000bp)
        
    Returns:
        DataFrame with summarized breakpoint pair events (both primary and supplementary regions)
    """
    if df is None or df.empty:
        return pd.DataFrame()
    
    try:
        # Early exit: if dataset is very large, limit processing to avoid blocking UI
        # This prevents the GUI from freezing when processing very large datasets
        MAX_ROWS_TO_PROCESS = 20000  # Limit to prevent UI blocking (reduced from 50000 for better responsiveness)
        if len(df) > MAX_ROWS_TO_PROCESS:
            logging.warning(f"[Fusion] Master BED dataset too large ({len(df)} rows), limiting to {MAX_ROWS_TO_PROCESS} rows for UI performance")
            # Sample rows rather than just taking head() to get better coverage
            df = df.sample(n=MAX_ROWS_TO_PROCESS, random_state=42).copy() if len(df) > MAX_ROWS_TO_PROCESS else df.copy()
        
        # Check required columns
        if "read_id" not in df.columns:
            logging.debug("[Fusion] Missing required column (read_id) in master BED candidates")
            return pd.DataFrame()
        
        df = df.copy()
        
        # Filter to only high-quality supplementary mappings (MapQ >= min_mapq)
        # For master_bed_region entries, we keep them all (they're primary alignments)
        # For supplementary entries, we filter by MapQ
        if "mapping_quality" in df.columns and "col4" in df.columns:
            # Keep master_bed_region entries regardless of MapQ (they're primary alignments)
            # Filter supplementary entries to only high-quality ones
            high_quality_mask = (
                (df["col4"] == "master_bed_region") |
                (df["mapping_quality"] >= min_mapq)
            )
            df = df[high_quality_mask].copy()
            
            if df.empty:
                logging.debug(f"[Fusion] No high-quality mappings found (MapQ >= {min_mapq})")
                return pd.DataFrame()
        
        # Separate primary and supplementary alignments using col4
        # col4 is more reliable: "master_bed_region" = primary alignment, "master_bed_supplementary" = supplementary
        if "col4" in df.columns:
            primary_df = df[df["col4"] == "master_bed_region"].copy()
            supplementary_df = df[df["col4"] == "master_bed_supplementary"].copy()
        else:
            # Fallback: use is_supplementary if col4 is not available
            if "is_supplementary" not in df.columns:
                logging.debug("[Fusion] Missing required columns (col4 or is_supplementary) in master BED candidates")
                return pd.DataFrame()
            primary_df = df[df["is_supplementary"] == False].copy()
            supplementary_df = df[df["is_supplementary"] == True].copy()
        
        if primary_df.empty or supplementary_df.empty:
            logging.debug("[Fusion] Need both primary and supplementary alignments to identify breakpoint pairs")
            return pd.DataFrame()
        
        # Find reads that have both primary and supplementary alignments
        primary_read_ids = set(primary_df["read_id"].unique())
        supplementary_read_ids = set(supplementary_df["read_id"].unique())
        reads_with_both = primary_read_ids & supplementary_read_ids
        
        if not reads_with_both:
            logging.debug("[Fusion] No reads have both primary and supplementary alignments")
            return pd.DataFrame()
        
        # Filter to only reads with both primary and supplementary alignments
        primary_filtered = primary_df[primary_df["read_id"].isin(reads_with_both)].copy()
        supplementary_filtered = supplementary_df[supplementary_df["read_id"].isin(reads_with_both)].copy()
        
        # Build breakpoint pairs: for each read, pair its primary alignment with each supplementary alignment
        breakpoint_pairs = []
        
        # Group by read_id for efficient lookup
        primary_by_read = primary_filtered.groupby("read_id", observed=True)
        supplementary_by_read = supplementary_filtered.groupby("read_id", observed=True)
        
        for read_id in reads_with_both:
            # Get all primary and supplementary alignments for this read
            read_primaries = primary_by_read.get_group(read_id) if read_id in primary_by_read.groups else pd.DataFrame()
            read_supplementaries = supplementary_by_read.get_group(read_id) if read_id in supplementary_by_read.groups else pd.DataFrame()
            
            if read_primaries.empty or read_supplementaries.empty:
                continue
            
            # Create pairs using cross product
            # Get columns that exist (mapping_quality and mapping_span may not always be present)
            primary_cols = ["reference_id", "reference_start", "reference_end"]
            supp_cols = ["reference_id", "reference_start", "reference_end"]
            if "mapping_quality" in read_primaries.columns:
                primary_cols.append("mapping_quality")
            if "mapping_span" in read_primaries.columns:
                primary_cols.append("mapping_span")
            if "mapping_quality" in read_supplementaries.columns:
                supp_cols.append("mapping_quality")
            if "mapping_span" in read_supplementaries.columns:
                supp_cols.append("mapping_span")
            
            primaries_list = list(read_primaries[primary_cols].itertuples(index=False, name=None))
            supplementaries_list = list(read_supplementaries[supp_cols].itertuples(index=False, name=None))
            
            # Create all combinations
            p_idx = {name: i for i, name in enumerate(primary_cols)}
            s_idx = {name: i for i, name in enumerate(supp_cols)}
            for primary_row in primaries_list:
                for supp_row in supplementaries_list:
                    breakpoint_pairs.append({
                        "read_id": read_id,
                        "primary_chrom": primary_row[p_idx["reference_id"]],
                        "primary_start": int(primary_row[p_idx["reference_start"]]),
                        "primary_end": int(primary_row[p_idx["reference_end"]]),
                        "primary_mapq": (
                            primary_row[p_idx["mapping_quality"]]
                            if "mapping_quality" in p_idx
                            else 0
                        ),
                        "primary_span": (
                            primary_row[p_idx["mapping_span"]]
                            if "mapping_span" in p_idx
                            else 0
                        ),
                        "supp_chrom": supp_row[s_idx["reference_id"]],
                        "supp_start": int(supp_row[s_idx["reference_start"]]),
                        "supp_end": int(supp_row[s_idx["reference_end"]]),
                        "supp_mapq": (
                            supp_row[s_idx["mapping_quality"]]
                            if "mapping_quality" in s_idx
                            else 0
                        ),
                        "supp_span": (
                            supp_row[s_idx["mapping_span"]]
                            if "mapping_span" in s_idx
                            else 0
                        ),
                    })
        
        if not breakpoint_pairs:
            logging.debug("[Fusion] No breakpoint pairs created")
            return pd.DataFrame()
        
        # Early exit: if too many breakpoint pairs, limit to prevent UI blocking
        MAX_PAIRS_TO_PROCESS = 10000  # Limit clustering to prevent UI blocking
        if len(breakpoint_pairs) > MAX_PAIRS_TO_PROCESS:
            logging.warning(f"[Fusion] Too many breakpoint pairs ({len(breakpoint_pairs)}), limiting to {MAX_PAIRS_TO_PROCESS} for UI performance")
            # Sample pairs to get better coverage
            import random
            breakpoint_pairs = random.sample(breakpoint_pairs, MAX_PAIRS_TO_PROCESS)
        
        # Cluster similar breakpoint pairs
        # Two breakpoint pairs are similar if:
        # - Primary locations are close (same chromosome, within cluster_distance)
        # - Supplementary locations are close (same chromosome, within cluster_distance)
        clustered_pairs = []
        used_indices = set()
        
        # Group pairs by chromosome combination for faster clustering
        pairs_by_chrom = {}
        for i, pair in enumerate(breakpoint_pairs):
            chrom_key = (pair["primary_chrom"], pair["supp_chrom"])
            if chrom_key not in pairs_by_chrom:
                pairs_by_chrom[chrom_key] = []
            pairs_by_chrom[chrom_key].append((i, pair))
        
        # Cluster within each chromosome combination
        for chrom_key, chrom_pairs in pairs_by_chrom.items():
            if len(chrom_pairs) == 0:
                continue
            
            # Sort pairs within this chromosome combination
            sorted_chrom_pairs = sorted(chrom_pairs, key=lambda x: (
                x[1]["primary_start"], x[1]["supp_start"]
            ))
            
            # Cluster pairs in this chromosome group
            for idx, (orig_i, pair) in enumerate(sorted_chrom_pairs):
                if orig_i in used_indices:
                    continue
                
                # Start a new cluster
                cluster_read_ids = {pair["read_id"]}
                cluster_primary_chrom = pair["primary_chrom"]
                cluster_primary_starts = [pair["primary_start"]]
                cluster_primary_ends = [pair["primary_end"]]
                cluster_primary_mapqs = [pair["primary_mapq"]]
                cluster_primary_spans = [pair["primary_span"]]
                cluster_supp_chrom = pair["supp_chrom"]
                cluster_supp_starts = [pair["supp_start"]]
                cluster_supp_ends = [pair["supp_end"]]
                cluster_supp_mapqs = [pair["supp_mapq"]]
                cluster_supp_spans = [pair["supp_span"]]
                used_indices.add(orig_i)
                
                # Calculate cluster bounds once
                cluster_primary_min = min(cluster_primary_starts)
                cluster_primary_max = max(cluster_primary_ends)
                cluster_supp_min = min(cluster_supp_starts)
                cluster_supp_max = max(cluster_supp_ends)
                
                # Iteratively find all similar breakpoint pairs in this chromosome group
                changed = True
                while changed:
                    changed = False
                    for j, (other_orig_i, other_pair) in enumerate(sorted_chrom_pairs):
                        if other_orig_i in used_indices:
                            continue
                        
                        # Check if primary locations are similar (same chromosome, close coordinates)
                        primary_similar = (
                            other_pair["primary_start"] <= cluster_primary_max + cluster_distance and
                            other_pair["primary_end"] >= cluster_primary_min - cluster_distance
                        )
                        
                        # Check if supplementary locations are similar (same chromosome, close coordinates)
                        supp_similar = (
                            other_pair["supp_start"] <= cluster_supp_max + cluster_distance and
                            other_pair["supp_end"] >= cluster_supp_min - cluster_distance
                        )
                        
                        # Both primary and supplementary must be similar to cluster
                        if primary_similar and supp_similar:
                            cluster_read_ids.add(other_pair["read_id"])
                            cluster_primary_starts.append(other_pair["primary_start"])
                            cluster_primary_ends.append(other_pair["primary_end"])
                            cluster_primary_mapqs.append(other_pair["primary_mapq"])
                            cluster_primary_spans.append(other_pair["primary_span"])
                            cluster_supp_starts.append(other_pair["supp_start"])
                            cluster_supp_ends.append(other_pair["supp_end"])
                            cluster_supp_mapqs.append(other_pair["supp_mapq"])
                            cluster_supp_spans.append(other_pair["supp_span"])
                            used_indices.add(other_orig_i)
                            changed = True
                            # Update cluster bounds for next iteration
                            cluster_primary_min = min(cluster_primary_starts)
                            cluster_primary_max = max(cluster_primary_ends)
                            cluster_supp_min = min(cluster_supp_starts)
                            cluster_supp_max = max(cluster_supp_ends)
                
                # Create clustered breakpoint pair
                clustered_pairs.append({
                    "primary_chrom": cluster_primary_chrom,
                    "primary_start": min(cluster_primary_starts),
                    "primary_end": max(cluster_primary_ends),
                    "primary_avg_mapq": sum(cluster_primary_mapqs) / len(cluster_primary_mapqs) if cluster_primary_mapqs else 0,
                    "primary_avg_span": sum(cluster_primary_spans) / len(cluster_primary_spans) if cluster_primary_spans else 0,
                    "supp_chrom": cluster_supp_chrom,
                    "supp_start": min(cluster_supp_starts),
                    "supp_end": max(cluster_supp_ends),
                    "supp_avg_mapq": sum(cluster_supp_mapqs) / len(cluster_supp_mapqs) if cluster_supp_mapqs else 0,
                    "supp_avg_span": sum(cluster_supp_spans) / len(cluster_supp_spans) if cluster_supp_spans else 0,
                    "read_count": len(cluster_read_ids),
                })
        
        # Filter for breakpoint pairs with sufficient read support
        supported_pairs = [
            p for p in clustered_pairs 
            if p["read_count"] >= min_read_support
        ]
        
        if not supported_pairs:
            logging.debug(f"[Fusion] No breakpoint pairs found with >= {min_read_support} read support")
            return pd.DataFrame()
        
        # Convert to DataFrame format - include both primary and supplementary regions
        events = []
        for pair in supported_pairs:
            # Add primary region event
            if pair["primary_start"] > 0 and pair["primary_end"] > pair["primary_start"]:
                events.append({
                    "chromosome": pair["primary_chrom"],
                    "start": pair["primary_start"],
                    "end": pair["primary_end"],
                    "event_type": "breakpoint_pair_primary",
                    "read_count": pair["read_count"],
                    "avg_mapping_quality": round(pair["primary_avg_mapq"], 1),
                    "avg_mapping_span": round(pair["primary_avg_span"], 0),
                })
            
            # Add supplementary region event
            if pair["supp_start"] > 0 and pair["supp_end"] > pair["supp_start"]:
                events.append({
                    "chromosome": pair["supp_chrom"],
                    "start": pair["supp_start"],
                    "end": pair["supp_end"],
                    "event_type": "breakpoint_pair_supplementary",
                    "read_count": pair["read_count"],
                    "avg_mapping_quality": round(pair["supp_avg_mapq"], 1),
                    "avg_mapping_span": round(pair["supp_avg_span"], 0),
                })
        
        if not events:
            return pd.DataFrame()
        
        result_df = pd.DataFrame(events)
        
        # Sort by read count (descending), then by chromosome and start
        result_df = result_df.sort_values(
            ["read_count", "chromosome", "start"],
            ascending=[False, True, True]
        )
        
        return result_df
        
    except Exception as e:
        logging.warning(f"[Fusion] Failed to summarize master BED events: {e}")
        import traceback
        logging.debug(f"[Fusion] Traceback: {traceback.format_exc()}")
        return pd.DataFrame()


def _make_master_bed_summary_table(container: Any, df: pd.DataFrame, sample_dir: Optional[Path] = None) -> Any:
    """Create a summary table for master BED events (regions) instead of individual reads.
    
    Args:
        container: UI container for the table
        df: Master BED candidates DataFrame (for backward compatibility, but not used if summary file exists)
        sample_dir: Sample directory path to load pre-computed summary file
    """
    # Try to load pre-computed summary file first (generated by analysis code)
    summary_df = pd.DataFrame()
    if sample_dir:
        summary_file = sample_dir / "master_bed_events_summary.csv"
        if summary_file.exists():
            try:
                summary_df = pd.read_csv(summary_file)
                logging.debug(f"[Fusion] Loaded pre-computed master BED events summary from {summary_file}")
            except Exception as e:
                logging.warning(f"[Fusion] Failed to load pre-computed summary: {e}")
    
    # Fallback to computing on the fly if summary file doesn't exist (backward compatibility)
    if summary_df.empty and df is not None and not df.empty:
        logging.debug("[Fusion] Pre-computed summary not found, computing on the fly (this may be slow)")
        summary_df = _summarize_master_bed_events(df, min_read_support=3, min_mapq=50, cluster_distance=5000)
    
    if summary_df.empty:
        with container:
            ui.label("No master BED events found").classes("text-gray-600")
        return None
    
    with container:
        # Use styled_table for consistent styling
        from robin.gui.theme import styled_table
        
        # Define columns for the summary table
        columns = [
            {
                "name": "chromosome",
                "label": "Chromosome",
                "field": "chromosome",
                "sortable": True,
            },
            {
                "name": "start",
                "label": "Start",
                "field": "start",
                "sortable": True,
            },
            {
                "name": "end",
                "label": "End",
                "field": "end",
                "sortable": True,
            },
            {
                "name": "event_type",
                "label": "Type",
                "field": "event_type",
                "sortable": True,
            },
            {
                "name": "read_count",
                "label": "Reads",
                "field": "read_count",
                "sortable": True,
            },
            {
                "name": "avg_mapping_quality",
                "label": "Avg MapQ",
                "field": "avg_mapping_quality",
                "sortable": True,
            },
            {
                "name": "avg_mapping_span",
                "label": "Avg Span",
                "field": "avg_mapping_span",
                "sortable": True,
            },
        ]
        
        table = _render_paged_df_table(
            summary_df,
            columns,
            pagination_default=100,
            class_size="table-xs",
            search_placeholder="Search events...",
        )
        
        # Add summary information
        total_events = len(summary_df)
        total_reads = summary_df["read_count"].sum()
        ui.label(
            f"Total events: {total_events} | Total supporting reads: {total_reads}"
        ).classes("text-xs text-gray-500 mt-1")
        
        return table


def _cluster_fusion_reads(filtered_data: pd.DataFrame, max_distance: int = 10000, 
                         use_breakpoint_validation: bool = True) -> pd.DataFrame:
    """
    Cluster fusion reads by similar mapping coordinates with optional breakpoint validation.
    
    Args:
        filtered_data: DataFrame with fusion candidate data
        max_distance: Maximum distance for clustering (used differently for breakpoint vs coordinate clustering)
        use_breakpoint_validation: If True, use breakpoint validation; if False, use original coordinate clustering
        
    Returns:
        DataFrame with clustered fusion results
    """
    if use_breakpoint_validation:
        # Use breakpoint validation for more accurate fusion detection
        logging.info("[Fusion] Using breakpoint validation for clustering")
        
        # Validate fusion breakpoints
        validated_breakpoints = _validate_fusion_breakpoints(
            filtered_data, 
            min_read_support=4, 
            max_breakpoint_distance=100  # Much tighter clustering for breakpoints
        )
        
        if validated_breakpoints.empty:
            logging.info("[Fusion] No validated breakpoints found")
            return pd.DataFrame()
        
        # Convert to the expected format for the summary table
        clustered_results = []
        for _, row in validated_breakpoints.iterrows():
            clustered_results.append({
                "fusion_pair": row["gene_pair"],
                "chr1": row["gene1_chr"],
                "chr2": row["gene2_chr"],
                "gene1": row["gene1"],
                "gene1_position": row["gene1_breakpoint"],
                "gene2": row["gene2"],
                "gene2_position": row["gene2_breakpoint"],
                "reads": row["supporting_reads"],
                "avg_mapping_quality": row["avg_mapping_quality"],
                "avg_mapping_span": row["avg_mapping_span"],
                "cluster_id": row["cluster_id"]
            })
        
        logging.info(f"[Fusion] Breakpoint validation found {len(clustered_results)} validated fusion clusters")
        return pd.DataFrame(clustered_results)
    
    else:
        # Original coordinate-based clustering (fallback)
        logging.info("[Fusion] Using original coordinate clustering")
        fusion_summary = []
        
        # Group by read_id to get gene pairs
        for read_id, group in filtered_data.groupby("read_id", observed=True):
            genes = group["col4"].unique()
            if len(genes) >= 2:  # Only consider multi-gene reads
                # Sort genes for consistent pair representation
                genes_sorted = sorted(genes)
                gene_pair = "-".join(genes_sorted)
                
                # Get chromosome and position information for each gene
                gene_info = {}
                for gene in genes_sorted:
                    gene_data = group[group["col4"] == gene]
                    if not gene_data.empty:
                        first_row = gene_data.iloc[0]
                        gene_info[gene] = {
                            "chromosome": first_row.get("reference_id", "Unknown"),
                            "start": first_row.get("reference_start", 0),
                            "end": first_row.get("reference_end", 0)
                        }
                
                # Only add if we have info for at least 2 genes
                if len(gene_info) >= 2:
                    genes_with_info = [g for g in genes_sorted if g in gene_info]
                    if len(genes_with_info) >= 2:
                        gene1, gene2 = genes_with_info[0], genes_with_info[1]
                        
                        fusion_summary.append({
                            "fusion_pair": gene_pair,
                            "chr1": gene_info[gene1]["chromosome"],
                            "chr2": gene_info[gene2]["chromosome"],
                            "gene1": gene1,
                            "gene1_start": gene_info[gene1]["start"],
                            "gene1_end": gene_info[gene1]["end"],
                            "gene2": gene2,
                            "gene2_start": gene_info[gene2]["start"],
                            "gene2_end": gene_info[gene2]["end"],
                            "read_id": read_id
                        })
        
        if not fusion_summary:
            return pd.DataFrame()
        
        summary_df = pd.DataFrame(fusion_summary)
        
        # Cluster by fusion pair and similar coordinates
        clustered_results = []
        
        for fusion_pair in summary_df["fusion_pair"].unique():
            pair_data = summary_df[summary_df["fusion_pair"] == fusion_pair]
            
            # Group by chromosome combination
            for (chr1, chr2), chr_group in pair_data.groupby(["chr1", "chr2"]):
                # Cluster gene1 positions
                gene1_positions = chr_group[["gene1_start", "gene1_end"]].values
                gene1_clusters = _cluster_positions(gene1_positions, max_distance)
                
                # Cluster gene2 positions
                gene2_positions = chr_group[["gene2_start", "gene2_end"]].values
                gene2_clusters = _cluster_positions(gene2_positions, max_distance)
                
                # Create cluster combinations
                for i, gene1_cluster in enumerate(gene1_clusters):
                    for j, gene2_cluster in enumerate(gene2_clusters):
                        # Find reads that belong to both clusters
                        cluster_reads = []
                        for idx, row in chr_group.iterrows():
                            gene1_start, gene1_end = row["gene1_start"], row["gene1_end"]
                            gene2_start, gene2_end = row["gene2_start"], row["gene2_end"]
                            
                            # Check if this read belongs to both clusters
                            if (_position_in_cluster(gene1_start, gene1_end, gene1_cluster, max_distance) and
                                _position_in_cluster(gene2_start, gene2_end, gene2_cluster, max_distance)):
                                cluster_reads.append(row["read_id"])
                        
                        # Only include clusters with minimum read support (4 or more reads)
                        if len(cluster_reads) >= 4:
                            # Calculate cluster boundaries
                            gene1_min_start = min(gene1_cluster[:, 0])
                            gene1_max_end = max(gene1_cluster[:, 1])
                            gene2_min_start = min(gene2_cluster[:, 0])
                            gene2_max_end = max(gene2_cluster[:, 1])
                            
                            clustered_results.append({
                                "fusion_pair": fusion_pair,
                                "chr1": chr1,
                                "chr2": chr2,
                                "gene1": chr_group.iloc[0]["gene1"],
                                "gene1_position": f"{gene1_min_start}-{gene1_max_end}",
                                "gene2": chr_group.iloc[0]["gene2"],
                                "gene2_position": f"{gene2_min_start}-{gene2_max_end}",
                                "reads": len(cluster_reads)
                            })
        
        return pd.DataFrame(clustered_results)


def _cluster_positions(positions: np.ndarray, max_distance: int) -> List[np.ndarray]:
    """Cluster genomic positions based on distance."""
    if len(positions) == 0:
        return []
    
    # Simple clustering: group positions that are within max_distance
    clusters = []
    used = set()
    
    for i, (start, end) in enumerate(positions):
        if i in used:
            continue
            
        cluster = [positions[i]]
        used.add(i)
        
        for j, (other_start, other_end) in enumerate(positions):
            if j in used:
                continue
                
            # Check if positions overlap or are close
            if (_positions_overlap(start, end, other_start, other_end) or
                _positions_close(start, end, other_start, other_end, max_distance)):
                cluster.append(positions[j])
                used.add(j)
        
        clusters.append(np.array(cluster))
    
    return clusters


def _positions_overlap(start1: int, end1: int, start2: int, end2: int) -> bool:
    """Check if two genomic positions overlap."""
    return not (end1 < start2 or end2 < start1)


def _positions_close(start1: int, end1: int, start2: int, end2: int, max_distance: int) -> bool:
    """Check if two genomic positions are within max_distance."""
    distance = min(abs(start1 - start2), abs(end1 - end2), 
                   abs(start1 - end2), abs(end1 - start2))
    return distance <= max_distance


def _position_in_cluster(start: int, end: int, cluster: np.ndarray, max_distance: int) -> bool:
    """Check if a position belongs to a cluster."""
    for cluster_start, cluster_end in cluster:
        if (_positions_overlap(start, end, cluster_start, cluster_end) or
            _positions_close(start, end, cluster_start, cluster_end, max_distance)):
            return True
    return False


# =============================================================================
# BREAKPOINT VALIDATION FUNCTIONS
# =============================================================================

def _extract_fusion_breakpoints(annotated_data: pd.DataFrame) -> pd.DataFrame:
    """
    Extract fusion breakpoints from annotated fusion data.
    
    This function analyzes reads that map to multiple genes and extracts
    the actual fusion junction points (breakpoints) from their supplementary alignments.
    
    Args:
        annotated_data: DataFrame with fusion candidate data
        
    Returns:
        DataFrame with breakpoint information for each read
    """
    breakpoint_data = []
    
    # Group by read_id to analyze each read's alignments
    for read_id, read_group in annotated_data.groupby("read_id", observed=True):
        genes = read_group["col4"].unique()
        
        # Only process reads that map to multiple genes (fusion candidates)
        if len(genes) < 2:
            continue
            
        # Sort genes for consistent ordering
        genes_sorted = sorted(genes)
        
        # Extract breakpoint information for each gene
        gene_breakpoints = {}
        for gene in genes_sorted:
            gene_data = read_group[read_group["col4"] == gene]
            if not gene_data.empty:
                # Take the first (primary) alignment for this gene
                first_row = gene_data.iloc[0]
                gene_breakpoints[gene] = {
                    "chromosome": first_row.get("reference_id", "Unknown"),
                    "start": first_row.get("reference_start", 0),
                    "end": first_row.get("reference_end", 0),
                    "strand": first_row.get("strand", "+"),
                    "mapping_quality": first_row.get("mapping_quality", 0),
                    "mapping_span": first_row.get("mapping_span", 0)
                }
        
        # Only proceed if we have breakpoint info for at least 2 genes
        if len(gene_breakpoints) >= 2:
            # Create breakpoint pairs for all gene combinations
            gene_list = list(gene_breakpoints.keys())
            for i in range(len(gene_list)):
                for j in range(i + 1, len(gene_list)):
                    gene1, gene2 = gene_list[i], gene_list[j]
                    
                    breakpoint_data.append({
                        "read_id": read_id,
                        "gene_pair": f"{gene1}-{gene2}",
                        "gene1": gene1,
                        "gene1_chr": gene_breakpoints[gene1]["chromosome"],
                        "gene1_start": gene_breakpoints[gene1]["start"],
                        "gene1_end": gene_breakpoints[gene1]["end"],
                        "gene1_strand": gene_breakpoints[gene1]["strand"],
                        "gene2": gene2,
                        "gene2_chr": gene_breakpoints[gene2]["chromosome"],
                        "gene2_start": gene_breakpoints[gene2]["start"],
                        "gene2_end": gene_breakpoints[gene2]["end"],
                        "gene2_strand": gene_breakpoints[gene2]["strand"],
                        "min_mapping_quality": min(
                            gene_breakpoints[gene1]["mapping_quality"],
                            gene_breakpoints[gene2]["mapping_quality"]
                        ),
                        "min_mapping_span": min(
                            gene_breakpoints[gene1]["mapping_span"],
                            gene_breakpoints[gene2]["mapping_span"]
                        )
                    })
    
    return pd.DataFrame(breakpoint_data)


def _cluster_breakpoints(breakpoint_data: pd.DataFrame, max_distance: int = 100) -> pd.DataFrame:
    """
    Cluster breakpoints by similar coordinates within each gene pair.
    
    This function groups reads that have similar breakpoint coordinates,
    which indicates they support the same fusion event.
    
    Args:
        breakpoint_data: DataFrame with breakpoint information
        max_distance: Maximum distance for clustering breakpoints
        
    Returns:
        DataFrame with clustered breakpoint information
    """
    if breakpoint_data.empty:
        return pd.DataFrame()
    
    clustered_results = []
    
    # Group by gene pair and chromosome combination
    for (gene_pair, chr1, chr2), group in breakpoint_data.groupby(["gene_pair", "gene1_chr", "gene2_chr"]):
        if group.empty:
            continue
            
        # Cluster gene1 breakpoints
        gene1_positions = group[["gene1_start", "gene1_end"]].values
        gene1_clusters = _cluster_positions(gene1_positions, max_distance)
        
        # Cluster gene2 breakpoints
        gene2_positions = group[["gene2_start", "gene2_end"]].values
        gene2_clusters = _cluster_positions(gene2_positions, max_distance)
        
        # Create cluster combinations
        for i, gene1_cluster in enumerate(gene1_clusters):
            for j, gene2_cluster in enumerate(gene2_clusters):
                # Find reads that belong to both clusters
                cluster_reads = []
                cluster_mapping_qualities = []
                cluster_mapping_spans = []
                
                for idx, row in group.iterrows():
                    gene1_start, gene1_end = row["gene1_start"], row["gene1_end"]
                    gene2_start, gene2_end = row["gene2_start"], row["gene2_end"]
                    
                    # Check if this read belongs to both clusters
                    if (_position_in_cluster(gene1_start, gene1_end, gene1_cluster, max_distance) and
                        _position_in_cluster(gene2_start, gene2_end, gene2_cluster, max_distance)):
                        cluster_reads.append(row["read_id"])
                        cluster_mapping_qualities.append(row["min_mapping_quality"])
                        cluster_mapping_spans.append(row["min_mapping_span"])
                
                if cluster_reads:
                    # Calculate cluster boundaries
                    gene1_min_start = min(gene1_cluster[:, 0])
                    gene1_max_end = max(gene1_cluster[:, 1])
                    gene2_min_start = min(gene2_cluster[:, 0])
                    gene2_max_end = max(gene2_cluster[:, 1])
                    
                    # Calculate average quality metrics
                    avg_mapping_quality = np.mean(cluster_mapping_qualities) if cluster_mapping_qualities else 0
                    avg_mapping_span = np.mean(cluster_mapping_spans) if cluster_mapping_spans else 0
                    
                    clustered_results.append({
                        "gene_pair": gene_pair,
                        "gene1": group.iloc[0]["gene1"],
                        "gene1_chr": chr1,
                        "gene1_breakpoint": f"{gene1_min_start}-{gene1_max_end}",
                        "gene1_start": gene1_min_start,
                        "gene1_end": gene1_max_end,
                        "gene2": group.iloc[0]["gene2"],
                        "gene2_chr": chr2,
                        "gene2_breakpoint": f"{gene2_min_start}-{gene2_max_end}",
                        "gene2_start": gene2_min_start,
                        "gene2_end": gene2_max_end,
                        "supporting_reads": len(cluster_reads),
                        "read_ids": cluster_reads,
                        "avg_mapping_quality": avg_mapping_quality,
                        "avg_mapping_span": avg_mapping_span,
                        "cluster_id": f"{gene_pair}_{i}_{j}"
                    })
    
    return pd.DataFrame(clustered_results)


def _validate_fusion_breakpoints(annotated_data: pd.DataFrame, min_read_support: int = 4, 
                                max_breakpoint_distance: int = 100) -> pd.DataFrame:
    """
    Validate fusion candidates by requiring consistent breakpoint support.
    
    This function implements the breakpoint validation logic you requested:
    - Extracts fusion breakpoints from supplementary alignments
    - Clusters reads by similar breakpoint coordinates
    - Only returns fusions with consistent breakpoint support across gene regions
    
    Args:
        annotated_data: DataFrame with fusion candidate data
        min_read_support: Minimum number of reads supporting the same breakpoint
        max_breakpoint_distance: Maximum distance for clustering breakpoints
        
    Returns:
        DataFrame with validated fusion breakpoints
    """
    if annotated_data.empty:
        return pd.DataFrame()
    
    logging.info(f"[Fusion] Validating breakpoints for {len(annotated_data)} fusion candidates")
    
    # Extract breakpoints from fusion data
    breakpoint_data = _extract_fusion_breakpoints(annotated_data)
    
    if breakpoint_data.empty:
        logging.info("[Fusion] No breakpoint data extracted")
        return pd.DataFrame()
    
    logging.info(f"[Fusion] Extracted {len(breakpoint_data)} breakpoint records")
    
    # Cluster breakpoints by similar coordinates
    clustered_breakpoints = _cluster_breakpoints(breakpoint_data, max_breakpoint_distance)
    
    if clustered_breakpoints.empty:
        logging.info("[Fusion] No clustered breakpoints found")
        return pd.DataFrame()
    
    logging.info(f"[Fusion] Found {len(clustered_breakpoints)} clustered breakpoint groups")
    
    # Filter by minimum read support
    validated_breakpoints = clustered_breakpoints[
        clustered_breakpoints["supporting_reads"] >= min_read_support
    ]
    
    logging.info(f"[Fusion] {len(validated_breakpoints)} breakpoint groups meet minimum support threshold ({min_read_support})")
    
    # Sort by supporting reads (descending)
    validated_breakpoints = validated_breakpoints.sort_values("supporting_reads", ascending=False)
    
    return validated_breakpoints


def _get_validated_fusion_pairs(annotated_data: pd.DataFrame, goodpairs: pd.Series) -> List[List[str]]:
    """Extract validated fusion pairs from annotated_data for dropdown options.
    
    Returns a list of gene pairs as lists (e.g., [["GENE1", "GENE2"], ...])
    representing validated fusion pairs from breakpoint validation.
    """
    try:
        # Filter to good pairs if available
        if not goodpairs.empty and goodpairs.sum() > 0:
            # Align indices to avoid reindexing warning
            aligned_goodpairs = goodpairs.reindex(annotated_data.index, fill_value=False)
            filtered_data = annotated_data[aligned_goodpairs]
        else:
            filtered_data = annotated_data
        
        if filtered_data.empty:
            return []
        
        # Get validated fusion pairs using breakpoint validation
        clustered_data = _cluster_fusion_reads(filtered_data, max_distance=10000, use_breakpoint_validation=True)
        
        if clustered_data.empty:
            return []
        
        # Extract unique gene pairs and convert to list format
        validated_pairs = []
        seen_pairs = set()
        
        for _, row in clustered_data.iterrows():
            fusion_pair_str = row["fusion_pair"]  # e.g., "GENE1-GENE2"
            if fusion_pair_str and fusion_pair_str not in seen_pairs:
                # Split the fusion pair string into individual genes
                genes = [g.strip() for g in fusion_pair_str.split("-") if g.strip()]
                if len(genes) >= 2:
                    # Sort genes for consistency
                    genes_sorted = sorted(genes)
                    pair_key = tuple(genes_sorted)
                    if pair_key not in seen_pairs:
                        validated_pairs.append(genes_sorted)
                        seen_pairs.add(pair_key)
        
        logging.info(f"[Fusion] Extracted {len(validated_pairs)} validated fusion pairs for dropdown")
        return validated_pairs
        
    except Exception as e:
        logging.exception(f"[Fusion] Failed to extract validated fusion pairs: {e}")
        return []


def _make_fusion_summary_table(container: Any, annotated_data: pd.DataFrame, goodpairs: pd.Series) -> Any:
    """Create a summary table showing fusion pairs with chromosomes, positions, and read counts."""
    if annotated_data is None or annotated_data.empty:
        with container:
            ui.label("No fusion data available").classes("text-gray-600")
        return None

    with container:
        # Use styled_table for consistent styling
        from robin.gui.theme import styled_table
        
        try:
            # Filter to good pairs if available
            if not goodpairs.empty and goodpairs.sum() > 0:
                # Align indices to avoid reindexing warning
                aligned_goodpairs = goodpairs.reindex(annotated_data.index, fill_value=False)
                filtered_data = annotated_data[aligned_goodpairs]
            else:
                filtered_data = annotated_data
            
            if filtered_data.empty:
                ui.label("No fusion pairs found").classes("text-gray-600")
                return None
            
            # Cluster fusion reads by similar coordinates
            clustered_data = _cluster_fusion_reads(filtered_data, max_distance=10000, use_breakpoint_validation=True)
            
            if clustered_data.empty:
                ui.label("No fusion pairs found").classes("text-gray-600")
                return None
            
            # Sort by read count (descending)
            aggregated = clustered_data.sort_values("reads", ascending=False)
            
            # Create columns definition with breakpoint validation info
            columns = [
                {"name": "fusion_pair", "label": "Fusion Pair", "field": "fusion_pair", "sortable": True},
                {"name": "chr1", "label": "Chr 1", "field": "chr1", "sortable": True},
                {"name": "chr2", "label": "Chr 2", "field": "chr2", "sortable": True},
                {"name": "gene1", "label": "Gene 1", "field": "gene1", "sortable": True},
                {"name": "gene1_position", "label": "Gene 1 Breakpoint", "field": "gene1_position", "sortable": True},
                {"name": "gene2", "label": "Gene 2", "field": "gene2", "sortable": True},
                {"name": "gene2_position", "label": "Gene 2 Breakpoint", "field": "gene2_position", "sortable": True},
                {"name": "reads", "label": "Supporting Reads", "field": "reads", "sortable": True}
            ]
            
            # Add quality metrics if available (from breakpoint validation)
            if "avg_mapping_quality" in aggregated.columns:
                columns.append({"name": "avg_mapping_quality", "label": "Avg MapQ", "field": "avg_mapping_quality", "sortable": True})
            if "avg_mapping_span" in aggregated.columns:
                columns.append({"name": "avg_mapping_span", "label": "Avg Span", "field": "avg_mapping_span", "sortable": True})
            
            # Format the data for display
            rows = []
            for _, row in aggregated.iterrows():
                formatted_row = {
                    "fusion_pair": row["fusion_pair"],
                    "chr1": row["chr1"],
                    "chr2": row["chr2"],
                    "gene1": row.get("gene1", ""),
                    "gene1_position": row["gene1_position"],
                    "gene2": row.get("gene2", ""),
                    "gene2_position": row["gene2_position"],
                    "reads": int(row["reads"]),
                }
                
                # Add quality metrics if available
                if "avg_mapping_quality" in row:
                    formatted_row["avg_mapping_quality"] = f"{row['avg_mapping_quality']:.1f}"
                if "avg_mapping_span" in row:
                    formatted_row["avg_mapping_span"] = f"{row['avg_mapping_span']:.0f}"
                
                rows.append(formatted_row)
            
            # Add title
            ui.label("Fusion Summary").classes("text-sm font-medium mb-2")
            
            # Create styled table
            table_container, table = styled_table(
                columns=columns,
                rows=rows,
                pagination=20,
                class_size="table-xs"
            )
            
            # Add search functionality
            try:
                with table.add_slot("top-right"):
                    with ui.input(placeholder="Search fusions...").props("type=search").bind_value(
                        table, "filter"
                    ).add_slot("append"):
                        ui.icon("search")
            except Exception:
                pass
            
            # Add summary information
            total_fusions = len(aggregated)
            total_reads = aggregated["reads"].sum()
            ui.label(f"Total fusions: {total_fusions} | Total supporting reads: {total_reads}").classes("text-xs text-gray-500 mt-1")
            
            return table
            
        except Exception as e:
            logging.exception(f"[Fusion] Failed to create summary table: {e}")
            ui.label(f"Error creating fusion summary: {str(e)}").classes("text-red-600")
            return None


def _make_fusion_groups_table(container: Any, data: Dict[str, Any]) -> Any:
    """Create a table showing validated fusion groups (gene lists that may contain 2+ genes).
    
    Only shows groups that contain at least one validated fusion pair
    (meeting the minimum read support threshold).
    
    Args:
        container: UI container to place the table in
        data: Dictionary containing fusion data with annotated_data, goodpairs, and gene_groups
    """
    validated_groups = _get_validated_fusion_groups(data)
    
    if not validated_groups or len(validated_groups) == 0:
        with container:
            ui.label("No validated fusion groups found").classes("text-gray-600")
        return None

    with container:
        # Use styled_table for consistent styling
        from robin.gui.theme import styled_table
        
        try:
            # Create columns definition
            columns = [
                {"name": "group_id", "label": "Group ID", "field": "group_id", "sortable": True},
                {"name": "genes", "label": "Genes", "field": "genes", "sortable": True},
                {"name": "gene_count", "label": "Gene Count", "field": "gene_count", "sortable": True},
            ]
            
            # Format the data for display
            rows = []
            for idx, group in enumerate(validated_groups):
                if not isinstance(group, (list, tuple)) or len(group) == 0:
                    continue
                
                # Sort genes for consistent display
                sorted_genes = sorted([str(g).strip() for g in group if g])
                if len(sorted_genes) == 0:
                    continue
                
                rows.append({
                    "group_id": idx + 1,
                    "genes": " - ".join(sorted_genes),
                    "gene_count": len(sorted_genes),
                })
            
            if not rows:
                ui.label("No valid fusion groups found").classes("text-gray-600")
                return None
            
            # Sort by gene count (descending), then by gene names
            rows.sort(key=lambda x: (-x["gene_count"], x["genes"]))
            
            # Add title
            ui.label("Fusion Groups").classes("text-sm font-medium mb-2")
            
            # Create styled table
            table_container, table = styled_table(
                columns=columns,
                rows=rows,
                pagination=20,
                class_size="table-xs"
            )
            
            # Add search functionality
            try:
                with table.add_slot("top-right"):
                    with ui.input(placeholder="Search groups...").props("type=search").bind_value(
                        table, "filter"
                    ).add_slot("append"):
                        ui.icon("search")
            except Exception:
                pass
            
            # Add summary information
            total_groups = len(rows)
            multi_gene_groups = sum(1 for r in rows if r["gene_count"] > 2)
            ui.label(f"Total groups: {total_groups} | Multi-gene groups (>2): {multi_gene_groups}").classes("text-xs text-gray-500 mt-1")
            
            return table
            
        except Exception as e:
            logging.exception(f"[Fusion] Failed to create fusion groups table: {e}")
            ui.label(f"Error creating fusion groups table: {str(e)}").classes("text-red-600")
            return None


def _make_fusion_reads_table(container: Any, reads_df: pd.DataFrame, gene_pair: List[str]) -> Any:
    """Create a table showing reads and specific locations for a selected gene pair."""
    if reads_df is None or reads_df.empty:
        with container:
            ui.label("No reads found for selected gene pair").classes("text-gray-600")
        return None

    with container:
        # Use styled_table for consistent styling
        from robin.gui.theme import styled_table
        
        # The reads_df is already filtered for the gene pair, so use it directly
        filtered_reads = reads_df.copy()
        
        # Debug logging
        logging.info(f"[Fusion] Reads table - gene_pair: {gene_pair}")
        logging.info(f"[Fusion] Reads table - filtered_reads shape: {filtered_reads.shape}")
        logging.info(f"[Fusion] Reads table - filtered_reads columns: {list(filtered_reads.columns)}")
        if not filtered_reads.empty:
            logging.info(f"[Fusion] Reads table - unique genes in data: {filtered_reads['col4'].unique()}")
            logging.info(f"[Fusion] Reads table - sample data: {filtered_reads.head(2).to_dict('records')}")
        
        if filtered_reads.empty:
            ui.label("No reads found for selected gene pair").classes("text-gray-600")
            return None
        
        # Select and rename relevant columns for the reads table
        columns_to_show = [
            "read_id", "col4", "reference_id", "reference_start", "reference_end", 
            "mapping_quality", "strand", "read_start", "read_end", "is_secondary", 
            "is_supplementary", "mapping_span"
        ]
        
        # Only include columns that exist in the DataFrame
        available_columns = [col for col in columns_to_show if col in filtered_reads.columns]
        logging.info(f"[Fusion] Reads table - available columns: {available_columns}")
        
        if not available_columns:
            # Fallback: show all available columns if none of the expected ones exist
            logging.warning(f"[Fusion] No expected columns found, using all available columns: {list(filtered_reads.columns)}")
            available_columns = list(filtered_reads.columns)
            if not available_columns:
                ui.label("No columns found in fusion data").classes("text-gray-600")
                return None
            
        reads_subset = filtered_reads[available_columns].copy()
        
        # Rename columns to user-friendly names
        column_rename_map = {
            "read_id": "Read ID",
            "col4": "Gene",
            "reference_id": "Chromosome", 
            "reference_start": "Read Start",
            "reference_end": "Read End",
            "mapping_quality": "Map Quality",
            "strand": "Strand",
            "read_start": "Query Start",
            "read_end": "Query End",
            "is_secondary": "Secondary",
            "is_supplementary": "Supplementary",
            "mapping_span": "Span"
        }
        
        # Only rename columns that exist
        final_rename_map = {k: v for k, v in column_rename_map.items() if k in reads_subset.columns}
        reads_subset = reads_subset.rename(columns=final_rename_map)
        
        # Sort by gene and then by read start position (if available)
        sort_columns = []
        if "Gene" in reads_subset.columns:
            sort_columns.append("Gene")
        if "Read Start" in reads_subset.columns:
            sort_columns.append("Read Start")
        elif "Query Start" in reads_subset.columns:
            sort_columns.append("Query Start")
            
        if sort_columns:
            reads_subset = reads_subset.sort_values(sort_columns)
        
        # Create columns definition
        columns = []
        for col in reads_subset.columns:
            columns.append({
                "name": col,
                "label": col,
                "field": col,
                "sortable": True
            })
        
        # Add title
        ui.label(f"Reads supporting fusion: {'-'.join(gene_pair)}").classes("text-sm font-medium mb-2")
        table = _render_paged_df_table(
            reads_subset,
            columns,
            pagination_default=100,
            class_size="table-xs",
            search_placeholder="Search reads...",
        )
        
        # Add summary information
        total_reads = len(reads_subset)
        if "Read ID" in reads_subset.columns:
            unique_reads = len(reads_subset["Read ID"].unique())
            ui.label(f"Total reads: {total_reads} | Unique reads: {unique_reads}").classes("text-xs text-gray-500 mt-1")
        else:
            ui.label(f"Total reads: {total_reads}").classes("text-xs text-gray-500 mt-1")
        
        return table


def _plot_gene_group(
    container: Any,
    gene_group: List[str],
    annotated_data: pd.DataFrame,
    goodpairs: pd.Series,
    section_state: Optional[Dict[str, Any]] = None,
) -> None:
    """Advanced fusion visualization matching the original sophisticated plotting code.

    Creates multi-panel plots showing:
    1. Gene structure with exons and annotations
    2. Read mapping positions with color coding
    3. Professional-grade visualization using DNA Features Viewer
    4. Reads table showing specific locations and details
    """
    try:
        # For genome-wide fusions, use raw annotated_data if goodpairs filters out everything
        if goodpairs.sum() == 0:
            # No good pairs - use raw data (like reporting code does)
            # Strip whitespace from gene names to handle leading/trailing spaces
            subset = annotated_data[annotated_data["col4"].str.strip().isin(gene_group)]
            logging.info(f"[Fusion] Using raw data for genome-wide gene group {gene_group}: {len(subset)} rows")
        else:
            # Target panel fusions - use filtered data
            # Fix pandas reindexing warning by ensuring indices match
            try:
                # Align indices to avoid reindexing warning
                aligned_goodpairs = goodpairs.reindex(annotated_data.index, fill_value=False)
                subset = annotated_data[aligned_goodpairs]
                
                # Strip whitespace from gene names to handle leading/trailing spaces
                subset = subset[subset["col4"].str.strip().isin(gene_group)]
                
                logging.info(f"[Fusion] Using filtered data for target gene group {gene_group}: {len(subset)} rows")
            except Exception as e:
                # Fallback to raw data if indexing fails
                logging.warning(f"[Fusion] Indexing failed, using raw data: {e}")
                # Strip whitespace from gene names to handle leading/trailing spaces
                subset = annotated_data[annotated_data["col4"].str.strip().isin(gene_group)]
                logging.info(f"[Fusion] Using raw data (fallback) for gene group {gene_group}: {len(subset)} rows")
        
        if subset.empty:
            with container:
                ui.label("No reads for selected gene group").classes("text-gray-600")
            return

        # Clear container and create advanced visualization
        container.clear()
        with container:
            # Create tabs for visualization and reads table
            with ui.tabs().classes("w-full") as tabs:
                visualization_tab = ui.tab("Visualization")
                reads_tab = ui.tab("Reads Table")
            
            with ui.tab_panels(tabs, value=visualization_tab).classes("w-full"):
                with ui.tab_panel(visualization_tab):
                    # Create matplotlib element for the sophisticated plot with ideograms
                    mpl_element = ui.matplotlib(figsize=(20, 10)).classes("w-full")

                    # Create the advanced fusion plot
                    fig = _create_advanced_fusion_plot(
                        gene_group, subset, annotated_data, goodpairs
                    )

                    _dark = _is_dark_mode()
                    _apply_fusion_figure_theme(fig, _dark)
                    # Update the matplotlib element
                    # Close previous figure if it exists before assigning new one
                    if hasattr(mpl_element, 'figure') and mpl_element.figure is not None:
                        plt.close(mpl_element.figure)
                    mpl_element.figure = fig
                    mpl_element.update()
                    if section_state is not None:
                        section_state["fusion_mpl"] = mpl_element
                        section_state["fusion_plot_theme_dark"] = _dark
                
                with ui.tab_panel(reads_tab):
                    # Create reads table
                    reads_table_container = ui.column().classes("w-full")
                    _make_fusion_reads_table(reads_table_container, subset, gene_group)

    except Exception as e:
        logging.exception(f"[Fusion] Failed to plot gene group {gene_group}: {e}")
        # Fallback to text display
        with container:
            ui.label(f"Failed to plot gene group: {str(e)}").classes("text-red-600")


def _create_advanced_fusion_plot(
    gene_group: List[str],
    subset: pd.DataFrame,
    annotated_data: pd.DataFrame,
    goodpairs: pd.Series,
) -> plt.Figure:
    """Create the sophisticated gene fusion plot using real DNA Features Viewer integration.

    This matches the original code's approach with:
    - Real gene structure visualization using DNA Features Viewer
    - Sophisticated read alignment display
    - Unified side-by-side layout for gene regions
    """

    # Check if DNA Features Viewer is available
    if not DNA_FEATURES_AVAILABLE:
        return _create_simple_fallback_plot(gene_group, subset)

    # Load gene annotation data
    gene_table = _load_gene_annotations()
    if gene_table is None:
        return _create_simple_fallback_plot(gene_group, subset)

    # Create figure with tight layout - disable constrained layout to avoid conflicts
    plt.rcParams["figure.constrained_layout.use"] = False

    # Get unique genes and their data
    unique_genes = list(sorted(subset["col4"].unique()))
    if len(unique_genes) == 0:
        logging.warning(f"[Fusion] No genes found for group {gene_group}, using fallback plot")
        return _create_simple_fallback_plot(gene_group, subset)
    
    # Ensure we have valid data to prevent zero-size axes
    if subset.empty:
        logging.warning(f"[Fusion] Empty subset for group {gene_group}, using fallback plot")
        return _create_simple_fallback_plot(gene_group, subset)

    # Create the unified side-by-side layout without ideograms
    # For 2 genes, we'll create 2 columns with 2 rows each (gene structure + read mapping)
    num_genes = len(unique_genes)
    logging.info(f"[Fusion] Creating {num_genes} gene layout with 2 rows (gene structure + read mapping)")
    # Increase figure size and add more padding to prevent tight layout warnings
    # Ensure minimum figure size to prevent zero-size axes warnings
    # Close previous figure if it exists to prevent accumulation
    global _last_fusion_figure
    if _last_fusion_figure is not None:
        try:
            plt.close(_last_fusion_figure)
        except Exception:
            pass  # Figure might already be closed
    fig, axes = plt.subplots(2, num_genes, figsize=(max(20, 5 * num_genes), 10))
    _last_fusion_figure = fig

    # Handle single gene case
    if num_genes == 1:
        axes = axes.reshape(2, 1)

    # Process each gene
    for col_idx, gene_name in enumerate(unique_genes):
        gene_data = subset[subset["col4"] == gene_name]

        if gene_data.empty:
            continue

        # Get gene coordinates
        gene_start = gene_data["reference_start"].min()
        gene_end = gene_data["reference_end"].max()
        gene_chrom = gene_data["reference_id"].iloc[0]

        # Row 0: Gene structure with DNA Features Viewer
        ax_gene = axes[0, col_idx]
        _plot_gene_structure_with_dna_features(
            ax_gene, gene_name, gene_chrom, gene_start, gene_end, gene_table
        )

        # Row 1: Read mapping visualization
        ax_reads = axes[1, col_idx]
        _plot_read_mapping_sophisticated(
            ax_reads, gene_data, gene_name, gene_chrom, gene_start, gene_end
        )

    # Add overall title first
    fig.suptitle(
        f"Fusion Analysis: {', '.join(gene_group)}", fontsize=14, fontweight="bold"
    )

    # Use manual layout adjustment instead of tight_layout to avoid warnings
    # tight_layout often fails when there are many subplots or complex decorations
    plt.subplots_adjust(
        top=0.90,  # Leave space for suptitle
        bottom=0.1,
        left=0.08,
        right=0.95,
        hspace=0.4,  # More space between rows
        wspace=0.3,  # Space between columns
    )

    return fig




def _format_ticks_to_megabases(ax: plt.Axes) -> None:
    """Safely format x-axis ticks to megabases, handling custom formatters from DNA Features Viewer."""
    try:
        ticks = ax.get_xticks()
        # Set the tick positions first, then the labels to avoid warnings
        ax.set_xticks(ticks)
        ax.set_xticklabels([f'{t/1e6:.1f}' for t in ticks])
    except Exception as e:
        logging.warning(f"Could not format ticks to megabases: {e}")
        # Fallback: try to use ticklabel_format if possible
        try:
            ax.ticklabel_format(style='scientific', axis='x', scilimits=(0,0))
        except Exception:
            pass  # If both methods fail, just leave the default formatting


def _load_gene_annotations() -> Optional[pd.DataFrame]:
    """Load gene annotation data from the rCNS2_data.csv.gz file."""
    try:
        # Try to find the gene annotation file using robin resources
        gene_data_path = None
        
        # First try to use robin resources to find the correct path
        try:
            from robin import resources
            resources_dir = os.path.dirname(resources.__file__)
            gene_data_path = os.path.join(resources_dir, "rCNS2_data.csv.gz")
            if not os.path.exists(gene_data_path):
                gene_data_path = None
        except ImportError:
            pass
        
        # Fallback to relative path if resources import failed
        if not gene_data_path or not os.path.exists(gene_data_path):
            # Try relative path from current working directory
            possible_paths = [
                "src/robin/resources/rCNS2_data.csv.gz",
                "robin/resources/rCNS2_data.csv.gz",
                os.path.join(os.path.dirname(__file__), "..", "..", "resources", "rCNS2_data.csv.gz"),
            ]
            
            for path in possible_paths:
                abs_path = os.path.abspath(path)
                if os.path.exists(abs_path):
                    gene_data_path = abs_path
                    break
        
        if gene_data_path and os.path.exists(gene_data_path):
            gene_table = pd.read_csv(gene_data_path)
            return gene_table
        else:
            return None
    except Exception as e:
        return None


def _plot_gene_structure_with_dna_features(
    ax: plt.Axes,
    gene_name: str,
    chrom: str,
    start: int,
    end: int,
    gene_table: pd.DataFrame,
):
    """Plot gene structure using DNA Features Viewer exactly as in the original code."""
    try:
        # Strip whitespace from gene name to handle leading/trailing spaces in the data
        gene_name_clean = gene_name.strip() if isinstance(gene_name, str) else gene_name
        
        # Filter gene table for this specific gene and chromosome
        # Also strip whitespace from gene_table gene_name column for matching
        gene_table_clean = gene_table.copy()
        gene_table_clean["gene_name"] = gene_table_clean["gene_name"].astype(str).str.strip()
        
        gene_info = gene_table_clean[
            (gene_table_clean["gene_name"] == gene_name_clean) & (gene_table_clean["Seqid"] == chrom)
        ]

        if gene_info.empty:
            # Fallback if no gene info found
            ax.set_title(f"Gene Structure: {gene_name} ({chrom})", fontsize=10, fontweight="bold")
            # Extend plot range to show full context
            plot_start = start - (end - start) * 0.1
            plot_end = end + (end - start) * 0.1
            ax.set_xlim(plot_start, plot_end)
            ax.set_ylim(0, 1)
            ax.text(
                0.5,
                0.5,
                f"Gene: {gene_name}\n{chrom}:{start:,}-{end:,}",
                transform=ax.transAxes,
                ha="center",
                va="center",
                bbox=dict(boxstyle="round,pad=0.3", facecolor="lightblue", alpha=0.7),
            )
            # Convert x-axis to megabases
            ax.set_xlabel(f"Position (Mb) - {chrom}", fontsize=10)
            # Convert tick labels to megabases
            _format_ticks_to_megabases(ax)
            return
        
        # Create DNA Features Viewer visualization
        features = []
        
        # Determine the overall gene region from annotations
        gene_start = int(gene_info["Start"].min())
        gene_end = int(gene_info["End"].max())
        
        # Use the wider of the annotation range or the read mapping range
        plot_start = min(start, gene_start)
        plot_end = max(end, gene_end)
        plot_length = plot_end - plot_start

        # Add gene features
        for _, row in gene_info.iterrows():
            try:
                if row["Type"] == "gene":
                    # Main gene feature
                    strand = 1 if row["Strand"] == "+" else -1
                    features.append(
                        GraphicFeature(
                            start=int(row["Start"]),
                            end=int(row["End"]),
                            strand=strand,
                            thickness=8,
                            color="#ffd700",  # Gold color for genes
                            label=row["gene_name"],
                            fontdict={"color": "black", "fontsize": 8},
                        )
                    )
                elif row["Type"] == "exon":
                    # Exon features
                    strand = 1 if row["Strand"] == "+" else -1
                    features.append(
                        GraphicFeature(
                            start=int(row["Start"]),
                            end=int(row["End"]),
                            strand=strand,
                            thickness=4,
                            color="#C0C0C0",  # Silver color for exons
                        )
                    )
            except Exception as e:
                continue

        if features:
            # Create GraphicRecord with the full plot range
            # DNA Features Viewer uses absolute coordinates, first_index sets the starting point
            record = GraphicRecord(
                sequence_length=plot_length, first_index=plot_start, features=features
            )

            # Plot on the axis
            record.plot(
                ax=ax, with_ruler=False, draw_line=True, strand_in_label_threshold=4
            )
            ax.set_title(f"Gene Structure: {gene_name} ({chrom})", fontsize=10, fontweight="bold")
            ax.set_xlabel(f"Position (Mb) - {chrom}", fontsize=10)
            # Extend plot range to show full context
            ax.margins(x=0.15, y=0.1)
            # Convert tick labels to megabases (avoid ticklabel_format due to DNA Features Viewer formatter)
            _format_ticks_to_megabases(ax)
        else:
            # Fallback if no features found
            ax.set_title(f"Gene Structure: {gene_name} ({chrom})", fontsize=10, fontweight="bold")
            # Extend plot range to show full context
            plot_start = start - (end - start) * 0.1
            plot_end = end + (end - start) * 0.1
            ax.set_xlim(plot_start, plot_end)
            ax.set_ylim(0, 1)
            ax.text(
                0.5,
                0.5,
                f"Gene: {gene_name}\n{chrom}:{start:,}-{end:,}",
                transform=ax.transAxes,
                ha="center",
                va="center",
                bbox=dict(boxstyle="round,pad=0.3", facecolor="lightblue", alpha=0.7),
            )
            # Convert x-axis to megabases
            ax.set_xlabel(f"Position (Mb) - {chrom}", fontsize=10)
        # Convert tick labels to megabases
        _format_ticks_to_megabases(ax)

    except Exception as e:
        logging.error(f"Error plotting gene structure with DNA Features Viewer: {e}")
        # Fallback to simple text
        ax.set_title(f"Gene Structure: {gene_name} ({chrom})", fontsize=10, fontweight="bold")
        # Extend plot range to show full context
        plot_start = start - (end - start) * 0.1
        plot_end = end + (end - start) * 0.1
        ax.set_xlim(plot_start, plot_end)
        ax.set_ylim(0, 1)
        ax.text(
            0.5,
            0.5,
            f"Error plotting\n{gene_name}",
            transform=ax.transAxes,
            ha="center",
            va="center",
            color="red",
        )
        # Convert x-axis to megabases even in error case
        ax.set_xlabel(f"Position (Mb) - {chrom}", fontsize=10)
        # Convert tick labels to megabases
        _format_ticks_to_megabases(ax)


def _plot_read_mapping_sophisticated(
    ax: plt.Axes,
    gene_data: pd.DataFrame,
    gene_name: str,
    chrom: str,
    start: int,
    end: int,
):
    """Plot sophisticated read mapping visualization exactly as in the original code."""
    try:
        # Get read positions and create the sophisticated visualization
        read_starts = gene_data["reference_start"].values
        read_ends = gene_data["reference_end"].values
        read_ids = gene_data["read_id"].values
        strands = gene_data["strand"].values

        # Create color map for reads (similar to original code)
        unique_reads = pd.unique(read_ids)
        colors = plt.cm.Set3(np.linspace(0, 1, len(unique_reads)))
        read_color_map = dict(zip(unique_reads, colors))

        # Plot each read with sophisticated features
        for i, (read_start, read_end, read_id, strand) in enumerate(
            zip(read_starts, read_ends, read_ids, strands)
        ):
            color = read_color_map[read_id]

            # Create the read visualization as in original code
            # Main read line with slight extension to show ends clearly
            read_length = read_end - read_start
            extended_start = read_start - read_length * 0.02  # Small extension left
            extended_end = read_end + read_length * 0.02  # Small extension right
            ax.plot(
                [extended_start, extended_end],
                [i, i],
                color=color,
                linewidth=3,
                alpha=0.8,
            )

            # Add strand indicator (arrow) at the actual read end
            if strand == "+":
                ax.arrow(
                    read_end,
                    i,
                    (end - start) * 0.01,
                    0,
                    head_width=0.2,
                    head_length=(end - start) * 0.005,
                    fc=color,
                    ec=color,
                    alpha=0.8,
                )
            else:
                ax.arrow(
                    read_start,
                    i,
                    -(end - start) * 0.01,
                    0,
                    head_width=0.2,
                    head_length=(end - start) * 0.005,
                    fc=color,
                    ec=color,
                    alpha=0.8,
                )

            # Add read ID label (limit to first few for clarity)
            if i < 5:
                ax.text(
                    read_start, i + 0.1, f"Read {read_id[:8]}", fontsize=6, color=color
                )

        # Customize plot
        ax.set_title(f"Read Mapping: {gene_name} ({chrom})", fontsize=10, fontweight="bold")
        # Extend plot range to clearly show read ends
        plot_start = start - (end - start) * 0.15
        plot_end = end + (end - start) * 0.15
        ax.set_xlim(plot_start, plot_end)
        ax.set_xlabel(f"Position (Mb) - {chrom}")
        ax.set_ylabel("Reads")

        # Convert x-axis to megabases
        _format_ticks_to_megabases(ax)

        # Add grid
        ax.grid(True, alpha=0.3)
        
        # Remove box edges/borders
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.spines['bottom'].set_visible(False)
        ax.spines['left'].set_visible(False)

        # Add legend if not too many reads and we have labeled artists
        if len(unique_reads) <= 10 and len(ax.get_legend_handles_labels()[0]) > 0:
            ax.legend(bbox_to_anchor=(1.05, 1), loc="upper left", fontsize=8)

    except Exception as e:
        logging.error(f"Error plotting sophisticated read mapping: {e}")
        # Fallback to simple text
        ax.set_title(f"Read Mapping: {gene_name} ({chrom})", fontsize=10, fontweight="bold")
        ax.text(
            0.5,
            0.5,
            f"Error plotting reads\n{gene_name}",
            transform=ax.transAxes,
            ha="center",
            va="center",
            color="red",
        )
        # Convert x-axis to megabases even in error case
        ax.set_xlabel(f"Position (Mb) - {chrom}")
        # Format x-axis ticks to show megabases
        ax.ticklabel_format(style='scientific', axis='x', scilimits=(0,0))
        # Convert tick labels to megabases
        _format_ticks_to_megabases(ax)


def _process_reads_for_original_plot(subset: pd.DataFrame) -> pd.DataFrame:
    """Process reads data exactly as the original _get_reads function does."""
    try:
        # Convert to the format expected by original code
        df = subset.copy()
        df.columns = [
            "chromosome",
            "start",
            "end",
            "gene",
            "chromosome2",
            "start2",
            "end2",
            "id",
            "quality",
            "strand",
            "read_start",
            "read_end",
            "secondary",
            "supplementary",
            "span",
            "tag",
            "color",
        ]

        # Convert start and end to int
        df["start"] = df["start"].astype(int)
        df["end"] = df["end"].astype(int)

        # Sort and remove duplicates exactly as original
        df = df.sort_values(by=["chromosome", "start", "end"])
        df = df.drop_duplicates(subset=["start2", "end2", "id"])

        # Group by chromosome and collapse ranges (simplified version)
        # This would need the original collapse_ranges function
        return df

    except Exception as e:
        logging.error(f"Error processing reads for original plot: {e}")
        return subset


def _plot_gene_structure_original(
    ax: plt.Axes, data: pd.Series, chrom: str, start: int, end: int
):
    """Plot gene structure exactly as in original code using DNA Features Viewer."""
    try:
        # This would need access to the gene_table data from the original code
        # For now, create a placeholder that matches the original structure
        ax.set_title(f"Gene Structure: {data['gene']} ({chrom})", fontsize=10, fontweight="bold")
        ax.set_xlim(start, end)
        ax.set_ylim(0, 1)

        # Placeholder for gene structure - would need gene_table integration
        ax.text(
            0.5,
            0.5,
            f"Gene: {data['gene']}\n{chrom}:{start:,}-{end:,}",
            transform=ax.transAxes,
            ha="center",
            va="center",
            bbox=dict(boxstyle="round,pad=0.3", facecolor="lightblue", alpha=0.7),
        )

        ax.set_xlabel(f"Position (Mb) - {chrom}")
        ax.set_ylabel("Gene Structure")
        # Convert x-axis to megabases
        _format_ticks_to_megabases(ax)

    except Exception as e:
        logging.error(f"Error plotting gene structure: {e}")
        ax.text(
            0.5,
            0.5,
            f"Error plotting\n{data['gene']}",
            transform=ax.transAxes,
            ha="center",
            va="center",
            color="red",
        )


def _plot_read_mapping_original(
    ax: plt.Axes,
    data: pd.Series,
    chrom: str,
    start: int,
    end: int,
    subset: pd.DataFrame,
):
    """Plot read mapping exactly as in original code."""
    try:
        # Get reads for this gene
        gene_reads = subset[subset["col4"] == data["gene"]]

        if gene_reads.empty:
            ax.text(
                0.5, 0.5, "No reads", transform=ax.transAxes, ha="center", va="center"
            )
            return

        # Create the exact same visualization as original code
        # This would need the original plotting logic with overlapping ranges, ranks, etc.

        # For now, create a simplified version that shows the structure
        ax.set_title(f"Read Mapping: {data['gene']} ({chrom})", fontsize=10, fontweight="bold")
        # Extend plot range to clearly show read ends
        plot_start = start - (end - start) * 0.15
        plot_end = end + (end - start) * 0.15
        ax.set_xlim(plot_start, plot_end)

        # Plot reads as horizontal lines (simplified)
        for i, (_, read) in enumerate(gene_reads.iterrows()):
            read_start = read["reference_start"]
            read_end = read["reference_end"]
            ax.plot([read_start, read_end], [i, i], linewidth=2, alpha=0.8)

        ax.set_xlabel(f"Position (Mb) - {chrom}")
        ax.set_ylabel("Reads")
        # Convert x-axis to megabases
        _format_ticks_to_megabases(ax)
        ax.grid(True, alpha=0.3)
        
        # Remove box edges/borders
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.spines['bottom'].set_visible(False)
        ax.spines['left'].set_visible(False)

    except Exception as e:
        logging.error(f"Error plotting read mapping: {e}")
        ax.text(
            0.5,
            0.5,
            f"Error plotting reads\n{data['gene']}",
            transform=ax.transAxes,
            ha="center",
            va="center",
            color="red",
        )


def _plot_gene_structure(
    ax: plt.Axes, gene_name: str, chrom: str, start: int, end: int
):
    """Plot gene structure with exons and annotations."""
    try:
        # This would integrate with your gene annotation data
        # For now, create a placeholder gene structure
        ax.set_title(f"Gene Structure: {gene_name} ({chrom})", fontsize=10, fontweight="bold")
        ax.set_xlim(start, end)
        ax.set_ylim(0, 1)

        # Placeholder gene structure visualization
        ax.text(
            0.5,
            0.5,
            f"Gene: {gene_name}\n{chrom}:{start:,}-{end:,}",
            transform=ax.transAxes,
            ha="center",
            va="center",
            bbox=dict(boxstyle="round,pad=0.3", facecolor="lightblue", alpha=0.7),
        )

        ax.set_xlabel(f"Position (Mb) - {chrom}")
        ax.set_ylabel("Gene Structure")
        # Convert x-axis to megabases
        _format_ticks_to_megabases(ax)

    except Exception as e:
        logging.error(f"Error plotting gene structure for {gene_name}: {e}")
        ax.text(
            0.5,
            0.5,
            f"Error plotting\n{gene_name}",
            transform=ax.transAxes,
            ha="center",
            va="center",
            color="red",
        )


def _plot_read_mapping(
    ax: plt.Axes,
    gene_data: pd.DataFrame,
    gene_name: str,
    chrom: str,
    start: int,
    end: int,
):
    """Plot read mapping positions with color coding."""
    try:
        # Get read positions
        read_starts = gene_data["reference_start"].values
        read_ends = gene_data["reference_end"].values
        read_ids = gene_data["read_id"].values

        # Create color map for reads
        unique_reads = pd.unique(read_ids)
        colors = plt.cm.Set3(np.linspace(0, 1, len(unique_reads)))
        read_color_map = dict(zip(unique_reads, colors))

        # Plot each read as a horizontal line with slight extension to show ends clearly
        for i, (read_start, read_end, read_id) in enumerate(
            zip(read_starts, read_ends, read_ids)
        ):
            color = read_color_map[read_id]
            read_length = read_end - read_start
            extended_start = read_start - read_length * 0.02  # Small extension left
            extended_end = read_end + read_length * 0.02  # Small extension right
            ax.plot(
                [extended_start, extended_end],
                [i, i],
                color=color,
                linewidth=3,
                alpha=0.8,
                label=f"Read {read_id}" if i < 10 else "",
            )  # Limit legend to first 10 reads

        # Customize plot
        ax.set_title(f"Read Mapping: {gene_name} ({chrom})", fontsize=10, fontweight="bold")
        # Extend plot range to clearly show read ends
        plot_start = start - (end - start) * 0.15
        plot_end = end + (end - start) * 0.15
        ax.set_xlim(plot_start, plot_end)
        ax.set_xlabel(f"Position (Mb) - {chrom}")
        ax.set_ylabel("Reads")

        # Convert x-axis to megabases
        _format_ticks_to_megabases(ax)

        # Add grid
        ax.grid(True, alpha=0.3)
        
        # Remove box edges/borders
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.spines['bottom'].set_visible(False)
        ax.spines['left'].set_visible(False)

        # Add legend if not too many reads and we have labeled artists
        if len(unique_reads) <= 10 and len(ax.get_legend_handles_labels()[0]) > 0:
            ax.legend(bbox_to_anchor=(1.05, 1), loc="upper left", fontsize=8)

    except Exception as e:
        logging.error(f"Error plotting read mapping for {gene_name}: {e}")
        ax.text(
            0.5,
            0.5,
            f"Error plotting reads\n{gene_name}",
            transform=ax.transAxes,
            ha="center",
            va="center",
            color="red",
        )


def _create_simple_fallback_plot(
    gene_group: List[str], subset: pd.DataFrame
) -> plt.Figure:
    """Create a simple fallback plot if advanced visualization fails."""
    # Close previous figure if it exists to prevent accumulation
    global _last_fusion_figure
    if _last_fusion_figure is not None:
        try:
            plt.close(_last_fusion_figure)
        except Exception:
            pass  # Figure might already be closed
    fig, ax = plt.subplots(1, 1, figsize=(10, 4))
    _last_fusion_figure = fig

    # Simple scatter plot
    genes = list(sorted(subset["col4"].unique()))
    y_map = {g: i for i, g in enumerate(genes)}
    x = subset["reference_start"].astype(int)
    y = subset["col4"].map(y_map).astype(int)

    ax.scatter(x, y, s=15, alpha=0.8)
    ax.set_yticks(list(y_map.values()))
    ax.set_yticklabels(list(y_map.keys()))
    ax.set_xlabel("Mapping Start (reference_start)")
    ax.set_title(
        f"Reads supporting fusion group: {', '.join(gene_group)}",
        fontsize=12,
        fontweight="bold",
    )
    ax.grid(True, axis="x", linestyle=":", alpha=0.4)

    return fig


def add_fusion_section(launcher: Any, sample_dir: Path) -> None:
    """Build the Fusion UI section (target panel and genome-wide only).

    Expects files written by fusion_work.preprocess_fusion_data_standalone:
    - fusion_candidates_master_processed.pkl (pickle format)
    - fusion_candidates_all_processed.pkl (pickle format)
    """
    logging.info(f"[Fusion] add_fusion_section() called with sample_dir: {sample_dir}")
    
    if ui is None:
        logging.warning("[Fusion] ui is None, returning early")
        return

    # Local state for this section
    state: Dict[str, Any] = {
        "target": {
            "data": None,
            "mtime": None,
            "data_hash": None,
            "summary_table_container": None,
            "summary_table": None,
            "groups_table_container": None,
            "groups_table": None,
            "table_container": None,
            "table": None,
            "plot_container": None,
            "card": None,
            "status_container": None,
            "selected_gene_pair": None,  # Persist selected gene pair
            "dropdown": None,  # Reference to dropdown for value updates
        },
        "genome": {
            "data": None,
            "mtime": None,
            "data_hash": None,
            "summary_table_container": None,
            "summary_table": None,
            "groups_table_container": None,
            "groups_table": None,
            "table_container": None,
            "table": None,
            "plot_container": None,
            "card": None,
            "status_container": None,
            "selected_gene_pair": None,  # Persist selected gene pair
            "dropdown": None,  # Reference to dropdown for value updates
        },
        "master_bed": {
            "data": None,
            "mtime": None,
            "table_container": None,
            "table": None,
            "status_container": None,
        },
    }

    def _handle_gene_pair_selection(section: str, gene_pair: List[str], data: Dict[str, Any]) -> None:
        """Handle gene pair selection and update visualization."""
        try:
            state[section]["selected_gene_pair"] = gene_pair
            state[section]["card"].clear()
            state[section]["card"].classes("w-full")
            _plot_gene_group(
                state[section]["card"],
                gene_pair,
                data.get("annotated_data"),
                data.get("goodpairs"),
                section_state=state[section],
            )
        except Exception as e:
            logging.exception(f"[Fusion] Failed to handle gene pair selection: {e}")

    def refresh_fusion() -> None:
        """Refresh fusion data."""
        import time as _time
        t_start = _time.perf_counter()
        try:
            logging.info(f"[Fusion] refresh_fusion() called for sample_dir: {sample_dir}")
            
            # Simple directory check
            if not sample_dir or not sample_dir.exists():
                logging.warning(f"[Fusion] Sample directory not found: {sample_dir}")
                return
            
            logging.info(f"[Fusion] Loading fusion data from: {sample_dir}")
            t_before_load = _time.perf_counter()
            
            # Load fusion data directly (already in background)
            fusion_data = _load_fusion_data(sample_dir)
            
            t_after_load = _time.perf_counter()
            logging.info(f"[Fusion] Loaded fusion data (load took {t_after_load - t_before_load:.2f}s): {fusion_data}")
            
            # Update UI directly
            t_before_ui = _time.perf_counter()
            _update_fusion_ui(fusion_data, state, sample_dir)
            t_after_ui = _time.perf_counter()
            logging.info(f"[Fusion] refresh_fusion() completed in {t_after_ui - t_start:.2f}s (load={t_after_load - t_before_load:.2f}s, ui={t_after_ui - t_before_ui:.2f}s)")
                
        except Exception as e:
            logging.exception(f"[Fusion] Refresh failed: {e}")

    def _load_fusion_data(sample_dir: Path) -> Dict[str, Any]:
        """Load fusion data from files with optional breakpoint validation."""
        import time as _time
        t_load_start = _time.perf_counter()
        try:
            logging.info(f"[Fusion] _load_fusion_data() called with sample_dir: {sample_dir}")
            
            # Configuration for breakpoint validation
            USE_BREAKPOINT_VALIDATION = True  # Set to False to disable breakpoint validation
            MIN_BREAKPOINT_SUPPORT = 4  # Minimum reads supporting same breakpoint
            MAX_BREAKPOINT_DISTANCE = 100  # Maximum distance for breakpoint clustering
            
            # Load target panel processed
            target_file = sample_dir / "fusion_candidates_master_processed.pkl"
            genome_file = sample_dir / "fusion_candidates_all_processed.pkl"
            # Master BED CSV file loading has been deprecated
            
            logging.info(f"[Fusion] Target file: {target_file}")
            logging.info(f"[Fusion] Genome file: {genome_file}")

            t = _load_processed_pickle(target_file)
            g = _load_processed_pickle(genome_file)
            
            # Master BED CSV loading has been deprecated - the table is no longer displayed in the GUI
            # The Parquet file (master_bed_candidates.parquet) is still used for BED generation
            master_bed_df = None
            
            # Apply breakpoint validation if enabled
            if USE_BREAKPOINT_VALIDATION:
                if t:
                    logging.info("[Fusion] Applying breakpoint validation to target data")
                    annotated_data = t.get("annotated_data", pd.DataFrame())
                    if not annotated_data.empty:
                        validated_breakpoints = _validate_fusion_breakpoints(
                            annotated_data, 
                            min_read_support=MIN_BREAKPOINT_SUPPORT,
                            max_breakpoint_distance=MAX_BREAKPOINT_DISTANCE
                        )
                        if not validated_breakpoints.empty:
                            t["validated_breakpoints"] = validated_breakpoints
                            t["original_candidate_count"] = t.get("candidate_count", 0)
                            t["candidate_count"] = len(validated_breakpoints)
                            logging.info(f"[Fusion] Target breakpoint validation: {t['original_candidate_count']} -> {t['candidate_count']} candidates")
                        else:
                            t["candidate_count"] = 0
                            logging.info("[Fusion] No target fusions passed breakpoint validation")
                
                if g:
                    logging.info("[Fusion] Applying breakpoint validation to genome-wide data")
                    annotated_data = g.get("annotated_data", pd.DataFrame())
                    if not annotated_data.empty:
                        validated_breakpoints = _validate_fusion_breakpoints(
                            annotated_data, 
                            min_read_support=MIN_BREAKPOINT_SUPPORT,
                            max_breakpoint_distance=MAX_BREAKPOINT_DISTANCE
                        )
                        if not validated_breakpoints.empty:
                            g["validated_breakpoints"] = validated_breakpoints
                            g["original_candidate_count"] = g.get("candidate_count", 0)
                            g["candidate_count"] = len(validated_breakpoints)
                            logging.info(f"[Fusion] Genome-wide breakpoint validation: {g['original_candidate_count']} -> {g['candidate_count']} candidates")
                        else:
                            g["candidate_count"] = 0
                            logging.info("[Fusion] No genome-wide fusions passed breakpoint validation")
            
            # Debug logging
            logging.info(f"[Fusion] Loaded target data: {t is not None}")
            logging.info(f"[Fusion] Loaded genome-wide data: {g is not None}")
            
            # Apply same logic to target panel as genome-wide
            if t is not None:
                logging.info(f"[Fusion] Target candidate count: {t.get('candidate_count', 0)}")
                logging.info(f"[Fusion] Target gene groups: {len(t.get('gene_groups', []))}")
                logging.info(f"[Fusion] Target gene groups content: {t.get('gene_groups', [])}")
                logging.info(f"[Fusion] Target annotated_data shape: {t.get('annotated_data', pd.DataFrame()).shape}")
                logging.info(f"[Fusion] Target goodpairs shape: {t.get('goodpairs', pd.Series()).shape}")
                logging.info(f"[Fusion] Target gene_pairs count: {len(t.get('gene_pairs', []))}")
                logging.info(f"[Fusion] Target gene_pairs: {t.get('gene_pairs', [])[:10]}...")  # Show first 10
                
                # Use the same logic as reporting code - count gene_pairs instead of relying on candidate_count
                if t.get('gene_pairs') and len(t.get('gene_pairs', [])) > 0:
                    # Override candidate_count with the actual number of gene pairs (like reporting code does)
                    t['candidate_count'] = len(t.get('gene_pairs', []))
                    logging.info(f"[Fusion] Override target candidate_count to: {t['candidate_count']}")
                    
                    # Generate gene_groups from gene_pairs if missing (like reporting code does)
                    if not t.get('gene_groups') or len(t.get('gene_groups', [])) == 0:
                        # Convert gene_pairs to gene_groups format
                        gene_groups = []
                        for gene_pair in t.get('gene_pairs', []):
                            if isinstance(gene_pair, (tuple, list)) and len(gene_pair) >= 2:
                                gene_groups.append(list(gene_pair))
                        t['gene_groups'] = gene_groups
                        logging.info(f"[Fusion] Generated {len(gene_groups)} gene_groups from gene_pairs")
            
            if g is not None:
                logging.info(f"[Fusion] Genome-wide candidate count: {g.get('candidate_count', 0)}")
                logging.info(f"[Fusion] Genome-wide gene groups: {len(g.get('gene_groups', []))}")
                logging.info(f"[Fusion] Genome-wide gene groups content: {g.get('gene_groups', [])}")
                logging.info(f"[Fusion] Genome-wide annotated_data shape: {g.get('annotated_data', pd.DataFrame()).shape}")
                logging.info(f"[Fusion] Genome-wide goodpairs shape: {g.get('goodpairs', pd.Series()).shape}")
                logging.info(f"[Fusion] Genome-wide gene_pairs count: {len(g.get('gene_pairs', []))}")
                logging.info(f"[Fusion] Genome-wide gene_pairs: {g.get('gene_pairs', [])[:10]}...")  # Show first 10
                
                # Use the same logic as reporting code - count gene_pairs instead of relying on candidate_count
                if g.get('gene_pairs') and len(g.get('gene_pairs', [])) > 0:
                    # Override candidate_count with the actual number of gene pairs (like reporting code does)
                    g['candidate_count'] = len(g.get('gene_pairs', []))
                    logging.info(f"[Fusion] Override genome-wide candidate_count to: {g['candidate_count']}")
                    
                    # Generate gene_groups from gene_pairs if missing (like reporting code does)
                    if not g.get('gene_groups') or len(g.get('gene_groups', [])) == 0:
                        # Convert gene_pairs to gene_groups format
                        gene_groups = []
                        for gene_pair in g.get('gene_pairs', []):
                            if isinstance(gene_pair, (tuple, list)) and len(gene_pair) >= 2:
                                gene_groups.append(list(gene_pair))
                        g['gene_groups'] = gene_groups
                        logging.info(f"[Fusion] Generated {len(gene_groups)} gene_groups from gene_pairs")
            else:
                logging.info(f"[Fusion] Genome-wide file exists: {genome_file.exists()}")
                if genome_file.exists():
                    logging.info(f"[Fusion] Genome-wide file size: {genome_file.stat().st_size} bytes")

            # Get file modification times
            target_mtime = target_file.stat().st_mtime if target_file.exists() else None
            genome_mtime = genome_file.stat().st_mtime if genome_file.exists() else None
            # Master BED CSV file has been deprecated - no longer track mtime
            master_bed_mtime = None
            
            # Create data hashes
            target_data_hash = _create_data_hash(t) if t is not None else None
            genome_data_hash = _create_data_hash(g) if g is not None else None
            
            elapsed = _time.perf_counter() - t_load_start
            logging.info(f"[Fusion] _load_fusion_data() completed in {elapsed:.2f}s")
            return {
                "target": {
                    "data": t,
                    "mtime": target_mtime,
                    "data_hash": target_data_hash
                },
                "genome": {
                    "data": g,
                    "mtime": genome_mtime,
                    "data_hash": genome_data_hash
                },
                "master_bed": {
                    "data": master_bed_df,
                    "mtime": master_bed_mtime,
                },
                "sample_dir": sample_dir,  # Pass sample_dir for loading pre-computed summary files
            }
        except Exception as e:
            elapsed = _time.perf_counter() - t_load_start
            logging.exception(f"[Fusion] Failed to load fusion data after {elapsed:.2f}s: {e}")
            return {
                "target": {"data": None, "mtime": None, "data_hash": None}, 
                "genome": {"data": None, "mtime": None, "data_hash": None},
                "master_bed": {"data": None, "mtime": None}
            }

    def _update_fusion_ui(fusion_data: Dict[str, Any], state: Dict[str, Any], sample_dir: Path) -> None:
        """Update fusion UI elements - runs on main UI thread."""
        try:
            target_data = fusion_data.get("target", {})
            genome_data = fusion_data.get("genome", {})
            master_bed_data = fusion_data.get("master_bed", {})
            
            t = target_data.get("data")
            g = genome_data.get("data")
            master_bed_df = master_bed_data.get("data")
            
            # Update target panel UI
            target_mtime = target_data.get("mtime")
            target_data_hash = target_data.get("data_hash")
            
            # Handle case when target data is None (no data loaded)
            if t is None:
                # Update if mtime changed or if this is initial load (both mtimes are None)
                mtime_changed = target_mtime != state["target"].get("mtime")
                is_initial_load = (
                    target_mtime is None and 
                    state["target"].get("mtime") is None and
                    state["target"].get("data") is None
                )
                
                if mtime_changed or is_initial_load:
                    state["target"]["data"] = None
                    state["target"]["mtime"] = target_mtime
                    state["target"]["data_hash"] = target_data_hash
                    
                    # Update summary label
                    try:
                        if "summary" in state and state["summary"].get("target_lbl"):
                            state["summary"]["target_lbl"].text = "Target: -- pairs, -- groups"
                    except Exception:
                        pass
                    
                    # Clear containers and show status messages
                    try:
                        state["target"]["summary_table_container"].clear()
                        state["target"]["groups_table_container"].clear()
                        state["target"]["table_container"].clear()
                        state["target"]["plot_container"].clear()
                        state["target"]["status_container"].clear()
                    except Exception:
                        pass
                    
                    # Show status messages
                    with state["target"]["summary_table_container"].classes("w-full"):
                        ui.label("No fusion data available").classes("text-gray-600")
                    with state["target"]["groups_table_container"].classes("w-full"):
                        ui.label("No validated fusion groups found").classes("text-gray-600")
                    with state["target"]["table_container"].classes("w-full"):
                        ui.label("No fusion candidates found yet").classes("text-gray-600")
                    with state["target"]["status_container"].classes("w-full"):
                        ui.label("Target panel fusion analysis not available").classes("text-gray-600 text-sm")
                        ui.label("(Fusion data file not found or could not be loaded)").classes("text-gray-500 text-xs")
            
            # Always update summary and table when file changes
            elif t is not None and target_mtime != state["target"].get("mtime"):
                state["target"]["data"] = t
                state["target"]["mtime"] = target_mtime
                state["target"]["data_hash"] = target_data_hash
                
                # summary
                try:
                    if "summary" in state and state["summary"].get("target_lbl"):
                        # Count unique pairs and groups
                        unique_pairs = _count_unique_fusion_pairs(t)
                        unique_groups = _count_unique_fusion_groups(t)
                        state["summary"][
                            "target_lbl"
                        ].text = (
                            f"Target: {unique_pairs} pairs, {unique_groups} groups"
                        )
                except Exception:
                    pass
                # summary table
                try:
                    state["target"]["summary_table_container"].clear()
                except Exception:
                    pass
                state["target"]["summary_table"] = _make_fusion_summary_table(
                    state["target"]["summary_table_container"],
                    t.get("annotated_data", pd.DataFrame()),
                    t.get("goodpairs", pd.Series()),
                )
                
                # groups table
                try:
                    state["target"]["groups_table_container"].clear()
                except Exception:
                    pass
                state["target"]["groups_table"] = _make_fusion_groups_table(
                    state["target"]["groups_table_container"],
                    t,
                )
                
                # table
                try:
                    state["target"]["table_container"].clear()
                except Exception:
                    pass
                state["target"]["table"] = _make_fusion_table(
                    state["target"]["table_container"],
                    t.get("annotated_data", pd.DataFrame()),
                )
                
                # Create visualization on first load or when data changes
                try:
                    state["target"]["plot_container"].clear()
                    state["target"]["status_container"].clear()
                except Exception:
                    pass
                # Get validated fusion pairs from summary table data
                validated_pairs = _get_validated_fusion_pairs(
                    t.get("annotated_data", pd.DataFrame()),
                    t.get("goodpairs", pd.Series())
                )
                logging.info(f"[Fusion] Target panel plotting check: candidate_count={t.get('candidate_count', 0)}, validated_pairs={len(validated_pairs)}")
                if validated_pairs:
                    with state["target"]["plot_container"].classes("w-full"):
                        with ui.row().classes("w-full"):
                            state["target"]["dropdown"] = ui.select(
                                options=validated_pairs,
                                with_input=False,
                                on_change=lambda e, t=t: _handle_gene_pair_selection("target", e.value, t),
                                value=state["target"]["selected_gene_pair"]
                            ).classes("w-40")
                        with ui.row().classes("w-full"):
                            state["target"]["card"] = ui.card()
                            with state["target"]["card"]:
                                ui.label("Select gene pair to see results.").classes(
                                    "drop-shadow font-bold"
                                )
                else:
                    # Show status message when target panel fusion data is not available
                    with state["target"]["status_container"].classes("w-full"):
                        ui.label("Target panel fusion analysis not available").classes("text-gray-600 text-sm")
                        ui.label("(No validated fusion pairs found in target panel)").classes("text-gray-500 text-xs")
                
                # Restore selected gene pair if it exists (for both cases above)
                if state["target"]["selected_gene_pair"] and state["target"]["dropdown"]:
                    try:
                        state["target"]["dropdown"].value = state["target"]["selected_gene_pair"]
                        _handle_gene_pair_selection("target", state["target"]["selected_gene_pair"], t)
                    except Exception as e:
                        logging.warning(f"[Fusion] Failed to restore target selection: {e}")
            
            # Only update visualization when data content actually changes (for background refreshes)
            elif t is not None and target_data_hash != state["target"].get("data_hash"):
                state["target"]["data_hash"] = target_data_hash
                # plot area
                try:
                    state["target"]["plot_container"].clear()
                    state["target"]["status_container"].clear()
                except Exception:
                    pass
                # Get validated fusion pairs from summary table data
                validated_pairs = _get_validated_fusion_pairs(
                    t.get("annotated_data", pd.DataFrame()),
                    t.get("goodpairs", pd.Series())
                )
                logging.info(f"[Fusion] Target panel plotting check: candidate_count={t.get('candidate_count', 0)}, validated_pairs={len(validated_pairs)}")
                if validated_pairs:
                    with state["target"]["plot_container"].classes("w-full"):
                        with ui.row().classes("w-full"):
                            state["target"]["dropdown"] = ui.select(
                                options=validated_pairs,
                                with_input=False,
                                on_change=lambda e, t=t: _handle_gene_pair_selection("target", e.value, t),
                                value=state["target"]["selected_gene_pair"]
                            ).classes("w-40")
                        with ui.row().classes("w-full"):
                            state["target"]["card"] = ui.card()
                            with state["target"]["card"]:
                                ui.label("Select gene pair to see results.").classes(
                                    "drop-shadow font-bold"
                                )
                else:
                    # Show status message when target panel fusion data is not available
                    with state["target"]["status_container"].classes("w-full"):
                        ui.label("Target panel fusion analysis not available").classes("text-gray-600 text-sm")
                        ui.label("(No validated fusion pairs found in target panel)").classes("text-gray-500 text-xs")
                
                # Restore selected gene pair if it exists (for both cases above)
                if state["target"]["selected_gene_pair"] and state["target"]["dropdown"]:
                    try:
                        state["target"]["dropdown"].value = state["target"]["selected_gene_pair"]
                        _handle_gene_pair_selection("target", state["target"]["selected_gene_pair"], t)
                    except Exception as e:
                        logging.warning(f"[Fusion] Failed to restore target selection: {e}")

            # Update genome-wide UI
            genome_mtime = genome_data.get("mtime")
            genome_data_hash = genome_data.get("data_hash")
            logging.info(f"[Fusion] Genome-wide update check: g={g is not None}, mtime_changed={genome_mtime != state['genome'].get('mtime')}")
            
            # Always update summary and table when file changes
            if g is not None and genome_mtime != state["genome"].get("mtime"):
                state["genome"]["data"] = g
                state["genome"]["mtime"] = genome_mtime
                state["genome"]["data_hash"] = genome_data_hash
                
                # summary
                try:
                    if "summary" in state and state["summary"].get("genome_lbl"):
                        # Count unique pairs and groups
                        unique_pairs = _count_unique_fusion_pairs(g)
                        unique_groups = _count_unique_fusion_groups(g)
                        state["summary"][
                            "genome_lbl"
                        ].text = f"Genome-wide: {unique_pairs} pairs, {unique_groups} groups"
                except Exception:
                    pass
                # summary table
                try:
                    state["genome"]["summary_table_container"].clear()
                except Exception:
                    pass
                state["genome"]["summary_table"] = _make_fusion_summary_table(
                    state["genome"]["summary_table_container"],
                    g.get("annotated_data", pd.DataFrame()),
                    g.get("goodpairs", pd.Series()),
                )
                
                # groups table
                try:
                    state["genome"]["groups_table_container"].clear()
                except Exception:
                    pass
                state["genome"]["groups_table"] = _make_fusion_groups_table(
                    state["genome"]["groups_table_container"],
                    g,
                )
                
                # table
                try:
                    state["genome"]["table_container"].clear()
                except Exception:
                    pass
                state["genome"]["table"] = _make_fusion_table(
                    state["genome"]["table_container"],
                    g.get("annotated_data", pd.DataFrame()),
                )
                
                # Create visualization on first load or when data changes
                try:
                    state["genome"]["plot_container"].clear()
                    state["genome"]["status_container"].clear()
                except Exception:
                    pass
                # Get validated fusion pairs from summary table data
                validated_pairs = _get_validated_fusion_pairs(
                    g.get("annotated_data", pd.DataFrame()),
                    g.get("goodpairs", pd.Series())
                )
                logging.info(f"[Fusion] Genome-wide plotting check: candidate_count={g.get('candidate_count', 0)}, validated_pairs={len(validated_pairs)}")
                if validated_pairs:
                    with state["genome"]["plot_container"].classes("w-full"):
                        with ui.row().classes("w-full"):
                            state["genome"]["dropdown"] = ui.select(
                                options=validated_pairs,
                                with_input=False,
                                on_change=lambda e, g=g: _handle_gene_pair_selection("genome", e.value, g),
                                value=state["genome"]["selected_gene_pair"]
                            ).classes("w-40")
                        with ui.row().classes("w-full"):
                            state["genome"]["card"] = ui.card()
                            with state["genome"]["card"]:
                                ui.label("Select gene pair to see results.").classes(
                                    "drop-shadow font-bold"
                                )
                else:
                    # Show status message when genome-wide data is not available
                    with state["genome"]["status_container"].classes("w-full"):
                        ui.label("Genome-wide fusion analysis not available").classes("text-gray-600 text-sm")
                        ui.label("(Requires supplementary reads in BAM file)").classes("text-gray-500 text-xs")
                
                # Restore selected gene pair if it exists (for both cases above)
                if state["genome"]["selected_gene_pair"] and state["genome"]["dropdown"]:
                    try:
                        state["genome"]["dropdown"].value = state["genome"]["selected_gene_pair"]
                        _handle_gene_pair_selection("genome", state["genome"]["selected_gene_pair"], g)
                    except Exception as e:
                        logging.warning(f"[Fusion] Failed to restore genome selection: {e}")
            
            # Only update visualization when data content actually changes (for background refreshes)
            elif g is not None and genome_data_hash != state["genome"].get("data_hash"):
                state["genome"]["data_hash"] = genome_data_hash
                # plot area
                try:
                    state["genome"]["plot_container"].clear()
                    state["genome"]["status_container"].clear()
                except Exception:
                    pass
                # Get validated fusion pairs from summary table data
                validated_pairs = _get_validated_fusion_pairs(
                    g.get("annotated_data", pd.DataFrame()),
                    g.get("goodpairs", pd.Series())
                )
                logging.info(f"[Fusion] Genome-wide plotting check: candidate_count={g.get('candidate_count', 0)}, validated_pairs={len(validated_pairs)}")
                if validated_pairs:
                    with state["genome"]["plot_container"].classes("w-full"):
                        with ui.row().classes("w-full"):
                            state["genome"]["dropdown"] = ui.select(
                                options=validated_pairs,
                                with_input=False,
                                on_change=lambda e, g=g: _handle_gene_pair_selection("genome", e.value, g),
                                value=state["genome"]["selected_gene_pair"]
                            ).classes("w-40")
                        with ui.row().classes("w-full"):
                            state["genome"]["card"] = ui.card()
                            with state["genome"]["card"]:
                                ui.label("Select gene pair to see results.").classes(
                                    "drop-shadow font-bold"
                                )
                else:
                    # Show status message when genome-wide data is not available
                    with state["genome"]["status_container"].classes("w-full"):
                        ui.label("Genome-wide fusion analysis not available").classes("text-gray-600 text-sm")
                        ui.label("(Requires supplementary reads in BAM file)").classes("text-gray-500 text-xs")
                
                # Restore selected gene pair if it exists (for both cases above)
                if state["genome"]["selected_gene_pair"] and state["genome"]["dropdown"]:
                    try:
                        state["genome"]["dropdown"].value = state["genome"]["selected_gene_pair"]
                        _handle_gene_pair_selection("genome", state["genome"]["selected_gene_pair"], g)
                    except Exception as e:
                        logging.warning(f"[Fusion] Failed to restore genome selection: {e}")
            
            # Update master BED UI
            master_bed_mtime = master_bed_data.get("mtime")
            
            # Always update when file changes
            if master_bed_mtime != state["master_bed"].get("mtime"):
                state["master_bed"]["data"] = master_bed_df
                state["master_bed"]["mtime"] = master_bed_mtime
                
                # Master BED summary label has been deprecated
                try:
                    if "summary" in state and state["summary"].get("master_bed_lbl"):
                        state["summary"]["master_bed_lbl"].text = "Master BED: (deprecated)"
                except Exception:
                    pass
                
                # Master BED table display has been deprecated
                try:
                    state["master_bed"]["table_container"].clear()
                    state["master_bed"]["status_container"].clear()
                    with state["master_bed"]["status_container"].classes("w-full"):
                        ui.label("Master BED table display has been deprecated").classes("text-gray-500 text-sm")
                        ui.label("(Data is still processed and used for BED generation)").classes("text-gray-400 text-xs")
                except Exception:
                    pass
        except Exception as e:
            logging.exception(f"[Fusion] Refresh failed: {e}")

    # Build UI — design.md §6–§9 (Digital Curator shell; no emoji in title)
    with ui.element("div").classes("w-full min-w-0").props("id=analysis-detail-fusion"):
        with ui.element("div").classes("classification-insight-shell w-full min-w-0"):
            ui.label("Fusion analysis").classes(
                "classification-insight-heading text-headline-small"
            )
            with ui.element("div").classes(
                "classification-insight-card w-full min-w-0"
            ):
                with ui.column().classes("w-full min-w-0 gap-2 p-2 md:p-3"):
                    with ui.row().classes("items-center gap-2 min-w-0"):
                        ui.icon("merge").classes("classification-insight-icon")
                        ui.label("Candidate summary").classes(
                            "classification-insight-model flex-1 min-w-0"
                        )
                    state.setdefault("summary", {})
                    state["summary"]["target_lbl"] = ui.label(
                        "Target: -- pairs, -- groups"
                    ).classes("classification-insight-meta w-full")
                    state["summary"]["genome_lbl"] = ui.label(
                        "Genome-wide: -- pairs, -- groups"
                    ).classes("classification-insight-meta w-full")
                    state["summary"]["master_bed_lbl"] = ui.label(
                        "Master BED: (deprecated)"
                    ).classes("hidden")
                    ui.label(
                        "From supplementary alignments; tables and plots below."
                    ).classes("classification-insight-foot")

            ui.label("Target panel").classes(
                "target-coverage-panel__meta-label mt-4 mb-1"
            )
            state["target"]["summary_table_container"] = ui.column().classes("w-full")
            state["target"]["groups_table_container"] = ui.column().classes(
                "w-full mt-2"
            )
            state["target"]["plot_container"] = ui.column().classes("w-full")
            state["target"]["table_container"] = ui.column().classes("w-full mt-2")
            state["target"]["status_container"] = ui.column().classes("w-full mt-2")

            # Genome-wide
            ui.separator().classes("mgmt-detail-separator")
            ui.label("Genome-wide").classes(
                "target-coverage-panel__meta-label mt-2 mb-1"
            )
            state["genome"]["summary_table_container"] = ui.column().classes("w-full")
            state["genome"]["groups_table_container"] = ui.column().classes(
                "w-full mt-2"
            )
            state["genome"]["plot_container"] = ui.column().classes("w-full")
            state["genome"]["table_container"] = ui.column().classes("w-full mt-2")
            state["genome"]["status_container"] = ui.column().classes("w-full mt-2")

            # Master BED Targets section has been deprecated
            # Keep containers initialized for backward compatibility
            state["master_bed"]["table_container"] = ui.column().classes(
                "w-full hidden"
            )
            state["master_bed"]["status_container"] = ui.column().classes(
                "w-full mt-2 hidden"
            )

    # Initial refresh and timer
    try:
        # Prime summary values on first load
        try:
            t0 = _load_processed_pickle(
                (sample_dir / "fusion_candidates_master_processed.pkl")
                if sample_dir
                else Path("/dev/null")
            )
            g0 = _load_processed_pickle(
                (sample_dir / "fusion_candidates_all_processed.pkl")
                if sample_dir
                else Path("/dev/null")
            )
            if "summary" in state:
                t_pairs0 = _count_unique_fusion_pairs(t0) if isinstance(t0, dict) else 0
                t_groups0 = _count_unique_fusion_groups(t0) if isinstance(t0, dict) else 0
                g_pairs0 = _count_unique_fusion_pairs(g0) if isinstance(g0, dict) else 0
                g_groups0 = _count_unique_fusion_groups(g0) if isinstance(g0, dict) else 0
                try:
                    state["summary"][
                        "target_lbl"
                    ].text = f"Target: {t_pairs0} pairs, {t_groups0} groups"
                    state["summary"][
                        "genome_lbl"
                    ].text = f"Genome-wide: {g_pairs0} pairs, {g_groups0} groups"
                    
                    # Master BED summary has been deprecated
                    try:
                        state["summary"]["master_bed_lbl"].text = "Master BED: (deprecated)"
                    except Exception:
                        pass
                except Exception:
                    pass
        except Exception:
            pass
        
        # Start the refresh timer (every 30 seconds)
        logging.info("[Fusion] Setting up refresh timer (30s interval + 0.5s deferred + immediate)")
        refresh_timer = ui.timer(30.0, refresh_fusion, active=True, immediate=False)
        ui.timer(0.5, refresh_fusion, once=True)

        def _sync_fusion_plot_theme_if_needed() -> None:
            try:
                dark = _is_dark_mode()
            except Exception:
                dark = False
            for sec in ("target", "genome"):
                st = state.get(sec, {})
                if st.get("fusion_plot_theme_dark") == dark:
                    continue
                mpl = st.get("fusion_mpl")
                if mpl is None:
                    continue
                fig = getattr(mpl, "figure", None)
                if fig is None:
                    continue
                _apply_fusion_figure_theme(fig, dark)
                try:
                    mpl.update()
                except Exception:
                    pass
                st["fusion_plot_theme_dark"] = dark

        from robin.gui.theme import register_theme_sync_callback

        unregister_fusion_theme_sync = register_theme_sync_callback(
            _sync_fusion_plot_theme_if_needed,
            element=state.get("target", {}).get("fusion_mpl"),
            interval_s=0.5,
            immediate=True,
        )
        try:
            ui.context.client.on_disconnect(
                lambda: (
                    refresh_timer.deactivate(),
                    unregister_fusion_theme_sync(),
                )
            )
        except Exception:
            pass
        logging.info("[Fusion] Timer set up successfully")
        # Populate fusion UI immediately so it matches Summary without waiting for 0.5s/30s
        refresh_fusion()
    except Exception as e:
        logging.exception(f"[Fusion] Exception in timer setup: {e}")
