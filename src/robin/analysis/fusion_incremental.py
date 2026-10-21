"""Incremental reporting over immutable fusion candidate Parquet parts.

Detection stays in fusion_work. State is disposable: changed/deleted parts,
threshold changes or corrupt state cause a rebuild from the evidence dataset.
Only reported rows are cached; read-to-part indexes recover historical evidence
when a pair first qualifies or a previously observed read acquires another gene.
"""

from __future__ import annotations

import os
import pickle
from itertools import combinations
from pathlib import Path

import pandas as pd

from robin.analysis.master_bed_generator import FileLock
from robin.classification_config import get_fusion_threshold

OUTPUTS = (
    (
        "target_candidates",
        "fusion_candidates_master.csv",
        "fusion_candidates_master_processed.pkl",
    ),
    (
        "genome_wide_candidates",
        "fusion_candidates_all.csv",
        "fusion_candidates_all_processed.pkl",
    ),
)
STATE_VERSION = 2


def _state_path(directory: Path, candidate_type: str) -> Path:
    return directory / f"{candidate_type}_incremental.pkl"


def update_gene_pair_calls(candidate_type, csv_name, pickle_name, work_dir, sample_id):
    """Update pair support from new parts; atomically commit replayable output state.

    The state contains the reported table as well as support indexes. Committing
    it before publishing outputs makes retries safe after partial publication.
    This lock serializes gene reporting only, not structural-variant processing.
    """
    from robin.analysis import fusion_work as fw

    directory = Path(work_dir) / sample_id
    directory.mkdir(parents=True, exist_ok=True)
    with FileLock(str(directory / "_locks" / "fusion_gene_calls.lock"), timeout=300):
        fw._migrate_legacy_parquet_to_dataset(work_dir, sample_id)
        parts = sorted(
            Path(fw._get_parquet_dataset_dir(work_dir, sample_id, candidate_type)).glob(
                "*.parquet"
            )
        )
        signatures = {str(p): (p.stat().st_size, p.stat().st_mtime_ns) for p in parts}
        support = get_fusion_threshold("read_support")
        state_path = _state_path(directory, candidate_type)
        try:
            with state_path.open("rb") as handle:
                state = pickle.load(handle)
            if state["version"] != STATE_VERSION or state["support"] != support:
                raise ValueError("Incompatible gene calling state")
            if any(signatures.get(p) != sig for p, sig in state["parts"].items()):
                raise ValueError("Previously processed evidence changed")
            for key in ("genes", "read_parts", "pairs", "rows"):
                if key not in state:
                    raise ValueError("Incomplete gene calling state")
        except (
            OSError,
            EOFError,
            pickle.UnpicklingError,
            AssertionError,
            KeyError,
            TypeError,
            ValueError,
            AttributeError,
        ):
            state = dict(
                version=STATE_VERSION,
                support=support,
                parts={},
                genes={},
                read_parts={},
                pairs={},
                rows=pd.DataFrame(),
            )

        new_parts = [str(p) for p in parts if str(p) not in state["parts"]]
        touched = set()
        previously_valid = {
            pair for pair, reads in state["pairs"].items() if len(reads) >= support
        }
        new_support = {}
        for part in new_parts:
            for batch in fw.pq.ParquetFile(part).iter_batches(
                columns=["read_id", "col4"], batch_size=200000
            ):
                frame = batch.to_pandas()
                for rid, group in frame.groupby("read_id", observed=True):
                    if pd.isna(rid):
                        continue
                    rid = str(rid)
                    touched.add(rid)
                    state["read_parts"].setdefault(rid, set()).add(part)
                    genes = state["genes"].setdefault(rid, set())
                    old_pairs = set(combinations(sorted(genes), 2))
                    genes.update(str(g) for g in group["col4"] if not pd.isna(g))
                    for pair in set(combinations(sorted(genes), 2)) - old_pairs:
                        state["pairs"].setdefault(pair, set()).add(rid)
                        new_support.setdefault(pair, set()).add(rid)

        affected = {
            pair
            for rid in touched
            for pair in combinations(sorted(state["genes"][rid]), 2)
            if len(state["pairs"][pair]) >= support
        }
        rows = state["rows"]
        if affected:
            # Established pairs only need new parts and historical evidence for
            # newly supporting reads. Newly qualifying pairs need all supporters.
            part_pairs = {}
            for pair in affected:
                relevant = set(new_parts)
                reads = (
                    new_support.get(pair, set())
                    if pair in previously_valid
                    else state["pairs"][pair]
                )
                for rid in reads:
                    relevant.update(state["read_parts"][rid])
                for part in relevant:
                    part_pairs.setdefault(part, set()).add(pair)
            updates = []
            for part, pairs in sorted(part_pairs.items()):
                if not rows.empty:
                    tags = {",".join(pair) for pair in pairs}
                    rows = rows[
                        ~((rows["__fusion_part"] == part) & rows["tag"].isin(tags))
                    ]
                for batch in fw.pq.ParquetFile(part).iter_batches(batch_size=200000):
                    frame = batch.to_pandas()
                    frame["read_id"] = frame["read_id"].astype(str)
                    tagged = []
                    for pair in sorted(pairs):
                        subset = frame[
                            frame["read_id"].isin(state["pairs"][pair])
                            & frame["col4"].astype(str).isin(pair)
                        ].copy()
                        if not subset.empty:
                            subset["tag"] = ",".join(pair)
                            subset["__fusion_part"] = part
                            tagged.append(subset)
                    if tagged:
                        # Match full caller's per-Parquet-batch deduplication.
                        updates.append(
                            pd.concat(tagged, ignore_index=True).drop_duplicates()
                        )
            rows = (
                pd.concat([rows, *updates], ignore_index=True)
                if not rows.empty
                else pd.concat(updates, ignore_index=True)
            )

        publish = bool(new_parts) or state.get("publish_pending", True)
        state["publish_pending"] = publish
        state["rows"] = rows
        state["parts"] = signatures
        temporary = state_path.with_suffix(".tmp")
        with temporary.open("wb") as handle:
            pickle.dump(state, handle, protocol=pickle.HIGHEST_PROTOCOL)
        os.replace(temporary, state_path)
        csv_path, processed_path = directory / csv_name, directory / pickle_name
        if rows.empty:
            csv_path.unlink(missing_ok=True)
            processed_path.unlink(missing_ok=True)
        elif publish or not csv_path.exists() or not processed_path.exists():
            temporary_csv = csv_path.with_suffix(".tmp")
            rows.drop(columns="__fusion_part").to_csv(temporary_csv, index=False)
            os.replace(temporary_csv, csv_path)
            fw.preprocess_fusion_data_standalone(
                pd.read_csv(csv_path), str(processed_path)
            )

        state["publish_pending"] = False
        with temporary.open("wb") as handle:
            pickle.dump(state, handle, protocol=pickle.HIGHEST_PROTOCOL)
        os.replace(temporary, state_path)


def reconcile_gene_pair_calls(work_dir, sample_id):
    """Authoritatively rebuild both gene reports, even with no pending BAMs.

    Discard incremental state so a subsequent live batch bootstraps from the
    reconciled evidence. No structural-variant or master-BED work is performed.
    """
    from robin.analysis import fusion_work as fw

    directory = Path(work_dir) / sample_id
    with FileLock(str(directory / "_locks" / "fusion_gene_calls.lock"), timeout=300):
        # A single-file job can leave gene evidence below its staging threshold.
        # Flush only gene parts here; leave SV staging and pending counters alone.
        fw._migrate_legacy_parquet_to_dataset(work_dir, sample_id)
        counts = fw._load_fusion_counts(work_dir, sample_id)
        for candidate_type, prefix in (
            ("target_candidates", "target"),
            ("genome_wide_candidates", "genome"),
        ):
            for part in sorted(
                (directory / "_fusion_staging").glob(f"{prefix}_*.parquet")
            ):
                batch_id = int(part.stem.split("_")[-1])
                fw._append_fusion_candidates_parquet_from_file(
                    str(part), candidate_type, work_dir, sample_id, batch_id
                )
            counts[candidate_type] = fw._count_parquet_rows(
                [
                    str(p)
                    for p in Path(
                        fw._get_parquet_dataset_dir(work_dir, sample_id, candidate_type)
                    ).glob("*.parquet")
                ]
            )
        fw._save_fusion_counts(work_dir, sample_id, counts)
        for candidate_type, csv_name, pickle_name in OUTPUTS:
            fw._rebuild_filtered_fusion_candidates(
                candidate_type, csv_name, pickle_name, work_dir, sample_id
            )
            _state_path(directory, candidate_type).unlink(missing_ok=True)
