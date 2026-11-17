"""Shared helpers for MNP-Flex hierarchical classifier summaries."""

from __future__ import annotations

from typing import Any, Dict, Iterator, List, Optional, Tuple


def format_mnpflex_hierarchy_score(value: Any) -> str:
    """Format hierarchy node scores (Epignostix uses two decimal places)."""
    try:
        return f"{float(value):.2f}"
    except (TypeError, ValueError):
        return "--"


def mnpflex_hierarchy_has_content(
    nodes: Optional[List[Dict[str, Any]]],
) -> bool:
    """Return True when summary_hierarchical contains renderable nodes."""
    for node in nodes or []:
        if (node.get("group") or "").strip():
            return True
        if mnpflex_hierarchy_has_content(node.get("members")):
            return True
    return False


def iter_mnpflex_hierarchy_nodes(
    nodes: Optional[List[Dict[str, Any]]],
    *,
    depth: int = 0,
    path: Optional[List[str]] = None,
) -> Iterator[Tuple[int, List[str], Dict[str, Any]]]:
    """Yield depth, breadcrumb path, and node for each entry in the tree."""
    current_path = path or []
    for node in nodes or []:
        group = (node.get("group") or "Unknown").strip() or "Unknown"
        next_path = current_path + [group]
        yield depth, next_path, node
        yield from iter_mnpflex_hierarchy_nodes(
            node.get("members"),
            depth=depth + 1,
            path=next_path,
        )


def flatten_mnpflex_hierarchy(
    nodes: Optional[List[Dict[str, Any]]],
    path: Optional[List[str]] = None,
) -> List[Tuple[Optional[float], List[str]]]:
    """Flatten hierarchy to leaf (score, path) pairs."""
    if path is None:
        path = []
    flat: List[Tuple[Optional[float], List[str]]] = []
    for node in nodes or []:
        group = node.get("group", "Unknown")
        score = node.get("score")
        current = path + [group]
        members = node.get("members") or []
        if members:
            flat.extend(flatten_mnpflex_hierarchy(members, current))
        else:
            flat.append((score, current))
    return flat


def best_mnpflex_hierarchy_path(
    nodes: Optional[List[Dict[str, Any]]],
) -> Tuple[Optional[float], List[str]]:
    """Return the highest-scoring leaf path in the hierarchy."""
    flat = flatten_mnpflex_hierarchy(nodes)
    if not flat:
        return None, []
    return max(flat, key=lambda item: item[0] or 0)


def mnpflex_hierarchy_export_rows(
    nodes: Optional[List[Dict[str, Any]]],
) -> List[Dict[str, Any]]:
    """Flatten hierarchy nodes for CSV/XLSX export."""
    rows: List[Dict[str, Any]] = []
    for depth, path, node in iter_mnpflex_hierarchy_nodes(nodes):
        rows.append(
            {
                "depth": depth,
                "path": " > ".join(path),
                "group": node.get("group", "Unknown"),
                "score": node.get("score"),
                "description": (node.get("description") or "").strip(),
            }
        )
    return rows
