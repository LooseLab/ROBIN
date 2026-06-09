"""
MNP-Flex reporting section.
"""

from __future__ import annotations

import json
import logging
import os
from typing import Any, Dict, List, Optional, Tuple

from reportlab.platypus import Paragraph, Spacer
from reportlab.lib.units import inch

from .base import ReportSection

logger = logging.getLogger(__name__)


class MNPFlexSection(ReportSection):
    """Report section for MNP-Flex results."""

    def _find_results_dir(self) -> Optional[str]:
        sample_id = self.report.sample_id
        base = self.report.output
        candidates = [
            os.path.join(base, f"mnpflex_results_{sample_id}"),
            os.path.join(base, "mnpflex_results"),
            base,
        ]
        for candidate in candidates:
            if os.path.exists(os.path.join(candidate, "bundle_summary.json")):
                return candidate
        return None

    def _load_bundle_summary(self, results_dir: str) -> Optional[Dict[str, Any]]:
        summary_path = os.path.join(results_dir, "bundle_summary.json")
        try:
            with open(summary_path, "r") as f:
                return json.load(f)
        except Exception as exc:
            logger.error(f"Failed to load MNP-Flex summary: {exc}")
            return None

    def _format_score(self, value: Optional[float]) -> str:
        try:
            return f"{float(value):.4f}"
        except Exception:
            return "--"

    def _format_mgmt_value(self, value: Any, decimals: int = 1) -> str:
        """Format MGMT/percentage values with sensible rounding."""
        if value is None:
            return "N/A"
        try:
            return f"{float(value):.{decimals}f}"
        except (TypeError, ValueError):
            return str(value)

    def _format_coverage_value(self, value: Any) -> str:
        """Format coverage values with 2 decimal places."""
        if value is None:
            return "N/A"
        try:
            return f"{float(value):.2f}"
        except (TypeError, ValueError):
            return str(value)

    def _collect_descriptions(self, node: Dict[str, Any]) -> str:
        descriptions: List[str] = []
        desc = (node.get("description") or "").strip()
        if desc:
            descriptions.append(desc)
        for member in node.get("members") or []:
            member_desc = self._collect_descriptions(member)
            if member_desc:
                descriptions.append(member_desc)
        seen = set()
        combined = []
        for text in descriptions:
            if text not in seen:
                seen.add(text)
                combined.append(text)
        return " ".join(combined).strip()

    def _flatten_hierarchy(
        self, nodes: List[Dict[str, Any]], path: Optional[List[str]] = None
    ) -> List[Tuple[Optional[float], List[str]]]:
        if path is None:
            path = []
        flat = []
        for node in nodes or []:
            group = node.get("group", "Unknown")
            score = node.get("score")
            current = path + [group]
            members = node.get("members") or []
            if members:
                flat.extend(self._flatten_hierarchy(members, current))
            else:
                flat.append((score, current))
        return flat

    def _sum_scores_by_field(
        self, scores: List[Dict[str, Any]], field_name: str, match_value: Optional[str]
    ) -> Optional[float]:
        if not match_value:
            return None
        total = 0.0
        for item in scores:
            ref = item.get("reference_group") or {}
            if ref.get(field_name) == match_value:
                total += float(item.get("score", 0) or 0)
        return total

    def add_content(self):
        """Add MNP-Flex content to the report if available."""
        results_dir = self._find_results_dir()
        if not results_dir:
            logger.info("No MNP-Flex results found for report.")
            return

        summary = self._load_bundle_summary(results_dir)
        if not summary:
            return

        qc = summary.get("qc", {}) or {}
        mgmt = summary.get("mgmt", {}) or {}
        classifier_summary = summary.get("classifier_summary", {}) or {}
        classifier = classifier_summary.get("classifier", {}) or {}
        hierarchy = classifier_summary.get("summary_hierarchical", []) or []
        has_hierarchical_summary = len(hierarchy) > 0

        # Summary card removed - MNP-Flex legend now appears as table legend
        # in ClassificationSection (with MNP-Flex Hierarchical Summary)

        # Main section header
        self.add_section_header("MNP-Flex")

        classifier_name = classifier.get("name", "Unknown")
        classifier_version = classifier.get("version", "Unknown")
        classifier_type = classifier.get("classifier_type", "Unknown")
        self.elements.append(
            Paragraph(
                f"<b>Classifier</b>: {classifier_name} (v{classifier_version}, {classifier_type})",
                self.styles.styles["Normal"],
            )
        )
        self.elements.append(Spacer(1, 4))

        # QC / MGMT summary
        summary_text = (
            f"<b>QC status</b>: {qc.get('status', 'Unknown')}<br/>"
            f"<b>Average coverage</b>: {self._format_coverage_value(qc.get('avg_coverage'))}<br/>"
            f"<b>Missing sites</b>: {qc.get('missing_site_count', 'N/A')}<br/>"
            f"<b>MGMT status</b>: {mgmt.get('status', 'Unknown')}<br/>"
            f"<b>MGMT average</b>: {self._format_mgmt_value(mgmt.get('average'))}"
        )
        self.elements.append(Paragraph(summary_text, self.styles.styles["Normal"]))
        self.elements.append(Spacer(1, 4))

        # Hierarchical summary
        flat = self._flatten_hierarchy(hierarchy)

        if not has_hierarchical_summary:
            self.elements.append(
                Paragraph(
                    "This is not a confirmed classification. No hierarchical summary "
                    "is available from the classifier.",
                    self.styles.styles["Warning"],
                )
            )
            self.elements.append(Spacer(1, 6))

        if flat:
            best_score, best_path = max(flat, key=lambda x: x[0] or 0)
            self.elements.append(
                Paragraph(
                    f"<b>Top path</b>: {' > '.join(best_path)} "
                    f"({self._format_score(best_score)})",
                    self.styles.styles["Normal"],
                )
            )
            self.elements.append(Spacer(1, 4))

        if hierarchy:
            hierarchy_rows = [
                ["Group", "Score", "Description"],
            ]
            for node in hierarchy:
                hierarchy_rows.append(
                    [
                        node.get("group", "Unknown"),
                        self._format_score(node.get("score")),
                        self._collect_descriptions(node),
                    ]
                )
            self.elements.append(self.create_table(hierarchy_rows))
            self.elements.append(Spacer(1, 6))

        # Export frames for CSV/XLSX/ZIP artifacts
        try:
            import pandas as pd

            top_path = ""
            top_score = None
            if flat:
                top_score, best_path = max(flat, key=lambda x: x[0] or 0)
                top_path = " > ".join(best_path)

            self.export_frames["mnpflex_summary"] = pd.DataFrame(
                [
                    {
                        "results_dir": results_dir,
                        "qc_status": qc.get("status", "Unknown"),
                        "qc_avg_coverage": qc.get("avg_coverage"),
                        "qc_missing_site_count": qc.get("missing_site_count"),
                        "mgmt_status": mgmt.get("status", "Unknown"),
                        "mgmt_average": mgmt.get("average"),
                        "mgmt_site_count": mgmt.get("site_count"),
                        "classifier_name": classifier.get("name", "Unknown"),
                        "classifier_version": classifier.get("version", "Unknown"),
                        "classifier_type": classifier.get("classifier_type", "Unknown"),
                        "has_hierarchical_summary": bool(has_hierarchical_summary),
                        "top_path": top_path,
                        "top_path_score": top_score,
                    }
                ]
            )

            if hierarchy:
                self.export_frames["mnpflex_hierarchy"] = pd.DataFrame(
                    [
                        {
                            "group": node.get("group", "Unknown"),
                            "score": node.get("score"),
                            "description": self._collect_descriptions(node),
                        }
                        for node in hierarchy
                    ]
                )

            scores = classifier_summary.get("scores") or []
            if scores:
                top_scores = sorted(
                    scores,
                    key=lambda item: float(item.get("score", 0) or 0),
                    reverse=True,
                )[:10]
                self.export_frames["mnpflex_top_scores"] = pd.DataFrame(
                    [
                        {
                            "score": item.get("score"),
                            "subclass": (item.get("reference_group") or {}).get("molecular_subclass")
                            or (item.get("reference_group") or {}).get("name"),
                            "class": (item.get("reference_group") or {}).get("molecular_class"),
                            "family": (item.get("reference_group") or {}).get("molecular_family"),
                            "superfamily": (item.get("reference_group") or {}).get("molecular_superfamily"),
                        }
                        for item in top_scores
                    ]
                )
        except Exception as ex:
            logger.error("Error building MNP-Flex export DataFrames: %s", ex, exc_info=True)

        # Top 10 classifier scores
        scores = classifier_summary.get("scores") or []
        if scores:
            top = sorted(
                scores,
                key=lambda item: float(item.get("score", 0) or 0),
                reverse=True,
            )[:10]

            top_rows = [
                ["Score", "Subclass", "Class", "Family", "Superfamily"],
            ]
            for item in top:
                ref = item.get("reference_group") or {}
                top_rows.append(
                    [
                        self._format_score(item.get("score")),
                        ref.get("molecular_subclass") or ref.get("name") or "",
                        ref.get("molecular_class") or "",
                        ref.get("molecular_family") or "",
                        ref.get("molecular_superfamily") or "",
                    ]
                )
            self.elements.append(Paragraph("Top 10 classifier scores", self.styles.styles["Heading3"]))
            self.elements.append(self.create_table(top_rows))
            self.elements.append(Spacer(1, 4))

            # Aggregate scores for top entry
            top_ref = (top[0].get("reference_group") or {})
            top_subclass = top_ref.get("molecular_subclass") or top_ref.get("name")
            top_class = top_ref.get("molecular_class")
            top_family = top_ref.get("molecular_family")
            top_superfamily = top_ref.get("molecular_superfamily")

            subclass_sum = self._sum_scores_by_field(scores, "molecular_subclass", top_subclass)
            class_sum = self._sum_scores_by_field(scores, "molecular_class", top_class)
            family_sum = self._sum_scores_by_field(scores, "molecular_family", top_family)
            superfamily_sum = self._sum_scores_by_field(scores, "molecular_superfamily", top_superfamily)

            agg_rows = [
                ["Level", "Name", "Aggregate Score"],
                ["Subclass", top_subclass or "N/A", self._format_score(subclass_sum)],
                ["Class", top_class or "N/A", self._format_score(class_sum)],
                ["Family", top_family or "N/A", self._format_score(family_sum)],
                ["Superfamily", top_superfamily or "N/A", self._format_score(superfamily_sum)],
            ]
            self.elements.append(Paragraph("Aggregate scores for top entry", self.styles.styles["Heading3"]))
            if not has_hierarchical_summary:
                self.elements.append(
                    Paragraph(
                        "For reference only — not a confirmed classification.",
                        self.styles.styles["Smaller"],
                    )
                )
                self.elements.append(Spacer(1, 4))
            self.elements.append(self.create_table(agg_rows))
            self.elements.append(Spacer(1, 4))

        # Plots (if available)
        plot_specs = [
            ("QC coverage plot", os.path.join(results_dir, "qc_coverage_plot.png")),
            ("QC methylation density plot", os.path.join(results_dir, "qc_methylation_density_plot.png")),
            ("MGMT region plot", os.path.join(results_dir, "mgmt_region_plot.png")),
        ]
        for title, path in plot_specs:
            if os.path.exists(path):
                self.add_figure(path, caption=title, width=6 * inch, height=3 * inch)
