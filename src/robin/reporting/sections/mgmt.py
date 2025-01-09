"""
mgmt.py

This module contains the MGMT methylation analysis section of the report.
"""

import os
import pandas as pd
import natsort
from reportlab.lib.units import inch
from reportlab.platypus import PageBreak, Paragraph, Image, Spacer, Table, TableStyle
from reportlab.lib.colors import white

from ..sections.base import ReportSection


class MGMTSection(ReportSection):
    """Section containing the MGMT methylation analysis."""

    def add_content(self):
        """Add the MGMT analysis content to the report."""
        # Find the latest MGMT results
        last_seen = 0
        mgmt_results = None
        plot_out = None

        for file in natsort.natsorted(os.listdir(self.report.output)):
            if file.endswith("_mgmt.csv"):
                count = int(file.split("_")[0])
                if count > last_seen:
                    mgmt_results = pd.read_csv(os.path.join(self.report.output, file))
                    plot_out = os.path.join(
                        self.report.output, file.replace(".csv", ".png")
                    )
                    last_seen = count

        if last_seen == 0 or mgmt_results is None:
            return

        try:
            # Add summary to elements_summary
            self.summary_elements.append(
                Paragraph("MGMT Promoter Methylation", self.styles.styles["Heading3"])
            )

            # Extract key metrics for summary
            methylation_status = (
                mgmt_results["status"].iloc[0]
                if "status" in mgmt_results.columns
                else "Unknown"
            )
            methylation_average = (
                mgmt_results["average"].iloc[0]
                if "average" in mgmt_results.columns
                else None
            )

            summary_text = f"Status: {methylation_status}"
            if methylation_average is not None:
                summary_text += (
                    f" (Average methylation: {round(methylation_average, 3)}%)"
                )

            self.summary_elements.append(
                Paragraph(summary_text, self.styles.styles["Normal"])
            )

            # Add detailed MGMT section to main report
            self.elements.append(PageBreak())
            self.add_section_header("MGMT Promoter Methylation Analysis", level=2)

            # Add the plot if it exists
            if os.path.exists(plot_out):
                self.elements.append(Image(plot_out, width=6 * inch, height=4 * inch))
                self.elements.append(Spacer(1, 12))

            # Add detailed results table
            data = [["Metric", "Value"]]
            data.append(["Status", methylation_status])
            data.append(["Average Methylation (%)", f"{round(methylation_average, 3)}"])

            if "pred" in mgmt_results.columns:
                data.append(
                    ["Prediction Score", f"{round(mgmt_results['pred'].iloc[0], 3)}"]
                )

            # Create table style with compact formatting
            style = TableStyle(
                [
                    ("BACKGROUND", (0, 0), (-1, 0), self.styles.COLORS["background"]),
                    ("TEXTCOLOR", (0, 0), (-1, 0), self.styles.COLORS["primary"]),
                    ("ALIGN", (0, 0), (-1, -1), "LEFT"),
                    ("FONTNAME", (0, 0), (-1, 0), "Helvetica-Bold"),
                    ("FONTSIZE", (0, 0), (-1, 0), 7),  # Smaller header font
                    ("BOTTOMPADDING", (0, 0), (-1, 0), 6),  # Reduced padding
                    ("TOPPADDING", (0, 0), (-1, 0), 6),  # Reduced padding
                    ("BACKGROUND", (0, 1), (-1, -1), white),
                    ("TEXTCOLOR", (0, 1), (-1, -1), self.styles.COLORS["text"]),
                    ("FONTNAME", (0, 1), (-1, -1), "Helvetica"),
                    ("FONTSIZE", (0, 1), (-1, -1), 7),  # Smaller body font
                    (
                        "GRID",
                        (0, 0),
                        (-1, -1),
                        0.5,
                        self.styles.COLORS["border"],
                    ),  # Thinner grid lines
                    # Align values to the right
                    ("ALIGN", (1, 0), (1, -1), "RIGHT"),
                    # Enable text wrapping for all cells
                    ("LEFTPADDING", (0, 0), (-1, -1), 4),  # Reduced padding
                    ("RIGHTPADDING", (0, 0), (-1, -1), 4),  # Reduced padding
                    ("TOPPADDING", (0, 1), (-1, -1), 4),  # Reduced padding
                    ("BOTTOMPADDING", (0, 1), (-1, -1), 4),  # Reduced padding
                ]
            )

            # Color code the methylation status
            if methylation_status.lower() == "methylated":
                style.add(
                    "TEXTCOLOR", (1, 1), (1, 1), self.styles.COLORS["success"]
                )  # Green for methylated
            else:
                style.add(
                    "TEXTCOLOR", (1, 1), (1, 1), self.styles.COLORS["error"]
                )  # Red for unmethylated

            # Create and add the table with adjusted column widths
            table = Table(
                data,
                colWidths=[inch * 2.5, inch * 2.0],  # Adjust column widths
                style=style,
            )
            self.elements.append(table)
            self.elements.append(Spacer(1, 12))

        except Exception as e:
            self.log_error("Error processing MGMT results", e)
