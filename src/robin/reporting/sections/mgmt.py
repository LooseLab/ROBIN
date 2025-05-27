"""
MGMT Analysis Section for ROBIN Reports.

This module handles the MGMT (O6-methylguanine-DNA methyltransferase) promoter methylation analysis section of the report.
"""

import os
import pandas as pd
import natsort
from reportlab.lib.units import inch
from reportlab.platypus import PageBreak, Paragraph, Image, Spacer, Table, TableStyle
from reportlab.lib.colors import HexColor
from reportlab.lib.styles import ParagraphStyle
from ..sections.base import ReportSection

import logging

logger = logging.getLogger(__name__)


class MGMTSection(ReportSection):
    """Section containing the MGMT methylation analysis."""

    def add_content(self):
        """Add the MGMT analysis content to the report."""
        # Find the latest MGMT results
        last_seen = 0
        mgmt_results = None
        plot_out = None
        specific_sites_file = None

        for file in natsort.natsorted(os.listdir(self.report.output)):
            if file.endswith("_mgmt.csv"):
                count = int(file.split("_")[0])
                if count > last_seen:
                    mgmt_results = pd.read_csv(os.path.join(self.report.output, file))
                    plot_out = os.path.join(
                        self.report.output, file.replace(".csv", ".png")
                    )
                    specific_sites_file = os.path.join(
                        self.report.output, f"{count}_specific_sites.csv"
                    )
                    last_seen = count

        if last_seen == 0 or mgmt_results is None:
            return

        try:
            self.elements.append(PageBreak())
            # Add MGMT section header
            self.elements.append(
                Paragraph("MGMT Promoter Methylation", self.styles.styles["Heading2"])
            )
            self.elements.append(Spacer(1, 12))

            # Extract key metrics
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
            prediction_score = (
                mgmt_results["pred"].iloc[0] if "pred" in mgmt_results.columns else None
            )

            # Create analysis summary table
            summary_data = []
            summary_data.append(
                [
                    Paragraph("Status:", self.styles.styles["Normal"]),
                    Paragraph(
                        methylation_status,
                        ParagraphStyle(
                            "StatusStyle",
                            parent=self.styles.styles["Normal"],
                            textColor=(
                                HexColor("#2563eb")
                                if methylation_status.lower() == "methylated"
                                else HexColor("#d97706")
                            ),
                            fontName="Helvetica-Bold",
                            fontSize=7,
                            leading=9,
                        ),
                    ),
                ]
            )

            if methylation_average is not None:
                summary_data.append(
                    [
                        Paragraph("Average:", self.styles.styles["Normal"]),
                        Paragraph(
                            f"{methylation_average:.1f}%",
                            ParagraphStyle(
                                "ValueStyle",
                                parent=self.styles.styles["Normal"],
                                fontSize=7,
                                leading=9,
                            ),
                        ),
                    ]
                )

            if prediction_score is not None:
                summary_data.append(
                    [
                        Paragraph("Score:", self.styles.styles["Normal"]),
                        Paragraph(
                            f"{prediction_score:.1f}%",
                            ParagraphStyle(
                                "ValueStyle",
                                parent=self.styles.styles["Normal"],
                                fontSize=7,
                                leading=9,
                            ),
                        ),
                    ]
                )

            # Create summary table with more compact styling
            summary_table = Table(summary_data, colWidths=[1 * inch, 1 * inch])
            summary_table.setStyle(
                TableStyle(
                    [
                        # Inherit modern table style
                        *self.MODERN_TABLE_STYLE._cmds,
                        # More compact styling
                        ("FONTSIZE", (0, 0), (-1, -1), 7),  # Smaller font
                        ("LEADING", (0, 0), (-1, -1), 9),  # Tighter line spacing
                        ("ALIGN", (0, 0), (0, -1), "LEFT"),  # Labels left-aligned
                        ("ALIGN", (1, 0), (1, -1), "RIGHT"),  # Values right-aligned
                        # Reduce padding
                        ("TOPPADDING", (0, 0), (-1, -1), 2),
                        ("BOTTOMPADDING", (0, 0), (-1, -1), 2),
                        ("LEFTPADDING", (0, 0), (-1, -1), 4),
                        ("RIGHTPADDING", (0, 0), (-1, -1), 4),
                    ]
                )
            )

            self.elements.append(summary_table)
            self.elements.append(Spacer(1, 6))  # Reduced spacing

            # Add explanation text in a more compact format
            self.elements.append(
                Paragraph(
                    "MGMT promoter methylation status is based on 137 predictive CpG sites (cutoff: 25%). "
                    "Analysis includes strand-specific methylation levels for key CpG pairs.",
                    ParagraphStyle(
                        "Explanation",
                        parent=self.styles.styles["Normal"],
                        fontSize=6,
                        leading=8,
                        textColor=HexColor("#4B5563"),
                    ),
                )
            )
            self.elements.append(Spacer(1, 6))  # Reduced spacing

            # Add specific CpG sites table if available
            if specific_sites_file and os.path.exists(specific_sites_file):
                specific_sites = pd.read_csv(specific_sites_file)
                if not specific_sites.empty:
                    self.elements.append(
                        Paragraph(
                            "Key CpG Sites Analysis", self.styles.styles["Heading3"]
                        )
                    )
                    self.elements.append(Spacer(1, 6))

                    # Create table header
                    cpg_data = [
                        [
                            "Site",  # Shortened header
                            "Position",
                            "Fwd\nCov",  # Shortened headers
                            "Rev\nCov",
                            "Total\nCov",
                            "Combined\nMeth %",
                            "Fwd\nMeth %",
                            "Rev\nMeth %",
                        ]
                    ]

                    # Add data rows
                    for _, row in specific_sites.iterrows():
                        cpg_data.append(
                            [
                                f"Site {row['Site_Label'].split(' ')[1]}",  # Just the number
                                row["Position"].split("/")[0],  # Just first position
                                str(row["Coverage_Forward"]),
                                str(row["Coverage_Reverse"]),
                                str(row["Total_Coverage"]),
                                f"{row['Methylation_Percentage']:.1f}",  # Removed % symbol
                                f"{row['Forward_Methylation']:.1f}",
                                f"{row['Reverse_Methylation']:.1f}",
                            ]
                        )

                    # Create the table with more compact column widths
                    cpg_table = Table(
                        cpg_data,
                        colWidths=[
                            0.5 * inch,  # Site (reduced)
                            0.8 * inch,  # Position (reduced)
                            0.5 * inch,  # Forward Coverage (reduced)
                            0.5 * inch,  # Reverse Coverage (reduced)
                            0.5 * inch,  # Total Coverage (reduced)
                            0.6 * inch,  # Combined Methylation (reduced)
                            0.6 * inch,  # Forward Methylation (reduced)
                            0.6 * inch,  # Reverse Methylation (reduced)
                        ],
                        repeatRows=1,
                    )

                    # Style the table with more compact settings
                    cpg_table.setStyle(
                        TableStyle(
                            [
                                # Inherit modern table style
                                *self.MODERN_TABLE_STYLE._cmds,
                                # More compact styling
                                ("FONTSIZE", (0, 0), (-1, -1), 6),  # Even smaller font
                                (
                                    "LEADING",
                                    (0, 0),
                                    (-1, -1),
                                    7,
                                ),  # Tighter line spacing
                                ("ALIGN", (0, 0), (-1, -1), "CENTER"),
                                ("VALIGN", (0, 0), (-1, -1), "MIDDLE"),
                                # Header styling
                                ("BACKGROUND", (0, 0), (-1, 0), HexColor("#F3F4F6")),
                                ("TEXTCOLOR", (0, 0), (-1, 0), HexColor("#374151")),
                                ("FONTNAME", (0, 0), (-1, 0), "Helvetica-Bold"),
                                # Add subtle gridlines
                                ("GRID", (0, 0), (-1, -1), 0.25, HexColor("#E5E7EB")),
                                # Reduce padding
                                ("TOPPADDING", (0, 0), (-1, -1), 2),
                                ("BOTTOMPADDING", (0, 0), (-1, -1), 2),
                                ("LEFTPADDING", (0, 0), (-1, -1), 2),
                                ("RIGHTPADDING", (0, 0), (-1, -1), 2),
                            ]
                        )
                    )

                    # Add a caption above the table
                    self.elements.append(
                        Paragraph(
                            "Key CpG Sites - Strand-specific methylation analysis",
                            ParagraphStyle(
                                "TableCaption",
                                parent=self.styles.styles["Caption"],
                                fontSize=7,
                                leading=8,
                                spaceBefore=6,
                                spaceAfter=2,
                            ),
                        )
                    )
                    self.elements.append(cpg_table)
                    self.elements.append(Spacer(1, 6))  # Reduced spacing

                    # Add average methylation in a more compact format
                    avg_methylation = specific_sites["Methylation_Percentage"].mean()
                    self.elements.append(
                        Paragraph(
                            f"Mean methylation: {avg_methylation:.1f}%",  # Shortened text
                            ParagraphStyle(
                                "SiteAverage",
                                parent=self.styles.styles["Normal"],
                                fontSize=7,
                                leading=8,
                                textColor=HexColor("#2563eb"),
                            ),
                        )
                    )
                    self.elements.append(Spacer(1, 6))  # Reduced spacing

            # Add file sources information in a more compact format
            file_sources = [
                ["Source", "Location"],  # Shorter headers
                [
                    "Results",  # Shorter labels
                    os.path.basename(
                        os.path.join(self.report.output, f"{last_seen}_mgmt.csv")
                    ),  # Just filename
                ],
                [
                    "Plot",
                    os.path.basename(
                        os.path.join(self.report.output, f"{last_seen}_mgmt.png")
                    ),
                ],
            ]

            if specific_sites_file and os.path.exists(specific_sites_file):
                file_sources.append(
                    [
                        "CpG Sites",
                        os.path.basename(specific_sites_file),
                    ]
                )

            # Create file sources table with compact styling
            sources_table = Table(file_sources, colWidths=[0.8 * inch, 2.5 * inch])
            sources_table.setStyle(
                TableStyle(
                    [
                        # Inherit modern table style
                        *self.MODERN_TABLE_STYLE._cmds,
                        # Compact styling
                        ("FONTSIZE", (0, 0), (-1, -1), 6),  # Smaller font
                        ("LEADING", (0, 0), (-1, -1), 8),  # Tighter line spacing
                        ("ALIGN", (0, 0), (0, -1), "LEFT"),
                        ("ALIGN", (1, 0), (1, -1), "LEFT"),
                        # Header styling
                        ("BACKGROUND", (0, 0), (-1, 0), HexColor("#F3F4F6")),
                        ("TEXTCOLOR", (0, 0), (-1, 0), HexColor("#374151")),
                        ("FONTNAME", (0, 0), (-1, 0), "Helvetica-Bold"),
                        # Reduce padding
                        ("TOPPADDING", (0, 0), (-1, -1), 2),
                        ("BOTTOMPADDING", (0, 0), (-1, -1), 2),
                        ("LEFTPADDING", (0, 0), (-1, -1), 4),
                        ("RIGHTPADDING", (0, 0), (-1, -1), 4),
                    ]
                )
            )

            self.elements.append(sources_table)
            self.elements.append(Spacer(1, 6))  # Reduced spacing

            # Add the methylation plot if it exists
            if plot_out and os.path.exists(plot_out):
                self.elements.append(Image(plot_out, width=6 * inch, height=4 * inch))
                self.elements.append(
                    Paragraph(
                        "MGMT promoter methylation plot showing methylation levels across CpG sites",
                        self.styles.styles["Caption"],
                    )
                )
            else:
                self.elements.append(
                    Paragraph(
                        "MGMT promoter methylation plot is not available due to insufficient coverage.",
                        ParagraphStyle(
                            "Warning",
                            parent=self.styles.styles["Normal"],
                            textColor=HexColor("#DC2626"),
                            fontSize=8,
                            leading=10,
                        ),
                    )
                )

            # Add summary to summary section
            self.summary_elements.append(
                Paragraph("MGMT Promoter Methylation", self.styles.styles["Heading3"])
            )

            summary_text = [f"Status: {methylation_status}"]
            if methylation_average is not None:
                summary_text.append(f"Average methylation: {methylation_average:.1f}%")
            if prediction_score is not None:
                summary_text.append(f"Prediction score: {prediction_score:.1f}%")

            self.summary_elements.append(
                Paragraph(" | ".join(summary_text), self.styles.styles["Normal"])
            )

        except Exception as e:
            logger.error("Error processing MGMT section: %s", str(e), exc_info=True)
            self.elements.append(
                Paragraph(
                    "Error processing MGMT analysis data", self.styles.styles["Normal"]
                )
            )
