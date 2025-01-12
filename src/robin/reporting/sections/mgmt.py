"""
MGMT Analysis Section for ROBIN Reports.

This module handles the MGMT (O6-methylguanine-DNA methyltransferase) promoter methylation analysis section of the report.
"""

import os
import pandas as pd
import natsort
from reportlab.lib.units import inch
from reportlab.platypus import PageBreak, Paragraph, Image, Spacer, Table, TableStyle
from reportlab.lib.colors import HexColor, white
from reportlab.lib.styles import ParagraphStyle
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
                mgmt_results["pred"].iloc[0]
                if "pred" in mgmt_results.columns
                else None
            )

            # Create analysis summary table
            summary_data = []
            summary_data.append([
                Paragraph("Methylation Status:", self.styles.styles["Normal"]),
                Paragraph(methylation_status, 
                         ParagraphStyle(
                             'StatusStyle',
                             parent=self.styles.styles["Normal"],
                             textColor=HexColor("#2563eb") if methylation_status.lower() == "methylated" 
                             else HexColor("#d97706"),
                             fontName="Helvetica-Bold",
                             fontSize=8,
                             leading=10
                         ))
            ])
            
            if methylation_average is not None:
                summary_data.append([
                    Paragraph("Average Methylation:", self.styles.styles["Normal"]),
                    Paragraph(f"{methylation_average:.1f}%", self.styles.styles["Normal"])
                ])
            
            if prediction_score is not None:
                summary_data.append([
                    Paragraph("Prediction Score:", self.styles.styles["Normal"]),
                    Paragraph(f"{prediction_score:.1f}%", self.styles.styles["Normal"])
                ])

            # Create summary table with styling
            summary_table = Table(summary_data, colWidths=[2 * inch, 1.5 * inch])
            summary_table.setStyle(TableStyle([
                # Inherit modern table style
                *self.MODERN_TABLE_STYLE._cmds,
                
                # Preserve specific alignments
                ('ALIGN', (0, 0), (0, -1), 'LEFT'),  # Labels left-aligned
                ('ALIGN', (1, 0), (1, -1), 'RIGHT'),  # Values right-aligned
            ]))

            self.elements.append(summary_table)
            self.elements.append(Spacer(1, 12))

            # Add file sources information
            file_sources = [
                ["Data Source", "File Location"],
                ["MGMT Results", os.path.join(self.report.output, f"{last_seen}_mgmt.csv")],
                ["MGMT Plot", os.path.join(self.report.output, f"{last_seen}_mgmt.png")]
            ]
            
            # Create file sources table
            sources_table = Table(file_sources, colWidths=[2 * inch, 4 * inch])
            sources_table.setStyle(TableStyle([
                # Inherit modern table style
                *self.MODERN_TABLE_STYLE._cmds,
                
                # Preserve specific alignments
                ('ALIGN', (0, 0), (0, -1), 'LEFT'),  # Labels left-aligned
                ('ALIGN', (1, 0), (1, -1), 'LEFT'),  # Values left-aligned
            ]))
            
            self.elements.append(sources_table)
            self.elements.append(Spacer(1, 12))

            # Add explanation text
            self.elements.append(
                Paragraph(
                    "MGMT promoter methylation status is assessed based on a logistic regression prediction model "
                    "that considers 137 most predictive CpG sites in the MGMT promoter region. "
                    "The methylation cutoff is assigned as 25%.",
                    ParagraphStyle(
                        'Explanation',
                        parent=self.styles.styles["Normal"],
                        fontSize=8,
                        leading=10,
                        textColor=HexColor("#4B5563")
                    )
                )
            )
            self.elements.append(Spacer(1, 12))

            # Add the methylation plot if it exists
            if plot_out and os.path.exists(plot_out):
                self.elements.append(Image(plot_out, width=6 * inch, height=4 * inch))
                self.elements.append(
                    Paragraph(
                        "MGMT promoter methylation plot showing methylation levels across CpG sites",
                        self.styles.styles["Caption"]
                    )
                )
            else:
                self.elements.append(
                    Paragraph(
                        "MGMT promoter methylation plot is not available due to insufficient coverage.",
                        ParagraphStyle(
                            'Warning',
                            parent=self.styles.styles["Normal"],
                            textColor=HexColor("#DC2626"),
                            fontSize=8,
                            leading=10
                        )
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
                Paragraph(
                    " | ".join(summary_text),
                    self.styles.styles["Normal"]
                )
            )

        except Exception as e:
            logger.error("Error processing MGMT section: %s", str(e), exc_info=True)
            self.elements.append(
                Paragraph(
                    "Error processing MGMT analysis data",
                    self.styles.styles["Normal"]
                )
            )
