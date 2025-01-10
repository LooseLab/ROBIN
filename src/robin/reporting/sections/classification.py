"""
Classification Section for ROBIN Reports.

This module handles the methylation-based classification section of the report,
including results from Sturgeon, Random Forest, NanoDX, and PannanoDX classifiers.
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

class ClassificationSection(ReportSection):
    """Section containing the methylation classification results."""

    def add_content(self):
        """Add the classification content to the report."""
        # Add section header
        self.elements.append(
            Paragraph("Methylation Classification", self.styles.styles["Heading2"])
        )
        self.elements.append(Spacer(1, 12))

        # Dictionary of classifiers and their corresponding files
        classifiers = {
            "Sturgeon": "sturgeon_scores.csv",
            "NanoDX": "nanodx_scores.csv",
            "PanNanoDX": "pannanodx_scores.csv",
            "Random Forest": "random_forest_scores.csv"
        }

        # Add summary table of all classifications
        summary_data = [
            [
                Paragraph("Classifier", self.styles.styles["Normal"]),
                Paragraph("Predicted Class", self.styles.styles["Normal"]),
                Paragraph("Confidence", self.styles.styles["Normal"]),
                Paragraph("Status", self.styles.styles["Normal"])
            ]
        ]

        # Process each classifier
        for name, filename in classifiers.items():
            try:
                file_path = None
                for f in os.listdir(self.report.output):
                    if f.lower() == filename.lower():
                        file_path = os.path.join(self.report.output, f)
                        break

                if file_path and os.path.exists(file_path):
                    df = pd.read_csv(file_path)
                    df = df.drop(columns=["timestamp"]) if "timestamp" in df.columns else df
                    df = df.drop(columns=["number_probes"]) if "number_probes" in df.columns else df
                    
                    # Get the last row and find top prediction
                    last_row = df.iloc[-1]
                    top_prediction = last_row.sort_values(ascending=False).head(1)
                    predicted_class = top_prediction.index[0]
                    raw_confidence = float(top_prediction.values[0])
                    
                    # Normalize confidence value
                    confidence_value = raw_confidence / 100.0 if name == "Random Forest" else raw_confidence
                    
                    # Determine confidence level
                    if confidence_value >= 0.75:
                        confidence_status = "High"
                        status_color = HexColor("#059669")  # Green
                    elif confidence_value >= 0.5:
                        confidence_status = "Medium"
                        status_color = HexColor("#D97706")  # Amber
                    else:
                        confidence_status = "Low"
                        status_color = HexColor("#DC2626")  # Red

                    # Add to summary table
                    summary_data.append([
                        Paragraph(name, self.styles.styles["Normal"]),
                        Paragraph(predicted_class, self.styles.styles["Normal"]),
                        Paragraph(f"{confidence_value:.1%}", self.styles.styles["Normal"]),
                        Paragraph(confidence_status, 
                                ParagraphStyle(
                                    'ConfidenceStatus',
                                    parent=self.styles.styles["Normal"],
                                    textColor=status_color,
                                    fontName="Helvetica-Bold"
                                ))
                    ])

            except Exception as e:
                logger.error(f"Error processing {name} classification: {str(e)}")
                continue

        # Create and style the summary table
        if len(summary_data) > 1:  # If we have any results
            summary_table = Table(summary_data, colWidths=[1.5*inch, 3*inch, inch, inch])
            summary_table.setStyle(TableStyle([
                # Header styling
                ('BACKGROUND', (0, 0), (-1, 0), HexColor("#F5F6FA")),
                ('TEXTCOLOR', (0, 0), (-1, 0), HexColor("#2C3E50")),
                ('FONTNAME', (0, 0), (-1, 0), 'Helvetica-Bold'),
                ('FONTSIZE', (0, 0), (-1, -1), 8),
                # Cell padding
                ('TOPPADDING', (0, 0), (-1, -1), 6),
                ('BOTTOMPADDING', (0, 0), (-1, -1), 6),
                ('LEFTPADDING', (0, 0), (-1, -1), 3),
                ('RIGHTPADDING', (0, 0), (-1, -1), 3),
                # Grid styling
                ('GRID', (0, 0), (-1, -1), 0.5, HexColor("#E2E8F0")),
                ('LINEBELOW', (0, 0), (-1, -1), 0.5, HexColor("#CBD5E1")),
                # Alignment
                ('ALIGN', (0, 0), (-1, -1), 'LEFT'),
                ('VALIGN', (0, 0), (-1, -1), 'MIDDLE'),
                # Alternating row colors
                ('ROWBACKGROUNDS', (0, 0), (-1, -1), [HexColor("#FFFFFF"), HexColor("#F8FAFC")])
            ]))
            
            self.elements.append(summary_table)
            self.elements.append(Spacer(1, 12))

            # Add detailed results for each classifier
            self.elements.append(PageBreak())
            self.elements.append(
                Paragraph("Detailed Classification Results", self.styles.styles["Heading2"])
            )
            self.elements.append(Spacer(1, 12))

            # Process each classifier for detailed view
            for name, filename in classifiers.items():
                try:
                    file_path = None
                    for f in os.listdir(self.report.output):
                        if f.lower() == filename.lower():
                            file_path = os.path.join(self.report.output, f)
                            break

                    if file_path and os.path.exists(file_path):
                        df = pd.read_csv(file_path)
                        df = df.drop(columns=["timestamp"]) if "timestamp" in df.columns else df
                        df = df.drop(columns=["number_probes"]) if "number_probes" in df.columns else df
                        
                        last_row = df.iloc[-1]
                        top_predictions = last_row.sort_values(ascending=False).head(10)
                        
                        # Add classifier header
                        self.elements.append(
                            Paragraph(f"{name} Classification", self.styles.styles["Heading3"])
                        )
                        self.elements.append(Spacer(1, 6))

                        # Create detailed table for top 10 predictions
                        detailed_data = [
                            [
                                Paragraph("Predicted Class", self.styles.styles["Normal"]),
                                Paragraph("Confidence Score", self.styles.styles["Normal"])
                            ]
                        ]

                        for class_name, score in top_predictions.items():
                            confidence = score / 100.0 if name == "Random Forest" else score
                            detailed_data.append([
                                Paragraph(class_name, self.styles.styles["Normal"]),
                                Paragraph(f"{confidence:.1%}", self.styles.styles["Normal"])
                            ])

                        detailed_table = Table(detailed_data, colWidths=[4*inch, 1.5*inch])
                        detailed_table.setStyle(TableStyle([
                            # Header styling
                            ('BACKGROUND', (0, 0), (-1, 0), HexColor("#F5F6FA")),
                            ('TEXTCOLOR', (0, 0), (-1, 0), HexColor("#2C3E50")),
                            ('FONTNAME', (0, 0), (-1, 0), 'Helvetica-Bold'),
                            ('FONTSIZE', (0, 0), (-1, -1), 8),
                            # Cell padding
                            ('TOPPADDING', (0, 0), (-1, -1), 6),
                            ('BOTTOMPADDING', (0, 0), (-1, -1), 6),
                            ('LEFTPADDING', (0, 0), (-1, -1), 3),
                            ('RIGHTPADDING', (0, 0), (-1, -1), 3),
                            # Grid styling
                            ('GRID', (0, 0), (-1, -1), 0.5, HexColor("#E2E8F0")),
                            # Alignment
                            ('ALIGN', (0, 0), (0, -1), 'LEFT'),
                            ('ALIGN', (1, 0), (1, -1), 'RIGHT'),
                            ('VALIGN', (0, 0), (-1, -1), 'MIDDLE'),
                            # Alternating row colors
                            ('ROWBACKGROUNDS', (0, 0), (-1, -1), [HexColor("#FFFFFF"), HexColor("#F8FAFC")])
                        ]))

                        self.elements.append(detailed_table)
                        self.elements.append(Spacer(1, 12))

                except Exception as e:
                    logger.error(f"Error processing detailed {name} classification: {str(e)}")
                    continue

            # Add explanation text
            self.elements.append(
                Paragraph(
                    "Note: Classification confidence levels are defined as High (≥75%), "
                    "Medium (≥50%), and Low (<50%). Multiple classifiers may provide different "
                    "results based on their training data and methodology.",
                    ParagraphStyle(
                        'Explanation',
                        parent=self.styles.styles["Normal"],
                        fontSize=8,
                        leading=10,
                        textColor=HexColor("#4B5563")
                    )
                )
            )

        else:
            self.elements.append(
                Paragraph(
                    "No classification results available.",
                    self.styles.styles["Normal"]
                )
            )

        # Add summary to summary section
        self.summary_elements.append(
            Paragraph("Methylation Classification", self.styles.styles["Heading3"])
        )
        if len(summary_data) > 1:
            summary_text = []
            for row in summary_data[1:]:  # Skip header row
                classifier = row[0].text
                predicted = row[1].text
                confidence = row[2].text
                summary_text.append(f"{classifier}: {predicted} ({confidence})")
            
            self.summary_elements.append(
                Paragraph(
                    " | ".join(summary_text),
                    self.styles.styles["Normal"]
                )
            )
