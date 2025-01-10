"""
Classification Section for ROBIN Reports.

This module handles the methylation-based classification section of the report,
including results from Sturgeon, Random Forest, NanoDX, and PannanoDX classifiers.
"""

import os
import pandas as pd
import natsort
from reportlab.platypus import PageBreak, Paragraph, Image, Spacer
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
            ["Classifier", "Predicted Class", "Confidence", "Status"]
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
                        status_color = "#059669"  # Green
                    elif confidence_value >= 0.5:
                        confidence_status = "Medium"
                        status_color = "#D97706"  # Amber
                    else:
                        confidence_status = "Low"
                        status_color = "#DC2626"  # Red

                    # Add to summary table with HTML-like color formatting
                    summary_data.append([
                        name,
                        predicted_class,
                        f"{confidence_value:.1%}",
                        f'<font color="{status_color}">{confidence_status}</font>'
                    ])

            except Exception as e:
                logger.error(f"Error processing {name} classification: {str(e)}")
                continue

        # Create and style the summary table
        if len(summary_data) > 1:  # If we have any results
            # Create table with automatic column width calculation
            summary_table = self.create_table(summary_data, repeat_rows=1)
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
                        detailed_data = [["Predicted Class", "Confidence Score"]]

                        for class_name, score in top_predictions.items():
                            confidence = score / 100.0 if name == "Random Forest" else score
                            detailed_data.append([
                                class_name,
                                f"{confidence:.1%}"
                            ])

                        # Create table with right-aligned confidence scores
                        detailed_table = self.create_table(
                            detailed_data,
                            repeat_rows=1,
                            auto_col_width=True
                        )
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
                classifier = row[0]
                predicted = row[1]
                confidence = row[2]
                summary_text.append(f"{classifier}: {predicted} ({confidence})")
            
            self.summary_elements.append(
                Paragraph(
                    " | ".join(summary_text),
                    self.styles.styles["Normal"]
                )
            )
