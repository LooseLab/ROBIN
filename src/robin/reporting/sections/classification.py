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
from reportlab.lib.units import inch
from ..sections.base import ReportSection
import logging
import matplotlib.pyplot as plt
import matplotlib
matplotlib.use('Agg')
import io

logger = logging.getLogger(__name__)

class ClassificationSection(ReportSection):
    """Section containing the methylation classification results."""

    def _create_time_plot(self, df, classifier_name):
        """Create a time series plot for classifier predictions.
        
        Args:
            df: DataFrame containing classification data
            classifier_name: Name of the classifier
            
        Returns:
            BytesIO object containing the plot image
        """
        # Drop non-classification columns
        df = df.drop(columns=["number_probes"]) if "number_probes" in df.columns else df
        
        # Convert confidence values to percentages if not already
        if classifier_name != "Random Forest":  # Random Forest is already in percentages
            df = df * 100
        
        # Get top classifications (those that exceed 5% at any point)
        threshold = 5
        top_classes = df.columns[df.max() > threshold]
        
        # Get current highest confidence for subtitle
        last_row = df.iloc[-1]
        top_prediction = last_row.sort_values(ascending=False).head(1)
        predicted_class = top_prediction.index[0]
        confidence_value = float(top_prediction.values[0])
        
        # Create the plot with white background and adjusted figure size for better spacing
        plt.figure(figsize=(8, 5), facecolor='white')  # Increased height for better spacing
        ax = plt.gca()
        ax.set_facecolor('white')
        
        # Website-like colors - extended for more classes
        colors = {
            # Sturgeon colors
            "Embryonal - HGNET - BCOR": "#34C759",  # Green
            "Ependymal - EPN - PF B": "#FF9500",    # Orange
            "Ependymal - EPN - RELA": "#FF2D55",    # Red
            # Additional colors for other classifications
            "color4": "#007AFF",  # Blue
            "color5": "#5856D6",  # Purple
            "color6": "#FF3B30",  # Red-Orange
            "color7": "#5AC8FA",  # Light Blue
            "color8": "#4CD964",  # Light Green
        }
        
        # Convert index to datetime
        df.index = pd.to_datetime(df.index, unit='ms')
        
        # Plot each classification
        for idx, column in enumerate(top_classes):
            # Get color from map or use from additional colors
            if column in colors:
                color = colors[column]
            else:
                color = colors[f"color{(idx % 5) + 4}"]  # Cycle through additional colors
                
            plt.plot(df.index, df[column], 
                    'o-',  # Line with circles
                    label=column, 
                    color=color,
                    linewidth=1,
                    markersize=3,
                    markeredgewidth=0)
        
        # Add title and subtitle with better spacing
        plt.suptitle(f'{classifier_name} Classification Confidence Over Time', 
                    y=0.95, fontsize=10)
        plt.title(f'Current highest confidence: {predicted_class} ({confidence_value:.1f}%)',
                 pad=10, fontsize=9)
        
        # Format x-axis to show HH:MM
        ax.xaxis.set_major_formatter(matplotlib.dates.DateFormatter('%H:%M'))
        plt.xlabel('Time', fontsize=9, labelpad=5)
        
        # Format y-axis with better spacing
        plt.ylabel('Confidence (%)', fontsize=9, labelpad=10)
        ax.yaxis.set_major_formatter(matplotlib.ticker.PercentFormatter(decimals=0))
        
        # Customize grid
        plt.grid(True, linestyle='--', alpha=0.2)
        
        # Set axis ranges
        plt.ylim(0, 100)
        
        # Customize legend with better positioning
        plt.legend(bbox_to_anchor=(0, 1.15, 1, 0),
                  loc='lower left',
                  mode='expand',
                  ncol=3,
                  fontsize=8,
                  frameon=False,
                  markerscale=2)
        
        # Remove spines
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.spines['left'].set_alpha(0.2)
        ax.spines['bottom'].set_alpha(0.2)
        
        # Adjust tick parameters
        ax.tick_params(axis='both', labelsize=8)
        plt.xticks(rotation=0)
        
        # Adjust layout to prevent label cutoff
        plt.subplots_adjust(top=0.85, bottom=0.15, left=0.1, right=0.95)
        
        # Save plot to bytes buffer with high DPI for crisp rendering
        buf = io.BytesIO()
        plt.savefig(buf, format='png', dpi=300, bbox_inches='tight', facecolor='white')
        plt.close()
        buf.seek(0)
        return buf

    def add_content(self):
        """Add the classification content to the report."""
        # Add summary to summary section
        self.summary_elements.append(
            Paragraph("Methylation Classification", self.styles.styles["Heading3"])
        )
        
        # Add section header
        self.elements.append(
            Paragraph("Methylation Classification", self.styles.styles["Heading2"])
        )
        self.elements.append(Spacer(1, 6))

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
                    
                    # Create and add time series plot for each classifier
                    try:
                        # Set timestamp as index for the plot
                        df.set_index("timestamp", inplace=True)
                        plot_buf = self._create_time_plot(df, name)
                        self.add_figure(
                            plot_buf,
                            caption=f"{name} Classification Confidence Over Time",
                            width=6*inch,
                            height=3*inch
                        )
                        self.elements.append(Spacer(1, 12))
                    except Exception as e:
                        logger.error(f"Error creating {name} time series plot: {str(e)}")
                    
                    # Continue with existing classification processing
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
            self.summary_elements.append(summary_table)
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

            Explanation_text = Paragraph(
                    "Note: Classification confidence levels are defined as High (>75%), "
                    "Medium (>50%), and Low (<50%). Multiple classifiers may provide different "
                    "results based on their training data and methodology.",
                    ParagraphStyle(
                        'Explanation',
                        parent=self.styles.styles["Normal"],
                        fontSize=8,
                        leading=10,
                        textColor=HexColor("#4B5563")
                    )
                )

            # Add explanation text
            self.elements.append(
                Explanation_text
            )

        else:
            self.elements.append(
                Paragraph(
                    "No classification results available.",
                    self.styles.styles["Normal"]
                )
            )

        
    
        # Add explanation text
        self.summary_elements.append(
            Explanation_text
        )
