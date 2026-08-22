"""
MGMT Analysis Section for ROBIN Reports.

This module handles the MGMT (O6-methylguanine-DNA methyltransferase) promoter methylation analysis section of the report.
"""

import io
import logging
import os

import natsort
import pandas as pd
from reportlab.lib.colors import HexColor
from reportlab.lib.styles import ParagraphStyle
from reportlab.lib.units import inch
from reportlab.platypus import Image, PageBreak, Paragraph, Spacer, Table, TableStyle

from ..sections.base import ReportSection

logger = logging.getLogger(__name__)


def _load_image_buffer(path: str) -> io.BytesIO | None:
    """Load an image file into memory for ReportLab (avoids late path/casing issues)."""
    resolved = os.path.realpath(os.path.expanduser(path))
    if not os.path.isfile(resolved):
        return None
    try:
        with open(resolved, "rb") as fh:
            buf = io.BytesIO(fh.read())
        buf.seek(0)
        return buf
    except OSError as exc:
        logger.warning("Failed to read image %s: %s", resolved, exc)
        return None


def _extract_mgmt_specific_sites_from_bed(bed_path: str) -> pd.DataFrame:
    """
    Extract MGMT-specific CpG site methylation data from a BED file.

    This function mirrors the logic from the GUI component to extract
    strand-specific methylation data for key CpG sites.

    Args:
        bed_path: Path to the BED file

    Returns:
        DataFrame with CpG site methylation data, or empty DataFrame if extraction fails
    """
    try:
        df = pd.read_csv(bed_path, sep="\t", header=None)

        # Check if column 10 contains space-separated values (old format)
        has_space_separated_col10 = False
        if df.shape[1] > 10 and len(df) > 0:
            sample_val = str(df.iloc[0, 9])
            has_space_separated_col10 = " " in sample_val or "\t" in sample_val

        # Check if this is the new bedmethyl format (separate columns) or old format
        if df.shape[1] >= 12 and not has_space_separated_col10:
            # New bedmethyl format with separate columns
            cols = [
                "Chromosome",
                "Start",
                "End",
                "Modified_Base_Code",
                "Score",
                "Strand",
                "Start2",
                "End2",
                "RGB",
                "Nvalid_cov",
                "Fraction_Modified",
                "Nmod",
            ]
            num_cols_to_read = min(len(cols), df.shape[1])
            df = df.iloc[:, :num_cols_to_read]
            df.columns = cols[:num_cols_to_read]

            df["Nvalid_cov"] = df["Nvalid_cov"].astype(float)
            df["Fraction_Modified"] = df["Fraction_Modified"].astype(float)
            if "Nmod" in df.columns:
                df["Nmod"] = df["Nmod"].astype(float)
            else:
                df["Nmod"] = df["Nvalid_cov"] * df["Fraction_Modified"]

            df["Start"] = df["Start"].astype(int)
            df["Coverage"] = df["Nvalid_cov"]
            df["Modified_Fraction"] = df["Fraction_Modified"] * 100.0
        elif df.shape[1] >= 10:
            # Old format with space-separated Coverage_Info in column 10
            cols = [
                "Chromosome",
                "Start",
                "End",
                "Name",
                "Score",
                "Strand",
                "Start2",
                "End2",
                "RGB",
                "Coverage_Info",
            ]
            df = df.iloc[:, : len(cols)]
            df.columns = cols

            cov_split = df["Coverage_Info"].astype(str).str.split()
            df["Coverage"] = cov_split.str[0].astype(float)
            fraction_val = cov_split.str[1].astype(float).fillna(0.0)
            is_percentage = (
                (fraction_val > 1.0).any() if len(fraction_val) > 0 else False
            )

            if is_percentage:
                df["Modified_Fraction"] = fraction_val
                df["Fraction_Modified"] = df["Modified_Fraction"] / 100.0
            else:
                df["Fraction_Modified"] = fraction_val
                df["Modified_Fraction"] = df["Fraction_Modified"] * 100.0

            df["Nvalid_cov"] = df["Coverage"]
            df["Nmod"] = df["Coverage"] * df["Fraction_Modified"]
            df["Start"] = df["Start"].astype(int)
        else:
            return pd.DataFrame()

        # Define CpG pairs and labels
        cpg_pairs = [
            (129467255, 129467256),
            (129467258, 129467259),
            (129467262, 129467263),
            (129467272, 129467273),
        ]
        label_map = {
            "129467255/129467256": "1",
            "129467258/129467259": "2",
            "129467262/129467263": "3",
            "129467272/129467273": "4",
        }

        rows = []
        for p1, p2 in cpg_pairs:
            pos_key = f"{p1}/{p2}"
            site_label = label_map.get(pos_key, "Unknown")

            # Check forward strand reads at position p1
            fwd_p1 = df[
                (df["Chromosome"] == "chr10")
                & (df["Start"] == p1 - 1)
                & (df["Strand"] == "+")
            ]

            # Check reverse strand reads at position p2
            rev_p2 = df[
                (df["Chromosome"] == "chr10")
                & (df["Start"] == p2 - 1)
                & (df["Strand"] == "-")
            ]

            # Get forward strand data
            if not fwd_p1.empty:
                cov_f = float(fwd_p1["Nvalid_cov"].iloc[0])
                mf = float(fwd_p1["Fraction_Modified"].iloc[0])
                nmod_f = float(fwd_p1["Nmod"].iloc[0])
                meth_fwd_count = int(round(nmod_f))
                meth_fwd_pct = float(fwd_p1["Modified_Fraction"].iloc[0])
            else:
                cov_f = 0.0
                mf = 0.0
                meth_fwd_count = 0
                meth_fwd_pct = 0.0

            # Get reverse strand data
            if not rev_p2.empty:
                cov_r = float(rev_p2["Nvalid_cov"].iloc[0])
                mr = float(rev_p2["Fraction_Modified"].iloc[0])
                nmod_r = float(rev_p2["Nmod"].iloc[0])
                meth_rev_count = int(round(nmod_r))
                meth_rev_pct = float(rev_p2["Modified_Fraction"].iloc[0])
            else:
                cov_r = 0.0
                mr = 0.0
                meth_rev_count = 0
                meth_rev_pct = 0.0

            # Only add row if we have data
            if cov_f > 0 or cov_r > 0:
                tot = cov_f + cov_r
                weighted = ((cov_f * mf) + (cov_r * mr)) / tot if tot > 0 else 0.0
                weighted_pct = weighted * 100.0

                rows.append(
                    {
                        "Site_Label": f"Site {site_label}",
                        "Position": pos_key,
                        "Coverage_Forward": int(cov_f),
                        "Coverage_Reverse": int(cov_r),
                        "Total_Coverage": int(tot),
                        "Methylation_Percentage": weighted_pct,
                        "Forward_Methylation": meth_fwd_pct,
                        "Reverse_Methylation": meth_rev_pct,
                        "Forward_Methylated_Count": meth_fwd_count,
                        "Reverse_Methylated_Count": meth_rev_count,
                    }
                )

        return pd.DataFrame(rows)
    except Exception as e:
        logger.warning(f"Failed to extract CpG sites from BED file {bed_path}: {e}")
        return pd.DataFrame()


class MGMTSection(ReportSection):
    """Section containing the MGMT methylation analysis."""

    def add_content(self):
        """Add the MGMT analysis content to the report."""
        # Find the latest MGMT results
        last_seen = 0
        mgmt_results = None
        plot_out = None
        specific_sites_file = None
        bed_file_used = None
        specific_sites = None

        # First, look for "final_mgmt.csv" files (highest priority)
        final_files = [
            f for f in os.listdir(self.report.output) if f == "final_mgmt.csv"
        ]
        if final_files:
            file = final_files[0]
            mgmt_results = pd.read_csv(os.path.join(self.report.output, file))
            plot_out = os.path.join(self.report.output, file.replace(".csv", ".png"))
            specific_sites_file = os.path.join(
                self.report.output, "final_specific_sites.csv"
            )
            last_seen = 999999  # High number to indicate final result
        else:
            # Fallback to numeric-prefixed files
            for file in natsort.natsorted(os.listdir(self.report.output)):
                if file.endswith("_mgmt.csv"):
                    try:
                        # Try to extract numeric prefix from filename
                        prefix = file.split("_")[0]
                        count = int(prefix)
                        if count > last_seen:
                            mgmt_results = pd.read_csv(
                                os.path.join(self.report.output, file)
                            )
                            plot_out = os.path.join(
                                self.report.output, file.replace(".csv", ".png")
                            )
                            specific_sites_file = os.path.join(
                                self.report.output, f"{count}_specific_sites.csv"
                            )
                            last_seen = count
                    except ValueError:
                        # Skip files that don't start with a numeric prefix
                        logger.debug(f"Skipping file with non-numeric prefix: {file}")
                        continue

        if last_seen == 0 or mgmt_results is None:
            return

        try:
            self.elements.append(PageBreak())
            # Add MGMT section header
            self.elements.append(
                Paragraph("MGMT Promoter Methylation", self.styles.styles["Heading2"])
            )
            self.elements.append(Spacer(1, 4))

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
            self.elements.append(Spacer(1, 4))  # Reduced spacing

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
            self.elements.append(Spacer(1, 4))  # Reduced spacing

            # Try to get CpG sites data from CSV file first, then fall back to BED file
            # First, try to load from CSV file if it exists
            if specific_sites_file and os.path.exists(specific_sites_file):
                try:
                    specific_sites = pd.read_csv(specific_sites_file)
                    logger.debug(f"Loaded CpG sites from CSV: {specific_sites_file}")
                except Exception as e:
                    logger.warning(f"Failed to load CpG sites CSV: {e}")

            # If CSV doesn't exist or is empty, try to extract from BED file
            if specific_sites is None or specific_sites.empty:
                # Determine BED file path
                if last_seen == 999999:  # Final file
                    bed_candidates = [
                        os.path.join(self.report.output, "final_mgmt.bed"),
                        os.path.join(self.report.output, "final_mgmt_mgmt.bed"),
                    ]
                else:
                    bed_candidates = [
                        os.path.join(self.report.output, f"{last_seen}_mgmt.bed"),
                        os.path.join(self.report.output, f"{last_seen}_mgmt_mgmt.bed"),
                    ]

                bed_path = None
                for candidate in bed_candidates:
                    if os.path.exists(candidate):
                        bed_path = candidate
                        break

                if bed_path:
                    try:
                        specific_sites = _extract_mgmt_specific_sites_from_bed(bed_path)
                        bed_file_used = bed_path
                        logger.debug(f"Extracted CpG sites from BED: {bed_path}")
                    except Exception as e:
                        logger.warning(f"Failed to extract CpG sites from BED: {e}")
                        specific_sites = None

            # Add specific CpG sites table if we have data
            if specific_sites is not None and not specific_sites.empty:
                self.elements.append(
                    Paragraph(
                        "MGMT CpG Site Methylation Data", self.styles.styles["Heading3"]
                    )
                )
                self.elements.append(Spacer(1, 4))

                # Create table header - include methylation counts if available
                has_meth_counts = "Forward_Methylated_Count" in specific_sites.columns
                cpg_data = [
                    [
                        "Site",
                        "Position",
                        "Fwd\nCov",
                        "Rev\nCov",
                        "Total\nCov",
                        "Meth %",
                        "Fwd\nMeth %",
                        "Rev\nMeth %",
                    ]
                ]

                # Add methylation count columns if available
                if has_meth_counts:
                    cpg_data[0].extend(["Fwd\nMeth\nCount", "Rev\nMeth\nCount"])

                # Add data rows
                for _, row in specific_sites.iterrows():
                    # Extract site number from Site_Label
                    site_label = str(row.get("Site_Label", ""))
                    site_num = (
                        site_label.split(" ")[-1] if " " in site_label else site_label
                    )

                    # Extract position
                    position = str(row.get("Position", ""))
                    pos_display = (
                        position.split("/")[0] if "/" in position else position
                    )

                    row_data = [
                        f"Site {site_num}",
                        pos_display,
                        str(int(row.get("Coverage_Forward", 0))),
                        str(int(row.get("Coverage_Reverse", 0))),
                        str(int(row.get("Total_Coverage", 0))),
                        f"{row.get('Methylation_Percentage', 0.0):.1f}",
                        f"{row.get('Forward_Methylation', 0.0):.1f}",
                        f"{row.get('Reverse_Methylation', 0.0):.1f}",
                    ]

                    # Add methylation counts if available
                    if has_meth_counts:
                        row_data.extend(
                            [
                                str(int(row.get("Forward_Methylated_Count", 0))),
                                str(int(row.get("Reverse_Methylated_Count", 0))),
                            ]
                        )

                    cpg_data.append(row_data)

                # Create the table with more compact column widths
                # Adjust column widths based on whether counts are included
                if has_meth_counts:
                    col_widths = [
                        0.5 * inch,  # Site
                        0.8 * inch,  # Position
                        0.5 * inch,  # Forward Coverage
                        0.5 * inch,  # Reverse Coverage
                        0.5 * inch,  # Total Coverage
                        0.5 * inch,  # Combined Methylation %
                        0.5 * inch,  # Forward Methylation %
                        0.5 * inch,  # Reverse Methylation %
                        0.5 * inch,  # Forward Methylated Count
                        0.5 * inch,  # Reverse Methylated Count
                    ]
                else:
                    col_widths = [
                        0.5 * inch,  # Site
                        0.8 * inch,  # Position
                        0.5 * inch,  # Forward Coverage
                        0.5 * inch,  # Reverse Coverage
                        0.5 * inch,  # Total Coverage
                        0.6 * inch,  # Combined Methylation %
                        0.6 * inch,  # Forward Methylation %
                        0.6 * inch,  # Reverse Methylation %
                    ]

                cpg_table = Table(
                    cpg_data,
                    colWidths=col_widths,
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
                        "MGMT CpG Site Methylation Data - Strand-specific methylation analysis",
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
                self.elements.append(Spacer(1, 4))  # Reduced spacing

                # Add average methylation in a more compact format
                avg_methylation = specific_sites["Methylation_Percentage"].mean()
                self.elements.append(
                    Paragraph(
                        f"Mean methylation: {avg_methylation:.1f}%",
                        ParagraphStyle(
                            "SiteAverage",
                            parent=self.styles.styles["Normal"],
                            fontSize=7,
                            leading=8,
                            textColor=HexColor("#2563eb"),
                        ),
                    )
                )
                self.elements.append(Spacer(1, 4))  # Reduced spacing

            # Add file sources information in a more compact format
            if last_seen == 999999:  # Final file
                results_file = "final_mgmt.csv"
                plot_file = "final_mgmt.png"
            else:
                results_file = f"{last_seen}_mgmt.csv"
                plot_file = f"{last_seen}_mgmt.png"

            file_sources = [
                ["Source", "Location"],  # Shorter headers
                [
                    "Results",  # Shorter labels
                    results_file,
                ],
                [
                    "Plot",
                    plot_file,
                ],
            ]

            if specific_sites_file and os.path.exists(specific_sites_file):
                file_sources.append(
                    [
                        "CpG Sites",
                        os.path.basename(specific_sites_file),
                    ]
                )
            elif bed_file_used:
                file_sources.append(
                    [
                        "CpG Sites",
                        os.path.basename(bed_file_used),
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
            self.elements.append(Spacer(1, 4))  # Reduced spacing

            # Try to add the methylation plot
            plot_added = False

            # First, try to load from existing PNG file (e.g. final_mgmt.png)
            if plot_out:
                plot_buf = _load_image_buffer(plot_out)
                if plot_buf is not None:
                    try:
                        self.elements.append(
                            Image(plot_buf, width=6 * inch, height=4 * inch)
                        )
                        self.elements.append(
                            Paragraph(
                                "MGMT promoter methylation plot showing methylation levels across CpG sites",
                                self.styles.styles["Caption"],
                            )
                        )
                        plot_added = True
                    except Exception as e:
                        logger.warning(f"Failed to load MGMT plot from file: {e}")

            # If not found, try to generate from BAM file using locus_figure
            if not plot_added:
                # Try sorted BAM first, then fall back to unsorted
                bam_candidates = [
                    os.path.join(self.report.output, "mgmt_sorted.bam"),
                    os.path.join(self.report.output, "mgmt.bam"),
                ]

                bam_path = None
                for candidate in bam_candidates:
                    if os.path.exists(candidate):
                        bam_path = candidate
                        break

                if bam_path:
                    try:
                        import warnings

                        import matplotlib.pyplot as plt

                        from robin.analysis.methylation_wrapper import locus_figure

                        logger.info(
                            f"Generating MGMT plot from BAM file for report: {os.path.basename(bam_path)}"
                        )
                        fig = locus_figure(
                            interval="chr10:129466536-129467536",
                            bam_path=bam_path,
                            motif="CG",
                            mods="m",
                        )

                        # Embed in memory (same pattern as CNV/coverage sections).
                        # Avoids ReportLab failing later if path casing differs on disk.
                        plot_buf = io.BytesIO()
                        with warnings.catch_warnings():
                            warnings.simplefilter("ignore", UserWarning)
                            fig.savefig(
                                plot_buf, format="png", dpi=150, bbox_inches="tight"
                            )
                        plot_buf.seek(0)
                        plt.close(fig)

                        if plot_buf.getbuffer().nbytes > 0:
                            self.elements.append(
                                Image(plot_buf, width=6 * inch, height=4 * inch)
                            )
                            self.elements.append(
                                Paragraph(
                                    "MGMT promoter methylation plot showing methylation levels across CpG sites",
                                    self.styles.styles["Caption"],
                                )
                            )
                            plot_added = True
                        else:
                            logger.warning(
                                "MGMT plot generation produced empty image data"
                            )

                    except Exception as e:
                        logger.warning(f"Failed to generate MGMT plot from BAM: {e}")

            # If still not added, show error message
            if not plot_added:
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

        # Build export DataFrames
        try:
            summary_row = {
                "Status": methylation_status,
                "AveragePercent": (
                    round(float(methylation_average), 2)
                    if methylation_average is not None
                    else None
                ),
                "PredictionScorePercent": (
                    round(float(prediction_score), 2)
                    if prediction_score is not None
                    else None
                ),
                "ResultsFile": (
                    "final_mgmt.csv"
                    if last_seen == 999999
                    else (f"{last_seen}_mgmt.csv" if last_seen else None)
                ),
                "PlotFile": (
                    "final_mgmt.png"
                    if last_seen == 999999
                    else (f"{last_seen}_mgmt.png" if last_seen else None)
                ),
                "CpGSitesFile": (
                    os.path.basename(specific_sites_file)
                    if specific_sites_file and os.path.exists(specific_sites_file)
                    else (os.path.basename(bed_file_used) if bed_file_used else None)
                ),
            }
            self.export_frames["mgmt_summary"] = pd.DataFrame([summary_row])

            # Export CpG sites data if available
            if specific_sites is not None and not specific_sites.empty:
                self.export_frames["mgmt_cpg_sites"] = specific_sites
        except Exception as ex:
            logger.error(
                "Error building MGMT export DataFrames: %s", str(ex), exc_info=True
            )
