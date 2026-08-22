"""
pdf_styles.py

This module contains PDF-specific styling utilities and helper functions
for implementing Material Design 3 and Apple HIG principles in PDF reports.
"""

from reportlab.lib import colors
from reportlab.lib.units import inch
from reportlab.platypus import Paragraph, Spacer, Table, TableStyle


class PDFStyleUtils:
    """Utility class for applying Material Design 3 styling to PDF elements."""

    @staticmethod
    def create_section_header(text, style_name="HeadlineLarge", styles_dict=None):
        """Create a section header with proper M3 spacing and typography."""
        if styles_dict and style_name in styles_dict:
            return Paragraph(text, styles_dict[style_name])
        else:
            # Fallback to default style
            from reportlab.lib.styles import getSampleStyleSheet

            fallback_styles = getSampleStyleSheet()
            return Paragraph(text, fallback_styles["Heading1"])

    @staticmethod
    def create_section_spacing():
        """Create consistent spacing between sections using M3 principles."""
        return Spacer(1, 16)  # 16pt spacing (M3 spacing unit)

    @staticmethod
    def create_subsection_header(text, style_name="HeadlineMedium", styles_dict=None):
        """Create a subsection header with M3 typography."""
        if styles_dict and style_name in styles_dict:
            return [Spacer(1, 8), Paragraph(text, styles_dict[style_name])]
        else:
            # Fallback to default style
            from reportlab.lib.styles import getSampleStyleSheet

            fallback_styles = getSampleStyleSheet()
            return [Spacer(1, 8), Paragraph(text, fallback_styles["Heading2"])]

    @staticmethod
    def create_body_text(text, style_name="Normal", styles_dict=None):
        """Create body text (10pt)."""
        if styles_dict and style_name in styles_dict:
            return Paragraph(text, styles_dict[style_name])
        else:
            # Fallback to default style
            from reportlab.lib.styles import getSampleStyleSheet

            fallback_styles = getSampleStyleSheet()
            return Paragraph(text, fallback_styles["Normal"])

    @staticmethod
    def create_caption(text, style_name="Caption", styles_dict=None):
        """Create caption text with M3 typography."""
        if styles_dict and style_name in styles_dict:
            return Paragraph(text, styles_dict[style_name])
        else:
            # Fallback to default style
            from reportlab.lib.styles import getSampleStyleSheet

            fallback_styles = getSampleStyleSheet()
            return Paragraph(text, fallback_styles["Normal"])

    @staticmethod
    def create_metric_display(
        value, label, unit="", style_name="Metric", styles_dict=None
    ):
        """Create a metric display with value, label, and optional unit."""
        if unit:
            metric_text = f"{value} {unit}"
        else:
            metric_text = str(value)

        # Create the metric value
        if styles_dict and style_name in styles_dict:
            metric_para = Paragraph(metric_text, styles_dict[style_name])
        else:
            # Fallback to default style
            from reportlab.lib.styles import getSampleStyleSheet

            fallback_styles = getSampleStyleSheet()
            metric_para = Paragraph(metric_text, fallback_styles["Normal"])

        # Create the label below
        if styles_dict and "LabelMedium" in styles_dict:
            label_para = Paragraph(label, styles_dict["LabelMedium"])
        else:
            # Fallback to default style
            from reportlab.lib.styles import getSampleStyleSheet

            fallback_styles = getSampleStyleSheet()
            label_para = Paragraph(label, fallback_styles["Normal"])

        return [metric_para, Spacer(1, 4), label_para]

    @staticmethod
    def create_info_card(title, content, card_style="SummaryCard", styles_dict=None):
        """Create an information card with title and content."""
        elements = []

        # Add title
        if title:
            if styles_dict and "TitleMedium" in styles_dict:
                elements.append(Paragraph(title, styles_dict["TitleMedium"]))
            else:
                # Fallback to default style
                from reportlab.lib.styles import getSampleStyleSheet

                fallback_styles = getSampleStyleSheet()
                elements.append(Paragraph(title, fallback_styles["Normal"]))
            elements.append(Spacer(1, 8))

        # Add content
        if isinstance(content, str):
            if styles_dict and "Normal" in styles_dict:
                elements.append(Paragraph(content, styles_dict["Normal"]))
            else:
                # Fallback to default style
                from reportlab.lib.styles import getSampleStyleSheet

                fallback_styles = getSampleStyleSheet()
                elements.append(Paragraph(content, fallback_styles["Normal"]))
        elif isinstance(content, list):
            for item in content:
                if isinstance(item, str):
                    if styles_dict and "Normal" in styles_dict:
                        elements.append(Paragraph(f"• {item}", styles_dict["Normal"]))
                    else:
                        # Fallback to default style
                        from reportlab.lib.styles import getSampleStyleSheet

                        fallback_styles = getSampleStyleSheet()
                        elements.append(
                            Paragraph(f"• {item}", fallback_styles["Normal"])
                        )
                else:
                    elements.append(item)

        # Wrap in card style if available
        if styles_dict and card_style in styles_dict:
            card_style_obj = styles_dict[card_style]
            card_elements = []
            for element in elements:
                if hasattr(element, "style") and hasattr(element, "text"):
                    # Clone the element with card styling
                    from reportlab.lib.styles import ParagraphStyle

                    new_style = ParagraphStyle(
                        name="CardStyle",
                        parent=element.style,
                        backColor=card_style_obj.backColor,
                        borderColor=card_style_obj.borderColor,
                        borderWidth=card_style_obj.borderWidth,
                        borderPadding=card_style_obj.borderPadding,
                        leftIndent=card_style_obj.leftIndent,
                        rightIndent=card_style_obj.rightIndent,
                        spaceBefore=card_style_obj.spaceBefore,
                        spaceAfter=card_style_obj.spaceAfter,
                    )
                    card_elements.append(Paragraph(element.text, new_style))
                else:
                    # For non-Paragraph elements (like Spacer), just add them as-is
                    card_elements.append(element)
            return card_elements
        else:
            return elements

    @staticmethod
    def create_status_indicator(status, text, style_map=None, styles_dict=None):
        """Create a status indicator with appropriate styling."""
        if style_map is None:
            style_map = {
                "success": "Success",
                "warning": "Warning",
                "error": "Error",
                "info": "Alert",
            }

        style_name = style_map.get(status.lower(), "Alert")

        if styles_dict and style_name in styles_dict:
            return Paragraph(text, styles_dict[style_name])
        else:
            # Fallback to default style
            from reportlab.lib.styles import getSampleStyleSheet

            fallback_styles = getSampleStyleSheet()
            return Paragraph(text, fallback_styles["Normal"])

    @staticmethod
    def create_data_table(
        data, headers, table_style="MODERN_TABLE_STYLE", styles_dict=None
    ):
        """Create a data table with M3 styling."""
        # Ensure data is properly formatted
        if not data:
            data = [["No data available"]]

        # Add headers if provided
        if headers:
            table_data = [headers] + data
        else:
            table_data = data

        # Create table

        table = Table(table_data)

        # Apply styling
        if styles_dict and hasattr(styles_dict, table_style):
            style = getattr(styles_dict, table_style)
            table.setStyle(style)
        else:
            # Fallback to basic table styling
            from reportlab.lib import colors
            from reportlab.platypus import TableStyle

            basic_style = TableStyle(
                [
                    ("BACKGROUND", (0, 0), (-1, 0), colors.grey),
                    ("TEXTCOLOR", (0, 0), (-1, 0), colors.whitesmoke),
                    ("ALIGN", (0, 0), (-1, -1), "CENTER"),
                    ("FONTNAME", (0, 0), (-1, 0), "Helvetica-Bold"),
                    ("FONTSIZE", (0, 0), (-1, 0), 14),
                    ("BOTTOMPADDING", (0, 0), (-1, 0), 12),
                    ("BACKGROUND", (0, 1), (-1, -1), colors.beige),
                    ("GRID", (0, 0), (-1, -1), 1, colors.black),
                ]
            )
            table.setStyle(basic_style)

        return table

    @staticmethod
    def create_compact_table(
        data, headers=None, table_style="COMPACT_TABLE_STYLE", styles_dict=None
    ):
        """Create a compact table for dense data display."""
        return PDFStyleUtils.create_data_table(data, headers, table_style, styles_dict)

    @staticmethod
    def create_summary_table(
        data, headers=None, table_style="SUMMARY_TABLE_STYLE", styles_dict=None
    ):
        """Create a summary table with highlighted styling."""
        return PDFStyleUtils.create_data_table(data, headers, table_style, styles_dict)

    @staticmethod
    def create_code_block(
        code_text, language="", style_name="CodeBlock", styles_dict=None
    ):
        """Create a code block with monospace font and background."""
        if language:
            header_text = f"```{language}"
            footer_text = "```"

            if styles_dict and "LabelMedium" in styles_dict:
                header_style = styles_dict["LabelMedium"]
                footer_style = styles_dict["LabelMedium"]
            else:
                # Fallback to default style
                from reportlab.lib.styles import getSampleStyleSheet

                fallback_styles = getSampleStyleSheet()
                header_style = fallback_styles["Normal"]
                footer_style = fallback_styles["Normal"]

            if styles_dict and style_name in styles_dict:
                code_style = styles_dict[style_name]
            else:
                # Fallback to default style
                from reportlab.lib.styles import getSampleStyleSheet

                fallback_styles = getSampleStyleSheet()
                code_style = fallback_styles["Normal"]

            elements = [
                Paragraph(header_text, header_style),
                Spacer(1, 4),
                Paragraph(code_text, code_style),
                Spacer(1, 4),
                Paragraph(footer_text, footer_style),
            ]
        else:
            if styles_dict and style_name in styles_dict:
                code_style = styles_dict[style_name]
            else:
                # Fallback to default style
                from reportlab.lib.styles import getSampleStyleSheet

                fallback_styles = getSampleStyleSheet()
                code_style = fallback_styles["Normal"]

            elements = [Paragraph(code_text, code_style)]

        return elements

    @staticmethod
    def create_bullet_list(items, style_name="Normal", styles_dict=None):
        """Create a bulleted list with consistent spacing."""
        elements = []
        for item in items:
            if isinstance(item, str):
                if styles_dict and style_name in styles_dict:
                    elements.append(Paragraph(f"• {item}", styles_dict[style_name]))
                else:
                    # Fallback to default style
                    from reportlab.lib.styles import getSampleStyleSheet

                    fallback_styles = getSampleStyleSheet()
                    elements.append(Paragraph(f"• {item}", fallback_styles["Normal"]))
            else:
                elements.append(item)
            elements.append(Spacer(1, 2))  # Small spacing between items

        return elements

    @staticmethod
    def create_numbered_list(items, style_name="Normal", styles_dict=None):
        """Create a numbered list with consistent spacing."""
        elements = []
        for i, item in enumerate(items, 1):
            if isinstance(item, str):
                if styles_dict and style_name in styles_dict:
                    elements.append(Paragraph(f"{i}. {item}", styles_dict[style_name]))
                else:
                    # Fallback to default style
                    from reportlab.lib.styles import getSampleStyleSheet

                    fallback_styles = getSampleStyleSheet()
                    elements.append(
                        Paragraph(f"{i}. {item}", fallback_styles["Normal"])
                    )
            else:
                elements.append(item)
            elements.append(Spacer(1, 2))  # Small spacing between items

        return elements

    @staticmethod
    def create_highlight_box(content, highlight_style="Success", styles_dict=None):
        """Create a highlighted box for important information."""
        if isinstance(content, str):
            if styles_dict and highlight_style in styles_dict:
                return Paragraph(content, styles_dict[highlight_style])
            else:
                # Fallback to default style
                from reportlab.lib.styles import getSampleStyleSheet

                fallback_styles = getSampleStyleSheet()
                return Paragraph(content, fallback_styles["Normal"])
        elif isinstance(content, list):
            elements = []
            for item in content:
                if isinstance(item, str):
                    if styles_dict and highlight_style in styles_dict:
                        elements.append(Paragraph(item, styles_dict[highlight_style]))
                    else:
                        # Fallback to default style
                        from reportlab.lib.styles import getSampleStyleSheet

                        fallback_styles = getSampleStyleSheet()
                        elements.append(Paragraph(item, fallback_styles["Normal"]))
                else:
                    elements.append(item)
            return elements
        else:
            return content

    @staticmethod
    def create_divider():
        """Create a visual divider between sections."""
        return Spacer(1, 16)  # Simple spacing divider

    @staticmethod
    def create_page_break():
        """Create a page break."""
        from reportlab.platypus import PageBreak

        return PageBreak()


class PDFLayoutUtils:
    """Utility class for PDF layout and positioning."""

    @staticmethod
    def calculate_column_widths(total_width, num_columns, margins=0.1):
        """Calculate column widths for tables with proper spacing."""
        available_width = total_width * (1 - margins)
        column_width = available_width / num_columns
        return [column_width] * num_columns

    @staticmethod
    def create_two_column_layout(
        left_content, right_content, left_width=0.48, right_width=0.48
    ):
        """Create a two-column layout with specified proportions."""

        # Ensure content is in list format
        if not isinstance(left_content, list):
            left_content = [left_content]
        if not isinstance(right_content, list):
            right_content = [right_content]

        # Pad shorter column with spacers
        max_height = max(len(left_content), len(right_content))

        while len(left_content) < max_height:
            left_content.append(Spacer(1, 1))
        while len(right_content) < max_height:
            right_content.append(Spacer(1, 1))

        # Create table for layout
        layout_data = [[left_content, right_content]]
        layout_table = Table(
            layout_data, colWidths=[left_width * 8.5 * inch, right_width * 8.5 * inch]
        )

        # Remove borders and spacing
        layout_style = TableStyle(
            [
                ("LEFTPADDING", (0, 0), (-1, -1), 0),
                ("RIGHTPADDING", (0, 0), (-1, -1), 0),
                ("TOPPADDING", (0, 0), (-1, -1), 0),
                ("BOTTOMPADDING", (0, 0), (-1, -1), 0),
                ("GRID", (0, 0), (-1, -1), 0, colors.white),
            ]
        )

        layout_table.setStyle(layout_style)
        return layout_table


class PDFColorUtils:
    """Utility class for PDF color management."""

    @staticmethod
    def get_m3_colors():
        """Get Material Design 3 color palette."""
        return {
            "primary": colors.HexColor("#4F9153"),
            "primary_container": colors.HexColor("#E8F5E8"),
            "on_primary": colors.HexColor("#FFFFFF"),
            "secondary": colors.HexColor("#6750A4"),
            "secondary_container": colors.HexColor("#EADDFF"),
            "on_secondary": colors.HexColor("#FFFFFF"),
            "surface": colors.HexColor("#FFFFFF"),
            "surface_variant": colors.HexColor("#F3F4F6"),
            "on_surface": colors.HexColor("#1C1B1F"),
            "on_surface_variant": colors.HexColor("#49454F"),
            "outline": colors.HexColor("#79747E"),
            "outline_variant": colors.HexColor("#CAC4D0"),
            "success": colors.HexColor("#2E7D32"),
            "warning": colors.HexColor("#F57C00"),
            "error": colors.HexColor("#C62828"),
            "info": colors.HexColor("#1976D2"),
        }

    @staticmethod
    def get_apple_hig_colors():
        """Get Apple Human Interface Guidelines color palette."""
        return {
            "system_blue": colors.HexColor("#007AFF"),
            "system_green": colors.HexColor("#34C759"),
            "system_indigo": colors.HexColor("#5856D6"),
            "system_orange": colors.HexColor("#FF9500"),
            "system_pink": colors.HexColor("#FF2D92"),
            "system_purple": colors.HexColor("#AF52DE"),
            "system_red": colors.HexColor("#FF3B30"),
            "system_teal": colors.HexColor("#5AC8FA"),
            "system_yellow": colors.HexColor("#FFCC02"),
            "system_gray": colors.HexColor("#8E8E93"),
            "system_gray2": colors.HexColor("#AEAEB2"),
            "system_gray3": colors.HexColor("#C7C7CC"),
            "system_gray4": colors.HexColor("#D1D1D6"),
            "system_gray5": colors.HexColor("#E5E5EA"),
            "system_gray6": colors.HexColor("#F2F2F7"),
        }
