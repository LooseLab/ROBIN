            "has been specified.",
            "Different classifiers may produce different results because "
            "their training sets and methods differ.",
        ]
        explanation_text = Paragraph(
            "<br/>".join(escape(line) for line in explanation_lines),
            ParagraphStyle(
                "ClassificationExplanation",
                parent=self.styles.styles["Normal"],
                fontSize=8, leading=10,
                textColor=HexColor("#4B5563"),
            ),
        )

        if len(summary_data) > 1:
            self.summary_elements.append(
                self.create_table(summary_data, repeat_rows=1)
            )
            self.elements.append(Spacer(1, 6))

        self.export_frames["classification_summary"] = pd.DataFrame(
            summary_export,
            columns=["Classifier", "PredictedClass",
                     "CalibratedScore", "Status"],
        )

        self.elements.append(
            Paragraph("Detailed Classification Results",
                      self.styles.styles["Heading2"])
        )
        self.elements.append(Spacer(1, 4))

        block_width = self.report.doc.width / 2
        col_widths = [block_width * 0.65, block_width * 0.35]
        classifier_blocks = []

        for name, predictions in valid_results.items():
            try:
                top_predictions = predictions.head(10)
                detailed_data = [
                    ["Predicted Class", "Calibrated Score"]
                ]
                for class_name, score in top_predictions.items():
                    detailed_data.append([
                        escape(str(class_name)),
                        format_calibrated_score(score),
                    ])

                detailed_table = self.create_table(
                    detailed_data,
                    col_widths=col_widths,
                    repeat_rows=1,
                    auto_col_width=False,
                    compact=True,
                )
                detailed_table.setStyle([
                    ("TOPPADDING", (0, 0), (-1, -1), 2),
                    ("BOTTOMPADDING", (0, 0), (-1, -1), 2),
                    ("LEFTPADDING", (0, 0), (-1, -1), 4),
                    ("RIGHTPADDING", (0, 0), (-1, -1), 4),
                ])
                header = Paragraph(
                    f"{escape(name)} Classification",
                    self.styles.styles["Heading3"],
                )
                block = Table(
                    [[header], [detailed_table]],
                    colWidths=[block_width],
                )
                block.setStyle([
                    ("TOPPADDING", (0, 0), (-1, -1), 0),
                    ("BOTTOMPADDING", (0, 0), (-1, 0), 2),
                    ("BOTTOMPADDING", (0, 1), (-1, -1), 0),
                ])
                classifier_blocks.append(block)

                self.export_frames[
                    f"classification_{name.lower().replace(' ', '_')}_top10"
                ] = pd.DataFrame([
                    {
                        "Classifier": name,
                        "PredictedClass": str(class_name),
                        "CalibratedScore": float(score),
                    }
                    for class_name, score in top_predictions.items()
                ])
            except Exception:
                logger.exception("Error building detailed table for %s", name)

        for start in range(0, len(classifier_blocks), 2):
            cells = classifier_blocks[start:start + 2]
            while len(cells) < 2:
                cells.append(Spacer(block_width, 1))
            grid_table = Table(
                [cells], colWidths=[block_width, block_width]
            )
            grid_table.setStyle([
                ("VALIGN", (0, 0), (-1, -1), "TOP"),
                ("ALIGN", (0, 0), (-1, -1), "LEFT"),
                ("LEFTPADDING", (1, 0), (1, -1), 8),
            ])
            self.elements.append(grid_table)
            self.elements.append(Spacer(1, 4))

        if len(summary_data) > 1:
            self.elements.append(explanation_text)
        else:
            self.elements.append(
                Paragraph("No classification results available.",
                          self.styles.styles["Normal"])
            )
        self.summary_elements.append(explanation_text)
