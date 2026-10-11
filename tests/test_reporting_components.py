from pathlib import Path

from reportlab.lib.pagesizes import A4
from reportlab.platypus import SimpleDocTemplate

from robin.gui import fonts
from robin.reporting.components import metadata_table, results_table
from robin.reporting.styling.styles import ReportStyles


def _styles():
    return ReportStyles(str(Path(fonts.__file__).parent))


def test_report_styles_resolve_all_component_fonts():
    styles = _styles()

    assert styles.fonts.regular
    assert styles.fonts.medium
    assert styles.fonts.bold
    assert styles.fonts.mono
    for role in (
        "ReportTitle",
        "SectionHeading",
        "SubsectionHeading",
        "TableHeader",
        "TableCell",
        "TableCellCompact",
        "Technical",
        "Warning",
    ):
        assert role in styles.styles


def test_metadata_and_results_tables_escape_values_and_use_repeat_rows(tmp_path):
    styles = _styles()
    metadata = metadata_table(
        [("Sample ID", "sample <long>&identifier")],
        styles,
        title="Sample Information",
        width=420,
    )
    results = results_table(
        ["Gene", "Finding"],
        [["A<B", "not detected & requires review"]],
        styles,
        column_proportions=(0.3, 0.7),
        width=420,
    )

    output = tmp_path / "components.pdf"
    SimpleDocTemplate(str(output), pagesize=A4).build([metadata, results])

    assert output.exists()
    assert results.repeatRows == 1
    assert "sample &lt;long&gt;&amp;identifier" in metadata._cellvalues[1][1][0].text
    assert "A&lt;B" in results._cellvalues[1][0][0].text
