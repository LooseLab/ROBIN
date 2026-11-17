from types import SimpleNamespace

import pandas as pd
from reportlab.lib import colors
from reportlab.lib.styles import getSampleStyleSheet
from reportlab.platypus import Paragraph

from robin.reporting.sections.run_data import RunDataSection


def _report(master_data):
    styles = getSampleStyleSheet()
    styles.add(styles["Normal"].clone("Warning"))
    report_styles = SimpleNamespace(
        styles=styles,
        COLORS={
            "background": colors.whitesmoke,
            "primary": colors.blue,
            "text": colors.black,
            "border": colors.grey,
        },
    )
    return SimpleNamespace(
        masterdf=pd.DataFrame([master_data]),
        sample_id="sample-1",
        sample_identifiers=None,
        styles=report_styles,
    )


def _paragraph_text(elements):
    return " ".join(
        element.getPlainText()
        for element in elements
        if isinstance(element, Paragraph)
    )


def test_report_includes_all_context_modbase_warning():
    section = RunDataSection(
        _report(
            {
                "run_time": "2024-12-17T16:46:52.622000+00:00",
                "devices": "p2soloTower",
                "flowcell_ids": "PAW67899",
                "basecall_models": "hac@v6.0.0",
                "modbase_models": (
                    "dna_r10.4.1_e8.2_400bps_hac@v6.0.0_5mC_5hmC@v1"
                ),
                "counter_bam_passed": 1,
                "counter_bam_failed": 0,
                "counter_bases_count": 100,
                "counter_mapped_reads_num": 1,
                "counter_unmapped_reads_num": 0,
            }
        )
    )

    section.add_content()

    assert "incorrect and slower than expected" in _paragraph_text(section.elements)
    assert "incorrect and slower than expected" in _paragraph_text(
        section.summary_elements
    )


def test_report_warns_when_modbase_model_is_unresolved_placeholder():
    section = RunDataSection(
        _report(
            {
                "run_time": "2024-12-17T16:46:52.622000+00:00",
                "devices": "p2soloTower",
                "flowcell_ids": "PAW67899",
                "basecall_models": "hac@v6.0.0",
                "modbase_models": "modbase_model_version_id",
                "counter_bam_passed": 1,
                "counter_bam_failed": 0,
                "counter_bases_count": 100,
                "counter_mapped_reads_num": 1,
                "counter_unmapped_reads_num": 0,
            }
        )
    )

    section.add_content()

    text = _paragraph_text(section.elements)
    assert "Methylation model note" in text
    assert "does not record which modbase model was used" in text
    assert "may be incorrect" not in text


def test_report_warns_when_modbase_model_is_missing():
    section = RunDataSection(
        _report(
            {
                "run_time": "2024-12-17T16:46:52.622000+00:00",
                "devices": "p2soloTower",
                "flowcell_ids": "PAW67899",
                "basecall_models": "",
                "counter_bam_passed": 1,
                "counter_bam_failed": 0,
                "counter_bases_count": 100,
                "counter_mapped_reads_num": 1,
                "counter_unmapped_reads_num": 0,
            }
        )
    )

    section.add_content()

    assert "does not report modbase_models" in _paragraph_text(section.elements)
