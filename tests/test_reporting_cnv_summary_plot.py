from __future__ import annotations

import numpy as np
import pandas as pd

from robin.reporting.plotting import create_CNV_plot


class _FakeResult:
    def __init__(self, cnv: dict) -> None:
        self.cnv = cnv


def _sample_cnv() -> dict:
    return {
        "chr1": np.array([2.0, 2.1, 1.9, 2.0]),
        "chr2": np.array([2.0, 2.2, 2.1, 2.0]),
    }


def test_create_cnv_plot_ploidy_mode_returns_jpeg() -> None:
    result = _FakeResult(_sample_cnv())
    buf = create_CNV_plot(result, {"bin_width": 1_000_000})
    data = buf.getvalue()
    assert data[:2] == b"\xff\xd8"


def test_log2_linear_axis_limits_uses_percentile_with_minimum_span() -> None:
    from robin.reporting.plotting import _log2_linear_axis_limits

    # Tight cluster should still show at least +/- 2.
    tight = np.array([0.0, 0.05, -0.05, 0.1, -0.1])
    y_min, y_max = _log2_linear_axis_limits(tight)
    assert y_min <= -2.0
    assert y_max >= 2.0

    # Outliers can expand up to the upper cap.
    wide = np.array([0.0] * 95 + [2.5, -2.5, 3.0, -3.0, 0.1])
    y_min, y_max = _log2_linear_axis_limits(wide)
    assert y_max <= 4.0 + 0.01
    assert y_min >= -4.0 - 0.01
    assert y_min < 0 < y_max


def test_report_log2_from_ploidy_matches_gui() -> None:
    """Report CNV log2 mode uses the same transform as the GUI."""
    from robin.analysis.cnv_analysis import compute_cnv_log2_from_ploidy

    cnv = {"chr5": np.array([3.5, 3.5, 3.5])}
    log2 = compute_cnv_log2_from_ploidy(cnv, "Male")
    np.testing.assert_allclose(log2["chr5"], np.log2(3.5 / 2.0), rtol=1e-6)


def test_create_cnv_plot_log2_ratio_mode_returns_jpeg() -> None:
    result = _FakeResult(_sample_cnv())
    log2_ratios = {
        "chr1": np.array([0.1, -0.2, 0.0, 0.3]),
        "chr2": np.array([-0.1, 0.4, -0.3, 0.0]),
    }
    buf = create_CNV_plot(
        result,
        {"bin_width": 1_000_000},
        normalized_cnv=log2_ratios,
        use_normalized_difference=True,
        plot_bin_width=1_000_000,
        sex_estimate="Female",
    )
    data = buf.getvalue()
    assert data[:2] == b"\xff\xd8"


def test_cnv_plot_point_state_uses_calling_thresholds() -> None:
    from robin.reporting.plotting import _cnv_plot_point_state

    assert _cnv_plot_point_state(0.5, "chr7", "Female") == "gain"
    assert _cnv_plot_point_state(-0.5, "chr7", "Female") == "loss"
    assert _cnv_plot_point_state(0.1, "chr7", "Female") == "neutral"


def test_collect_genome_significant_panel_points_offsets_by_chromosome() -> None:
    from robin.reporting.plotting import _collect_genome_significant_panel_points

    panel_genes = pd.DataFrame(
        {
            "chrom": ["chr1", "chr2"],
            "start_pos": [1_000_000, 1_000_000],
            "end_pos": [2_000_000, 2_000_000],
            "gene": ["GENE1", "GENE2"],
        }
    )
    chr1_vals = np.zeros(100, dtype=float)
    chr1_vals[1] = 1.0
    chr2_vals = np.zeros(80, dtype=float)
    chr2_vals[1] = -1.0
    cnv_source = {
        "chr1": chr1_vals,
        "chr2": chr2_vals,
    }
    chrom_offsets = {"chr1": 0.0, "chr2": 100_000_000.0}
    target_coverage = pd.DataFrame(
        {
            "chrom": ["chr1", "chr2"],
            "startpos": [1_000_000, 1_000_000],
            "endpos": [2_000_000, 2_000_000],
            "name": ["GENE1", "GENE2"],
            "length": [1_000_000, 1_000_000],
            "coverage": [30.0, 25.0],
            "bases": [30_000_000, 25_000_000],
        }
    )
    significant_regions = {}

    points = _collect_genome_significant_panel_points(
        panel_genes,
        cnv_source,
        ["chr1", "chr2"],
        chrom_offsets,
        1_000_000,
        significant_regions,
        target_coverage,
        use_log2=True,
    )

    assert len(points) == 2
    by_label = {point["label"]: point for point in points}
    assert by_label["GENE1"]["direction"] == "gain"
    assert by_label["GENE2"]["direction"] == "loss"
    assert by_label["GENE1"]["position_bp"] == 1_500_000
    assert by_label["GENE2"]["position_bp"] == 101_500_000


def test_layout_panel_coverage_point_labels_staggers_overlapping_genes() -> None:
    from robin.reporting.plotting import _layout_panel_coverage_point_labels

    points = [
        {"label": "A", "position_bp": 1_000_000.0, "coverage_val": 30.0},
        {"label": "B", "position_bp": 1_100_000.0, "coverage_val": 30.0},
    ]
    layouts = _layout_panel_coverage_point_labels(
        points,
        cov_ylim=50.0,
        x_max=250_000_000.0,
        x_key="position_bp",
        min_x_spacing=1_800_000.0,
    )
    y_a = layouts[("A", 1_000_000.0)]
    y_b = layouts[("B", 1_100_000.0)]
    assert y_b > y_a


def test_should_label_panel_gene_ignores_called_regions_without_outlier() -> None:
    from robin.reporting.plotting import _should_label_panel_gene

    point = {"cnv_val": 0.1, "mid_mb": 1.5}
    assert _should_label_panel_gene(point, mean_cnv=0.0, std_cnv=0.1) is False


def test_is_gene_cnv_outlier_uses_three_standard_deviations() -> None:
    from robin.reporting.plotting import _is_gene_cnv_outlier

    mean_cnv = 0.0
    std_cnv = 0.1
    assert _is_gene_cnv_outlier(0.25, mean_cnv, std_cnv) is False
    assert _is_gene_cnv_outlier(0.35, mean_cnv, std_cnv) is True


def test_mean_target_coverage_uses_all_targets() -> None:
    from robin.reporting.plotting import _mean_target_coverage

    df = pd.DataFrame(
        {
            "chrom": ["chr1", "chr1", "chr2"],
            "coverage": [10.0, 20.0, 30.0],
        }
    )
    assert _mean_target_coverage(df) == 20.0
    assert _mean_target_coverage(None) is None
    assert _mean_target_coverage(pd.DataFrame()) is None


def test_add_panel_coverage_points_chromosome_axis_uses_mid_mb() -> None:
    from unittest.mock import MagicMock

    from robin.reporting.plotting import _add_panel_coverage_points

    ax_cnv = MagicMock()
    ax_cov = MagicMock()
    ax_cnv.twinx.return_value = ax_cov

    panel_points = [
        {"mid_mb": 12.5, "coverage_val": 30.0, "direction": "gain", "label": "GENE1"},
    ]

    assert (
        _add_panel_coverage_points(
            ax_cnv,
            panel_points,
            120.0,
            x_key="mid_mb",
            min_x_spacing=1.8,
            mean_cov=20.0,
        )
        is True
    )
    scatter_x = ax_cov.scatter.call_args[0][0]
    assert scatter_x == [12.5]
    assert ax_cov.axhline.call_args[0][0] == 20.0


def test_add_genome_panel_coverage_points_draws_scatter_and_mean_line() -> None:
    from unittest.mock import MagicMock

    from robin.reporting.plotting import _add_genome_panel_coverage_points

    ax_cnv = MagicMock()
    ax_cov = MagicMock()
    ax_cnv.twinx.return_value = ax_cov

    # Plotted points are outliers only; all-target mean is passed separately.
    panel_points = [
        {"position_bp": 1_000_000.0, "coverage_val": 20.0, "direction": "gain", "label": "GENE1"},
        {"position_bp": 2_000_000.0, "coverage_val": 40.0, "direction": "loss", "label": "GENE2"},
        {"position_bp": 3_000_000.0, "coverage_val": 30.0, "direction": "gain", "label": "GENE3"},
    ]

    assert (
        _add_genome_panel_coverage_points(
            ax_cnv, panel_points, 250_000_000.0, mean_cov=15.0
        )
        is True
    )
    ax_cnv.twinx.assert_called_once()
    ax_cov.axhline.assert_called_once()
    assert ax_cov.axhline.call_args[0][0] == 15.0  # all-target mean, not outlier mean
    assert ax_cov.scatter.call_count == 2
    assert ax_cov.text.call_count == 3


def test_add_genome_panel_coverage_points_falls_back_to_outlier_mean() -> None:
    from unittest.mock import MagicMock

    from robin.reporting.plotting import _add_genome_panel_coverage_points

    ax_cnv = MagicMock()
    ax_cov = MagicMock()
    ax_cnv.twinx.return_value = ax_cov

    panel_points = [
        {"position_bp": 1_000_000.0, "coverage_val": 20.0, "direction": "gain", "label": "GENE1"},
        {"position_bp": 2_000_000.0, "coverage_val": 40.0, "direction": "loss", "label": "GENE2"},
        {"position_bp": 3_000_000.0, "coverage_val": 30.0, "direction": "gain", "label": "GENE3"},
    ]

    assert _add_genome_panel_coverage_points(ax_cnv, panel_points, 250_000_000.0) is True
    assert ax_cov.axhline.call_args[0][0] == 30.0  # fallback: mean of 20, 40, 30


def test_downsample_cnv_for_plot_groups_values() -> None:
    from robin.analysis.cnv_analysis import downsample_cnv_for_plot

    values = np.array([1.0, 3.0, 5.0, 7.0], dtype=float)
    x_bp, out = downsample_cnv_for_plot(values, analysis_bin_width=12_000, plot_bin_width=24_000)
    assert len(out) == 2
    assert out[0] == 2.0
    assert out[1] == 6.0
    assert x_bp[0] == 12_000
    assert x_bp[1] == 36_000


def test_downsample_cnv_chromosome_track_keeps_full_x_axis() -> None:
    from robin.analysis.cnv_analysis import downsample_cnv_chromosome_track

    analysis_bw = 12_000
    n_bins = 1000
    values = np.linspace(0.0, 1.0, n_bins)
    x_mb, out, x_max_mb = downsample_cnv_chromosome_track(
        values, analysis_bw, plot_bin_width=500_000,
    )
    assert x_max_mb == n_bins * analysis_bw / 1_000_000
    assert len(out) < n_bins
    assert len(x_mb) == len(out)
    assert float(x_mb[-1]) < x_max_mb


def test_collect_chromosome_significant_panel_points_requires_cnv_outlier() -> None:
    from robin.reporting.plotting import _collect_chromosome_significant_panel_points

    panel_genes = pd.DataFrame(
        {
            "chrom": ["chr1", "chr1"],
            "start_pos": [1_000_000, 50_000_000],
            "end_pos": [2_000_000, 51_000_000],
            "gene": ["OUTLIER", "NORMAL"],
        }
    )
    values = np.zeros(100, dtype=float)
    values[1] = 1.0
    target_coverage = pd.DataFrame(
        {
            "chrom": ["chr1", "chr1"],
            "startpos": [1_000_000, 50_000_000],
            "endpos": [2_000_000, 51_000_000],
            "name": ["OUTLIER", "NORMAL"],
            "length": [1_000_000, 1_000_000],
            "coverage": [40.0, 20.0],
            "bases": [40_000_000, 20_000_000],
        }
    )
    mean_cnv = float(np.mean(values))
    std_cnv = float(np.std(values))

    points = _collect_chromosome_significant_panel_points(
        panel_genes,
        "chr1",
        values,
        1_000_000,
        mean_cnv,
        std_cnv,
        [],
        target_coverage,
        use_log2=True,
    )

    assert len(points) == 1
    assert points[0]["label"] == "OUTLIER"


def test_add_chromosome_panel_coverage_overlay_uses_right_axis() -> None:
    from unittest.mock import MagicMock

    from robin.reporting.plotting import _add_chromosome_panel_coverage_overlay

    ax_cnv = MagicMock()
    ax_cnv.get_xlim.return_value = (0.0, 120.0)
    ax_cov = MagicMock()
    ax_cnv.twinx.return_value = ax_cov

    panel_points = [
        {"mid_mb": 12.5, "coverage_val": 30.0, "direction": "gain", "label": "GENE1"},
    ]

    assert (
        _add_chromosome_panel_coverage_overlay(
            ax_cnv,
            panel_points,
            120.0,
            target_coverage_df=pd.DataFrame({"coverage": [20.0]}),
        )
        is True
    )
    ax_cnv.twinx.assert_called_once()
    ax_cov.set_ylabel.assert_called_once()


def test_create_cnv_plot_per_chromosome_log2_mode_returns_jpeg() -> None:
    from robin.reporting.plotting import create_CNV_plot_per_chromosome

    chr1_vals = np.zeros(40, dtype=float)
    chr1_vals[5] = 1.0
    sample = {
        "chr1": chr1_vals,
        "chr2": np.linspace(0.1, -0.1, 30),
    }
    log2 = {
        "chr1": chr1_vals,
        "chr2": np.linspace(0.1, -0.1, 30),
    }
    panel_genes = pd.DataFrame(
        {
            "chrom": ["chr1"],
            "start_pos": [5_000_000],
            "end_pos": [6_000_000],
            "gene": ["TEST1"],
        }
    )
    target_coverage = pd.DataFrame(
        {
            "chrom": ["chr1"],
            "startpos": [5_000_000],
            "endpos": [6_000_000],
            "name": ["TEST1"],
            "length": [1_000_000],
            "coverage": [25.0],
            "bases": [25_000_000],
        }
    )

    class _FakeResult:
        cnv = sample

    plots = create_CNV_plot_per_chromosome(
        _FakeResult(),
        {"bin_width": 1_000_000},
        chromosomes=["chr1"],
        panel_genes_df=panel_genes,
        normalized_cnv=log2,
        target_coverage_df=target_coverage,
        use_log2_ratio=True,
        sex_estimate="Female",
    )
    assert len(plots) == 1
    assert plots[0][0] == "chr1"
    assert plots[0][1].getvalue()[:2] == b"\xff\xd8"


def test_per_chromosome_plot_uses_analysis_bin_width_by_default() -> None:
    from unittest.mock import patch

    from robin.reporting.plotting import create_CNV_plot_per_chromosome

    n_bins = 80
    chr1_vals = np.linspace(-0.2, 0.2, n_bins)

    class _FakeResult:
        cnv = {"chr1": chr1_vals}

    with patch(
        "robin.analysis.cnv_analysis.downsample_cnv_chromosome_track",
        wraps=__import__(
            "robin.analysis.cnv_analysis", fromlist=["downsample_cnv_chromosome_track"]
        ).downsample_cnv_chromosome_track,
    ) as mock_downsample:
        plots = create_CNV_plot_per_chromosome(
            _FakeResult(),
            {"bin_width": 100_000},
            chromosomes=["chr1"],
            normalized_cnv={"chr1": chr1_vals},
            use_log2_ratio=True,
            sex_estimate="Female",
        )

    assert len(plots) == 1
    mock_downsample.assert_called_once()
    assert mock_downsample.call_args.args[1] == 100_000
    assert mock_downsample.call_args.args[2] == 100_000


def test_cnv_chromosome_fig_height_fits_four_per_page() -> None:
    from robin.reporting.plotting import (
        CNV_CHROMOSOME_PLOTS_PER_PAGE,
        CNV_CHROMOSOME_PLOT_SPACER_PT,
        CNV_REPORT_FRAME_PADDING_PT,
        cnv_chromosome_fig_height_for_page,
    )

    page_height = 9.34
    plot_height = cnv_chromosome_fig_height_for_page(page_height)
    spacer_inch = (CNV_CHROMOSOME_PLOTS_PER_PAGE - 1) * CNV_CHROMOSOME_PLOT_SPACER_PT / 72.0
    frame_inch = page_height - CNV_REPORT_FRAME_PADDING_PT / 72.0
    assert (
        CNV_CHROMOSOME_PLOTS_PER_PAGE * plot_height + spacer_inch
        <= frame_inch + 1e-6
    )
    assert plot_height < 2.5


def test_twelve_chromosome_pdf_images_fit_three_pages() -> None:
    import io

    from PIL import Image as PILImage
    from reportlab.lib.pagesizes import A4
    from reportlab.lib.units import inch
    from reportlab.platypus import Image, SimpleDocTemplate, Spacer

    from robin.reporting.plotting import (
        cnv_chromosome_fig_height_for_page,
        cnv_chromosome_fig_width_for_page,
    )

    doc = SimpleDocTemplate(
        "/tmp/chrom_layout_test.pdf",
        pagesize=A4,
        rightMargin=1.0 * inch,
        leftMargin=1.0 * inch,
        topMargin=1.35 * inch,
        bottomMargin=1.0 * inch,
    )
    width_inch = cnv_chromosome_fig_width_for_page(doc.width / inch)
    height_inch = cnv_chromosome_fig_height_for_page(doc.height / inch)
    width = width_inch * inch
    height = height_inch * inch

    elements = []
    for plot_idx in range(12):
        buf = io.BytesIO()
        PILImage.new("RGB", (int(width_inch * 100), int(height_inch * 100)), "white").save(
            buf,
            format="JPEG",
        )
        buf.seek(0)
        elements.append(Image(buf, width=width, height=height))
        last_on_page = (plot_idx + 1) % 4 == 0
        if plot_idx < 11 and not last_on_page:
            elements.append(Spacer(1, 6))

    doc.build(elements)
    from PyPDF2 import PdfReader

    page_count = len(PdfReader("/tmp/chrom_layout_test.pdf").pages)
    assert page_count == 3
