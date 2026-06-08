from pathlib import Path
from unittest.mock import patch

from robin.analysis.fusion_work import (
    FusionMetadata,
    _generate_fusion_breakpoint_bed,
)
from robin.analysis.master_bed_generator import generate_master_bed


def test_master_bed_regenerates_when_same_counter_source_changes(tmp_path: Path):
    sample_id = "S1"
    bed_dir = tmp_path / sample_id / "bed_files"
    bed_dir.mkdir(parents=True)
    fusion_source = bed_dir / "fusion_breakpoints_001.bed"
    fusion_source.write_text("chr1\t100\t200\tfusion-a\n")

    master_path = generate_master_bed(
        sample_id=sample_id,
        work_dir=str(tmp_path),
        analysis_counter=1,
    )
    assert master_path is not None
    first_content = Path(master_path).read_text()
    assert "100" in first_content

    with patch(
        "robin.analysis.master_bed_generator._build_master_bed_data"
    ) as build_master_bed:
        unchanged_path = generate_master_bed(
            sample_id=sample_id,
            work_dir=str(tmp_path),
            analysis_counter=1,
        )
    assert unchanged_path == master_path
    build_master_bed.assert_not_called()

    fusion_source.write_text("chr1\t300\t400\tfusion-b\n")
    regenerated_path = generate_master_bed(
        sample_id=sample_id,
        work_dir=str(tmp_path),
        analysis_counter=1,
    )

    assert regenerated_path == master_path
    regenerated_content = Path(master_path).read_text()
    assert regenerated_content != first_content
    assert "300" in regenerated_content
    assert "100" not in regenerated_content
    assert Path(f"{master_path}.sources.json").exists()


def test_fusion_breakpoint_bed_is_not_rewritten_when_content_is_unchanged(
    tmp_path: Path,
):
    sample_id = "S1"
    sample_dir = tmp_path / sample_id
    sample_dir.mkdir(parents=True)
    (sample_dir / "cnv_analysis_counter.txt").write_text("1")
    (sample_dir / "fusion_candidates_master.csv").write_text(
        "read_id,reference_id,reference_start,reference_end,col4\n"
        "read-a,chr2,100,200,GENE1\n"
    )

    metadata = FusionMetadata(
        sample_id=sample_id,
        file_path="test.bam",
        analysis_timestamp=0.0,
    )

    assert _generate_fusion_breakpoint_bed(sample_id, metadata, str(tmp_path))
    fusion_bed = sample_dir / "bed_files" / "fusion_breakpoints_001.bed"
    first_content = fusion_bed.read_text()

    assert not _generate_fusion_breakpoint_bed(sample_id, metadata, str(tmp_path))
    assert fusion_bed.read_text() == first_content
