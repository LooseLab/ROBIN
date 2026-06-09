import numpy as np
import pandas as pd

from robin.analysis.fusion_work import _optimize_fusion_dataframe


def _candidate_dataframe():
    return pd.DataFrame(
        {
            "col1": ["chr1", "chr2"],
            "col2": [100, 300],
            "col3": [200, 400],
            "col4": ["GENE1", "GENE2"],
            "reference_id": ["chr1", "chr2"],
            "reference_start": [110, 310],
            "reference_end": [190, 390],
            "read_id": ["read-1", "read-1"],
            "mapping_quality": [60, 30],
            "strand": ["+", "-"],
            "read_start": [0, 100],
            "read_end": [80, 180],
            "is_secondary": [False, False],
            "is_supplementary": [False, True],
            "mapping_span": [80, 80],
        }
    )


def test_fusion_dataframe_uses_native_fixed_width_dtypes():
    result = _optimize_fusion_dataframe(_candidate_dataframe())

    assert result["read_id"].dtype == object
    assert result["col4"].dtype == object
    assert result["reference_start"].dtype == np.dtype("int32")
    assert result["mapping_quality"].dtype == np.dtype("uint8")
    assert result["is_supplementary"].dtype == np.dtype("bool")
    assert not any(str(dtype).startswith("Int") for dtype in result.dtypes)


def test_fusion_dataframe_with_missing_integer_avoids_nullable_extension_dtype():
    candidates = _candidate_dataframe()
    candidates.loc[1, "reference_start"] = None

    result = _optimize_fusion_dataframe(candidates)

    assert result["reference_start"].dtype == np.dtype("float64")
    assert pd.isna(result.loc[1, "reference_start"])
    assert str(result["reference_start"].dtype) != "Int64"


def test_native_dtypes_round_trip_through_parquet(tmp_path):
    result = _optimize_fusion_dataframe(_candidate_dataframe())
    parquet_path = tmp_path / "candidates.parquet"

    result.to_parquet(parquet_path, index=False, engine="pyarrow", compression="snappy")
    restored = pd.read_parquet(parquet_path, engine="pyarrow")

    assert restored["reference_start"].tolist() == [110, 310]
    assert restored["mapping_quality"].tolist() == [60, 30]
    assert restored["read_id"].tolist() == ["read-1", "read-1"]
