from pathlib import Path

from robin.gui.components.summary import _extract_classification_data


def test_classification_data_marks_unused_models_explicitly(tmp_path: Path) -> None:
    data = _extract_classification_data(tmp_path, active_classifiers={"sturgeon"})

    assert data["sturgeon"]["classification"] == "Not available"
    assert data["marlin"]["classification"] == "Classifier not used"
    assert data["marlin"]["confidence"] == "Classifier not used"
    assert data["marlin"]["confidence_level"] == "Classifier not used"
    assert data["marlin"]["features"] == "Classifier not used"


def test_classification_data_keeps_missing_active_model_distinct(tmp_path: Path) -> None:
    data = _extract_classification_data(tmp_path, active_classifiers={"marlin"})

    assert data["marlin"]["classification"] == "Not available"
    assert data["marlin"]["confidence"] == 0.0
    assert data["sturgeon"]["classification"] == "Classifier not used"
