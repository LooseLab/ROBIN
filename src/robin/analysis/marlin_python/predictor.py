from __future__ import annotations

import csv
import gzip
import re
import xml.etree.ElementTree as ET
from collections import OrderedDict
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Iterable, Mapping, MutableMapping, Optional, Sequence
from zipfile import ZipFile


_CG_PATTERN = re.compile(rb"cg\d{8}")
_XML_NS = {"main": "http://schemas.openxmlformats.org/spreadsheetml/2006/main"}


def extract_probe_names_from_rdata(path: str | Path) -> list[str]:
    """Extract ordered CpG probe IDs from the bundled RData file.

    The repository stores probe names as a gzipped `.RData`. The probe IDs are
    plain ASCII strings inside the payload, so we can recover the ordered list
    without depending on an R runtime.
    """

    with gzip.open(path, "rb") as handle:
        payload = handle.read()

    probes = [match.decode("ascii") for match in _CG_PATTERN.findall(payload)]
    if not probes:
        raise ValueError(f"No probe IDs were found in {path!s}.")

    return probes


def read_class_names_from_xlsx(path: str | Path) -> list[str]:
    """Read class labels from the MARLIN annotation workbook."""

    with ZipFile(path) as archive:
        shared_strings = _read_shared_strings(archive)
        sheet_xml = archive.read("xl/worksheets/sheet1.xml")

    root = ET.fromstring(sheet_xml)
    rows = []
    for row in root.findall(".//main:sheetData/main:row", _XML_NS):
        rows.append(_parse_xlsx_row(row, shared_strings))

    if not rows:
        raise ValueError(f"No rows were found in {path!s}.")

    header = rows[0]
    body = rows[1:]
    model_id_index = header.index("model_id")
    class_name_index = header.index("class_name_current")

    body.sort(key=lambda row: int(row[model_id_index]))
    return [row[class_name_index] for row in body]


def _read_shared_strings(archive: ZipFile) -> list[str]:
    xml_bytes = archive.read("xl/sharedStrings.xml")
    root = ET.fromstring(xml_bytes)
    values: list[str] = []
    for item in root.findall("main:si", _XML_NS):
        text = "".join(item.itertext())
        values.append(text)
    return values


def _parse_xlsx_row(row: ET.Element, shared_strings: Sequence[str]) -> list[str]:
    values: list[str] = []
    for cell in row.findall("main:c", _XML_NS):
        cell_type = cell.attrib.get("t")
        value_node = cell.find("main:v", _XML_NS)
        if value_node is None or value_node.text is None:
            values.append("")
            continue

        value = value_node.text
        if cell_type == "s":
            values.append(shared_strings[int(value)])
        else:
            values.append(value)
    return values


@dataclass(frozen=True)
class MARLINPrediction:
    scores: OrderedDict[str, float]
    covered_cpgs: int
    timestamp: datetime

    @property
    def top_class(self) -> str:
        return max(self.scores, key=self.scores.__getitem__)

    @property
    def top_score(self) -> float:
        return self.scores[self.top_class]

    def as_dict(self) -> dict[str, object]:
        return {
            "scores": dict(self.scores),
            "covered_cpgs": self.covered_cpgs,
            "timestamp": self.timestamp.isoformat(),
            "top_class": self.top_class,
            "top_score": self.top_score,
        }


class MARLINPredictor:
    """Python port of the R prediction preprocessing.

    Input expectations match `MARLIN_prediction.R`:
    chromosome, start, end, methylation call, probe name
    """

    def __init__(
        self,
        model,
        probe_names: Sequence[str],
        class_names: Optional[Sequence[str]] = None,
    ) -> None:
        self.model = model
        self.probe_names = list(probe_names)
        self.class_names = list(class_names) if class_names is not None else None
        self._probe_index = {probe: idx for idx, probe in enumerate(self.probe_names)}

    @classmethod
    def from_paths(
        cls,
        model_path: str | Path,
        feature_path: str | Path,
        annotation_path: Optional[str | Path] = None,
        *,
        custom_objects: Optional[Mapping[str, object]] = None,
    ) -> "MARLINPredictor":
        model = _load_tensorflow_model(model_path, custom_objects=custom_objects)
        probe_names = extract_probe_names_from_rdata(feature_path)
        class_names = (
            read_class_names_from_xlsx(annotation_path) if annotation_path is not None else None
        )
        return cls(model=model, probe_names=probe_names, class_names=class_names)

    def predict_bed(self, bed_path: str | Path) -> MARLINPrediction:
        probe_values = self._read_bed_values(bed_path)
        return self.predict_from_probe_values(probe_values)

    def predict_from_probe_values(
        self,
        probe_values: Mapping[str, float | int | str | None],
    ) -> MARLINPrediction:
        feature_vector = self.build_feature_vector(probe_values)
        return self._predict_from_feature_vector(feature_vector)

    def build_feature_vector(
        self,
        probe_values: Mapping[str, float | int | str | None],
    ) -> list[int]:
        vector = [0] * len(self.probe_names)
        for probe_name, raw_value in probe_values.items():
            index = self._probe_index.get(probe_name)
            if index is None:
                continue

            beta = _coerce_optional_float(raw_value)
            if beta is None:
                vector[index] = 0
            else:
                vector[index] = 1 if beta >= 0.5 else -1
        return vector

    def _predict_from_feature_vector(self, feature_vector: Sequence[int]) -> MARLINPrediction:
        # Keras 3 / TF 2.16+ rejects bare Python lists; use a float32 batch array.
        # Prefer model(batch, training=False) over model.predict(...): on macOS,
        # predict() can hang indefinitely in the TF data-adapter path.
        try:
            import numpy as np
        except ModuleNotFoundError as exc:  # pragma: no cover
            raise ModuleNotFoundError(
                "NumPy is required for MARLIN prediction."
            ) from exc

        batch = np.asarray([list(feature_vector)], dtype=np.float32)
        try:
            raw_predictions = self.model(batch, training=False)
        except TypeError:
            raw_predictions = self.model.predict(batch, verbose=0)
        row = _coerce_prediction_row(raw_predictions)
        class_names = self.class_names or [f"class_{idx + 1}" for idx in range(len(row))]
        scores = OrderedDict((name, float(score)) for name, score in zip(class_names, row))
        covered_cpgs = sum(1 for value in feature_vector if value != 0)
        return MARLINPrediction(
            scores=scores,
            covered_cpgs=covered_cpgs,
            timestamp=datetime.now(timezone.utc),
        )

    @staticmethod
    def _read_bed_values(path: str | Path) -> dict[str, Optional[float]]:
        probe_values: dict[str, Optional[float]] = {}
        with open(path, "r", encoding="utf-8", newline="") as handle:
            reader = csv.reader(handle, delimiter="\t")
            for row in reader:
                if not row or len(row) < 5:
                    continue
                value = _coerce_optional_float(row[3])
                probe_values[row[4]] = value
        return probe_values


class LiveMARLINPredictor:
    """Stateful wrapper for asking for a prediction at any moment in time."""

    def __init__(self, predictor: MARLINPredictor) -> None:
        self.predictor = predictor
        self._current_probe_values: MutableMapping[str, Optional[float]] = {}

    def update_from_probe_values(
        self,
        probe_values: Mapping[str, float | int | str | None],
    ) -> None:
        for probe_name, value in probe_values.items():
            self._current_probe_values[probe_name] = _coerce_optional_float(value)

    def update_from_bed(self, bed_path: str | Path) -> None:
        self.update_from_probe_values(self.predictor._read_bed_values(bed_path))

    def predict_current(self) -> MARLINPrediction:
        return self.predictor.predict_from_probe_values(self._current_probe_values)

    @property
    def covered_probe_count(self) -> int:
        return sum(
            1
            for probe_name in self.predictor.probe_names
            if probe_name in self._current_probe_values
            and self._current_probe_values[probe_name] is not None
        )


def _load_tensorflow_model(
    model_path: str | Path,
    *,
    custom_objects: Optional[Mapping[str, object]] = None,
):
    """Load the MARLIN HDF5 checkpoint.

    MARLIN was trained/saved with Keras 2. TensorFlow 2.16+ defaults to Keras 3,
    which can hang or fail on this legacy HDF5. Prefer ``tf_keras`` (Keras 2 API
    shipped for TF 2.16+) and fall back to ``tensorflow.keras``.
    """
    keras = None
    backend_name = None

    try:
        import tf_keras as keras  # type: ignore

        backend_name = "tf_keras"
    except ModuleNotFoundError:
        try:
            # Must be set before importing tensorflow.keras on TF 2.16+.
            import os

            os.environ.setdefault("TF_USE_LEGACY_KERAS", "1")
            from tensorflow import keras  # type: ignore

            backend_name = "tensorflow.keras"
        except ModuleNotFoundError as exc:
            raise ModuleNotFoundError(
                "TensorFlow is required to load a MARLIN model. "
                "Install with: pip install 'robin[marlin]' "
                "(provides tensorflow and tf-keras)."
            ) from exc

    path = str(model_path)
    load_kwargs: dict = {"custom_objects": custom_objects, "compile": False}
    try:
        model = keras.models.load_model(path, **load_kwargs)
    except TypeError:
        # Older APIs may not accept compile=
        model = keras.models.load_model(path, custom_objects=custom_objects)

    # Attach for debugging / logs without changing predict API.
    try:
        model._robin_marlin_keras_backend = backend_name  # type: ignore[attr-defined]
    except Exception:
        pass
    return model


def _coerce_optional_float(value: object) -> Optional[float]:
    if value is None:
        return None

    if isinstance(value, str):
        stripped = value.strip()
        if not stripped or stripped.upper() == "NA":
            return None
        value = stripped

    return float(value)


def _coerce_prediction_row(raw_predictions) -> list[float]:
    if hasattr(raw_predictions, "tolist"):
        raw_predictions = raw_predictions.tolist()

    if not isinstance(raw_predictions, Iterable):
        raise TypeError("Model.predict(...) did not return an iterable result.")

    raw_predictions = list(raw_predictions)
    if not raw_predictions:
        raise ValueError("Model.predict(...) returned an empty result.")

    first_row = raw_predictions[0]
    if hasattr(first_row, "tolist"):
        first_row = first_row.tolist()

    if not isinstance(first_row, Iterable):
        raise TypeError("Prediction row is not iterable.")

    return [float(value) for value in first_row]
