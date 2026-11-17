from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence, Tuple
import json
import logging

import numpy as np
import pandas as pd

try:
    import polars as pl
except ImportError:  # pragma: no cover
    pl = None

from robin.analysis.variant_classification import classify_clinvar_significance


logger = logging.getLogger(__name__)


def _process_annotations(record: Dict[str, Any]) -> Tuple[Dict[int, Dict[str, Any]], Dict[str, Any]]:
    """
    Expand annotation information from a VCF record.

    Returns:
        Tuple of (annotation entries, record-level fields)
    """
    if "INFO" not in record:
        return {}, {}

    annotations = record["INFO"]
    significance = classify_clinvar_significance(annotations)
    rec_dict: Dict[str, Any] = {
        "is_pathogenic": significance.is_pathogenic,
        "is_clinvar_significant": significance.is_clinvar_significant,
        "is_oncogenic": significance.is_oncogenic,
        "is_vus": significance.is_vus,
        "is_somatic_significant": significance.is_somatic_significant,
    }
    ann_dict: Dict[int, Dict[str, Any]] = {}

    for ann in annotations.split(";"):
        if "=" not in ann:
            continue
        mykey, myvalue = ann.split("=", 1)
        if mykey != "ANN":
            rec_dict[mykey] = myvalue
            continue

        chunks = myvalue.split(",")
        for idx, chunk in enumerate(chunks):
            parts = chunk.split("|")
            if len(parts) < 16:
                continue
            ann_entry = {
                "Allele": parts[0],
                "Annotation": parts[1],
                "Annotation_Impact": parts[2],
                "Gene_Name": parts[3],
                "Gene_ID": parts[4],
                "Feature_Type": parts[5],
                "Feature_ID": parts[6],
                "Transcript_BioType": parts[7],
                "Rank": parts[8],
                "HGVS.c": parts[9],
                "HGVS.p": parts[10],
                "cDNA.pos / cDNA.length": parts[11],
                "CDS.pos / CDS.length": parts[12],
                "AA.pos / AA.length": parts[13],
                "Distance": parts[14],
                "ERRORS / WARNINGS / INFO": parts[15],
            }
            ann_dict[idx] = ann_entry

    return ann_dict, rec_dict


def parse_vcf(vcf_path: Path) -> Optional[pd.DataFrame]:
    """
    Parse a VCF file and return an exploded DataFrame with annotations.
    """
    try:
        if not vcf_path.exists():
            logger.warning(f"VCF file not found: {vcf_path}")
            return None

        header = "CHROM POS ID REF ALT QUAL FILTER INFO FORMAT GT".split()
        # Do not use pandas `comment="#"` here: INFO values can legally contain '#'
        # (e.g. ClinVar IDs like UniProtKB:...#VAR_...), which truncates lines.
        vcf = pd.read_csv(vcf_path, delimiter="\t", names=header, dtype=str)
        vcf = vcf[~vcf["CHROM"].str.startswith("#", na=False)]

        if len(vcf) == 0:
            logger.info(f"VCF file is empty: {vcf_path}")
            return pd.DataFrame()

        exploded_records: List[Dict[str, Any]] = []
        for record in vcf.to_dict("records"):
            annotations, record_fields = _process_annotations(record)
            if annotations:
                for ann in annotations.values():
                    exploded_records.append({**record, **ann, **record_fields})
            elif record_fields:
                exploded_records.append({**record, **record_fields})
            else:
                exploded_records.append(record)

        vcf_df = pd.DataFrame.from_records(exploded_records)
        if "INFO" in vcf_df.columns:
            vcf_df = vcf_df.drop(columns=["INFO"]).drop_duplicates()
        else:
            vcf_df = vcf_df.drop_duplicates()

        if vcf_df.empty:
            return pd.DataFrame()

        shared_columns = [
            "CHROM",
            "POS",
            "ID",
            "REF",
            "ALT",
            "QUAL",
            "FILTER",
            "FORMAT",
            "GT",
        ]
        if "Allele" in vcf_df.columns:
            shared_columns.append("Allele")

        non_shared_columns = [col for col in vcf_df.columns if col not in shared_columns]
        vcf_df = vcf_df.replace({np.nan: None})

        aggregated = (
            vcf_df.groupby(shared_columns)[non_shared_columns]
            .agg(lambda series: ", ".join(sorted({str(item) for item in series.dropna()})) or None)
            .reset_index()
        )
        return aggregated

    except Exception as exc:
        logger.error(f"Error parsing VCF file {vcf_path}: {exc}")
        logger.debug("VCF parsing failure", exc_info=True)
        return None


def build_snp_display_data(vcf_path: Path) -> Optional[Dict[str, Any]]:
    """
    Build pre-formatted SNP table data for the GUI from a VCF file.

    Returns:
        Dict with columns, rows (all and pathogenic-only), summary, and IGV regions.
    """
    vcf_df = parse_vcf(vcf_path)
    if vcf_df is None:
        return None

    columns: List[Dict[str, Any]] = []
    preferred_order = [
        "CHROM",
        "POS",
        "ID",
        "REF",
        "ALT",
        "Allele",
        "QUAL",
        "FILTER",
        "FORMAT",
        "GT",
        "Gene_Name",
        "Gene_ID",
        "Annotation",
        "Annotation_Impact",
        "Transcript_BioType",
        "Rank",
        "HGVS.c",
        "HGVS.p",
        "cDNA.pos / cDNA.length",
        "CDS.pos / CDS.length",
        "AA.pos / AA.length",
        "Distance",
        "CLNSIG",
        "CLNREVSTAT",
        "CLNDN",
        "ONC",
        "ONCDN",
        "SCI",
        "SCIDN",
        "is_pathogenic",
        "is_clinvar_significant",
        "is_vus",
    ]

    added_fields: set[str] = set()

    def add_column(field_name: str) -> None:
        if field_name in added_fields:
            return
        added_fields.add(field_name)
        label_overrides = {
            "is_clinvar_significant": "ClinVar significant",
            "is_pathogenic": "Germline pathogenic",
            "is_vus": "VUS",
            "ONCDN": "Oncogenic disease",
            "SCIDN": "Somatic disease",
            "CLNDN": "Germline disease",
            "CLNSIG": "CLNSIG (germline)",
            "ONC": "ONC (oncogenic)",
            "SCI": "SCI (somatic tier)",
        }
        columns.append(
            {
                "name": field_name,
                "label": label_overrides.get(field_name, field_name.replace("_", " ")),
                "field": field_name,
                "sortable": True,
            }
        )

    for col in preferred_order:
        if col in vcf_df.columns:
            add_column(col)

    remaining_columns = [
        col
        for col in vcf_df.columns
        if col not in added_fields and col != "INFO"
    ]
    for col in remaining_columns:
        add_column(col)

    if "is_clinvar_significant" in vcf_df.columns or "is_pathogenic" in vcf_df.columns:
        columns.append(
            {
                "name": "action",
                "label": "View in IGV",
                "field": "action",
                "sortable": False,
            }
        )

    rows_all: List[Dict[str, Any]] = []
    rows_pathogenic: List[Dict[str, Any]] = []
    rows_clinvar_significant: List[Dict[str, Any]] = []
    snp_regions_map: Dict[str, Dict[str, int]] = {}

    def _as_bool_flag(value: Any) -> bool:
        if isinstance(value, bool):
            return value
        if isinstance(value, str):
            return value.strip().upper() in {"YES", "TRUE", "1", "PATHOGENIC"}
        return False

    for _, variant in vcf_df.iterrows():
        row_dict: Dict[str, Any] = {}
        for col in columns:
            field = col["field"]
            if field == "action":
                continue
            value = variant.get(field)
            if value is None or (isinstance(value, float) and pd.isna(value)):
                row_dict[field] = ""
            elif isinstance(value, bool):
                row_dict[field] = "Yes" if value else "No"
            else:
                row_dict[field] = str(value)

        is_pathogenic_value = _as_bool_flag(variant.get("is_pathogenic", False))
        is_significant_value = _as_bool_flag(
            variant.get("is_clinvar_significant", is_pathogenic_value)
        )

        row_dict["is_pathogenic"] = "Yes" if is_pathogenic_value else "No"
        row_dict["is_clinvar_significant"] = "Yes" if is_significant_value else "No"
        row_dict["action"] = "🔍" if is_significant_value else ""
        rows_all.append(row_dict)

        if is_pathogenic_value:
            rows_pathogenic.append(dict(row_dict))

        if is_significant_value:
            rows_clinvar_significant.append(dict(row_dict))
            chrom = row_dict.get("CHROM")
            pos_str = row_dict.get("POS", "")
            try:
                pos = int(str(pos_str).replace(",", ""))
            except (TypeError, ValueError):
                continue

            if chrom:
                snp_key = f"{chrom}:{pos}"
                snp_regions_map[snp_key] = {"chrom": str(chrom), "pos": pos}

    summary = {
        "total_variants": len(rows_all),
        "pathogenic_variants": len(rows_pathogenic),
        "clinvar_significant_variants": len(rows_clinvar_significant),
    }

    return {
        "columns": columns,
        "rows_all": rows_all,
        "rows_pathogenic": rows_pathogenic,
        "rows_clinvar_significant": rows_clinvar_significant,
        "summary": summary,
        "snp_regions_map": snp_regions_map,
    }


SNP_DISPLAY_JSON = "snpsift_output_display.json"
SNP_DISPLAY_PARQUET = "snpsift_output_display.parquet"
INDEL_DISPLAY_JSON = "snpsift_indel_output_display.json"
INDEL_DISPLAY_PARQUET = "snpsift_indel_output_display.parquet"
VARIANT_SIDECAR_FORMAT = "sidecar_v1"
_PARQUET_WRITE_KWARGS = {"index": False, "engine": "pyarrow", "compression": "snappy"}
_TRUTHY_FLAGS = ("YES", "TRUE", "1", "PATHOGENIC")


def sidecar_from_display(display: Dict[str, Any], parquet_name: str) -> Dict[str, Any]:
    """Small JSON sidecar: schema/summary/IGV map, no variant rows."""
    return {
        "format": VARIANT_SIDECAR_FORMAT,
        "parquet": parquet_name,
        "columns": display.get("columns", []),
        "summary": display.get("summary", {}) or {},
        "snp_regions_map": display.get("snp_regions_map", {}) or {},
    }


def write_variant_display_bundle(
    display: Dict[str, Any],
    json_path: Path,
    parquet_path: Path,
) -> None:
    """Write Parquet rows plus a sidecar JSON that does not embed ``rows_all``."""
    rows = display.get("rows_all") or []
    fields = [
        str(col.get("field"))
        for col in display.get("columns", [])
        if isinstance(col, dict) and col.get("field")
    ]
    if rows:
        df = pd.DataFrame(rows)
    else:
        df = pd.DataFrame(columns=fields)
    if "__row_idx" in df.columns:
        df = df.drop(columns=["__row_idx"])
    df.insert(0, "__row_idx", range(len(df)))
    parquet_path.parent.mkdir(parents=True, exist_ok=True)
    df.to_parquet(parquet_path, **_PARQUET_WRITE_KWARGS)
    sidecar = sidecar_from_display(display, parquet_path.name)
    json_path.parent.mkdir(parents=True, exist_ok=True)
    with json_path.open("w", encoding="utf-8") as handle:
        json.dump(sidecar, handle)


def write_clair_variant_display_files(
    clair_dir: str | Path,
    *,
    summary_extra: Optional[Dict[str, Any]] = None,
) -> Dict[str, Optional[Path]]:
    """Rebuild SNP and INDEL display Parquet+sidecar under ``clair_dir``."""
    clair_dir = Path(clair_dir)
    written: Dict[str, Optional[Path]] = {"snp_json": None, "snp_parquet": None, "indel_json": None, "indel_parquet": None}

    def _write(vcf_name: str, json_name: str, parquet_name: str, key: str) -> None:
        vcf_path = clair_dir / vcf_name
        if not vcf_path.exists():
            logger.warning("Skipping %s display bundle; %s not found", key, vcf_path)
            return
        display = build_snp_display_data(vcf_path)
        if display is None:
            logger.warning("Could not build %s display data from %s", key, vcf_path)
            return
        if summary_extra:
            display.setdefault("summary", {}).update(summary_extra)
        json_path = clair_dir / json_name
        parquet_path = clair_dir / parquet_name
        write_variant_display_bundle(display, json_path, parquet_path)
        written[f"{key}_json"] = json_path
        written[f"{key}_parquet"] = parquet_path
        logger.info("Wrote %s variant display bundle: %s + %s", key, json_path, parquet_path)

    _write("snpsift_output.vcf", SNP_DISPLAY_JSON, SNP_DISPLAY_PARQUET, "snp")
    _write("snpsift_indel_output.vcf", INDEL_DISPLAY_JSON, INDEL_DISPLAY_PARQUET, "indel")
    return written


@dataclass
class VariantTableFilters:
    pass_only: bool = False
    significant_only: bool = False
    min_qual: Optional[float] = None
    min_dp: Optional[float] = None
    search_text: str = ""
    search_fields: Sequence[str] = field(default_factory=tuple)

    @property
    def is_active(self) -> bool:
        """True when any user filter would change the unfiltered row set."""
        return bool(
            self.pass_only
            or self.significant_only
            or self.min_qual is not None
            or self.min_dp is not None
            or str(self.search_text or "").strip()
        )


def _is_truthy_flag(value: Any) -> bool:
    if isinstance(value, bool):
        return value
    return str(value).strip().upper() in _TRUTHY_FLAGS


def _row_matches_filters(row: Dict[str, Any], filters: VariantTableFilters) -> bool:
    if filters.pass_only and str(row.get("FILTER", "")).strip().upper() != "PASS":
        return False
    if filters.significant_only and not _is_truthy_flag(
        row.get("is_clinvar_significant", row.get("is_pathogenic", ""))
    ):
        return False
    if filters.min_qual is not None:
        try:
            qual = float(str(row.get("QUAL", "")).replace(",", "").strip() or "nan")
        except (TypeError, ValueError):
            qual = float("nan")
        if qual != qual or qual < filters.min_qual:
            return False
    if filters.min_dp is not None:
        try:
            dp = float(str(row.get("DP", "")).replace(",", "").strip() or "nan")
        except (TypeError, ValueError):
            dp = float("nan")
        if dp != dp or dp < filters.min_dp:
            return False
    search = str(filters.search_text or "").strip().lower()
    if search:
        fields = filters.search_fields or [
            key for key in row.keys() if key not in {"action", "details", "__row_idx"}
        ]
        haystack = " ".join("" if row.get(field) is None else str(row.get(field)) for field in fields)
        if search not in haystack.lower():
            return False
    return True


def _apply_parquet_filters(lf: Any, filters: VariantTableFilters, column_names: Sequence[str]) -> Any:
    names = set(column_names)
    if filters.pass_only and "FILTER" in names:
        lf = lf.filter(pl.col("FILTER").cast(pl.Utf8, strict=False).str.strip_chars().str.to_uppercase() == "PASS")
    if filters.significant_only:
        flag_col = "is_clinvar_significant" if "is_clinvar_significant" in names else "is_pathogenic"
        if flag_col in names:
            lf = lf.filter(
                pl.col(flag_col)
                .cast(pl.Utf8, strict=False)
                .str.strip_chars()
                .str.to_uppercase()
                .is_in(list(_TRUTHY_FLAGS))
            )
    if filters.min_qual is not None and "QUAL" in names:
        lf = lf.filter(pl.col("QUAL").cast(pl.Float64, strict=False) >= float(filters.min_qual))
    if filters.min_dp is not None and "DP" in names:
        lf = lf.filter(pl.col("DP").cast(pl.Float64, strict=False) >= float(filters.min_dp))
    search = str(filters.search_text or "").strip().lower()
    if search:
        search_fields = [
            field
            for field in (filters.search_fields or column_names)
            if field in names and field not in {"action", "details", "__row_idx"}
        ]
        if search_fields:
            hay = pl.concat_str(
                [pl.col(field).fill_null("").cast(pl.Utf8, strict=False) for field in search_fields],
                separator=" ",
            )
            lf = lf.filter(hay.str.to_lowercase().str.contains(search, literal=True))
    return lf


class VariantDisplayStore:
    """Page variant rows from Parquet, with an in-memory JSON fallback for old samples."""

    def __init__(
        self,
        *,
        columns: List[Dict[str, Any]],
        summary: Dict[str, Any],
        regions_map: Optional[Dict[str, Dict[str, int]]] = None,
        parquet_path: Optional[Path] = None,
        rows_all: Optional[List[Dict[str, Any]]] = None,
    ) -> None:
        self.columns = columns
        self.summary = summary or {}
        self.regions_map = regions_map or {}
        self.parquet_path = Path(parquet_path) if parquet_path else None
        self._rows_all = rows_all
        self._parquet_columns: Optional[set[str]] = None

    @property
    def uses_parquet(self) -> bool:
        return self.parquet_path is not None and self.parquet_path.is_file()

    @property
    def total_variants(self) -> int:
        total = self.summary.get("total_variants")
        if total is not None:
            try:
                return int(total)
            except (TypeError, ValueError):
                pass
        if self._rows_all is not None:
            return len(self._rows_all)
        if self.uses_parquet:
            _rows, total = self.page(VariantTableFilters(), offset=0, limit=0)
            return total
        return 0

    def has_column(self, field: str) -> bool:
        if any(isinstance(col, dict) and col.get("field") == field for col in self.columns):
            return True
        return field in self._parquet_column_names()

    def _parquet_column_names(self) -> set[str]:
        if self._parquet_columns is not None:
            return self._parquet_columns
        names: set[str] = set()
        if self.uses_parquet:
            try:
                import pyarrow.parquet as pq

                names = set(pq.read_schema(self.parquet_path).names)
            except Exception:
                names = set()
        self._parquet_columns = names
        return names

    def page(
        self,
        filters: VariantTableFilters,
        *,
        offset: int = 0,
        limit: int = 100,
    ) -> Tuple[List[Dict[str, Any]], int]:
        offset = max(0, int(offset))
        limit = max(0, int(limit))
        if self.uses_parquet and pl is not None:
            return self._page_parquet(filters, offset=offset, limit=limit)
        return self._page_memory(filters, offset=offset, limit=limit)

    def fetch_row(self, row_idx: int) -> Optional[Dict[str, Any]]:
        idx = int(row_idx)
        if self.uses_parquet and pl is not None:
            try:
                lf = pl.scan_parquet(str(self.parquet_path))
                names = self._parquet_column_names()
                if "__row_idx" in names:
                    lf = lf.filter(pl.col("__row_idx") == idx)
                else:
                    lf = lf.slice(idx, 1)
                frame = lf.collect()
                if frame.height == 0:
                    return None
                row = frame.to_dicts()[0]
                row.setdefault("__row_idx", idx)
                return row
            except Exception as exc:
                logger.warning("Failed to fetch variant row %s from %s: %s", idx, self.parquet_path, exc)
                return None
        if self._rows_all is None or idx < 0 or idx >= len(self._rows_all):
            return None
        row = dict(self._rows_all[idx])
        row["__row_idx"] = idx
        return row

    def significant_rows(self, *, limit: int = 10_000) -> List[Dict[str, Any]]:
        rows, _total = self.page(
            VariantTableFilters(significant_only=True),
            offset=0,
            limit=limit,
        )
        return rows

    def _page_parquet(
        self,
        filters: VariantTableFilters,
        *,
        offset: int,
        limit: int,
    ) -> Tuple[List[Dict[str, Any]], int]:
        lf = pl.scan_parquet(str(self.parquet_path))
        lf = _apply_parquet_filters(lf, filters, list(self._parquet_column_names()))
        known_total: Optional[int] = None
        if not filters.is_active:
            raw_total = self.summary.get("total_variants")
            if raw_total is not None:
                try:
                    known_total = int(raw_total)
                except (TypeError, ValueError):
                    known_total = None
        if known_total is not None:
            total = known_total
        else:
            total = int(lf.select(pl.len()).collect().item() or 0)
        if limit <= 0 or total == 0:
            return [], total
        frame = lf.slice(offset, limit).collect()
        rows = frame.to_dicts()
        for i, row in enumerate(rows):
            row.setdefault("__row_idx", offset + i)
        return rows, total

    def _page_memory(
        self,
        filters: VariantTableFilters,
        *,
        offset: int,
        limit: int,
    ) -> Tuple[List[Dict[str, Any]], int]:
        source = self._rows_all or []
        matched: List[Dict[str, Any]] = []
        for idx, raw in enumerate(source):
            if not _row_matches_filters(raw, filters):
                continue
            row = dict(raw)
            row["__row_idx"] = idx
            matched.append(row)
        total = len(matched)
        if limit <= 0:
            return [], total
        return matched[offset : offset + limit], total

    @classmethod
    def open_json(
        cls,
        json_path: Path,
        parquet_path: Optional[Path] = None,
        *,
        vcf_fallback: Optional[Path] = None,
    ) -> Optional["VariantDisplayStore"]:
        sidecar: Dict[str, Any] = {}
        if json_path.is_file():
            try:
                with json_path.open("r", encoding="utf-8") as handle:
                    sidecar = json.load(handle)
            except Exception as exc:
                logger.error("Failed to load variant display sidecar %s: %s", json_path, exc)
                sidecar = {}

        resolved_parquet = parquet_path
        if resolved_parquet is None and sidecar.get("parquet"):
            resolved_parquet = json_path.parent / str(sidecar["parquet"])
        if resolved_parquet is not None and not resolved_parquet.is_file():
            resolved_parquet = None

        rows_all = sidecar.get("rows_all") if isinstance(sidecar.get("rows_all"), list) else None
        columns = sidecar.get("columns") if isinstance(sidecar.get("columns"), list) else []
        summary = sidecar.get("summary") if isinstance(sidecar.get("summary"), dict) else {}
        regions = sidecar.get("snp_regions_map") if isinstance(sidecar.get("snp_regions_map"), dict) else {}

        if resolved_parquet is not None:
            if not columns:
                columns = [
                    {"name": name, "label": name, "field": name, "sortable": True}
                    for name in VariantDisplayStore._schema_names(resolved_parquet)
                    if name != "__row_idx"
                ]
            if not summary:
                summary = {"total_variants": None}
            return cls(
                columns=columns,
                summary=summary,
                regions_map=regions,
                parquet_path=resolved_parquet,
            )

        if rows_all:
            if not summary:
                summary = {"total_variants": len(rows_all)}
            return cls(
                columns=columns,
                summary=summary,
                regions_map=regions,
                rows_all=rows_all,
            )

        if vcf_fallback is not None and vcf_fallback.is_file():
            display = build_snp_display_data(vcf_fallback)
            if display is None:
                return None
            return cls(
                columns=display.get("columns", []),
                summary=display.get("summary", {}),
                regions_map=display.get("snp_regions_map", {}),
                rows_all=display.get("rows_all", []),
            )
        return None

    @staticmethod
    def _schema_names(parquet_path: Path) -> List[str]:
        try:
            import pyarrow.parquet as pq

            return list(pq.read_schema(parquet_path).names)
        except Exception:
            return []

    @classmethod
    def open_snp(cls, clair_dir: str | Path) -> Optional["VariantDisplayStore"]:
        clair_dir = Path(clair_dir)
        return cls.open_json(
            clair_dir / SNP_DISPLAY_JSON,
            clair_dir / SNP_DISPLAY_PARQUET,
        )

    @classmethod
    def open_indel(
        cls, clair_dir: str | Path, *, vcf_fallback: bool = True
    ) -> Optional["VariantDisplayStore"]:
        clair_dir = Path(clair_dir)
        return cls.open_json(
            clair_dir / INDEL_DISPLAY_JSON,
            clair_dir / INDEL_DISPLAY_PARQUET,
            vcf_fallback=(clair_dir / "snpsift_indel_output.vcf") if vcf_fallback else None,
        )

