from __future__ import annotations

import logging
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

import numpy as np
import pandas as pd

from robin.analysis.variant_classification import classify_clinvar_significance

logger = logging.getLogger(__name__)


def _process_annotations(
    record: Dict[str, Any],
) -> Tuple[Dict[int, Dict[str, Any]], Dict[str, Any]]:
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

        non_shared_columns = [
            col for col in vcf_df.columns if col not in shared_columns
        ]
        vcf_df = vcf_df.replace({np.nan: None})

        aggregated = (
            vcf_df.groupby(shared_columns)[non_shared_columns]
            .agg(
                lambda series: ", ".join(
                    sorted({str(item) for item in series.dropna()})
                )
                or None
            )
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
        col for col in vcf_df.columns if col not in added_fields and col != "INFO"
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
