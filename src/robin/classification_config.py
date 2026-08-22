"""
Centralized configuration for classification confidence thresholds and CNV analysis rules.

This module provides a single source of truth for confidence thresholds
and CNV classification rules used across the robin application (GUI, reporting, etc.).
"""

from typing import Any, Dict, Tuple

# Confidence thresholds for different classifiers
CLASSIFIER_CONFIDENCE_THRESHOLDS: Dict[str, Dict[str, float]] = {
    "sturgeon": {
        "high": 95.0,
        "medium": 80.0,
        "low": 0.0,
    },
    "nanodx": {
        "high": 50.0,
        "medium": 25.0,
        "low": 0.0,
    },
    "pannanodx": {
        "high": 50.0,
        "medium": 25.0,
        "low": 0.0,
    },
    "random_forest": {
        "high": 85.0,
        "medium": 65.0,
        "low": 0.0,
    },
    "marlin": {
        "high": 80.0,
        "medium": 50.0,
        "low": 0.0,
    },
    "lamprey": {
        # Same tiers as Sturgeon (research/evaluation classifier)
        "high": 95.0,
        "medium": 80.0,
        "low": 0.0,
    },
    "tucan": {
        # Paper confident-call threshold (CFT) is ≥ 0.7
        "high": 70.0,
        "medium": 50.0,
        "low": 0.0,
    },
}

# Default thresholds (fallback)
DEFAULT_CONFIDENCE_THRESHOLDS = {
    "high": 80.0,
    "medium": 50.0,
    "low": 20.0,
}


def get_confidence_level(classifier: str, confidence: float) -> str:
    """
    Get confidence level based on classifier-specific thresholds.

    Args:
        classifier: Name of the classifier (e.g., 'sturgeon', 'nanodx')
        confidence: Confidence value as percentage (0-100)

    Returns:
        Confidence level string ('High confidence', 'Medium confidence', etc.)
    """
    thresholds = CLASSIFIER_CONFIDENCE_THRESHOLDS.get(
        classifier, DEFAULT_CONFIDENCE_THRESHOLDS
    )

    if confidence >= thresholds["high"]:
        return "High confidence"
    elif confidence >= thresholds["medium"]:
        return "Medium confidence"
    elif confidence >= thresholds["low"]:
        return "Low confidence"
    else:
        return "Very low confidence"


def get_confidence_ui_tier(classifier: str, confidence: float) -> str:
    """
    Map confidence percentage (0–100) to 'high', 'medium', or 'low' for UI bars and badges.

    Uses the same thresholds as get_confidence_level (CLASSIFIER_CONFIDENCE_THRESHOLDS /
    DEFAULT_CONFIDENCE_THRESHOLDS). Bands below ``medium`` map to ``low`` (including
    “very low” verbal labels).
    """
    thresholds = CLASSIFIER_CONFIDENCE_THRESHOLDS.get(
        classifier, DEFAULT_CONFIDENCE_THRESHOLDS
    )
    if confidence >= thresholds["high"]:
        return "high"
    if confidence >= thresholds["medium"]:
        return "medium"
    return "low"


def get_confidence_status(classifier: str, confidence: float) -> tuple[str, str]:
    """
    Get confidence status and color for reporting.

    Args:
        classifier: Name of the classifier (e.g., 'sturgeon', 'nanodx')
        confidence: Confidence value as decimal (0-1)

    Returns:
        Tuple of (status, color_hex) where status is 'High', 'Medium', or 'Low'
    """
    # Convert to percentage for threshold comparison
    confidence_percent = confidence * 100

    thresholds = CLASSIFIER_CONFIDENCE_THRESHOLDS.get(
        classifier, DEFAULT_CONFIDENCE_THRESHOLDS
    )

    if confidence_percent >= thresholds["high"]:
        return "High", "#059669"  # Green
    elif confidence_percent >= thresholds["medium"]:
        return "Medium", "#D97706"  # Amber
    else:
        return "Low", "#DC2626"  # Red


# CNV Analysis Configuration
# Thresholds apply to log2(ploidy / expected copy number); 0 = normal.
# Call thresholds are set for shifts visible on the genome-wide log2 plot
# (moderate aneuploidy / subclonal). Full single-copy anchors: gain log2(3/2)≈0.58,
# loss log2(1/2)=−1.0.
CNV_THRESHOLDS = {
    "autosomes": {
        "gain": 0.3,
        "loss": -0.3,
    },
    "chrX": {
        "male": {
            "gain": 0.3,
            "loss": -0.3,
        },
        "female": {
            "gain": 0.3,
            "loss": -0.3,
        },
    },
    "chrY": {
        "male": {
            "gain": 0.3,
            "loss": -0.3,
        },
        "female": {
            "gain": 0.3,
            "loss": -0.3,
        },
    },
}

# CNV Event Detection Rules
CNV_EVENT_RULES = {
    "whole_chromosome": {
        "min_proportion_affected": 0.7,  # 70% of bins must support the event
        "min_arm_proportion": 0.4,  # Each arm must show at least 40% effect
        "single_arm_multiplier": 1.5,  # For single-arm chromosomes, use 1.5x threshold
    },
    "arm_specific": {
        "min_proportion_affected": 0.4,  # 40% of arm must be affected (reduced from 70%)
    },
    "resolution": {
        "max_bin_width": 10_000_000,  # 10Mb - resolution too low for accurate CNV calling
    },
}


def get_cnv_thresholds(chromosome: str, sex_estimate: str) -> Tuple[float, float]:
    """
    Get CNV gain/loss thresholds for a specific chromosome and sex.

    Args:
        chromosome: Chromosome name (e.g., 'chr1', 'chrX', 'chrY')
        sex_estimate: Sex estimate ('XY', 'XX', 'Male', 'Female')

    Returns:
        Tuple of (gain_threshold, loss_threshold)
    """
    # Normalize sex estimate
    sex_key = "male" if sex_estimate.upper() in ("XY", "MALE") else "female"

    if chromosome == "chrX":
        thresholds = CNV_THRESHOLDS["chrX"][sex_key]
    elif chromosome == "chrY":
        thresholds = CNV_THRESHOLDS["chrY"][sex_key]
    else:
        # Autosomes
        thresholds = CNV_THRESHOLDS["autosomes"]

    return thresholds["gain"], thresholds["loss"]


def is_whole_chromosome_event(
    p_arm_mean: float,
    q_arm_mean: float,
    p_arm_proportion_gain: float,
    p_arm_proportion_loss: float,
    q_arm_proportion_gain: float,
    q_arm_proportion_loss: float,
    gain_threshold: float,
    loss_threshold: float,
) -> Tuple[bool, str]:
    """
    Determine if a chromosome shows a whole chromosome event.

    ``p_arm_mean`` / ``q_arm_mean`` are arm-level location statistics (median
    of finite bins). Proportion checks are direction-specific (gain bins vs
    loss bins) so a uniform log2 shift across an arm counts toward the call.
    """
    rules = CNV_EVENT_RULES["whole_chromosome"]

    both_arms_gained = (
        p_arm_mean > gain_threshold
        and q_arm_mean > gain_threshold
        and (
            p_arm_proportion_gain > rules["min_proportion_affected"]
            or q_arm_proportion_gain > rules["min_proportion_affected"]
        )
        and (
            p_arm_proportion_gain > rules["min_arm_proportion"]
            and q_arm_proportion_gain > rules["min_arm_proportion"]
        )
    )

    both_arms_lost = (
        p_arm_mean < loss_threshold
        and q_arm_mean < loss_threshold
        and (
            p_arm_proportion_loss > rules["min_proportion_affected"]
            or q_arm_proportion_loss > rules["min_proportion_affected"]
        )
        and (
            p_arm_proportion_loss > rules["min_arm_proportion"]
            and q_arm_proportion_loss > rules["min_arm_proportion"]
        )
    )

    if both_arms_gained:
        return True, "GAIN"
    if both_arms_lost:
        return True, "LOSS"
    return False, "NORMAL"


def is_arm_event(
    arm_mean: float,
    arm_proportion_gain: float,
    arm_proportion_loss: float,
    gain_threshold: float,
    loss_threshold: float,
) -> Tuple[bool, str]:
    """
    Determine if a chromosome arm shows a significant event.

    ``arm_mean`` is the arm-level location statistic (median of finite bins).
    """
    rules = CNV_EVENT_RULES["arm_specific"]

    if (
        arm_mean > gain_threshold
        and arm_proportion_gain > rules["min_proportion_affected"]
    ):
        return True, "GAIN"
    if (
        arm_mean < loss_threshold
        and arm_proportion_loss > rules["min_proportion_affected"]
    ):
        return True, "LOSS"

    return False, "NORMAL"


def is_resolution_sufficient(bin_width: int) -> bool:
    """
    Check if the resolution is sufficient for CNV calling.

    Args:
        bin_width: Bin width in base pairs

    Returns:
        True if resolution is sufficient for CNV calling
    """
    return bin_width <= CNV_EVENT_RULES["resolution"]["max_bin_width"]


# Fusion Detection Configuration
FUSION_DETECTION_THRESHOLDS = {
    "mapping_quality": {
        "min_threshold": 49,  # Minimum mapping quality score
        "description": "Minimum mapping quality to consider read alignment reliable",
    },
    "mapping_span": {
        "min_threshold": 249,  # Minimum mapping span in bases
        "description": "Minimum read mapping length for reliable fusion detection",
    },
    "gene_overlap": {
        "min_threshold": 99,  # Minimum overlap with gene region
        "description": "Minimum overlap between read and gene region (0 = any overlap)",
    },
    "read_support": {
        "min_threshold": 3,  # Minimum supporting reads per gene pair
        "description": "Minimum number of supporting reads required for reliable fusion detection",
    },
    "read_coordinate_overlap": {
        "max_threshold": 100,  # Maximum allowed overlap between read alignments
        "description": "Maximum allowed overlap between read alignments to prevent false positives",
    },
    "coordinate_similarity": {
        "max_threshold": 50,  # Maximum allowed coordinate difference for similar alignments
        "description": "Maximum coordinate difference (start or end) to consider alignments as similar and filter duplicates",
    },
}

# Fusion Detection Rules
FUSION_DETECTION_RULES = {
    "supplementary_alignment": {
        "required": True,
        "description": "Only process reads with supplementary alignments (SA tag)",
    },
    "multi_gene_mapping": {
        "required": True,
        "description": "Reads must map to more than 1 gene to be considered fusion candidates",
    },
    "overlapping_gene_filter": {
        "enabled": True,
        "description": "Filter out reads where same genomic alignment is annotated with multiple overlapping genes",
    },
    "exact_coordinate_matching": {
        "enabled": True,
        "description": "Use exact genomic coordinate matching to identify identical alignments",
    },
    "coordinate_similarity_filter": {
        "enabled": True,
        "description": "Filter out reads with very similar (but not identical) alignments to reduce mapping artifacts",
    },
}


def get_fusion_threshold(threshold_name: str) -> int:
    """
    Get fusion detection threshold by name.

    Args:
        threshold_name: Name of the threshold ('mapping_quality', 'mapping_span', etc.)

    Returns:
        Threshold value as integer
    """
    if threshold_name not in FUSION_DETECTION_THRESHOLDS:
        raise ValueError(f"Unknown fusion threshold: {threshold_name}")

    # Handle both min_threshold and max_threshold keys
    threshold_config = FUSION_DETECTION_THRESHOLDS[threshold_name]
    if "min_threshold" in threshold_config:
        return threshold_config["min_threshold"]
    elif "max_threshold" in threshold_config:
        return threshold_config["max_threshold"]
    else:
        raise ValueError(
            f"Threshold {threshold_name} has no min_threshold or max_threshold defined"
        )


def get_fusion_rule(rule_name: str) -> bool:
    """
    Get fusion detection rule setting by name.

    Args:
        rule_name: Name of the rule ('supplementary_alignment', 'multi_gene_mapping', etc.)

    Returns:
        Rule setting as boolean
    """
    if rule_name not in FUSION_DETECTION_RULES:
        raise ValueError(f"Unknown fusion rule: {rule_name}")

    return (
        FUSION_DETECTION_RULES[rule_name]["required"]
        if "required" in FUSION_DETECTION_RULES[rule_name]
        else FUSION_DETECTION_RULES[rule_name]["enabled"]
    )


def validate_fusion_candidate(
    mapping_quality: int,
    mapping_span: int,
    gene_overlap: int,
    supporting_reads: int,
    has_supplementary: bool = True,
    maps_multiple_genes: bool = True,
) -> Tuple[bool, str]:
    """
    Validate a fusion candidate against all detection rules and thresholds.

    Args:
        mapping_quality: Read mapping quality score
        mapping_span: Read mapping span in bases
        gene_overlap: Overlap with gene region in bases
        supporting_reads: Number of supporting reads for the gene pair
        has_supplementary: Whether read has supplementary alignments
        maps_multiple_genes: Whether read maps to multiple genes

    Returns:
        Tuple of (is_valid, reason) where reason explains why validation failed
    """
    # Check supplementary alignment requirement
    if get_fusion_rule("supplementary_alignment") and not has_supplementary:
        return False, "Missing supplementary alignment (SA tag)"

    # Check multi-gene mapping requirement
    if get_fusion_rule("multi_gene_mapping") and not maps_multiple_genes:
        return False, "Read does not map to multiple genes"

    # Check mapping quality threshold
    min_mq = get_fusion_threshold("mapping_quality")
    if mapping_quality <= min_mq:
        return False, f"Mapping quality {mapping_quality} below threshold {min_mq}"

    # Check mapping span threshold
    min_span = get_fusion_threshold("mapping_span")
    if mapping_span <= min_span:
        return False, f"Mapping span {mapping_span} below threshold {min_span}"

    # Check gene overlap threshold
    min_overlap = get_fusion_threshold("gene_overlap")
    if gene_overlap < min_overlap:
        return False, f"Gene overlap {gene_overlap} below threshold {min_overlap}"

    # Check read support threshold
    min_support = get_fusion_threshold("read_support")
    if supporting_reads < min_support:
        return (
            False,
            f"Supporting reads {supporting_reads} below threshold {min_support}",
        )

    return True, "All validation criteria passed"


def are_coordinates_similar(
    coord1_start: int,
    coord1_end: int,
    coord2_start: int,
    coord2_end: int,
    max_difference: int = None,
) -> bool:
    """
    Check if two coordinate ranges are similar (within threshold).

    Args:
        coord1_start: Start of first coordinate range
        coord1_end: End of first coordinate range
        coord2_start: Start of second coordinate range
        coord2_end: End of second coordinate range
        max_difference: Maximum allowed difference in coordinates (uses config if None)

    Returns:
        True if coordinates are similar (within threshold)
    """
    if max_difference is None:
        max_difference = FUSION_DETECTION_THRESHOLDS["coordinate_similarity"][
            "max_threshold"
        ]

    # Check if start and end coordinates are within the threshold
    start_diff = abs(coord1_start - coord2_start)
    end_diff = abs(coord1_end - coord2_end)

    return start_diff <= max_difference and end_diff <= max_difference


def get_fusion_config_summary() -> Dict[str, Any]:
    """
    Get a summary of all fusion detection configuration.

    Returns:
        Dictionary containing all thresholds and rules
    """
    return {
        "thresholds": FUSION_DETECTION_THRESHOLDS,
        "rules": FUSION_DETECTION_RULES,
        "description": "Centralized fusion detection configuration for consistent analysis across the application",
    }
