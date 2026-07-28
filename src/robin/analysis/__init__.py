"""Analysis package for robin.

This package consolidates all analytical modules under one namespace.
It re-exports the primary handler entry points for convenience.
"""

__all__ = []

# Important: keep imports lazy/optional.
# Several analysis handlers depend on optional third-party packages (e.g. sturgeon).
# Importing them unconditionally makes *any* import of `robin.analysis` fail, even for
# unrelated commands like `robin utils update-models`.

try:
    from .bam_preprocessor import bam_preprocessing_handler  # noqa: F401

    __all__.append("bam_preprocessing_handler")
except Exception:
    pass

try:
    from .cnv_analysis import cnv_handler  # noqa: F401

    __all__.append("cnv_handler")
except Exception:
    pass

try:
    from .mgmt_analysis import mgmt_handler  # noqa: F401
    from .mgmt_analysis import extract_mgmt_site_rows_from_bed  # noqa: F401

    __all__.extend(["mgmt_handler", "extract_mgmt_site_rows_from_bed"])
except Exception:
    pass

try:
    from .fusion_analysis import fusion_handler  # noqa: F401

    __all__.append("fusion_handler")
except Exception:
    pass

try:
    from .itd_analysis import itd_handler  # noqa: F401

    __all__.append("itd_handler")
except Exception:
    pass

try:
    from .target_analysis import target_handler  # noqa: F401

    __all__.append("target_handler")
except Exception:
    pass

try:
    from .sturgeon_analysis import sturgeon_handler  # noqa: F401

    __all__.append("sturgeon_handler")
except Exception:
    pass

try:
    from .random_forest_analysis import random_forest_handler  # noqa: F401

    __all__.append("random_forest_handler")
except Exception:
    pass

try:
    from .nanodx_analysis import nanodx_handler, pannanodx_handler  # noqa: F401

    __all__.extend(["nanodx_handler", "pannanodx_handler"])
except Exception:
    pass

try:
    from .bed_conversion import bed_conversion_handler  # noqa: F401

    __all__.append("bed_conversion_handler")
except Exception:
    pass

try:
    from .master_csv_manager import MasterCSVManager  # noqa: F401

    __all__.append("MasterCSVManager")
except Exception:
    pass

try:
    from .temp_utilities import merge_modkit_files  # noqa: F401

    __all__.append("merge_modkit_files")
except Exception:
    pass

try:
    from .methylation_wrapper import (  # noqa: F401
        figure_is_renderable,
        locus_figure,
        load_figure_pickle,
        save_figure_pickle,
        try_load_figure_pickle,
    )

    __all__.extend(
        [
            "figure_is_renderable",
            "locus_figure",
            "load_figure_pickle",
            "save_figure_pickle",
            "try_load_figure_pickle",
        ]
    )
except Exception:
    pass
