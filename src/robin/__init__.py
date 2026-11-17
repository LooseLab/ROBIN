"""ROBIN and robin - A Python CLI tool with file watching capabilities."""

# Suppress pkg_resources deprecation warnings from sorted_nearest
import warnings
warnings.filterwarnings(
    "ignore", message="pkg_resources is deprecated", category=UserWarning
)
# Suppress matplotlib tight_layout warnings
warnings.filterwarnings(
    "ignore", message="The figure layout has changed to tight", category=UserWarning
)

__version__ = "0.6.0"
__author__ = "Matt Loose"
__email__ = "matt.loose@nottingham.ac.uk"

# Keep the analysis namespace lazy. Importing ``robin`` is part of every CLI
# invocation, while importing ``robin.analysis`` pulls in plotting libraries
# and other optional dependencies. Loading those modules here makes even
# ``robin --help`` build Matplotlib's font cache.
def __getattr__(name):
    if name == "analysis":
        import importlib

        return importlib.import_module("robin.analysis")
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")


__all__ = ["analysis"]
