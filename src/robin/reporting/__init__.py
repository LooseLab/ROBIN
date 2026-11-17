"""ROBIN report generation package.

The report implementation is intentionally loaded lazily. Importing a small
reporting constant (for example the research-use disclaimer) should not load
Matplotlib, Seaborn, and the complete PDF rendering stack.
"""

__all__ = ["create_pdf", "RobinReport"]


def __getattr__(name: str):
    if name in __all__:
        from .report import RobinReport, create_pdf

        return {"create_pdf": create_pdf, "RobinReport": RobinReport}[name]
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")
