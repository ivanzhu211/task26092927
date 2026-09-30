"""Cross-media exposure frequency bound computation.

Public API:
    Summary:       aggregated input data
    BoundsResult:  interval result (or infeasibility)
    bounds:        compute exact integer min/max for each frequency bucket
"""

from .types import Summary, BoundsResult, MetricBounds, bounds, is_feasible

__all__ = ["Summary", "BoundsResult", "MetricBounds", "bounds", "is_feasible"]
