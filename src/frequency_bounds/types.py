"""Input/output data types for frequency bound computation."""

from __future__ import annotations

from dataclasses import dataclass

_CHANNEL_NAMES = ("tv", "digital", "social", "print")
_PAIR_FIELDS = (
    ("tv", "digital"),
    ("tv", "social"),
    ("tv", "print"),
    ("digital", "social"),
    ("digital", "print"),
    ("social", "print"),
)


@dataclass(frozen=True)
class Summary:
    """Aggregated 4-channel exposure data.

    All values are non-negative integers. ``n`` includes the completely
    unexposed population. A pair overlap counts people covered by both of
    the two channels, including people also covered by a third/fourth one.

    Pair fields (in canonical order):
        tv_digital, tv_social, tv_print,
        digital_social, digital_print, social_print.
    """

    n: int
    tv: int
    digital: int
    social: int
    print: int
    tv_digital: int
    tv_social: int
    tv_print: int
    digital_social: int
    digital_print: int
    social_print: int

    @property
    def singles(self) -> tuple[int, int, int, int]:
        return (self.tv, self.digital, self.social, self.print)

    @property
    def pair_fields(self) -> tuple[tuple[str, str], ...]:
        return _PAIR_FIELDS

    def pair_overlap(self, i: int, j: int) -> int:
        ni, nj = _CHANNEL_NAMES[i], _CHANNEL_NAMES[j]
        if i > j:
            ni, nj = nj, ni
        return getattr(self, f"{ni}_{nj}")

    def validate(self) -> None:
        """Raise ValueError on structurally invalid input.

        Non-negativity/type and obvious marginal bounds are checked here;
        joint consistency (whether any integer population exists) is
        decided later by the exact solver in :func:`bounds`.
        """
        values = {
            "n": self.n,
            "tv": self.tv,
            "digital": self.digital,
            "social": self.social,
            "print": self.print,
            "tv_digital": self.tv_digital,
            "tv_social": self.tv_social,
            "tv_print": self.tv_print,
            "digital_social": self.digital_social,
            "digital_print": self.digital_print,
            "social_print": self.social_print,
        }
        for name, value in values.items():
            if isinstance(value, bool) or not isinstance(value, int):
                raise ValueError(f"{name} must be an integer, got {value!r}")
            if value < 0:
                raise ValueError(f"{name} must be non-negative, got {value}")

        singles = self.singles
        for count in singles:
            if count > self.n:
                raise ValueError("channel reach cannot exceed total population N")

        for ni, nj in _PAIR_FIELDS:
            overlap = getattr(self, f"{ni}_{nj}")
            cap = min(getattr(self, ni), getattr(self, nj))
            if overlap > cap:
                raise ValueError(
                    f"overlap {ni}_{nj}={overlap} exceeds {cap}=min({ni},{nj})"
                )


@dataclass(frozen=True)
class MetricBounds:
    """Attainable integer interval [minimum, maximum]."""

    minimum: int
    maximum: int


@dataclass(frozen=True)
class BoundsResult:
    """Result of :func:`bounds`.

    When ``feasible`` is False the metric fields are None: no integer
    population matches the summary.
    """

    feasible: bool
    reach: MetricBounds | None = None
    frequency_1: MetricBounds | None = None
    frequency_2: MetricBounds | None = None
    frequency_3: MetricBounds | None = None
    frequency_4: MetricBounds | None = None


def bounds(summary: Summary) -> BoundsResult:
    """Compute exact attainable integer bounds for all frequency buckets.

    Each bound is attained independently by some integer population
    consistent with every aggregate in ``summary``.
    """
    summary.validate()
    from .core import compute_bounds

    return compute_bounds(summary)


def is_feasible(summary: Summary) -> bool:
    """Return whether any integer population matches the summary."""
    try:
        summary.validate()
    except ValueError:
        return False
    from .core import compute_bounds

    return compute_bounds(summary).feasible
