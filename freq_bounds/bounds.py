"""Exact frequency-histogram bounds from aggregate channel summaries.

Only aggregate data is consumed: total population ``N``, four marginal
reaches, and six pairwise overlaps.  Pairwise overlaps alone do not determine
a unique histogram, so for every metric we return the exact minimum and
maximum attainable by an integer population consistent with *all* inputs.

Integer Venn parameterization
-----------------------------
For four channels, define five non-negative integer variables:

* ``q``      -- people covered by all four channels (frequency 4);
* ``y_i``    -- people covered by the other three channels but not channel
  ``i`` (triple-exclusive cell; each contributes frequency 3).

The remaining ten Venn cells are then fixed by the summaries:

* pair-exclusive cell for pair (i, j)::

      e_ij = overlap(i,j) - q - y_a - y_b

  where {a, b} are the two channels other than i and j (non-negativity gives
  ``q + y_a + y_b <= overlap(i,j)``);

* single cell of channel ``i``::

      s_i = reach_i - sum(3 incident overlaps)
            + (y_j + y_k + y_l) + 2*q

  where {j, k, l} are the other channels (non-negativity gives
  ``-(y_j + y_k + y_l) - 2*q <= reach_i - sum(incident overlaps)``).

Feasibility also requires the union not to exceed ``N``::

      U = sum(reaches) - sum(overlaps) + (y_0+...+y_3) + 3*q <= N

These inequalities are necessary and sufficient: any non-negative integer
(q, y_0..y_3) satisfying them makes all 15 Venn cells non-negative integers
and hence defines a real population, and every real population induces such a
point.  Exact integer optimization therefore yields true extrema, not merely
conservative bounds.
"""

from __future__ import annotations

from dataclasses import dataclass

from .errors import InfeasibleError
from .solver import solve_ilp
from .summary import CHANNELS, PAIRS, Summary

_Q = 4  # variable index of q


@dataclass(frozen=True)
class MetricBounds:
    minimum: int
    maximum: int


@dataclass(frozen=True)
class BoundsResult:
    reach: MetricBounds
    frequencies: dict[int, MetricBounds]


def _omitted_channels(pair_index: int) -> tuple[int, int]:
    """Indices of the two channels not belonging to the given pair."""
    a_name, b_name = PAIRS[pair_index]
    return tuple(k for k, ch in enumerate(CHANNELS) if ch not in (a_name, b_name))


def _build_model(summary: Summary):
    reaches = [summary.reaches[ch] for ch in CHANNELS]
    overlaps = [summary.overlaps[pair] for pair in PAIRS]

    incident = [[] for _ in range(4)]
    for p, (a_name, b_name) in enumerate(PAIRS):
        incident[CHANNELS.index(a_name)].append(p)
        incident[CHANNELS.index(b_name)].append(p)

    constraints = []

    # Pair-exclusive cells non-negative: q + y_a + y_b <= overlap(p).
    for p, count in enumerate(overlaps):
        a, b = _omitted_channels(p)
        row = [0] * 5
        row[a] = row[b] = row[_Q] = 1
        constraints.append((row, count))

    # Single cells non-negative:
    #   -(y_j+y_k+y_l) - 2*q <= reach_i - sum(incident overlaps).
    for i in range(4):
        row = [0] * 5
        for j in range(4):
            if j != i:
                row[j] = -1
        row[_Q] = -2
        rhs = reaches[i] - sum(overlaps[p] for p in incident[i])
        constraints.append((row, rhs))

    # Union <= N:
    #   S + 3*q <= N - sum(reaches) + sum(overlaps)
    constraints.append(([1, 1, 1, 1, 3],
                        summary.n - sum(reaches) + sum(overlaps)))

    return reaches, overlaps, constraints


def _objectives(reaches, overlaps):
    """Return ``(constant, coefficients)`` for each metric.

    Derived from the Venn parameterization:

    * reach  = sum(r) - sum(o) + S + 3*q
    * f4     = q
    * f3     = S
    * f2     = sum(o) - 3*S - 6*q
    * f1     = sum(r) - 2*sum(o) + 3*S + 8*q
    """
    sum_r, sum_o = sum(reaches), sum(overlaps)
    return {
        "reach": (sum_r - sum_o, [1, 1, 1, 1, 3]),
        4: (0, [0, 0, 0, 0, 1]),
        3: (0, [1, 1, 1, 1, 0]),
        2: (sum_o, [-3, -3, -3, -3, -6]),
        1: (sum_r - 2 * sum_o, [3, 3, 3, 3, 8]),
    }


def bounds(summary: Summary) -> BoundsResult:
    """Compute exact bounds on deduplicated reach and frequency counts.

    Raises :class:`InfeasibleError` when no integer population is consistent
    with the supplied summaries.
    """
    if not isinstance(summary, Summary):
        raise TypeError("summary must be a Summary; use Summary.from_mapping()")

    reaches, overlaps, constraints = _build_model(summary)

    feasibility = solve_ilp([0] * 5, constraints)
    if feasibility is None:
        raise InfeasibleError(
            "no integer population is consistent with the supplied summaries"
        )

    results = {}
    for key, (constant, objective) in _objectives(reaches, overlaps).items():
        _, x_lo = solve_ilp(objective, constraints)
        _, x_hi = solve_ilp(objective, constraints, maximize=True)
        lo = constant + sum(objective[j] * x_lo[j] for j in range(5))
        hi = constant + sum(objective[j] * x_hi[j] for j in range(5))
        results[key] = MetricBounds(minimum=lo, maximum=hi)

    return BoundsResult(
        reach=results["reach"],
        frequencies={f: results[f] for f in (1, 2, 3, 4)},
    )


def _extremal_point(objective, constraints, maximize: bool):
    _, point = solve_ilp(objective, constraints, maximize=maximize)
    return point


def realize(summary: Summary, point) -> dict:
    """Expand a feasible ``(y0..y3, q)`` point into the 15 Venn cells.

    Keys are frozensets of channels; values are integer people counts.
    ``frozenset()`` represents people covered by no channel.
    """
    reaches = [summary.reaches[ch] for ch in CHANNELS]
    overlaps = [summary.overlaps[pair] for pair in PAIRS]
    y = point[:4]
    q = point[4]

    cells = {}
    cells[frozenset(CHANNELS)] = q
    for i, ch in enumerate(CHANNELS):
        cells[frozenset(ch for k, ch in enumerate(CHANNELS) if k != i)] = y[i]

    incident = [[] for _ in range(4)]
    for p, (a_name, b_name) in enumerate(PAIRS):
        incident[CHANNELS.index(a_name)].append(p)
        incident[CHANNELS.index(b_name)].append(p)

    for p, pair in enumerate(PAIRS):
        a, b = _omitted_channels(p)
        cells[frozenset(pair)] = overlaps[p] - q - y[a] - y[b]
    for i, ch in enumerate(CHANNELS):
        others_sum = sum(y[j] for j in range(4) if j != i)
        cells[frozenset({ch})] = (
            reaches[i]
            - sum(overlaps[k] for k in incident[i])
            + others_sum
            + 2 * q
        )

    covered = sum(cells.values())
    cells[frozenset()] = summary.n - covered
    return cells
