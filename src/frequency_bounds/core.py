"""Aggregate reduction and exact bound computation.

Every 4-channel population is fully described, for frequency purposes, by
the counts of the 15 non-empty channel subsets.  This module collapses
those 15 counts into five non-negative integer aggregate variables:

    D   = people covered by all four channels
    c_i = people covered by the triple-coverage atom that *misses*
          channel i (there is one such atom for each i)

All remaining atoms (single- and double-coverage people) are then linear
functions of (D, c_0..c_3) determined by the input marginals:

    a_i = people covered only by channel i
        = S1_i - P2_i + P3_i - P4
          - C + 2*D + c_i
    b_ij = people covered exactly by {i, j}
         = p_ij - (c_k + c_l) - D     (k, l are the other two channels)

so a feasible population exists iff these ten numbers are non-negative
integers, the five variables are non-negative integers, and the total
does not exceed N.

Variables in the ILP are ordered  x = (D, c_0, c_1, c_2, c_3).
"""

from __future__ import annotations

from itertools import combinations

from .ilp import solve_ilp
from .types import BoundsResult, MetricBounds, Summary

_NVARS = 5
_D = 0
_C = (1, 2, 3, 4)


def _complements(i: int, j: int) -> tuple[int, int]:
    return tuple(k for k in range(4) if k != i and k != j)  # type: ignore[return-value]


def _build_model(summary: Summary):
    s1 = list(summary.singles)
    p2 = [[0] * 4 for _ in range(4)]
    for i, j in combinations(range(4), 2):
        p2[i][j] = p2[j][i] = summary.pair_overlap(i, j)

    p3_i = [
        sum(p2[i][j] for j in range(4) if j != i) for i in range(4)
    ]

    rows: list[list[int]] = []
    rhs: list[int] = []

    # (1) b_ij = p_ij - c_k - c_l - D >= 0  ->  D + c_k + c_l <= p_ij.
    for i, j in combinations(range(4), 2):
        k, l = _complements(i, j)
        row = [0] * _NVARS
        row[_D] = 1
        row[_C[k]] = 1
        row[_C[l]] = 1
        rows.append(row)
        rhs.append(p2[i][j])

    # (2) a_i = s_i - p3_i + C - c_i + 2D >= 0
    #     ->  -C + c_i - 2D <= s_i - p3_i.
    for i in range(4):
        row = [0] * _NVARS
        row[_D] = -2
        for t in range(4):
            row[_C[t]] = -1
        row[_C[i]] = 0
        rows.append(row)
        rhs.append(s1[i] - p3_i[i])

    # (3) Population cap: R + U = N,  R = S1 - P2 + C + 3D, U >= 0
    #     ->  C + 3D <= N - S1 + P2.
    total_p2 = sum(p2[i][j] for i, j in combinations(range(4), 2))
    row = [0] * _NVARS
    row[_D] = 3
    for t in range(4):
        row[_C[t]] = 1
    rows.append(row)
    rhs.append(summary.n - sum(s1) + total_p2)

    # Simple variable upper bounds (keep every LP bounded and prune search).
    bound = summary.n
    for col in range(_NVARS):
        row = [0] * _NVARS
        row[col] = 1
        rows.append(row)
        rhs.append(bound)

    # Objective coefficients for (D, c0, c1, c2, c3).
    # Reach R   = S1 - P2 + C + 3D
    obj_reach = [3] + [1] * 4
    # F4 = D
    obj_f4 = [1, 0, 0, 0, 0]
    # F3 = C
    obj_f3 = [0, 1, 1, 1, 1]
    # F2 = P2 - 3C - 6D
    obj_f2 = [-6, -3, -3, -3, -3]
    # F1 = S1 - 2P2 + 3C + 8D
    obj_f1 = [8, 3, 3, 3, 3]

    objectives = {
        "reach": obj_reach,
        "f1": obj_f1,
        "f2": obj_f2,
        "f3": obj_f3,
        "f4": obj_f4,
    }
    constants = {
        "reach": sum(s1) - total_p2,
        "f1": sum(s1) - 2 * total_p2,
        "f2": total_p2,
        "f3": 0,
        "f4": 0,
    }
    return rows, rhs, objectives, constants


def _extremum(rows, rhs, objective, constant, maximize: bool):
    direction = objective if maximize else [-v for v in objective]
    result = solve_ilp(direction, rows, rhs)
    if result is None:
        return None
    value, _x = result
    signed = value if maximize else -value
    return constant + signed


def compute_bounds(summary: Summary) -> BoundsResult:
    rows, rhs, objectives, constants = _build_model(summary)

    found: dict[str, int] = {}
    for name, obj in objectives.items():
        hi = _extremum(rows, rhs, obj, constants[name], True)
        if hi is None:
            return BoundsResult(feasible=False)
        found[name + "_max"] = hi
        lo = _extremum(rows, rhs, obj, constants[name], False)
        assert lo is not None  # feasibility already established
        found[name + "_min"] = lo

    metrics = {
        "reach": MetricBounds(found["reach_min"], found["reach_max"]),
        "f1": MetricBounds(found["f1_min"], found["f1_max"]),
        "f2": MetricBounds(found["f2_min"], found["f2_max"]),
        "f3": MetricBounds(found["f3_min"], found["f3_max"]),
        "f4": MetricBounds(found["f4_min"], found["f4_max"]),
    }
    for mb in metrics.values():
        assert 0 <= mb.minimum <= mb.maximum <= summary.n, (name, mb)

    return BoundsResult(
        feasible=True,
        reach=metrics["reach"],
        frequency_1=metrics["f1"],
        frequency_2=metrics["f2"],
        frequency_3=metrics["f3"],
        frequency_4=metrics["f4"],
    )
