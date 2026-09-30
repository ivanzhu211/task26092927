"""Exhaustive cross-check for small N.

Every Venn population on four channels is an assignment of non-negative
integers to 16 cells (15 non-empty subsets plus the "exposed to nothing"
cell) whose sum is N.  We enumerate them exactly with a stars-and-bars
generator and compare library bounds against true extrema.
"""

import unittest
from itertools import combinations

from freq_bounds import CHANNELS, PAIRS, Summary, bounds, InfeasibleError

CELLS = [frozenset()]
for size in range(1, 5):
    CELLS.extend(frozenset(c) for c in combinations(CHANNELS, size))


def weak_compositions(total, parts):
    """Yield all non-negative integer tuples of length ``parts`` summing
    to ``total`` (stars and bars), via multiset positions."""
    for cuts in combinations_with_replacement(range(total + 1), parts - 1):
        vals = []
        prev = 0
        for cut in cuts:
            vals.append(cut - prev)
            prev = cut
        vals.append(total - prev)
        yield tuple(vals)


def combinations_with_replacement(iterable, r):
    pool = tuple(iterable)
    n = len(pool)
    if not n and r:
        return
    indices = [0] * r
    yield tuple(pool[i] for i in indices)
    while True:
        for i in reversed(range(r)):
            if indices[i] != n - 1:
                break
        else:
            return
        indices[i:] = [indices[i] + 1] * (r - i)
        yield tuple(pool[k] for k in indices)


def all_populations(n):
    for values in weak_compositions(n, len(CELLS)):
        yield dict(zip(CELLS, values))


def summary_of(cells):
    reaches = {ch: 0 for ch in CHANNELS}
    overlaps = {pair: 0 for pair in PAIRS}
    for cell, count in cells.items():
        for ch in cell:
            reaches[ch] += count
        for pair in PAIRS:
            if set(pair) <= cell:
                overlaps[pair] += count
    return reaches, overlaps


def histogram(cells):
    h = {f: 0 for f in (1, 2, 3, 4)}
    for cell, count in cells.items():
        if 1 <= len(cell) <= 4:
            h[len(cell)] += count
    return h


class BruteForceTest(unittest.TestCase):
    def test_required_infeasible_scenario(self):
        target = Summary.from_mapping({
            "n": 4,
            "reaches": {c: 2 for c in CHANNELS},
            "overlaps": {p: 1 for p in PAIRS},
        })
        with self.assertRaises(InfeasibleError):
            bounds(target)

    def test_exhaustive_bounds_n_up_to_3(self):
        for n in range(0, 4):
            groups = {}
            for cells in all_populations(n):
                reaches, overlaps = summary_of(cells)
                key = (
                    tuple(reaches[ch] for ch in CHANNELS),
                    tuple(overlaps[p] for p in PAIRS),
                )
                h = histogram(cells)
                reach_total = sum(h.values())
                groups.setdefault(key, []).append((reach_total, h))

            for (reach_vals, overlap_vals), points in groups.items():
                summary = Summary.from_mapping({
                    "n": n,
                    "reaches": dict(zip(CHANNELS, reach_vals)),
                    "overlaps": dict(zip(PAIRS, overlap_vals)),
                })
                result = bounds(summary)
                reach_values = [p[0] for p in points]
                self.assertEqual(result.reach.minimum, min(reach_values))
                self.assertEqual(result.reach.maximum, max(reach_values))
                for f in (1, 2, 3, 4):
                    vals = [p[1][f] for p in points]
                    self.assertEqual(result.frequencies[f].minimum, min(vals))
                    self.assertEqual(result.frequencies[f].maximum, max(vals))


if __name__ == "__main__":
    unittest.main()
