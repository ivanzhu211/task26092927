"""Brute-force oracle cross-check for small N.

Enumerates *every* allocation of up to ``n_max`` people across the 15
non-empty channel-subset atoms (budget-pruned DFS), derives the implied
single/pair marginals, and compares the exact solver bounds with the
exhaustive population range for every distinct summary.
"""

import os
import sys
import unittest
from itertools import combinations

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))

from frequency_bounds import Summary, bounds  # noqa: E402


def _all_nonempty_subsets():
    subsets = []
    for size in range(1, 5):
        for combo in combinations(range(4), size):
            subsets.append(combo)
    return subsets


SUBSETS = _all_nonempty_subsets()


def _enumerate_populations(n_max):
    seen: dict[tuple, dict] = {}

    def record(counts):
        s1 = [0] * 4
        p2 = [[0] * 4 for _ in range(4)]
        buckets = [0] * 5
        used = sum(counts)
        for count, subset in zip(counts, SUBSETS):
            k = len(subset)
            buckets[k] += count
            for i in subset:
                s1[i] += count
            for i in subset:
                for j in subset:
                    if i < j:
                        p2[i][j] += count
        key = (
            n_max,
            s1[0], s1[1], s1[2], s1[3],
            p2[0][1], p2[0][2], p2[0][3],
            p2[1][2], p2[1][3], p2[2][3],
        )
        metrics = (used, buckets[1], buckets[2], buckets[3], buckets[4])
        entry = seen.get(key)
        if entry is None:
            seen[key] = {
                name: [v, v]
                for name, v in zip(
                    ("reach", "f1", "f2", "f3", "f4"), metrics
                )
            }
        else:
            for name, val in zip(("reach", "f1", "f2", "f3", "f4"), metrics):
                entry[name][0] = min(entry[name][0], val)
                entry[name][1] = max(entry[name][1], val)

    def dfs(idx, remaining, counts):
        if idx == len(SUBSETS):
            record(counts)
            return
        for value in range(remaining + 1):
            counts.append(value)
            dfs(idx + 1, remaining - value, counts)
            counts.pop()

    dfs(0, n_max, [])
    return seen


class BruteForceOracleTest(unittest.TestCase):
    def test_n3_exhaustive(self):
        self._run_exhaustive(3)

    @unittest.skipUnless(
        os.environ.get("RUN_SLOW_TESTS") == "1",
        "set RUN_SLOW_TESTS=1 to run the full N=4 exhaustive cross-check",
    )
    def test_n4_exhaustive(self):
        self._run_exhaustive(4)

    def _run_exhaustive(self, n_max):
        oracle = _enumerate_populations(n_max)
        self.assertGreater(len(oracle), 200)
        for key, ext in oracle.items():
            r = bounds(Summary(*key))
            self.assertTrue(r.feasible, f"solver declared feasible pop infeasible: {key}")
            for got, name in (
                (r.reach, "reach"),
                (r.frequency_1, "f1"),
                (r.frequency_2, "f2"),
                (r.frequency_3, "f3"),
                (r.frequency_4, "f4"),
            ):
                lo, hi = ext[name]
                self.assertEqual(
                    (got.minimum, got.maximum),
                    (lo, hi),
                    f"{key} metric {name}: got "
                    f"{(got.minimum, got.maximum)} expected {(lo, hi)}",
                )


if __name__ == "__main__":
    unittest.main()
