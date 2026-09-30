import unittest
from itertools import combinations

from freq_bounds import CHANNELS, PAIRS, Summary, InfeasibleError, bounds
from freq_bounds.bounds import _build_model, _objectives, realize, _extremal_point


def make_summary(n, reaches, overlap_values):
    data = {"n": n, "reaches": dict(zip(CHANNELS, reaches)), "overlaps": {}}
    for pair, value in zip(PAIRS, overlap_values):
        data["overlaps"][pair] = value
    return Summary.from_mapping(data)


SCENARIO_1 = (4, [2, 2, 2, 0], [1, 1, 0, 1, 0, 0])
SCENARIO_2 = (4, [2, 2, 2, 2], [1, 1, 1, 1, 1, 1])


def histogram_from_cells(cells):
    hist = {k: 0 for k in range(5)}
    for cell, count in cells.items():
        hist[len(cell)] += count
    return hist


def cells_match_summary(summary, cells):
    for i, ch in enumerate(CHANNELS):
        total = sum(v for cell, v in cells.items() if ch in cell)
        if total != summary.reaches[ch]:
            return False
    for pair in PAIRS:
        total = sum(v for cell, v in cells.items() if set(pair) <= cell)
        if total != summary.overlaps[pair]:
            return False
    if sum(cells.values()) != summary.n:
        return False
    if any(v < 0 for v in cells.values()):
        return False
    return True


class BoundsCertifiedTest(unittest.TestCase):
    def assert_certified(self, summary):
        result = bounds(summary)
        reaches, overlaps, constraints = _build_model(summary)
        objectives = _objectives(reaches, overlaps)

        checks = [
            ("reach", result.reach, objectives["reach"]),
        ] + [(f, result.frequencies[f], objectives[f]) for f in (1, 2, 3, 4)]

        for label, metric, (constant, objective) in checks:
            for want_min, point in (
                (True, _extremal_point(objective, constraints, False)),
                (False, _extremal_point(objective, constraints, True)),
            ):
                cells = realize(summary, point)
                self.assertTrue(
                    cells_match_summary(summary, cells),
                    f"{label} {'min' if want_min else 'max'} "
                    f"extremum is not a valid population",
                )
                hist = histogram_from_cells(cells)
                if label == "reach":
                    value = sum(hist[k] for k in (1, 2, 3, 4))
                else:
                    value = hist[label]
                expected = metric.minimum if want_min else metric.maximum
                self.assertEqual(value, expected)
        return result

    def test_scenario_1_required(self):
        result = self.assert_certified(make_summary(*SCENARIO_1))
        self.assertEqual((result.reach.minimum, result.reach.maximum), (3, 4))
        self.assertEqual(
            (result.frequencies[3].minimum, result.frequencies[3].maximum),
            (0, 1),
        )

    def test_scenario_2_integer_infeasible(self):
        summary = make_summary(*SCENARIO_2)
        with self.assertRaises(InfeasibleError):
            bounds(summary)

    def test_no_exposure(self):
        result = self.assert_certified(make_summary(5, [0, 0, 0, 0], [0] * 6))
        self.assertEqual((result.reach.minimum, result.reach.maximum), (0, 0))
        for f in (1, 2, 3, 4):
            self.assertEqual(result.frequencies[f].minimum, 0)
            self.assertEqual(result.frequencies[f].maximum, 0)

    def test_everyone_all_four(self):
        result = self.assert_certified(make_summary(3, [3] * 4, [3] * 6))
        self.assertEqual((result.reach.minimum, result.reach.maximum), (3, 3))
        self.assertEqual((result.frequencies[4].minimum,
                          result.frequencies[4].maximum), (3, 3))
        for f in (1, 2, 3):
            self.assertEqual(result.frequencies[f].maximum, 0)

    def test_disjoint_channels(self):
        result = self.assert_certified(
            make_summary(6, [1, 2, 3, 0], [0] * 6)
        )
        self.assertEqual((result.reach.minimum, result.reach.maximum), (6, 6))
        self.assertEqual((result.frequencies[1].minimum,
                          result.frequencies[1].maximum), (6, 6))
        for f in (2, 3, 4):
            self.assertEqual(result.frequencies[f].maximum, 0)

    def test_larger_random_style(self):
        summary = make_summary(
            100, [50, 45, 40, 30],
            [20, 15, 8, 12, 6, 5],
        )
        result = self.assert_certified(summary)
        lo, hi = result.reach.minimum, result.reach.maximum
        self.assertGreaterEqual(lo, max(summary.reaches.values()))
        self.assertLessEqual(hi, summary.n)
        for f in (1, 2, 3, 4):
            mb = result.frequencies[f]
            self.assertLessEqual(mb.minimum, mb.maximum)
            self.assertGreaterEqual(mb.minimum, 0)

    def test_frequency_counts_sum_to_reach_at_extrema(self):
        summary = make_summary(
            50, [30, 20, 15, 10], [10, 8, 4, 6, 3, 2]
        )
        result = self.assert_certified(summary)
        for f in (1, 2, 3, 4):
            self.assertGreaterEqual(result.frequencies[f].minimum, 0)

    def test_locally_plausible_but_jointly_infeasible(self):
        # marginals and overlap are individually <= reaches, but their
        # inclusion-exclusion union is forced above N:
        # TV U Digital alone already needs 6 + 6 - 3 = 9 > 8 people.
        summary = make_summary(8, [6, 6, 0, 0], [3, 0, 0, 0, 0, 0])
        with self.assertRaises(InfeasibleError):
            bounds(summary)


class SummaryValidationTest(unittest.TestCase):
    def test_reach_exceeds_n(self):
        with self.assertRaises(ValueError):
            make_summary(2, [3, 0, 0, 0], [0] * 6)

    def test_overlap_exceeds_reach(self):
        with self.assertRaises(ValueError):
            make_summary(5, [1, 2, 0, 0], [2, 0, 0, 0, 0, 0])

    def test_negative_rejected(self):
        with self.assertRaises(ValueError):
            make_summary(5, [-1, 0, 0, 0], [0] * 6)

    def test_missing_pair(self):
        with self.assertRaises(ValueError):
            Summary.from_mapping({
                "n": 3,
                "reaches": {c: 0 for c in CHANNELS},
                "overlaps": {},
            })

    def test_string_pair_keys_accepted(self):
        summary = Summary.from_mapping({
            "n": 3,
            "reaches": {c: 1 for c in CHANNELS},
            "overlaps": {f"{a}|{b}": 1 for a, b in PAIRS},
        })
        result = bounds(summary)
        self.assertEqual(result.reach.minimum, 1)

    def test_wrong_type(self):
        with self.assertRaises(TypeError):
            bounds({"n": 1})


if __name__ == "__main__":
    unittest.main()
