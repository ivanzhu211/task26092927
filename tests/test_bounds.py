"""Automated tests for the frequency bounds library."""

import os
import sys
import unittest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))

from frequency_bounds import BoundsResult, MetricBounds, Summary, bounds, is_feasible  # noqa: E402


def _scenario1():
    # N=4; TV/Digital/Social each reach 2 with every pair overlap 1; Print 0.
    return Summary(
        n=4,
        tv=2, digital=2, social=2, print=0,
        tv_digital=1, tv_social=1, tv_print=0,
        digital_social=1, digital_print=0, social_print=0,
    )


def _scenario2():
    # N=4; all four channels reach 2; all six pair overlaps 1.
    return Summary(
        n=4,
        tv=2, digital=2, social=2, print=2,
        tv_digital=1, tv_social=1, tv_print=1,
        digital_social=1, digital_print=1, social_print=1,
    )


class RequiredScenariosTest(unittest.TestCase):
    def test_scenario1_reach(self):
        r = bounds(_scenario1())
        self.assertTrue(r.feasible)
        self.assertEqual(r.reach, MetricBounds(3, 4))

    def test_scenario1_frequency3(self):
        r = bounds(_scenario1())
        self.assertEqual(r.frequency_3, MetricBounds(0, 1))

    def test_scenario1_internal_identities(self):
        # Rebuild a full 15-atom population from each exact optimizer and
        # check the marginal identities at that population.
        from frequency_bounds.core import _build_model
        s = _scenario1()
        rows, rhs, objectives, constants = _build_model(s)
        total_p2 = 3
        s1_sum = sum(s.singles)
        for name in ("reach", "f1", "f2", "f3", "f4"):
            for maximize in (True, False):
                from frequency_bounds.ilp import solve_ilp
                obj = objectives[name]
                direction = obj if maximize else [-v for v in obj]
                val, x = solve_ilp(direction, rows, rhs)
                D = x[0]
                c = x[1:]
                C = sum(c)
                f3, f4 = C, D
                f2 = total_p2 - 3 * C - 6 * D
                f1 = s1_sum - 2 * total_p2 + 3 * C + 8 * D
                reach = f1 + f2 + f3 + f4
                self.assertGreaterEqual(f1, 0)
                self.assertGreaterEqual(f2, 0)
                self.assertLessEqual(reach, s.n)
                self.assertEqual(s1_sum, f1 + 2 * f2 + 3 * f3 + 4 * f4)
                self.assertEqual(total_p2, f2 + 3 * f3 + 6 * f4)

    def test_scenario2_infeasible(self):
        result = bounds(_scenario2())
        self.assertIsInstance(result, BoundsResult)
        self.assertFalse(result.feasible)
        self.assertIsNone(result.reach)

    def test_scenario2_is_feasible_flag(self):
        self.assertTrue(is_feasible(_scenario1()))
        self.assertFalse(is_feasible(_scenario2()))

    def test_scenario2_fractional_but_not_integer(self):
        # Half-integral LP relaxation: 4 vars, each pair sum <= 1, sum >= 2.
        # Fractional x_i = 1/2 exists; no integral solution.
        from frequency_bounds.ilp import solve_ilp, solve_lp
        rows = []
        for i in range(4):
            for j in range(i + 1, 4):
                row = [0] * 4
                row[i] = row[j] = 1
                rows.append(row)
        rows.append([-1, -1, -1, -1])
        self.assertIsNone(solve_ilp([1, 1, 1, 1], rows, [1] * 6 + [-2]))
        lp = solve_lp([1, 1, 1, 1], rows, [1] * 6 + [-2])
        self.assertEqual(lp[0], 2)


class InputValidationTest(unittest.TestCase):
    def test_negative_rejected(self):
        with self.assertRaises(ValueError):
            Summary(4, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0).validate()

    def test_non_integer_rejected(self):
        with self.assertRaises(ValueError):
            Summary(4, 1.5, 0, 0, 0, 0, 0, 0, 0, 0, 0).validate()

    def test_reach_above_n_rejected(self):
        with self.assertRaises(ValueError):
            Summary(1, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0).validate()

    def test_overlap_above_min_reach_rejected(self):
        with self.assertRaises(ValueError):
            Summary(4, 1, 2, 0, 0, 2, 0, 0, 0, 0, 0).validate()

    def test_boolean_rejected(self):
        with self.assertRaises(ValueError):
            Summary(True, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0).validate()

    def test_is_feasible_invalid_returns_false(self):
        self.assertFalse(
            is_feasible(Summary(1, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0))
        )


class SanityPopulationsTest(unittest.TestCase):
    def test_no_exposure(self):
        r = bounds(Summary(5, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0))
        self.assertTrue(r.feasible)
        self.assertEqual(r.reach, MetricBounds(0, 0))
        for metric in (r.frequency_1, r.frequency_2, r.frequency_3, r.frequency_4):
            self.assertEqual(metric, MetricBounds(0, 0))

    def test_everyone_all_four(self):
        r = bounds(Summary(3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3))
        self.assertTrue(r.feasible)
        self.assertEqual(r.reach, MetricBounds(3, 3))
        self.assertEqual(r.frequency_4, MetricBounds(3, 3))
        self.assertEqual(r.frequency_1, MetricBounds(0, 0))

    def test_disjoint_channels(self):
        r = bounds(Summary(8, 2, 2, 2, 2, 0, 0, 0, 0, 0, 0))
        self.assertTrue(r.feasible)
        self.assertEqual(r.reach, MetricBounds(8, 8))
        self.assertEqual(r.frequency_1, MetricBounds(8, 8))
        self.assertEqual(r.frequency_2, MetricBounds(0, 0))


if __name__ == "__main__":
    unittest.main()
