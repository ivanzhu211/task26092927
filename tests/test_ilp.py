import os
import sys
import unittest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))

from fractions import Fraction  # noqa: E402

from frequency_bounds.ilp import (  # noqa: E402
    INFEASIBLE,
    UNBOUNDED,
    solve_ilp,
    solve_lp,
)


class TestLP(unittest.TestCase):
    def test_simple_max(self):
        val, x = solve_lp([5, 4], [[6, 4], [1, 2]], [24, 6])
        self.assertEqual(val, 21)
        self.assertEqual(x, [3, Fraction(3, 2)])

    def test_infeasible(self):
        self.assertIs(solve_lp([1], [[1], [-1]], [1, -2]), INFEASIBLE)

    def test_unbounded(self):
        self.assertIs(solve_lp([1], [], []), UNBOUNDED)

    def test_negative_rhs(self):
        val, x = solve_lp([1], [[-1], [1]], [-3, 5])
        self.assertEqual(val, 5)
        self.assertEqual(x, [5])

    def test_redundant_rows(self):
        val, x = solve_lp([1, 1], [[1, 1], [-1, -1]], [4, -1])
        self.assertEqual(val, 4)

    def test_fractional_gap(self):
        # K4 vertex packing with all pair caps 1: LP 2, integer 1.
        rows = []
        for i in range(4):
            for j in range(i + 1, 4):
                row = [0] * 4
                row[i] = row[j] = 1
                rows.append(row)
        val, x = solve_lp([1] * 4, rows, [1] * 6)
        self.assertEqual(val, 2)


class TestILP(unittest.TestCase):
    def test_integer_gap(self):
        rows = []
        for i in range(4):
            for j in range(i + 1, 4):
                row = [0] * 4
                row[i] = row[j] = 1
                rows.append(row)
        result = solve_ilp([1] * 4, rows, [1] * 6)
        self.assertEqual(result[0], 1)

    def test_infeasible(self):
        self.assertIsNone(solve_ilp([1], [[1], [-1]], [1, -2]))

    def test_simple(self):
        val, x = solve_ilp([5, 4], [[6, 4], [1, 2]], [24, 6])
        self.assertEqual(val, 20)
        self.assertEqual(x, [4, 0])

    def test_feasible_with_lower_bounds(self):
        val, x = solve_ilp([1, 1], [[1, 1]], [5])
        self.assertEqual(val, 5)


if __name__ == "__main__":
    unittest.main()
