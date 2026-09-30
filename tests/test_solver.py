import unittest
from fractions import Fraction

from freq_bounds.solver import solve_lp, solve_ilp, UnboundedError


class LpTest(unittest.TestCase):
    def test_min_with_negated_constraint(self):
        v, x = solve_lp([1], [([-1], -1)])
        self.assertEqual(v, 1)
        self.assertEqual(x, [Fraction(1)])

    def test_infeasible(self):
        self.assertIsNone(solve_lp([1], [([-1], -1), ([1], 0)])[0])

    def test_fractional_vertex(self):
        v, x = solve_lp([3, 2], [([-2, -1], -6), ([-1, -2], -8)])
        self.assertEqual(v, Fraction(32, 3))
        self.assertEqual(x, [Fraction(4, 3), Fraction(10, 3)])

    def test_equality_via_two_rows(self):
        v, x = solve_lp([1, 1], [([1, 1], 2), ([-1, -1], -2)])
        self.assertEqual(v, 2)


class IlpTest(unittest.TestCase):
    def test_integral_case(self):
        v, x = solve_ilp([3, 2], [([-2, -1], -6), ([-1, -2], -8)])
        self.assertEqual((v, x), (11, [1, 4]))

    def test_fractional_relaxation_rounds_to_integer_optimum(self):
        v, x = solve_ilp([1, 2], [([2, 1], 6), ([1, 2], 8)], maximize=True)
        self.assertEqual(v, 7)
        self.assertTrue(all(isinstance(v, int) for v in x))

    def test_integer_infeasible(self):
        self.assertIsNone(solve_ilp([1], [([-1], -1), ([1], 0)]))

    def test_maximize(self):
        v, _ = solve_ilp([1, 1], [([1, 1], 2)], maximize=True)
        self.assertEqual(v, 2)

    def test_all_nonnegative(self):
        v, x = solve_ilp([1, 1], [([1, 1], 5)])
        self.assertEqual((v, x), (0, [0, 0]))
        self.assertTrue(all(c >= 0 for c in x))


if __name__ == "__main__":
    unittest.main()
