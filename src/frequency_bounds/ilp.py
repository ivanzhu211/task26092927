"""Exact integer linear programming (standard library only).

A Fraction-exact two-phase simplex (Bland's rule, no cycling) solves the
LP relaxation; branch-and-bound finds the exact integer optimum.

All problem variables are non-negative.  Constraints are of the form
``A x <= b``.  Coefficients may be negative.
"""

from __future__ import annotations

from fractions import Fraction


class ILPError(Exception):
    pass


INFEASIBLE = "infeasible"
UNBOUNDED = "unbounded"


def _as_fraction_matrix(objective, rows, rhs):
    n = len(objective)
    obj = [Fraction(v) for v in objective]
    a = [[Fraction(v) for v in row] for row in rows]
    b = [Fraction(v) for v in rhs]
    for i, row in enumerate(a):
        if len(row) != n:
            raise ILPError(f"row {i} has {len(row)} coefficients, expected {n}")
    if len(b) != len(a):
        raise ILPError("rhs length does not match number of rows")
    return obj, a, b


def solve_lp(objective, rows, rhs):
    """Maximize ``objective . x`` subject to ``rows @ x <= rhs``, ``x >= 0``.

    Returns ``(value, x)`` with exact Fractions; ``x`` has integral values
    converted to ``int``.  Returns ``INFEASIBLE`` if no feasible point
    exists and ``UNBOUNDED`` if the objective is unbounded above.
    """
    c, a, b = _as_fraction_matrix(objective, rows, rhs)
    n = len(c)
    m = len(a)

    # Column layout: [0,n) originals, [n,n+m) slacks, then artificials.
    art_cols: dict[int, int] = {}
    ncols = n + m
    tableau = [row[:] + [Fraction(0)] * m for row in a]
    trhs = b[:]
    for i in range(m):
        tableau[i][n + i] = Fraction(1)  # slack s_i with coefficient +1
        if trhs[i] < 0:
            # Flip sign first (over original + slack columns) so RHS >= 0;
            # the appended artificial then gets coefficient +1.
            for j in range(n + m):
                tableau[i][j] = -tableau[i][j]
            art = ncols
            ncols += 1
            art_cols[i] = art
            tableau[i].append(Fraction(1))
            trhs[i] = -trhs[i]
    for i in range(m):
        while len(tableau[i]) < ncols:
            tableau[i].append(Fraction(0))

    basis = {}  # row -> basic column
    for i in range(m):
        basis[i] = art_cols.get(i, n + i)

    obj = [Fraction(0)] * ncols

    def pivot(col, pr):
        pivot_val = tableau[pr][col]
        for j in range(ncols):
            tableau[pr][j] /= pivot_val
        trhs[pr] /= pivot_val
        for r in range(m):
            if r == pr:
                continue
            factor = tableau[r][col]
            if factor == 0:
                continue
            for j in range(ncols):
                tableau[r][j] -= factor * tableau[pr][j]
            trhs[r] -= factor * trhs[pr]
        factor = obj[col]
        if factor != 0:
            for j in range(ncols):
                obj[j] -= factor * tableau[pr][j]
            obj_rhs[0] -= factor * trhs[pr]
        basis[pr] = col

    obj_rhs = [Fraction(0)]

    # ---------- Phase I: maximize -(sum of artificials) ----------
    if art_cols:
        # Stored objective row is  z - c.x ; c = -1 for artificials.
        for art in art_cols.values():
            obj[art] = Fraction(1)
        for i, art in art_cols.items():
            factor = obj[art]  # = 1; subtract factor * row to zero it
            for j in range(ncols):
                obj[j] -= factor * tableau[i][j]
            obj_rhs[0] -= factor * trhs[i]

        while True:
            entering = next((j for j in range(ncols) if obj[j] < 0), None)
            if entering is None:
                break
            pr = _ratio_row(tableau, trhs, m, entering)
            if pr is None:
                raise ILPError("phase I reported unbounded (internal error)")
            pivot(entering, pr)

        if obj_rhs[0] != 0:
            return INFEASIBLE

        # An artificial that stays basic has value zero. Pivot it out on any
        # nonzero original/slack coefficient of that row (ratio 0 is valid);
        # a fully zero row is genuinely redundant and may stay as-is.
        art_set = set(art_cols.values())
        for i in range(m):
            if basis[i] in art_set:
                swap = next(
                    (j for j in range(n + m) if tableau[i][j] != 0), None
                )
                if swap is not None:
                    pivot(swap, i)

    # ---------- Phase II ----------
    obj = [Fraction(0)] * ncols
    for j in range(n):
        obj[j] = -c[j]
    obj_rhs[0] = Fraction(0)
    for i in range(m):
        bc = basis[i]
        factor = obj[bc]
        if factor != 0:
            for j in range(ncols):
                obj[j] -= factor * tableau[i][j]
            obj_rhs[0] -= factor * trhs[i]

    while True:
        entering = next((j for j in range(n + m) if obj[j] < 0), None)
        if entering is None:
            break
        pr = _ratio_row(tableau, trhs, m, entering)
        if pr is None:
            return UNBOUNDED
        pivot(entering, pr)

    x = [Fraction(0)] * n
    for i, col in basis.items():
        if col < n:
            x[col] = trhs[i]
    value = obj_rhs[0]
    return value, [int(v) if v.denominator == 1 else v for v in x]


def _ratio_row(tableau, trhs, m, col):
    best = None
    best_ratio = None
    for i in range(m):
        a = tableau[i][col]
        if a <= 0:
            continue
        ratio = trhs[i] / a
        if best_ratio is None or ratio < best_ratio:
            best_ratio = ratio
            best = i
    return best


def solve_ilp(objective, rows, rhs):
    """Integer maximization counterpart of :func:`solve_lp`.

    Returns ``(value, x)`` with Python ints, or ``None`` when infeasible.
    Objective coefficients must be integral so that LP bounds can be
    floored exactly.
    """
    for v in objective:
        if isinstance(v, Fraction):
            if v.denominator != 1:
                raise ILPError("integer objective coefficients required")
        elif int(v) != v:
            raise ILPError("integer objective coefficients required")

    best = _branch(objective, [list(r) for r in rows], list(rhs), None, 0)
    if best is None:
        return None
    value, x = best
    return int(value), [int(v) for v in x]


def _branch(objective, rows, rhs, incumbent, depth):
    lp = solve_lp(objective, rows, rhs)
    if lp is INFEASIBLE:
        return incumbent
    value, x = lp

    if incumbent is not None and value <= incumbent[0]:
        return incumbent  # cannot improve (strictly better needed)

    if all(isinstance(v, int) for v in x):
        candidate = (value, x)
        if incumbent is None or value > incumbent[0]:
            return candidate
        return incumbent

    # Bound for any integer solution in this node.
    upper = value.numerator // value.denominator  # floor(Fraction)
    if incumbent is not None and upper <= incumbent[0]:
        return incumbent

    frac_idx = _most_fractional(x)
    fv = x[frac_idx]
    floor_v = fv.numerator // fv.denominator

    upper_row = [0] * len(objective)
    upper_row[frac_idx] = 1
    lower_row = [0] * len(objective)
    lower_row[frac_idx] = -1

    for extra, extra_rhs in (
        (upper_row, floor_v),             # x_k <= floor(v)
        (lower_row, -(floor_v + 1)),      # x_k >= ceil(v)
    ):
        incumbent = _branch(
            objective,
            rows + [extra],
            rhs + [Fraction(extra_rhs) if not isinstance(extra_rhs, int) else extra_rhs],
            incumbent,
            depth + 1,
        )
    return incumbent


def _most_fractional(x):
    best_idx = None
    best_dist = None
    for i, v in enumerate(x):
        if isinstance(v, int):
            continue
        frac = v - Fraction(v.numerator // v.denominator)
        dist = abs(frac - Fraction(1, 2))
        if best_dist is None or dist < best_dist:
            best_dist = dist
            best_idx = i
    return best_idx
