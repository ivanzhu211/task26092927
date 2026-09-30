"""Exact pure-Python integer linear programming (standard library only).

Minuscule solver used by :mod:`freq_bounds`: a two-phase simplex with exact
:class:`fractions.Fraction` arithmetic for LP relaxations, wrapped in
branch-and-bound.  No floating point is ever used, so returned integer
optima are provably exact.
"""

from __future__ import annotations

from fractions import Fraction


class UnboundedError(RuntimeError):
    """The relaxation is unbounded (never expected for the frequency model)."""


def _f(value) -> Fraction:
    return value if isinstance(value, Fraction) else Fraction(value)



def _eliminate(tableau, row, pivot_row, col):
    factor = row[col]
    if factor != 0:
        pivot = tableau[pivot_row]
        for j in range(len(row)):
            row[j] -= factor * pivot[j]


def _reduce_row(row, pivot_row, col):
    factor = row[col]
    if factor != 0:
        for j in range(len(row)):
            row[j] -= factor * pivot_row[j]


def _simplex(tableau, basis, n_struct, var_order):
    """Minimize using Bland's rule. Mutates tableau/basis in place.

    ``var_order`` maps current columns to the *global* variable index, used
    by Bland's rule to choose the smallest-index entering variable.
    """
    n_cols = len(tableau[0])
    rhs = n_cols - 1
    n_rows = len(tableau) - 1
    while True:
        obj = tableau[-1]
        enter = None
        for j in range(n_cols - 1):
            if obj[j] < 0:
                if enter is None or var_order[j] < var_order[enter]:
                    enter = j
        if enter is None:
            return True  # optimum reached
        best_ratio = None
        leave = None
        for i in range(n_rows):
            aij = tableau[i][enter]
            if aij > 0:
                ratio = tableau[i][rhs] / aij
                if ratio < 0:
                    raise AssertionError("negative RHS encountered")
                if (
                    best_ratio is None
                    or ratio < best_ratio
                    or (ratio == best_ratio and basis[i] < basis[leave])
                ):
                    best_ratio = ratio
                    leave = i
        if leave is None:
            return False  # unbounded
        pivot = tableau[leave][enter]
        for j in range(n_cols):
            tableau[leave][j] /= pivot
        basis[leave] = var_order[enter]
        for i, row in enumerate(tableau):
            if i != leave:
                _eliminate(tableau, row, leave, enter)


def solve_lp(objective, constraints):
    """Solve an LP relaxation exactly (all variables non-negative).

    ``constraints`` is a list of ``(row, rhs)`` with semantics
    ``sum(row[j]) * x[j] <= rhs``.  Returns ``(value, x)`` with exact
    :class:`Fraction` values, or ``(None, None)`` when infeasible.
    """
    n_vars = len(objective)
    n_cons = len(constraints)
    parsed = [([_f(v) for v in row], _f(rhs_in)) for row, rhs_in in constraints]
    art_rows = [k for k, (_, rhs) in enumerate(parsed) if rhs < 0]
    n_art = len(art_rows)
    total_cols = n_vars + n_cons + n_art

    rows = []
    basis = []
    art_col_of = {}
    art_cursor = 0
    for k, (row_in, rhs) in enumerate(parsed):
        row = [Fraction(0)] * (total_cols + 1)
        if rhs < 0:
            for j, value in enumerate(row_in):
                row[j] = -value
            row[n_vars + k] = Fraction(-1)
            art_col = n_vars + n_cons + art_cursor
            row[art_col] = Fraction(1)
            row[-1] = -rhs
            basis.append(art_col)
            art_col_of[k] = art_col
            art_cursor += 1
        else:
            for j, value in enumerate(row_in):
                row[j] = value
            row[n_vars + k] = Fraction(1)
            row[-1] = rhs
            basis.append(n_vars + k)
        rows.append(row)

    var_order = list(range(total_cols))
    art_cols = set(art_col_of.values())

    if n_art:
        # Phase I: minimize the sum of artificial variables.  Build the
        # cost vector (1 per artificial, 0 otherwise) and eliminate the
        # basic columns, exactly as in phase II.
        phase_obj = [Fraction(0)] * total_cols + [Fraction(0)]
        for art_col in art_col_of.values():
            phase_obj[art_col] = Fraction(1)
        for i, b in enumerate(basis):
            _reduce_row(phase_obj, rows[i], b)
        tableau = rows + [phase_obj]
        if not _simplex(tableau, basis, n_vars, var_order):
            raise UnboundedError("phase I unbounded (unexpected)")
        if tableau[-1][-1] != 0:
            return None, None  # min of artificials > 0

        body = tableau[:-1]
        for idx in range(len(body)):
            if basis[idx] not in art_cols:
                continue
            pivot_col = None
            for j in range(n_vars + n_cons):
                if body[idx][j] != 0:
                    pivot_col = j
                    break
            if pivot_col is None:
                body[idx] = None  # redundant constraint
                continue
            pivot = body[idx][pivot_col]
            for j in range(total_cols + 1):
                body[idx][j] /= pivot
            basis[idx] = var_order[pivot_col]
            for i2, other in enumerate(body):
                if i2 != idx and other is not None:
                    _reduce_row(other, body[idx], pivot_col)

        rows = [r for r in body if r is not None]
        basis = [b for r, b in zip(body, basis) if r is not None]

    # Phase II objective: c @ x, expressed in the current basis.
    obj = [_f(c) for c in objective] + [Fraction(0)] * n_cons + [Fraction(0)]
    width = n_vars + n_cons + 1
    rows = [r[: n_vars + n_cons] + [r[-1]] for r in rows]
    for i, b in enumerate(basis):
        _reduce_row(obj, rows[i], b)
    tableau = rows + [obj]
    var_order2 = list(range(n_vars + n_cons))
    if not _simplex(tableau, basis, n_vars, var_order2):
        raise UnboundedError("LP relaxation unbounded")

    rhs_col = width - 1
    x = [Fraction(0)] * n_vars
    for i, b in enumerate(basis):
        if b < n_vars:
            x[b] = tableau[i][rhs_col]
    value = sum(_f(c) * x[j] for j, c in enumerate(objective))
    return value, x


def _is_integer(value: Fraction) -> bool:
    return value.denominator == 1


def solve_ilp(objective, constraints, maximize: bool = False):
    """Exactly optimize an integer linear program (variables non-negative).

    ``objective`` is a coefficient list; ``constraints`` is a list of
    ``(row, rhs)`` meaning ``row @ x <= rhs``.  Variables are non-negative
    integers.  Returns ``(value, x)`` with Python ints, or ``None`` when no
    integer feasible point exists.  Uses branch-and-bound on top of the exact
    simplex relaxation, so the returned optimum is provably exact.
    """
    obj = [-_f(c) for c in objective] if maximize else [_f(c) for c in objective]
    n_vars = len(obj)
    best = None  # (value, solution)

    stack = [(list(constraints), None)]
    while stack:
        base_cons, lp_bound = stack.pop()
        if best is not None and lp_bound is not None:
            if lp_bound > best[0]:
                continue
        value, x = solve_lp(obj, base_cons)
        if value is None:
            continue
        if best is not None and value > best[0]:
            continue
        frac = None
        frac_dist = None
        for j, coord in enumerate(x):
            if not _is_integer(coord):
                dist = min(coord - coord.__floor__(), coord.__ceil__() - coord)
                if frac_dist is None or dist > frac_dist:
                    frac = (j, coord)
                    frac_dist = dist
        if frac is None:
            candidate_x = [int(coord) for coord in x]
            candidate_v = sum(objective[j] * candidate_x[j] for j in range(n_vars))
            if best is None or candidate_v < best[0] or (
                maximize and candidate_v > best[0]
            ):
                best = (candidate_v, candidate_x)
            continue
        j, coord = frac
        lower = int(coord.__floor__())
        upper = lower + 1

        lower_row = [Fraction(0)] * n_vars
        lower_row[j] = Fraction(1)
        upper_row = [Fraction(0)] * n_vars
        upper_row[j] = Fraction(-1)
        # x_j <= floor ; x_j >= ceil  (== -x_j <= -ceil)
        stack.append((base_cons + [(lower_row, lower)], value))
        stack.append((base_cons + [(upper_row, -upper)], value))

    if best is None:
        return None
    return best[0], best[1]
