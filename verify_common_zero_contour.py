"""Exact rational-interval checks for the common-zero contour proof.

The proof and certificate specification appear in the accompanying paper,
"Divisible formal weight enumerators: Chebyshev determinants, extremal
bounds, and zeros of Duursma zeta polynomials."

This script checks the finite algebraic inequalities in the finite
certificate. Only ``fractions.Fraction`` arithmetic is used for the
certificate inequalities. SymPy is used separately to count real algebraic
roots in the stated rational isolating intervals.
"""

from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction as F
from math import comb

import sympy as sp

if not __debug__:
    raise RuntimeError(
        "This verifier must be run without Python's -O optimization option."
    )

@dataclass(frozen=True)
class I:
    lo: F
    hi: F

    @staticmethod
    def point(x: int | F) -> "I":
        x = F(x)
        return I(x, x)

    def __add__(self, other: "I | int | F") -> "I":
        other = as_i(other)
        return I(self.lo + other.lo, self.hi + other.hi)

    __radd__ = __add__

    def __neg__(self) -> "I":
        return I(-self.hi, -self.lo)

    def __sub__(self, other: "I | int | F") -> "I":
        return self + (-as_i(other))

    def __rsub__(self, other: "I | int | F") -> "I":
        return as_i(other) - self

    def __mul__(self, other: "I | int | F") -> "I":
        other = as_i(other)
        vals = (
            self.lo * other.lo,
            self.lo * other.hi,
            self.hi * other.lo,
            self.hi * other.hi,
        )
        return I(min(vals), max(vals))

    __rmul__ = __mul__

    def inv(self) -> "I":
        assert not (self.lo <= 0 <= self.hi), self
        vals = (1 / self.lo, 1 / self.hi)
        return I(min(vals), max(vals))

    def __truediv__(self, other: "I | int | F") -> "I":
        return self * as_i(other).inv()

    def __rtruediv__(self, other: "I | int | F") -> "I":
        return as_i(other) / self

    def __pow__(self, n: int) -> "I":
        assert n >= 0
        if n == 0:
            return I.point(1)
        if n % 2 == 0 and self.lo < 0 < self.hi:
            return I(F(0), max(abs(self.lo), abs(self.hi)) ** n)
        vals = (self.lo**n, self.hi**n)
        return I(min(vals), max(vals))

    def strictly_positive(self) -> bool:
        return self.lo > 0

    def strictly_below_one(self) -> bool:
        return self.hi < 1

    def outward_decimal(self, digits: int = 10) -> str:
        """Print a rigorous decimal enclosure, rounded outwards."""

        scale = 10**digits

        def floor_scaled(x: F) -> int:
            return (x.numerator * scale) // x.denominator

        def ceil_scaled(x: F) -> int:
            return -((-x.numerator * scale) // x.denominator)

        def render(k: int) -> str:
            sign = "-" if k < 0 else ""
            k = abs(k)
            return f"{sign}{k // scale}.{k % scale:0{digits}d}"

        return f"[{render(floor_scaled(self.lo))}, {render(ceil_scaled(self.hi))}]"


def as_i(x: I | int | F) -> I:
    return x if isinstance(x, I) else I.point(x)


# Sparse bivariate polynomials in (x,t), with interval coefficients.
Poly2 = dict[tuple[int, int], I]


def p_const(x: I | int | F) -> Poly2:
    return {(0, 0): as_i(x)}


def p_var_x() -> Poly2:
    return {(1, 0): I.point(1)}


def p_var_t() -> Poly2:
    return {(0, 1): I.point(1)}


def p_add(a: Poly2, b: Poly2) -> Poly2:
    out = dict(a)
    for key, value in b.items():
        out[key] = out.get(key, I.point(0)) + value
    return out


def p_neg(a: Poly2) -> Poly2:
    return {key: -value for key, value in a.items()}


def p_sub(a: Poly2, b: Poly2) -> Poly2:
    return p_add(a, p_neg(b))


def p_mul(a: Poly2, b: Poly2) -> Poly2:
    out: Poly2 = {}
    for (i, j), av in a.items():
        for (k, ell), bv in b.items():
            key = (i + k, j + ell)
            out[key] = out.get(key, I.point(0)) + av * bv
    return out


def p_scale(a: Poly2, x: I | int | F) -> Poly2:
    x = as_i(x)
    return {key: value * x for key, value in a.items()}


def p_pow(a: Poly2, n: int) -> Poly2:
    out = p_const(1)
    base = a
    while n:
        if n & 1:
            out = p_mul(out, base)
        base = p_mul(base, base)
        n //= 2
    return out


def bernstein_coefficients(poly: Poly2, nx: int, nt: int) -> list[I]:
    """Power basis on [0,1]^2 -> tensor Bernstein basis of degrees nx,nt."""
    out: list[I] = []
    for i in range(nx + 1):
        for j in range(nt + 1):
            total = I.point(0)
            for (k, ell), value in poly.items():
                if k <= i and ell <= j:
                    total += value * F(comb(i, k), comb(nx, k)) * F(
                        comb(j, ell), comb(nt, ell)
                    )
            out.append(total)
    return out


DEN = 10**12

RINGS = {
    4: {
        "rho": I(F(2613125929752, DEN), F(2613125929753, DEN)),
        "y": I(F(82392200292, DEN), F(82392200293, DEN)),
        "R": F(1, 8),
        "rho_poly": [8, 0, -8, 0, 1],
        "y_poly": [1, -12, -2, 4, 1],
        "conv_degree": 4,
        "conv_floor": F(18, 10),
        "q_floor": F(7, 1000),
        "ceil_A": 2,
        "ceil_Y": 1,
        "bad_ratio_ceil": F(72, 1000),
        "derivative_floor": F(39, 1000),
        "base_good_ceil": F(169, 1000),
        "base_bad_ceil": F(24, 1000),
    },
    5: {
        "rho": I(F(1701301616704, DEN), F(1701301616705, DEN)),
        "y": I(F(115073031670, DEN), F(115073031671, DEN)),
        "R": F(9, 40),
        "rho_poly": [16, 0, -20, 0, 5],
        "y_poly": [5, -40, -30, 0, 1],
        "conv_degree": 6,
        "conv_floor": F(120),
        "q_floor": F(2, 10000),
        "ceil_A": 2,
        "ceil_Y": 2,
        "bad_ratio_ceil": F(644, 1000),
        "derivative_floor": F(3, 1000),
        "base_good_ceil": F(3089, 10000),
        "base_bad_ceil": F(234, 1000),
    },
    6: {
        "rho": I(F(3863703305156, DEN), F(3863703305157, DEN)),
        "y": I(F(35276180410, DEN), F(35276180411, DEN)),
        "R": F(1, 16),
        "rho_poly": [16, 0, -16, 0, 1],
        "y_poly": [1, -28, -10, 4, 1],
        "conv_degree": 4,
        "conv_floor": F(3, 10),
        "q_floor": F(2, 1000),
        "ceil_A": 3,
        "ceil_Y": 2,
        "bad_ratio_ceil": F(10, 1000),
        "derivative_floor": F(77, 100),
        "base_good_ceil": F(253, 1000000),
        "base_bad_ceil": F(10, 1000000),
    },
}


def isolate_roots() -> None:
    x = sp.symbols("x")
    for data in RINGS.values():
        for name in ("rho", "y"):
            iv: I = data[name]
            coeffs = data[f"{name}_poly"]
            poly = sp.Poly(sum(sp.Integer(a) * x**k for k, a in enumerate(coeffs)), x)
            assert poly.count_roots(sp.Rational(iv.lo.numerator, iv.lo.denominator),
                                    sp.Rational(iv.hi.numerator, iv.hi.denominator)) == 1
            # The trigonometric values used in the proof are identified by
            # their root index: rho is the unique root in (3/2,+infinity),
            # and y_r is the unique root in (0,1).
            if name == "rho":
                assert poly.count_roots(sp.Rational(3, 2), sp.oo) == 1
            else:
                assert poly.count_roots(0, 1) == 1


def finite_checks(r: int) -> dict[str, I]:
    data = RINGS[r]
    rho: I = data["rho"]
    y: I = data["y"]
    R = I.point(data["R"])
    q = rho**2
    one = I.point(1)

    assert (q * R).hi < 1
    assert R.hi < (one / rho).lo

    h = rho - 1
    S0 = 4 * rho / (h**2)

    # Convexity of f_z(s).  c=-1/2+3x/2, s=S0*t.
    x_poly = p_var_x()
    t_poly = p_var_t()
    c_poly = p_add(p_const(F(-1, 2)), p_scale(x_poly, F(3, 2)))
    s_poly = p_scale(t_poly, S0)
    aa = q * R
    bb = one / R
    u_poly = p_scale(
        p_add(p_scale(c_poly, aa + bb), p_const(-2 * rho)), one / (h**2)
    )
    v_poly = p_scale(
        p_add(
            p_add(p_const((aa - bb) ** 2 + 4 * q), p_scale(c_poly, -4 * rho * (aa + bb))),
            p_scale(p_pow(c_poly, 2), 4 * q),
        ),
        one / (h**4),
    )
    modulus_sq = p_add(p_add(p_pow(s_poly, 2), p_scale(p_mul(u_poly, s_poly), 2)), v_poly)
    imag_sq = p_sub(v_poly, p_pow(u_poly, 2))
    curvature_num = p_add(
        p_pow(modulus_sq, 2),
        p_mul(
            p_sub(imag_sq, p_pow(p_add(s_poly, u_poly), 2)),
            p_pow(p_add(s_poly, p_const(1)), 2),
        ),
    )
    degree = data["conv_degree"]
    curvature_bern = bernstein_coefficients(curvature_num, degree, degree)
    curvature_min = min(curvature_bern, key=lambda value: value.lo)
    assert curvature_min.lo > data["conv_floor"], curvature_min

    # The n-step ratio on the good arc is <1.  Prove the quartic
    # R^4(q-1)^2-A(c)^2(B0^2-B1^2 c^2) is positive in Bernstein form.
    A0 = y**2 + (one - y) ** 2 * R**2
    A1 = 2 * y * (one - y) * R
    B0 = one + q * R**2
    B1 = 2 * rho * R
    A_poly = p_add(p_const(A0), p_scale(c_poly, A1))
    Q_poly = p_sub(
        p_const(R**4 * (q - 1) ** 2),
        p_mul(p_pow(A_poly, 2), p_sub(p_const(B0**2), p_scale(p_pow(c_poly, 2), B1**2))),
    )
    q_bern = bernstein_coefficients(Q_poly, 4, 0)
    q_min = min(q_bern, key=lambda value: value.lo)
    assert q_min.lo > data["q_floor"], q_min

    # d-step ratio on the good arc: U_-U_+>1 and U_-/U_+>1.
    product_at_one = ((one - q * R**2) / (R * (q - 1))) ** 2
    assert product_at_one.lo > 1
    assert (one - q * R).lo > 0  # equivalent to U_-(1)>U_+(1)

    # The n-step ratio on the bad arc is <1.
    c_bad = F(-1, 2)
    A_bad = A0 + A1 * c_bad
    z_star = (one + rho * R) ** 2 / (R * (rho - 1) ** 2)
    w_star = (S0 + z_star) / (S0 + 1)
    bad_n_ratio_sq = (A_bad / R) ** 2 * z_star * w_star
    assert bad_n_ratio_sq.strictly_below_one(), bad_n_ratio_sq
    assert bad_n_ratio_sq.hi < data["bad_ratio_ceil"], bad_n_ratio_sq
    assert z_star.lo > 1 and w_star.lo > 1

    # Base bounds n=r,d=2.  A+Y=r-2.  Integer ceilings suffice for
    # r=4,6.  For r=5 use A<5/3 and a tight maximizer bracket.
    C = (one + q * R**2) / ((one - R) * (one - q * R))
    D = F(r - 2)
    root_sum = (q - 2) / rho
    Aexp = (2 * D + root_sum) / 4
    Yexp = D - Aexp
    assert Aexp.lo > 0 and Yexp.lo > 0
    assert Aexp.hi < data["ceil_A"], Aexp
    assert Yexp.hi < data["ceil_Y"], Yexp

    derivative_left: I
    derivative_right: I
    if r in (4, 6):
        c_good = F(1)
        derivative = (
            r * A1 / (A0 + A1 * c_good)
            - Aexp * B1 / (B0 - B1 * c_good)
            + Yexp * B1 / (B0 + B1 * c_good)
        )
        assert derivative.lo > data["derivative_floor"]
        derivative_left = derivative
        derivative_right = derivative
        Uminus = (B0 - B1) / (R * (rho - 1) ** 2)
        Uplus = (B0 + B1) / (R * (rho + 1) ** 2)
        assert Uminus.lo > 1 and Uplus.lo > 1
        base_good = C * ((A0 + A1) / R) ** r * Uminus ** data["ceil_A"] * Uplus ** data["ceil_Y"]
        assert base_good.strictly_below_one(), base_good
        assert base_good.hi < data["base_good_ceil"], base_good
    else:
        c_lo, c_hi = F(943, 1000), F(118, 125)
        deriv_lo = (
            r * A1 / (A0 + A1 * c_lo)
            - Aexp * B1 / (B0 - B1 * c_lo)
            + Yexp * B1 / (B0 + B1 * c_lo)
        )
        deriv_hi = (
            r * A1 / (A0 + A1 * c_hi)
            - Aexp * B1 / (B0 - B1 * c_hi)
            + Yexp * B1 / (B0 + B1 * c_hi)
        )
        assert deriv_lo.lo > data["derivative_floor"]
        assert deriv_hi.hi < -data["derivative_floor"]
        derivative_left = deriv_lo
        derivative_right = deriv_hi
        assert Aexp.hi < F(5, 3)
        Uminus = (B0 - B1 * c_lo) / (R * (rho - 1) ** 2)
        Uplus = (B0 + B1 * c_hi) / (R * (rho + 1) ** 2)
        assert Uminus.lo > Uplus.hi > 1
        # [C A(c)^5]^3 U_-^5 U_+^4 <1 implies the desired base bound.
        base_good = (C * ((A0 + A1 * c_hi) / R) ** 5) ** 3 * Uminus**5 * Uplus**4
        assert base_good.strictly_below_one(), base_good
        assert base_good.hi < data["base_good_ceil"], base_good

    base_bad = C * (A_bad / R) ** r * z_star ** data["ceil_A"] * w_star ** data["ceil_Y"]
    assert base_bad.strictly_below_one(), base_bad
    assert base_bad.hi < data["base_bad_ceil"], base_bad

    return {
        "curvature_min": curvature_min,
        "q_min": q_min,
        "bad_n_ratio_sq": bad_n_ratio_sq,
        "base_derivative_left": derivative_left,
        "base_derivative_right": derivative_right,
        "base_good_certificate": base_good,
        "base_bad": base_bad,
    }


def main() -> None:
    isolate_roots()
    for r in (4, 5, 6):
        result = finite_checks(r)
        print(f"r={r}")
        for name, value in result.items():
            print(f"  {name:24s} {value.outward_decimal(12)}")
    print("All certificate checks passed.")
    

if __name__ == "__main__":
    main()
