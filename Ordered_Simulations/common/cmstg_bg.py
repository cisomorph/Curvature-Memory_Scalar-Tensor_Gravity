"""
cmstg_bg.py — shared background-equation helpers for the scalar-tensor sims.

Added in the 2026-10-22 audit (see AUDIT_NOTES.md, bugs B9/B10).

Action family (signature -+++, units 8πG = 1, so GR is F = 1/2):
    S = ∫√-g [ F(Ψ) R − ½(∂Ψ)² − V(Ψ) ] + S_m

Friedmann (lapse variation), with ρ in GR units (3H² = ρ) and ' = d/dN:
    H² (6F + 6F'Ψ' − ½Ψ'²) = ρ + V
The sims express densities as ω (GR: H² = Σω), so ρ = 3ω:
    H² (2F + 2F'Ψ' − Ψ'²/6) = Σω + V/3

Scalar EOM:  Ψ'' + (3 − ε)Ψ' = (F'R − V')/H²,   R/H² = 6(2 − ε)

Raychaudhuri (BEP/DEF form with F_BEP = 2F, Z = 1), closed form:
    ε ≡ −Ḣ/H² = [(ρ+p)/H² + (1+2F'')Ψ'² + 24F'² − 2F'V'/H² − 8F'Ψ'] / (4F + 12F'²)
GR checks: dust → 3/2, radiation → 2.

The published sims used H²(3F + 6F'Ψ' − ½Ψ'²) = Σω (+…), which breaks the GR limit
(H² = ω/1.5) and doubles the F'Ψ' and kinetic terms; and they took ε from ∂_N ln H²
at fixed (Ψ, Ψ'), dropping the field's own evolution.
"""

import numpy as np


def eps_closed_form(F, Fp, Fpp, y, H2, rho_plus_p, Vp=0.0):
    """Raychaudhuri ε for the F(Ψ)R + canonical scalar family (all in ρ-units, 3H²=ρ in GR)."""
    num = rho_plus_p / H2 + (1.0 + 2.0 * Fpp) * y**2 + 24.0 * Fp**2 - 2.0 * Fp * Vp / H2 - 8.0 * Fp * y
    return num / (4.0 * F + 12.0 * Fp**2)


def eps_from_constraint(H2_func, Psi, y, N, A, B, h=1e-6):
    """
    ε = −½ d ln H²/dN along the trajectory, for any constraint H²(Ψ, Ψ', N)
    and a scalar EOM of the form Ψ'' = A + B ε (true whenever the source is linear
    in R/H² = 6(2 − ε)).  Partial derivatives by central differences.
    """
    H2 = H2_func(Psi, y, N)
    if not np.isfinite(H2) or H2 <= 0:
        return np.nan
    lnH2 = lambda p, yy, n: np.log(H2_func(p, yy, n))
    dP = h * max(1.0, abs(Psi))
    dY = h * max(1.0, abs(y))
    dN = 1e-5
    aN = (lnH2(Psi, y, N + dN) - lnH2(Psi, y, N - dN)) / (2 * dN)
    aP = (lnH2(Psi + dP, y, N) - lnH2(Psi - dP, y, N)) / (2 * dP)
    aY = (lnH2(Psi, y + dY, N) - lnH2(Psi, y - dY, N)) / (2 * dY)
    return -0.5 * (aN + aP * y + aY * A) / (1.0 + 0.5 * aY * B)


def secant_match(fun, x0, x1, tol=1e-6, maxit=30):
    """Root of fun(x) = 0 by the secant method (fun returns nan on failure)."""
    f0, f1 = fun(x0), fun(x1)
    for _ in range(maxit):
        if not (np.isfinite(f0) and np.isfinite(f1)) or f1 == f0:
            return np.nan
        x2 = x1 - f1 * (x1 - x0) / (f1 - f0)
        if abs(x2 - x1) < tol:
            return x2
        x0, f0, x1, f1 = x1, f1, x2, fun(x2)
    return np.nan
