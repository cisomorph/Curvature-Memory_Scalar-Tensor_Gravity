"""
SIM157: background + linear growth for "dark curvature memory" (DARK_COUPLING_NOTES.md §2) and the reference models.

Units: H in 100 km/s/Mpc; densities in "ω units" (3H² = ρ ⇔ H² = Σω a^.. + Λ_b + V/3 + H²φ'²/6); N = ln a.
Dark matter couples through A(φ) = (1 − 2Λ0φ²)^{-1/2} (CMSTG family A around its attractor, δφ ≡ φ); baryons,
photons and neutrinos do not. Λ0 = 0 is minimally coupled thawing quintessence (TQ); v1 = 0 is ΛCDM.
"""
import os, sys
import numpy as np
from scipy.integrate import solve_ivp
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', 'common'))
sys.path.insert(0, os.path.join(HERE, '..', 'SIM156'))
import cmstg_obs as obs
import bglik as L                                     # SIM156's validated background likelihood

OMEGA_NU = L.OMEGA_NU
OMEGA_R = obs.OMEGA_R
Z_INI, Z_GROWTH = 3000.0, 100.0
SIGMA8_REF = 0.8076                                   # SIM155 hi_class ΛCDM D1 best fit (growth.json)
REF = dict(omega_b=0.02251, omega_cdm=0.11790, h=0.68300)   # SIM155 D1 ΛCDM best fit (σ8 normalisation)


class DarkBG:
    def __init__(self, omega_b, omega_c, L0=0.0, v1=0.0):
        self.ob, self.oc, self.L0, self.v1 = omega_b, omega_c, L0, v1

    def A(self, p):
        return (1.0 - 2.0 * self.L0 * p * p) ** -0.5

    def beta(self, p):
        return 2.0 * self.L0 * p / (1.0 - 2.0 * self.L0 * p * p)

    def H2(self, N, p, y, Lb):
        a = np.exp(N)
        r = (self.ob + OMEGA_NU + self.oc * self.A(p)) / a**3 + OMEGA_R / a**4 + Lb + self.v1 * p
        c = 1.0 - y * y / 6.0
        if np.ndim(c):                                    # vectorised use (H_func)
            return np.where(c > 1e-3, r / np.where(c > 1e-3, c, 1.0), np.nan)
        return r / c if c > 1e-3 else np.nan

    def rhs(self, N, s, Lb, growth):
        p, y = s[0], s[1]
        a = np.exp(N)
        H2 = self.H2(N, p, y, Lb)
        Ac = self.oc * self.A(p) / a**3
        rho_c = 3.0 * Ac                                 # 3H² = ρ units
        eps = (3.0 * (self.ob + OMEGA_NU) / a**3 + rho_c + 4.0 * OMEGA_R / a**4 + H2 * y * y) / (2.0 * H2)
        b = self.beta(p)
        dy = -(3.0 - eps) * y - (3.0 * self.v1 + b * rho_c) / H2
        out = [y, dy]
        if growth:
            dc, dcp, db, dbp = s[2:6]
            Om_c, Om_b = Ac / H2, self.ob / a**3 / H2
            out += [dcp, -(2.0 - eps + b * y) * dcp + 1.5 * ((1 + 2 * b * b) * Om_c * dc + Om_b * db),
                    dbp, -(2.0 - eps) * dbp + 1.5 * (Om_c * dc + Om_b * db)]
        return out

    def _wall(self, N, s, Lb, growth):
        """Terminal event: F = ½ − Λ0φ² → 0 makes A = (2F)^{-1/2} (the dark-matter mass) diverge; stop at 2Λ0φ² = 0.95."""
        return 0.95 - 2.0 * self.L0 * s[0] * s[0]

    def solve(self, h, n_iter=10, tol=1e-11):
        """Closure: iterate Λ_b so that H(z=0) = 100h. Returns dict with background and growth tables, or None if the
        field reaches the singular point of A (unphysical: diverging dark-matter mass)."""
        Lb = h * h - (self.ob + OMEGA_NU + self.oc + OMEGA_R)
        Ni, Ng = -np.log1p(Z_INI), -np.log1p(Z_GROWTH)
        wall = lambda N, s, Lb_, g_: self._wall(N, s, Lb_, g_)
        wall.terminal = True
        for _ in range(n_iter):
            s1 = solve_ivp(self.rhs, (Ni, Ng), [0.0, 0.0], args=(Lb, False), rtol=1e-9, atol=1e-12, max_step=0.1,
                           dense_output=True, events=wall)
            if not s1.success or s1.status == 1:
                return None
            ag = np.exp(Ng)
            s2 = solve_ivp(self.rhs, (Ng, 0.0), [s1.y[0, -1], s1.y[1, -1], ag, ag, ag, ag], args=(Lb, True),
                           rtol=1e-9, atol=1e-12, max_step=0.02, dense_output=True, events=wall)
            if not s2.success or s2.status == 1:
                return None
            p0, y0 = s2.y[0, -1], s2.y[1, -1]
            Lb_new = h * h * (1.0 - y0 * y0 / 6.0) - (self.ob + OMEGA_NU + self.oc * self.A(p0) + OMEGA_R) - self.v1 * p0
            if abs(Lb_new - Lb) < tol:
                Lb = Lb_new
                break
            Lb = Lb_new
        N = np.concatenate([s1.t[:-1], s2.t])
        P = np.concatenate([s1.y[0, :-1], s2.y[0]]); Y = np.concatenate([s1.y[1, :-1], s2.y[1]])
        H = 100.0 * np.sqrt([self.H2(n, p, y, Lb) for n, p, y in zip(N, P, Y)])
        if not np.all(np.isfinite(H)):
            return None
        z = np.exp(-N) - 1.0
        o = np.argsort(z)
        return dict(z=z[o], H=H[o], phi=P[o], y=Y[o], Lb=Lb, phi0=float(p0), y0=float(y0), A0=float(self.A(p0)),
                    growth=s2, early=s1, beta_hi=float(np.max(np.abs(self.beta(s1.y[0])[np.exp(-s1.t) - 1 >= 1000.0]))
                                         if np.any(np.exp(-s1.t) - 1 >= 1000.0) else 0.0))

    def H_func(self, bg):
        """H(z) evaluated from the constraint on the dense ODE solution (no table interpolation): needed because
        the compressed CMB prior constrains D_M(z*) to 3e-4. Above z = 3000 the field is frozen at φ = 0, so H is exact."""
        Lb, s1, s2 = bg['Lb'], bg['early'], bg['growth']
        Ng, Ni = -np.log1p(Z_GROWTH), -np.log1p(Z_INI)
        om = self.ob + OMEGA_NU + self.oc

        def H(z):
            z = np.atleast_1d(np.asarray(z, dtype=float))
            N = -np.log1p(z)
            out = np.empty_like(z)
            hi = N < Ni
            out[hi] = 100.0 * np.sqrt(om * (1 + z[hi])**3 + OMEGA_R * (1 + z[hi])**4 + Lb)
            mid = (~hi) & (N < Ng)
            lo = N >= Ng
            for m, sol in ((mid, s1), (lo, s2)):
                if np.any(m):
                    st = sol.sol(N[m])
                    out[m] = 100.0 * np.sqrt(self.H2(N[m], st[0], st[1], Lb))
            return out
        return H

    def growth_table(self, bg, zs):
        """δ_m(z) and dδ_m/dN at the requested redshifts (total matter, c + b)."""
        out_d, out_dp = [], []
        for z in zs:
            N = -np.log1p(z)
            s = bg['growth'].sol(N)
            p, dc, dcp, db, dbp = s[0], s[2], s[3], s[4], s[5]
            a = np.exp(N)
            wc, wb = self.oc * self.A(p) / a**3, self.ob / a**3
            dm = (wc * dc + wb * db) / (wc + wb)
            # dδ_m/dN including the change of the dark-matter weight (A changes with φ)
            y = s[1]
            dwc = wc * (self.beta(p) * y - 3.0); dwb = -3.0 * wb
            ddm = ((dwc * dc + wc * dcp + dwb * db + wb * dbp) * (wc + wb) - (wc * dc + wb * db) * (dwc + dwb)) / (wc + wb)**2
            out_d.append(dm); out_dp.append(ddm)
        return np.array(out_d), np.array(out_dp)


def growth_lcdm_like(H_of_z, omega_b, omega_c, zs):
    """GR growth with smooth dark energy (used for w0wa): δ'' + (2 − ε)δ' = 1.5 Ω_m δ from z = 100, where
    Ω_m counts baryons + CDM (neutrinos smooth, as in DarkBG) and ε = −d ln H/dN by central differences."""
    lnH = lambda N: np.log(H_of_z(np.array([np.exp(-N) - 1.0]))[0] / 100.0)

    def rhs(N, s):
        a = np.exp(N)
        d = 1e-4
        eps = -(lnH(N + d) - lnH(N - d)) / (2 * d)
        Om = (omega_b + omega_c) / a**3 / np.exp(2 * lnH(N))
        return [s[1], -(2.0 - eps) * s[1] + 1.5 * Om * s[0]]
    Ng = -np.log1p(Z_GROWTH); ag = np.exp(Ng)
    sol = solve_ivp(rhs, (Ng, 0.0), [ag, ag], rtol=1e-9, atol=1e-12, max_step=0.02, dense_output=True)
    d = np.array([sol.sol(-np.log1p(z))[0] for z in zs]); dp = np.array([sol.sol(-np.log1p(z))[1] for z in zs])
    return d, dp


_REF_D0 = None


def ref_delta0():
    """δ_m(z=0) of the reference ΛCDM (σ8 normalisation), same z = 100 start."""
    global _REF_D0
    if _REF_D0 is None:
        om = REF['omega_b'] + REF['omega_cdm'] + OMEGA_NU
        H = obs.H_lcdm(REF['h'], om, OMEGA_R)
        _REF_D0 = growth_lcdm_like(H, REF['omega_b'], REF['omega_cdm'], [0.0])[0][0]
    return _REF_D0


def w_eff(bg, ob, oc):
    """Apparent dark-energy equation of state for an observer assuming uncoupled DM with today's density."""
    z, H = bg['z'], bg['H'] / 100.0
    m = z < 3.0
    z, H = z[m], H[m]
    a = 1 / (1 + z)
    rde = H**2 - (ob + OMEGA_NU + oc * bg['A0']) / a**3 - OMEGA_R / a**4
    lnr = np.log(np.clip(rde, 1e-12, None)); lna = np.log(a)
    o = np.argsort(lna)                                   # differentiate in ln a (ascending) ...
    w = -1.0 - np.gradient(lnr[o], lna[o]) / 3.0
    zo = z[o]
    return zo[::-1], w[::-1]                              # ... return with z ascending (np.interp needs it)
