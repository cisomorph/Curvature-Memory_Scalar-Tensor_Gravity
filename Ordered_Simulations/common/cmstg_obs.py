"""
cmstg_obs.py — shared, validated observables for the CMSTG simulations.

Added in the 2026-10-22 audit (see AUDIT_NOTES.md at the repo root). It replaces
per-sim copies of code that had these bugs:

  * a hardcoded "DESI H(z)" table that is not DESI data
        -> official DESI DR1 (2024) / DR2 (2025) Gaussian BAO likelihood
           (D_M/r_d, D_H/r_d, D_V/r_d with full covariance; CobayaSampler/bao_data)
  * sound-horizon integrals truncated at z = 1e4 or 1e5
        -> r_s integrated to z = infinity (integration variable ln a, from a ~ 0)
  * theta* = r_s(z_drag) / D_C(z*)
        -> theta* = r_s(z*) / D_M(z*)
  * "tension" = sqrt(chi2 / N), reported as a number of sigma
        -> chi2, dof, p-value and Gaussian-equivalent significance

Conventions
  H_of_z : vectorised callable, z (array) -> H(z) in km/s/Mpc. It must be valid
           from z = 0 up to whatever redshift the observable needs (z ~ 1e9 for
           r_s). Use `H_from_table` to wrap a numerically solved background and
           extend it above its last tabulated redshift.
  omegas : physical densities omega_i = Omega_i h^2, in the same convention as
           the sims (omega_m includes all matter; omega_r = photons + 3.046
           massless neutrinos).
"""

import os
import numpy as np
from scipy import stats
from scipy.special import gammaln
from scipy.integrate import simpson, cumulative_trapezoid
from scipy.optimize import brentq

C_KMS = 299792.458
_DATA = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'data')

# ── Planck 2018 (VI, Table 2, TT,TE,EE+lowE+lensing) ─────────────────────────
THETA_STAR_100     = 1.04110      # 100 theta_*  (not theta_MC)
THETA_STAR_100_ERR = 0.00031
Z_STAR   = 1089.92
Z_DRAG   = 1059.94
R_STAR_PLANCK  = 144.43           # Mpc, used only for the calibration below
R_DRAG_PLANCK  = 147.09           # Mpc
H0_PLANCK      = 67.36

# Background densities used throughout the sims
OMEGA_M     = 0.1430
OMEGA_B     = 0.02237
OMEGA_R     = 4.18e-5             # photons + 3.046 massless neutrinos
OMEGA_GAMMA = 2.4728e-5           # photons only, T_cmb = 2.7255 K


# ─────────────────────────────────────────────────────────────────────────────
# Backgrounds
# ─────────────────────────────────────────────────────────────────────────────

def H_lcdm(h, omega_m=OMEGA_M, omega_r=OMEGA_R, w0=-1.0, wa=0.0):
    """Flat (w0, wa)CDM; H(z=0) = 100 h exactly."""
    omega_de = h**2 - omega_m - omega_r

    def H(z):
        z = np.asarray(z, dtype=float)
        a = 1.0 / (1.0 + z)
        f_de = a**(-3.0 * (1.0 + w0 + wa)) * np.exp(-3.0 * wa * (1.0 - a))
        return 100.0 * np.sqrt(omega_m / a**3 + omega_r / a**4 + omega_de * f_de)
    return H


def H_from_table(z_tab, H_tab, omega_m=OMEGA_M, omega_r=OMEGA_R):
    """
    Wrap a tabulated H(z) (km/s/Mpc) as a vectorised callable.

    Inside the table: interpolation of ln H in ln(1+z).
    Above the table : H^2 follows matter+radiation scaling, matched to the table's
                      last point, i.e. the modification of gravity is assumed frozen
                      beyond it. That holds when the scalar is Hubble-frozen in
                      radiation domination (R ~ 0 there, so curvature sourcing
                      switches off).
    Below z = 0 is never needed.
    """
    z_tab = np.asarray(z_tab, dtype=float)
    H_tab = np.asarray(H_tab, dtype=float)
    order = np.argsort(z_tab)
    z_tab, H_tab = z_tab[order], H_tab[order]
    keep = np.concatenate([[True], np.diff(z_tab) > 0])
    z_tab, H_tab = z_tab[keep], H_tab[keep]
    x_tab = np.log1p(z_tab)
    lnH_tab = np.log(H_tab)
    z_top = z_tab[-1]
    E2_top = omega_m * (1 + z_top)**3 + omega_r * (1 + z_top)**4
    H_top = H_tab[-1]

    def H(z):
        z = np.asarray(z, dtype=float)
        out = np.exp(np.interp(np.log1p(z), x_tab, lnH_tab))
        hi = z > z_top
        if np.any(hi):
            zh = z[hi]
            E2 = omega_m * (1 + zh)**3 + omega_r * (1 + zh)**4
            out = np.array(out, dtype=float, copy=True)
            out[hi] = H_top * np.sqrt(E2 / E2_top)
        return out
    return H


# ─────────────────────────────────────────────────────────────────────────────
# Distances
# ─────────────────────────────────────────────────────────────────────────────

def comoving_distance(H_of_z, z, n=4000):
    """D_M(z) = int_0^z c/H dz' for a single z (Mpc), Simpson in ln(1+z)."""
    if z <= 0:
        return 0.0
    x = np.linspace(0.0, np.log1p(z), n + 1)
    zz = np.expm1(x)
    return float(simpson(C_KMS * (1.0 + zz) / H_of_z(zz), x=x))


def _sound_horizon_raw(H_of_z, z, omega_b=OMEGA_B, a_min=1e-9, n=4000):
    """Uncalibrated r_s(z) = int_z^inf c_s/H dz', integrated in ln a from a_min."""
    a_z = 1.0 / (1.0 + z)
    lna = np.linspace(np.log(a_min), np.log(a_z), n + 1)
    a = np.exp(lna)
    R = 3.0 * omega_b * a / (4.0 * OMEGA_GAMMA)
    cs = C_KMS / np.sqrt(3.0 * (1.0 + R))
    Hz = H_of_z(1.0 / a - 1.0)
    return float(simpson(cs / (a * Hz), x=lna))


# Calibration of the simple r_s integral against Planck's Boltzmann-code values.
# The simple integral ignores the finite width of the visibility function and the
# massive-neutrino transition; it lands within ~0.1% of Planck. One fixed
# multiplicative factor per epoch, computed for the Planck 2018 reference ΛCDM
# in exactly the simplified background the sims use, removes that offset.
# Model-dependent changes to early-time H(z) are still captured by the integral.
def _calibration():
    H_ref = H_lcdm(H0_PLANCK / 100.0)
    k_star = R_STAR_PLANCK / _sound_horizon_raw(H_ref, Z_STAR)
    k_drag = R_DRAG_PLANCK / _sound_horizon_raw(H_ref, Z_DRAG)
    return k_star, k_drag

K_STAR, K_DRAG = _calibration()


def sound_horizon(H_of_z, z, omega_b=OMEGA_B, calibrated=True):
    """r_s(z) in Mpc, integrated to z = infinity."""
    rs = _sound_horizon_raw(H_of_z, z, omega_b)
    if not calibrated:
        return rs
    if abs(z - Z_STAR) < 1.0:
        return rs * K_STAR
    if abs(z - Z_DRAG) < 1.0:
        return rs * K_DRAG
    return rs


def r_drag(H_of_z, omega_b=OMEGA_B):
    return sound_horizon(H_of_z, Z_DRAG, omega_b)


def theta_star_100(H_of_z, omega_b=OMEGA_B):
    """100 theta_* = 100 r_s(z*) / D_M(z*)."""
    return 100.0 * sound_horizon(H_of_z, Z_STAR, omega_b) / comoving_distance(H_of_z, Z_STAR)


def chi2_theta(theta100):
    return float(((theta100 - THETA_STAR_100) / THETA_STAR_100_ERR)**2)


def h_from_theta(make_H, h_lo=0.40, h_hi=1.00, omega_b=OMEGA_B, xtol=1e-7):
    """
    Solve 100 theta_*(h) = THETA_STAR_100 for h, where make_H(h) returns H_of_z.
    Returns nan if there is no root in [h_lo, h_hi].
    """
    f = lambda h: theta_star_100(make_H(h), omega_b) - THETA_STAR_100
    try:
        f_lo, f_hi = f(h_lo), f(h_hi)
    except (ValueError, FloatingPointError):
        return np.nan
    if not (np.isfinite(f_lo) and np.isfinite(f_hi)) or f_lo * f_hi > 0:
        return np.nan
    return brentq(f, h_lo, h_hi, xtol=xtol)


# ─────────────────────────────────────────────────────────────────────────────
# DESI BAO likelihood
# ─────────────────────────────────────────────────────────────────────────────

class DesiBAO:
    """Official DESI Gaussian BAO likelihood (DR1 = 2024, DR2 = 2025)."""

    def __init__(self, release='DR1'):
        tag = {'DR1': 'desi_dr1_bao', 'DR2': 'desi_dr2_bao'}[release.upper()]
        self.release = release.upper()
        rows = [ln.split() for ln in open(os.path.join(_DATA, tag + '_mean.txt'))
                if ln.strip() and not ln.startswith('#')]
        self.z = np.array([float(r[0]) for r in rows])
        self.value = np.array([float(r[1]) for r in rows])
        self.kind = [r[2] for r in rows]
        self.cov = np.loadtxt(os.path.join(_DATA, tag + '_cov.txt'))
        self.icov = np.linalg.inv(self.cov)
        self.sigma = np.sqrt(np.diag(self.cov))
        self.n = len(self.z)
        self._zgrid = np.linspace(0.0, max(self.z) * 1.02, 3001)

    def model(self, H_of_z, rd):
        """Model vector in the data's order: D_M/r_d, D_H/r_d, D_V/r_d."""
        Hg = H_of_z(self._zgrid)
        DMg = np.concatenate([[0.0], cumulative_trapezoid(C_KMS / Hg, self._zgrid)])
        DM = np.interp(self.z, self._zgrid, DMg)
        DH = C_KMS / H_of_z(self.z)
        out = np.empty(self.n)
        for i, k in enumerate(self.kind):
            if k == 'DM_over_rs':
                out[i] = DM[i] / rd
            elif k == 'DH_over_rs':
                out[i] = DH[i] / rd
            elif k == 'DV_over_rs':
                out[i] = (self.z[i] * DM[i]**2 * DH[i])**(1.0 / 3.0) / rd
            else:
                raise ValueError(k)
        return out

    def chi2(self, H_of_z, rd):
        d = self.model(H_of_z, rd) - self.value
        return float(d @ self.icov @ d)

    def pulls(self, H_of_z, rd):
        """Per-point (model - data)/sigma; diagonal only, for display."""
        return (self.model(H_of_z, rd) - self.value) / self.sigma

    def table(self, H_of_z, rd):
        m = self.model(H_of_z, rd)
        return [dict(z=float(z), quantity=k, data=float(v), sigma=float(s),
                     model=float(mm), pull=float((mm - v) / s))
                for z, k, v, s, mm in zip(self.z, self.kind, self.value, self.sigma, m)]


# ─────────────────────────────────────────────────────────────────────────────
# Statistics
# ─────────────────────────────────────────────────────────────────────────────

def p_to_sigma(p):
    """Two-sided Gaussian-equivalent significance of a p-value (handles p -> 0)."""
    if p <= 0:
        return np.inf
    if p > 1e-300:
        return float(stats.norm.isf(p / 2.0))
    return np.inf


def chi2_to_sigma(chi2, dof):
    """Goodness of fit: (p-value, Gaussian-equivalent sigma) for chi2 with dof."""
    if not np.isfinite(chi2):
        return np.nan, np.nan
    dof = max(int(dof), 1)
    logp = float(stats.chi2.logsf(chi2, dof))
    if not np.isfinite(logp):
        # large-chi2 asymptote of the chi2 survival function
        logp = (dof / 2.0 - 1.0) * np.log(chi2 / 2.0) - chi2 / 2.0 - gammaln(dof / 2.0)
    if logp > -690:
        p = float(np.exp(logp))
        return p, p_to_sigma(p)
    # asymptotic inversion of log(p/2) = log Q(s), Q(s) ~ phi(s)/s
    target = logp - np.log(2.0)
    s = np.sqrt(-2.0 * target)
    for _ in range(50):
        s = np.sqrt(max(-2.0 * (target + np.log(s * np.sqrt(2.0 * np.pi))), 1e-12))
    return 0.0, float(s)


def delta_chi2_to_sigma(dchi2, ddof):
    """
    Significance of a chi2 improvement dchi2 (> 0 means the alternative fits
    better) from ddof extra parameters (likelihood-ratio test, Wilks).
    Returns 0 if the alternative is not better.
    """
    if not np.isfinite(dchi2) or dchi2 <= 0:
        return 0.0
    return chi2_to_sigma(dchi2, ddof)[1]


def summarize_fit(chi2_desi, n_data, n_params, chi2_th=0.0, label=''):
    """Standard reporting block used by every audited sim."""
    dof = n_data - n_params
    p, sig = chi2_to_sigma(chi2_desi, max(dof, 1))
    return dict(label=label, chi2_DESI=float(chi2_desi), N_DESI=int(n_data),
                n_params=int(n_params), dof=int(dof), p_value=float(p),
                gof_sigma=float(sig), chi2_theta=float(chi2_th),
                rms_pull_old_metric=float(np.sqrt(chi2_desi / n_data)))


# ─────────────────────────────────────────────────────────────────────────────
# ΛCDM reference in the same pipeline (omega_m, omega_b fixed; h from theta*)
# ─────────────────────────────────────────────────────────────────────────────

_REF_CACHE = {}

def lcdm_reference(release='DR1', omega_m=OMEGA_M, omega_r=OMEGA_R, omega_b=OMEGA_B):
    key = (release, omega_m, omega_r, omega_b)
    if key in _REF_CACHE:
        return _REF_CACHE[key]
    lik = DesiBAO(release)
    h = h_from_theta(lambda hh: H_lcdm(hh, omega_m, omega_r), omega_b=omega_b)
    Hf = H_lcdm(h, omega_m, omega_r)
    rd = r_drag(Hf, omega_b)
    c2 = lik.chi2(Hf, rd)
    out = summarize_fit(c2, lik.n, 0, 0.0, label=f'LCDM (Planck omegas, h from theta*), DESI {release}')
    out.update(h=float(h), H0=float(100 * h), r_drag=float(rd),
               theta_star_100=float(theta_star_100(Hf, omega_b)))
    _REF_CACHE[key] = out
    return out
