#!/usr/bin/env python3
# Sim 87 — RIFT BAO Full Covariance Refit
"""
Referee vulnerability fix (Phase 4):

The BAO chi-squared in sim13_4 used diagonal errors only:
  chi2_diag = sum_i  (d_i - m_i)^2 / sigma_i^2

BOSS DR12 and eBOSS DR16 publish off-diagonal DM-DH correlation
coefficients (rho ~ -0.5). These are physically significant: DM and DH
are anti-correlated because both are derived from the same BAO peak
position, so a positive shift in DM is compensated by a negative shift
in DH.  Ignoring this correlation inflates chi2.

This simulation:
  1. Assembles the full 12x12 covariance matrix C from published rho values
  2. Inverts C to form C^{-1}
  3. Computes chi2_full = (d - m)^T C^{-1} (d - m)
  4. Optimises over (H0, Omega_m, r_d, Lambda0) under both chi2 definitions
  5. Reports Delta_chi2 = chi2_diag_best - chi2_full_best

RIFT background model (FLRW, from section3_lagrangian.tex, verified SIM85-86):
  3H^2 = 8pi G_eff [rho_m + rho_r + (1/2)(Psi')^2 H^2
                     + (1/2) m^2 Psi^2 - 6H Lambda'(Psi) Psi' H]
  Psi'' + 3 Psi' + m^2 Psi/H^2 = Lambda'(Psi) R / H^2
  Lambda(Psi) = Lambda0 * Psi^2   [SIM84-verified]
  G_eff(Psi)  = G / (1 + 16pi G Lambda0 Psi^2)

  Coordinates: ln(a). Pi = dPsi/d(lna), Pi' = d^2Psi/d(lna)^2.
  H is solved analytically from the Friedmann constraint at each step
  (see H_from_constraint below).

Data vector ordering (12 elements):
  0:  z=0.15  DV/rd     (6dFGS+MGS)
  1:  z=0.38  DH/rd     (BOSS DR12)      <-- correlated pair
  2:  z=0.38  DM/rd     (BOSS DR12)      <--
  3:  z=0.51  DH/rd     (BOSS DR12)      <-- correlated pair
  4:  z=0.51  DM/rd     (BOSS DR12)      <--
  5:  z=0.70  DH/rd     (eBOSS DR16 LRG) <-- correlated pair
  6:  z=0.70  DM/rd     (eBOSS DR16 LRG) <--
  7:  z=0.85  DV/rd     (eBOSS DR16 QSO)
  8:  z=1.48  DH/rd     (eBOSS DR16 QSO) <-- correlated pair
  9:  z=1.48  DM/rd     (eBOSS DR16 QSO) <--
  10: z=2.33  DH/rd     (Ly-alpha DR16)  <-- correlated pair
  11: z=2.33  DM/rd     (Ly-alpha DR16)  <--

Outputs:
  Outputs/sim87_chi2_surface.png     -- 2D chi2 surface (Omega_m vs Lambda0)
  Outputs/sim87_bao_residuals.png    -- best-fit predictions vs data
  Outputs/sim87_covariance.png       -- covariance matrix visualisation
  Outputs/sim87_diagnostics.json
"""

import os, json, math, warnings
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
from scipy.optimize import minimize
from scipy.linalg import cho_factor, cho_solve

warnings.filterwarnings("ignore", category=RuntimeWarning)

BASE    = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
INPUTS  = os.path.join(BASE, "Inputs")
OUTPUTS = os.path.join(BASE, "Outputs")
PARAMS  = os.path.join(INPUTS, "sim87_params.json")
os.makedirs(OUTPUTS, exist_ok=True)

with open(PARAMS) as f:
    P = json.load(f)

c_km_s   = 299792.458   # speed of light, km/s
G_nat    = 1.0          # G = 1 in natural units (H0 absorbed into rho normalisation)

# ── BAO data vector (12 elements, fixed ordering) ───────────────────────────
# Ordering: see module docstring
BAO_DATA = [
    # (z, quantity, value, sigma, survey)
    (0.15,  "DV_over_rd",  4.47,                    0.17,                   "6dFGS+MGS"),
    (0.38,  "DH_over_rd",  25.0,                    0.76,                   "BOSS DR12"),
    (0.38,  "DM_over_rd",  10.23,                   0.17,                   "BOSS DR12"),
    (0.51,  "DH_over_rd",  22.33,                   0.58,                   "BOSS DR12"),
    (0.51,  "DM_over_rd",  13.36,                   0.21,                   "BOSS DR12"),
    (0.70,  "DH_over_rd",  19.33,                   0.53,                   "eBOSS DR16 LRG"),
    (0.70,  "DM_over_rd",  17.86,                   0.33,                   "eBOSS DR16 LRG"),
    (0.85,  "DV_over_rd",  18.33,                   0.595,                  "eBOSS DR16 ELG"),
    (1.48,  "DH_over_rd",  13.26,                   0.55,                   "eBOSS DR16 QSO"),
    (1.48,  "DM_over_rd",  30.69,                   0.80,                   "eBOSS DR16 QSO"),
    (2.33,  "DH_over_rd",  8.990618556701030,        0.21614046597277392,   "Ly-alpha DR16"),
    (2.33,  "DM_over_rd",  37.433384615384625,       1.26691023299267,      "Ly-alpha DR16"),
]
N_DATA = len(BAO_DATA)  # 12
z_obs   = np.array([d[0] for d in BAO_DATA])
d_obs   = np.array([d[2] for d in BAO_DATA])
sig_obs = np.array([d[3] for d in BAO_DATA])
qty_obs = [d[1] for d in BAO_DATA]

# ── Full covariance matrix ───────────────────────────────────────────────────
# Published DM-DH correlation coefficients at paired redshifts.
# All other pairs assumed uncorrelated (cross-survey, different z).
#
# Data vector indices of paired (DH, DM) at each z:
#   z=0.38: (1, 2), rho=-0.52  Alam+2017 Table 4
#   z=0.51: (3, 4), rho=-0.47  Alam+2017 Table 4
#   z=0.70: (5, 6), rho=-0.48  Alam+2021 Table 3
#   z=1.48: (8, 9), rho=-0.46  Alam+2021 Table 3
#   z=2.33: (10,11),rho=-0.43  Bourboux+2020 Table 3

RHO_DM_DH = {
    0.38: -0.52,
    0.51: -0.47,
    0.70: -0.48,
    1.48: -0.46,
    2.33: -0.43,
}

# Paired index map: z -> (i_DH, i_DM) in BAO_DATA list
PAIR_IDX = {}
for i, (z, qty, *_) in enumerate(BAO_DATA):
    z = float(z)
    if z not in PAIR_IDX:
        PAIR_IDX[z] = {}
    if "DH" in qty:
        PAIR_IDX[z]["DH"] = i
    elif "DM" in qty:
        PAIR_IDX[z]["DM"] = i

def build_covariance_matrix():
    """
    Assemble the full 12x12 covariance matrix C.
    Diagonal: sigma_i^2.
    Off-diagonal: C[i_DH, i_DM] = C[i_DM, i_DH] = rho * sigma_DH * sigma_DM
    for paired (DH, DM) measurements at the same redshift.
    """
    C = np.diag(sig_obs**2)
    for z_key, rho in RHO_DM_DH.items():
        pair = PAIR_IDX.get(z_key, {})
        if "DH" in pair and "DM" in pair:
            i_dh = pair["DH"]
            i_dm = pair["DM"]
            cov_off = rho * sig_obs[i_dh] * sig_obs[i_dm]
            C[i_dh, i_dm] = cov_off
            C[i_dm, i_dh] = cov_off
    return C

C_full = build_covariance_matrix()
# Cholesky decomposition for efficient repeated solve
C_cho  = cho_factor(C_full)
C_inv  = np.linalg.inv(C_full)   # also keep explicit inverse for diagnostics

def chi2_diagonal(residuals):
    """chi2 = sum_i (d_i - m_i)^2 / sigma_i^2"""
    return float(np.sum((residuals / sig_obs)**2))

def chi2_full_cov(residuals):
    """chi2 = (d - m)^T C^{-1} (d - m)"""
    v = cho_solve(C_cho, residuals)
    return float(np.dot(residuals, v))

# ── RIFT background integrator ───────────────────────────────────────────────
# Units: G = 1, H0 in km/s/Mpc folded into rho_crit via natural units
# We work in dimensionless units: define rho in units of rho_crit0 = 3H0^2/(8piG)
# => E^2(a) = H^2(a)/H0^2 is dimensionless
# Density parameters: rho_x(a) / rho_crit0 = Omega_x * a^{-3(1+w_x)}
#
# RIFT field: Psi dimensionless, Lambda0 dimensionless,
# m0 in units of H0 (so m_phys = m0 * H0)

def integrate_rift_background(H0, Omega_m, Lambda0, Psi_ini=0.01,
                               m0=1.0, alpha=0.1, beta=0.05,
                               Omega_b=0.049, Omega_r=9.2e-5,
                               a_ini=1e-5, a_fin=1.0, N=5000):
    """
    Integrate the RIFT FLRW system in ln(a) coordinates.

    State: y = [Psi, Pi]  where Pi = dPsi/d(lna)

    Friedmann constraint (Eq. friedmann from section3_lagrangian.tex):
      3H^2 = 8pi G_eff [ rho_m + rho_r
                         + (1/2)(Pi*H)^2
                         + (1/2) m_eff^2 Psi^2
                         - 6H * Lambda'(Psi) * Pi*H ]
    with G_eff = G / (1 + 16pi G Lambda0 Psi^2)
    and Lambda'(Psi) = 2 Lambda0 Psi.

    Solving for H^2 analytically:
      Let A = 8pi G_eff (dimensionless prefactor)
          B = A * (rho_m + rho_r + (1/2) m_eff^2 Psi^2)
          Q = A * ((1/2) Pi^2 - 12 Lambda0 Psi Pi)
      Then: H^2 (3 - Q) = B => H = sqrt(B / (3 - Q))

    Scalar field equation in ln(a):
      Pi' = -3 Pi - m_eff^2 Psi / H^2 + Lambda'(Psi) R / H^2
    where R = 6 (H dH/dlna + 2H^2) = 6 H^2 (dH/(H dlna) + 2).

    Returns: (a_arr, H_arr, E_arr) with E = H/H0
    """
    # Omega_cdm = Omega_m - Omega_b  (Omega_m includes baryons)
    Omega_L = 1.0 - Omega_m - Omega_r   # flat universe, RIFT field is extra

    # rho_crit0 in natural units: 3*H0^2 / (8*pi*G) = 3/(8*pi) in G=H0=1 units
    # Since we work with E^2 = H^2/H0^2, we set H0=1 and restore via km/s/Mpc at end.
    # All densities normalised to rho_crit0 = 3/(8*pi) in these units.
    rho_crit0_nat = 3.0 / (8.0 * math.pi)   # = 3H0^2/(8piG) with G=H0=1

    m0_phys = m0   # in units of H0

    def m_eff_sq(Psi):
        return m0_phys**2 * (1.0 + alpha * Psi**2 * math.exp(-beta * Psi**2))

    def H_from_friedmann(lna, Psi, Pi):
        """Solve Friedmann constraint analytically for E = H/H0 (H0=1 natural units).

        Derivation (section3_lagrangian.tex Eq. friedmann):
          3H^2 = 8pi G_eff [rho_m + rho_r + rho_L
                             + (1/2)(Pi H)^2 + (1/2)m^2 Psi^2
                             - 6H * Lambda'(Psi) * Pi * H]

        Physical densities (G=H0=1, rho_crit0 = 3/(8pi)):
          rho_x_phys = Omega_x * 3/(8pi) / a^n

        Substituting and collecting H^2 terms:
          3H^2 = G_eff * [3*Omega_bg + 4pi*Pi^2*H^2 + 4pi*m^2*Psi^2 - 48pi*H^2*Lambda'*Pi]
          H^2 * [3 - G_eff*(4pi*Pi^2 - 48pi*Lambda'*Pi)]
              = G_eff * [3*Omega_bg + 4pi*m^2*Psi^2]

        At Psi=0, Pi=0: H^2 = Omega_bg = Omega_m/a^3 + Omega_r/a^4 + Omega_L  (LCDM) ✓
        """
        a   = math.exp(lna)
        Lam = Lambda0 * Psi**2
        dLam_dPsi = 2.0 * Lambda0 * Psi
        Geff = 1.0 / (1.0 + 16.0 * math.pi * Lam)   # G_eff/G, G=1

        Omega_bg = Omega_m / a**3 + Omega_r / a**4 + Omega_L
        m2 = m_eff_sq(Psi)

        numerator = Geff * (3.0 * Omega_bg + 4.0 * math.pi * m2 * Psi**2)
        denom     = 3.0 - Geff * (4.0 * math.pi * Pi**2 - 48.0 * math.pi * dLam_dPsi * Pi)

        if denom <= 1e-10:
            return 1e-30
        H2 = numerator / denom
        if H2 <= 0:
            return 1e-30
        return math.sqrt(H2)   # E = H/H0 (H0=1 natural units)

    def rhs(lna, y):
        Psi, Pi = float(y[0]), float(y[1])
        H = H_from_friedmann(lna, Psi, Pi)
        if H < 1e-30:
            return [Pi, 0.0]

        m2 = m_eff_sq(Psi)
        dLam_dPsi = 2.0 * Lambda0 * Psi

        # Need dH/dlna for R; approximate using current H and radiation-era slope
        # dH/dlna ~ -2H in radiation era, ~ -3/2 H in matter era
        # We use the exact expression: R/H^2 = 6*(dH/(H*dlna) + 2)
        # For simplicity, compute R = 6*(2H^2 + H*dH/dlna) using
        # the equation of state: dH/dlna = H * (d lnH / d lna)
        # dE^2/dlna = -3 Omega_m a^{-3} - 4 Omega_r a^{-4}  [dominant at high z]
        # => dE/dlna = (dE^2/dlna) / (2E)
        a = math.exp(lna)
        dE2_dlna = -3.0 * Omega_m / a**3 - 4.0 * Omega_r / a**4
        dH_dlna  = dE2_dlna / (2.0 * max(H, 1e-30))

        R = 6.0 * (H * dH_dlna + 2.0 * H**2)

        dPsi_dlna  = Pi
        dPi_dlna   = (-3.0 * Pi
                      - m2 * Psi / H**2
                      + dLam_dPsi * R / H**2)
        return [dPsi_dlna, dPi_dlna]

    # Initial conditions at a_ini
    lna_ini = math.log(a_ini)
    lna_fin = math.log(a_fin)
    lna_span = (lna_ini, lna_fin)

    H_ini = H_from_friedmann(lna_ini, Psi_ini, 0.0)
    y0 = [Psi_ini, 0.0]

    lna_eval = np.linspace(lna_ini, lna_fin, N)

    sol = solve_ivp(rhs, lna_span, y0, method="RK45",
                    t_eval=lna_eval, rtol=1e-9, atol=1e-12,
                    dense_output=False)

    a_arr = np.exp(sol.t)
    Psi_arr = sol.y[0]

    # Recompute H at each stored point
    H_arr = np.array([
        H_from_friedmann(float(sol.t[i]), float(sol.y[0, i]), float(sol.y[1, i]))
        for i in range(len(sol.t))
    ])

    # Convert E -> H in km/s/Mpc
    H_arr_kms = H_arr * H0   # E * H0

    return a_arr, H_arr_kms, Psi_arr


def compute_distances(a_arr, H_arr_kms, rd_Mpc, z_query):
    """
    Compute D_M/rd, D_H/rd, D_V/rd at requested redshifts.

    H_arr_kms: H(a) in km/s/Mpc
    Uses trapezoidal integration for D_C(z) = c * integral_0^z dz'/H(z').
    """
    z_arr = 1.0 / a_arr - 1.0
    # Sort in ascending z
    idx   = np.argsort(z_arr)
    z_s   = z_arr[idx]
    H_s   = H_arr_kms[idx]

    out = {}
    for z in z_query:
        z = float(z)
        # D_H = c / H(z)
        H_z = float(np.interp(z, z_s, H_s))
        D_H = c_km_s / max(H_z, 1e-30)

        # D_C = c * integral_0^z dz'/H(z')
        m = z_s <= z
        z_use = z_s[m]
        H_use = H_s[m]
        if z_use[-1] < z:
            # append endpoint by linear interpolation
            z1, z2 = z_s[m][-1], z_s[~m][0] if np.any(~m) else z
            H1 = float(np.interp(z1, z_s, H_s))
            H2 = float(np.interp(z2, z_s, H_s))
            t  = (z - z1) / max(z2 - z1, 1e-30)
            H_end = H1 + t * (H2 - H1)
            z_use = np.append(z_use, z)
            H_use = np.append(H_use, H_end)

        integrand = 1.0 / np.maximum(H_use, 1e-30)
        D_C = c_km_s * float(np.trapezoid(integrand, z_use))

        D_M = D_C   # flat universe: D_M = D_C
        D_V = (z * D_H * D_M**2)**(1.0 / 3.0) if z > 0 else 0.0

        out[z] = {
            "DM_over_rd": D_M / rd_Mpc,
            "DH_over_rd": D_H / rd_Mpc,
            "DV_over_rd": D_V / rd_Mpc,
        }
    return out


def predict_bao(H0, Omega_m, r_d, Lambda0,
                Psi_ini=0.01, m0=1.0, alpha=0.1, beta=0.05,
                Omega_b=0.049, Omega_r=9.2e-5):
    """
    Run the RIFT background and return model predictions for the BAO data vector.
    Returns a 12-element numpy array aligned with BAO_DATA ordering.
    Returns None if integration failed.
    """
    try:
        a_arr, H_arr, _ = integrate_rift_background(
            H0=H0, Omega_m=Omega_m, Lambda0=Lambda0,
            Psi_ini=Psi_ini, m0=m0, alpha=alpha, beta=beta,
            Omega_b=Omega_b, Omega_r=Omega_r, a_ini=1e-5, a_fin=1.0, N=3000
        )
        z_unique = sorted(set(float(d[0]) for d in BAO_DATA))
        dist = compute_distances(a_arr, H_arr, r_d, z_unique)

        pred = np.zeros(N_DATA)
        for i, (z, qty, *_) in enumerate(BAO_DATA):
            z = float(z)
            pred[i] = dist[z][qty]
        return pred
    except Exception:
        return None


def _in_bounds(H0, Om, rd, L0):
    """
    Physically motivated bounds for BAO fit.
    H0: Planck ±10 km/s/Mpc; rd: BBN+CMB range ±5 Mpc; Omega_m: broad flat prior.
    BAO fits H0*rd as a combined parameter; independent H0 and rd are weakly constrained
    without CMB, so we impose narrow rd prior (r_d = 147.1 ± 5 Mpc from Planck+BBN)
    and moderate H0 range to avoid the H0-rd degeneracy collapsing the optimizer.
    """
    return (64.0 <= H0 <= 72.0 and
            0.26 <= Om <= 0.42 and
            142.0 <= rd <= 153.0 and
            0.0 <= L0 <= 0.025)


def objective_diagonal(params):
    H0, Om, rd, L0 = params
    if not _in_bounds(H0, Om, rd, L0):
        return 1e6
    pred = predict_bao(H0, Om, rd, L0)
    if pred is None:
        return 1e6
    return chi2_diagonal(d_obs - pred)


def objective_full(params):
    H0, Om, rd, L0 = params
    if not _in_bounds(H0, Om, rd, L0):
        return 1e6
    pred = predict_bao(H0, Om, rd, L0)
    if pred is None:
        return 1e6
    return chi2_full_cov(d_obs - pred)


# ── Parameter optimization ───────────────────────────────────────────────────
opt_cfg = P["optimization"]
x0      = opt_cfg["x0"]   # [H0, Omega_m, r_d, Lambda0]

print("  Optimizing chi2_diagonal ...", flush=True)
res_diag = minimize(objective_diagonal, x0, method="Nelder-Mead",
                    options={"maxiter": 3000, "xatol": 1e-5, "fatol": 1e-5})
best_diag = res_diag.x
chi2_diag_best = float(res_diag.fun)
print(f"    Best: H0={best_diag[0]:.2f}, Om={best_diag[1]:.4f}, "
      f"rd={best_diag[2]:.2f}, L0={best_diag[3]:.5f}  chi2={chi2_diag_best:.4f}")

print("  Optimizing chi2_full_cov ...", flush=True)
res_full = minimize(objective_full, x0, method="Nelder-Mead",
                    options={"maxiter": 3000, "xatol": 1e-5, "fatol": 1e-5})
best_full = res_full.x
chi2_full_best = float(res_full.fun)
print(f"    Best: H0={best_full[0]:.2f}, Om={best_full[1]:.4f}, "
      f"rd={best_full[2]:.2f}, L0={best_full[3]:.5f}  chi2={chi2_full_best:.4f}")

# Also evaluate old sim13_4 best-fit (diagonal, fixed params) for comparison
om023_params = [67.4, 0.265 + 0.1, 147.1, 0.003]  # Omega_m = cdm+b = 0.265+0.049≈0.314 but sim used Omega_psi=0.1
pred_fiducial = predict_bao(*[67.4, 0.314, 147.1, 0.003])
chi2_fiducial_diag = chi2_diagonal(d_obs - pred_fiducial) if pred_fiducial is not None else float("nan")
chi2_fiducial_full = chi2_full_cov(d_obs - pred_fiducial) if pred_fiducial is not None else float("nan")

n_dof = N_DATA - 4   # 12 - 4 free params = 8
print(f"\n  ── Summary ──")
print(f"  Old fiducial (sim13_4-like): chi2_diag={chi2_fiducial_diag:.3f}  "
      f"chi2_full={chi2_fiducial_full:.3f}  (n_data={N_DATA}, dof={n_dof})")
print(f"  Best-fit diagonal:           chi2={chi2_diag_best:.3f}  chi2/dof={chi2_diag_best/n_dof:.3f}")
print(f"  Best-fit full-cov:           chi2={chi2_full_best:.3f}  chi2/dof={chi2_full_best/n_dof:.3f}")
print(f"  Delta_chi2 (diag-full):      {chi2_diag_best - chi2_full_best:.3f}")

# ── 2D chi2 surface ──────────────────────────────────────────────────────────
scan_cfg = P["scan"]
Om_grid  = np.linspace(scan_cfg["Omega_m_grid"]["min"],
                       scan_cfg["Omega_m_grid"]["max"],
                       scan_cfg["Omega_m_grid"]["N"])
L0_grid  = np.linspace(scan_cfg["Lambda0_grid"]["min"],
                       scan_cfg["Lambda0_grid"]["max"],
                       scan_cfg["Lambda0_grid"]["N"])

# Fix H0 and r_d at best-fit full-cov values for the surface
H0_fix = float(best_full[0])
rd_fix = float(best_full[2])

print(f"\n  Computing 2D chi2 surface ({len(Om_grid)}x{len(L0_grid)} grid) ...", flush=True)
surf_diag = np.full((len(L0_grid), len(Om_grid)), np.nan)
surf_full = np.full((len(L0_grid), len(Om_grid)), np.nan)

for j, Om in enumerate(Om_grid):
    for i, L0 in enumerate(L0_grid):
        pred = predict_bao(H0_fix, Om, rd_fix, L0)
        if pred is not None:
            r = d_obs - pred
            surf_diag[i, j] = chi2_diagonal(r)
            surf_full[i, j] = chi2_full_cov(r)
    if j % 5 == 0:
        print(f"    Om column {j+1}/{len(Om_grid)}", flush=True)

# ── Best-fit residual plot ───────────────────────────────────────────────────
pred_bf = predict_bao(*best_full)
resid_bf = d_obs - pred_bf if pred_bf is not None else np.zeros(N_DATA)

# ── Plots ────────────────────────────────────────────────────────────────────

# Plot 1: Covariance matrix
fig, ax = plt.subplots(figsize=(6, 5), dpi=140)
labels = [f"z={d[0]}\n{d[1].replace('_over_rd','')}" for d in BAO_DATA]
norm_cov = C_inv * np.outer(sig_obs, sig_obs)   # C^{-1} normalised
im = ax.imshow(C_full / np.outer(sig_obs, sig_obs),
               cmap="RdBu_r", vmin=-1, vmax=1, aspect="auto")
ax.set_xticks(range(N_DATA)); ax.set_yticks(range(N_DATA))
ax.set_xticklabels(labels, rotation=90, fontsize=5)
ax.set_yticklabels(labels, fontsize=5)
ax.set_title("BAO covariance matrix (normalised by σ_i σ_j)", fontsize=9)
plt.colorbar(im, ax=ax, label=r"$C_{ij}/(\sigma_i\sigma_j)$")
plt.tight_layout()
p_cov = os.path.join(OUTPUTS, "sim87_covariance.png")
fig.savefig(p_cov); plt.close(fig)

# Plot 2: 2D chi2 surface
fig, axes = plt.subplots(1, 2, figsize=(12, 5), dpi=130)
for ax, surf, label in zip(axes, [surf_diag, surf_full],
                           ["diagonal $\\chi^2$", "full-cov $\\chi^2$"]):
    vmin = np.nanmin(surf)
    vmax = min(np.nanmax(surf), vmin + 20)
    im = ax.contourf(Om_grid, L0_grid, surf, levels=20,
                     cmap="plasma_r", vmin=vmin, vmax=vmax)
    ax.contour(Om_grid, L0_grid, surf, levels=[vmin + 2.3, vmin + 6.17],
               colors=["white", "cyan"], linewidths=[1.2, 0.8])
    plt.colorbar(im, ax=ax, label=label)
    ax.set_xlabel(r"$\Omega_m$", fontsize=11)
    ax.set_ylabel(r"$\Lambda_0$", fontsize=11)
    ax.set_title(f"Sim 87 — {label} surface\n(H0={H0_fix:.1f}, rd={rd_fix:.1f} Mpc)", fontsize=9)
    # Mark best-fit
    ax.axvline(best_full[1], color="lime", lw=0.8, ls="--")
    ax.axhline(best_full[3], color="lime", lw=0.8, ls="--")
plt.tight_layout()
p_surf = os.path.join(OUTPUTS, "sim87_chi2_surface.png")
fig.savefig(p_surf); plt.close(fig)

# Plot 3: BAO residuals
fig, axes = plt.subplots(1, 2, figsize=(13, 5), dpi=130)
for ax, pred_arr, title in zip(axes,
        [predict_bao(*best_diag), predict_bao(*best_full)],
        ["Best-fit diagonal $\\chi^2$", "Best-fit full-cov $\\chi^2$"]):
    if pred_arr is None:
        continue
    r = (d_obs - pred_arr) / sig_obs
    colors = plt.cm.tab10(np.linspace(0, 1, N_DATA))
    ax.bar(range(N_DATA), r, color=colors, edgecolor="k", linewidth=0.5)
    ax.axhline(0, color="k", lw=0.8)
    ax.axhline(+1, color="gray", lw=0.5, ls="--")
    ax.axhline(-1, color="gray", lw=0.5, ls="--")
    ax.set_xticks(range(N_DATA))
    ax.set_xticklabels([f"z={d[0]}\n{d[1].replace('_over_rd','')}" for d in BAO_DATA],
                       rotation=90, fontsize=6)
    ax.set_ylabel(r"$(d_i - m_i)/\sigma_i$", fontsize=10)
    ax.set_title(f"Sim 87 — {title}", fontsize=9)
plt.tight_layout()
p_resid = os.path.join(OUTPUTS, "sim87_bao_residuals.png")
fig.savefig(p_resid); plt.close(fig)

# ── Acceptance check ─────────────────────────────────────────────────────────
chi2_dof_full = chi2_full_best / n_dof
chi2_dof_max  = float(P["acceptance"]["chi2_full_per_dof_max"])
passed = bool(chi2_dof_full < chi2_dof_max)
if not passed:
    raise RuntimeError(
        f"BAO full-cov chi2/dof = {chi2_dof_full:.3f} exceeds threshold {chi2_dof_max}"
    )

# ── Diagnostics ──────────────────────────────────────────────────────────────
diag_out = {
    "description": (
        "BAO full-covariance refit. Fixes referee vulnerability: prior chi2 "
        "was diagonal-only. Off-diagonal DM-DH correlations from BOSS DR12 "
        "(Alam+2017) and eBOSS DR16 (Alam+2021, Bourboux+2020) now included."
    ),
    "n_data": N_DATA,
    "n_dof": n_dof,
    "covariance_matrix": {
        "type": "12x12 block-diagonal by survey",
        "off_diagonal_pairs": [
            {"z": z, "rho_DM_DH": rho, "i_DH": PAIR_IDX[z]["DH"], "i_DM": PAIR_IDX[z]["DM"]}
            for z, rho in RHO_DM_DH.items()
        ],
        "C_inv_condition_number": float(np.linalg.cond(C_full)),
    },
    "fiducial_sim13_4": {
        "params": {"H0": 67.4, "Omega_m": 0.314, "r_d": 147.1, "Lambda0": 0.003},
        "chi2_diagonal": float(chi2_fiducial_diag),
        "chi2_full_cov": float(chi2_fiducial_full),
        "delta_chi2": float(chi2_fiducial_full - chi2_fiducial_diag),
    },
    "best_fit_diagonal": {
        "params": {"H0": float(best_diag[0]), "Omega_m": float(best_diag[1]),
                   "r_d": float(best_diag[2]), "Lambda0": float(best_diag[3])},
        "chi2": chi2_diag_best,
        "chi2_per_dof": float(chi2_diag_best / n_dof),
        "converged": bool(res_diag.success),
    },
    "best_fit_full_cov": {
        "params": {"H0": float(best_full[0]), "Omega_m": float(best_full[1]),
                   "r_d": float(best_full[2]), "Lambda0": float(best_full[3])},
        "chi2": chi2_full_best,
        "chi2_per_dof": float(chi2_full_best / n_dof),
        "converged": bool(res_full.success),
    },
    "comparison": {
        "delta_chi2_best": float(chi2_diag_best - chi2_full_best),
        "interpretation": (
            "Positive delta_chi2 means diagonal chi2 overestimates goodness-of-fit "
            "when DM-DH residuals have opposite sign (expected for typical RIFT vs LCDM). "
            "Full-cov chi2 is the statistically correct quantity for referee submission."
        ),
    },
    "passed": passed,
    "chi2_full_per_dof": float(chi2_dof_full),
    "artifacts": {
        "covariance_matrix": os.path.relpath(p_cov, BASE),
        "chi2_surface":      os.path.relpath(p_surf, BASE),
        "bao_residuals":     os.path.relpath(p_resid, BASE),
    },
}

diag_path = os.path.join(OUTPUTS, "sim87_diagnostics.json")
with open(diag_path, "w") as f:
    json.dump(diag_out, f, indent=2)

print(f"\n  chi2_full/dof = {chi2_dof_full:.3f}  PASS" if passed else
      f"\n  chi2_full/dof = {chi2_dof_full:.3f}  FAIL")
print(f"  Delta_chi2 (diag_best - full_best) = {chi2_diag_best - chi2_full_best:+.3f}")
print(f"\nWrote diagnostics to {diag_path}")
