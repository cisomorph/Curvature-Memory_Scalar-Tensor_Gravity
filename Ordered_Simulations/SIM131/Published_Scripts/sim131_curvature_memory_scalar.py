"""
SIM131 — CMSTG Phase 3: Curvature-Memory Scalar (Conformal Coupling)
====================================================================
Phase 3 opens a first-principles redesign of the CMSTG action. Direction D:
replace V(Ψ) with a linear curvature-memory coupling ξΨR.

Motivation:
  - Phase 1/2 DESI tension (2.77σ) is structural — all parameter variation fails
  - Root cause: Ψ frozen at hilltop; SIM123/130 confirm kinematic freeze
  - Linear ξΨR term sources Ψ from R even at Ψ=0 (breaks the frozen attractor)
  - Ψ = □⁻¹R: causal integral of curvature over past light cone — spacetime memory
  - R ∝ T_trace (trace of energy-momentum) — field tied to energy history
  - ξ = 1/6 fixed by conformal symmetry in 4D: no new free parameter
  - Removes V(Ψ) (eliminates λ, v): net -1 parameter vs Phase 1

Phase 3 action (Jordan frame):
  S = ∫d⁴x√(-g) [ F_eff(Ψ)·R − ½(∂Ψ)² ] + S_SM

  F_eff(Ψ) = (1 + 2Λ₀Ψ²)/2 + ξΨ

  (The ξΨR and (1+2Λ₀Ψ²)/2 · R terms combine into F_eff·R)

Equations of motion (FRW, dN = d ln a):

  Field:
    Ψ'' + (3−ε_H)Ψ' = (2Λ₀Ψ + ξ) × R/H²
    R/H² = 6(2−ε_H)

  Friedmann:
    H²[3F_eff + 6(2Λ₀Ψ+ξ)Ψ' − ½Ψ'²] = ω_m a⁻³ + ω_r a⁻⁴ + Λ_bare

  where F_eff(Ψ) = (1+2Λ₀Ψ²)/2 + ξΨ

Key IC change: Ψ_ini = 0 (not tuned from Planck epoch)
  — ξ sources Ψ from curvature; no tuning required.

SIM131 structure:
  A. Reference run: ξ=1/6, Λ₀=0.003, Ψ_ini=0 → trajectory + DESI test
  B. Λ₀ scan: ξ=1/6 fixed, Λ₀ ∈ [0.001, 0.003, 0.005, 0.010]
  C. ξ scan: Λ₀=0.003 fixed, ξ ∈ [0.05, 1/12, 1/6, 1/4, 1/3]
  D. Best joint DESI+CMB fit — verdict

Baseline: SIM121C (frozen Phase 1), DESI tension = 2.77σ
PASS: DESI tension < 2.0σ AND |100θ* − 1.04101| < 2×0.00029
PARTIAL: tension < 2.77σ OR significant improvement in Δχ²
"""

import numpy as np
from scipy.integrate import quad, solve_ivp
import json, os, warnings
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
warnings.filterwarnings('ignore')

OUT = os.path.join(os.path.dirname(__file__), '..', 'Outputs')
os.makedirs(OUT, exist_ok=True)

plt.rcParams.update({'font.family': 'serif', 'font.size': 11,
                     'axes.labelsize': 12, 'legend.fontsize': 10})

# ─────────────────────────────────────────────────────────────────────────────
# CONSTANTS (same units as SIM130: H in H100, densities in H100² = ω units)
# ─────────────────────────────────────────────────────────────────────────────
c_kms    = 2.998e5
H100     = 100.0
omh2_m   = 0.1430       # Planck Ω_m h²
omh2_r   = 4.18e-5      # Planck Ω_r h²
h_target = 0.674        # H₀/H100

# CMB
theta_obs     = 1.04101
theta_obs_err = 0.00029
z_drag        = 1059.6
z_star        = 1089.8
omh2_b        = 0.02237

# DESI H(z) data [km/s/Mpc]
DESI_z = np.array([0.295, 0.510, 0.706, 0.930, 1.317, 2.330])
DESI_H = np.array([ 81.7,  97.9, 110.7, 128.1, 156.4, 240.8])
DESI_s = np.array([  4.5,   4.4,   6.2,   5.6,   8.6,  11.0])

# Baseline
SIM121C_chi2    = 18.26
SIM121C_tension = 2.77
PSI0_phase1     = 2.62
F0_phase1       = 0.5 + 0.003 * PSI0_phase1**2   # Phase 1 frozen F₀

# Phase 3 fixed parameter
XI_CONFORMAL = 1.0 / 6.0   # conformal coupling — fixed by symmetry


# ─────────────────────────────────────────────────────────────────────────────
# CORE ODE SYSTEM
# ─────────────────────────────────────────────────────────────────────────────

def F_eff(Psi, xi, Lambda0):
    """Effective gravitational coupling: (1+2Λ₀Ψ²)/2 + ξΨ"""
    return 0.5 + Lambda0 * Psi**2 + xi * Psi

def F_eff_prime(Psi, xi, Lambda0):
    """dF_eff/dΨ = 2Λ₀Ψ + ξ"""
    return 2.0 * Lambda0 * Psi + xi

def get_H2(Psi, y, N, Lambda_bare, xi, Lambda0):
    """
    Modified Friedmann. Returns H² in H100² units.
    H²[3F_eff + 6(2Λ₀Ψ+ξ)Ψ' − ½Ψ'²] = ω_m a⁻³ + ω_r a⁻⁴ + Λ_bare
    """
    a    = np.exp(N)
    rhs  = omh2_m / a**3 + omh2_r / a**4 + Lambda_bare
    Fv   = F_eff(Psi, xi, Lambda0)
    Fp   = F_eff_prime(Psi, xi, Lambda0)
    coef = 3.0 * Fv + 6.0 * Fp * y - 0.5 * y**2
    if coef < 0.3:
        coef = 3.0 * Fv   # stability fallback
    return rhs / coef

def get_eps_H(Psi, y, N, Lambda_bare, xi, Lambda0):
    """ε_H = −½ d ln H²/dN"""
    eps = 5e-4
    H2p = get_H2(Psi, y, N + eps, Lambda_bare, xi, Lambda0)
    H2m = get_H2(Psi, y, N - eps, Lambda_bare, xi, Lambda0)
    H2  = get_H2(Psi, y, N, Lambda_bare, xi, Lambda0)
    if H2 < 1e-40:
        return 0.0
    return -0.5 * (H2p - H2m) / (2.0 * eps * H2)

def phase3_ode(N, state, Lambda_bare, xi, Lambda0):
    """
    Phase 3 CMSTG ODE. State: [Ψ, y=dΨ/dN]
    Source = (2Λ₀Ψ + ξ) × R/H²  where R/H² = 6(2−ε_H)
    """
    Psi, y   = state
    H2       = get_H2(Psi, y, N, Lambda_bare, xi, Lambda0)
    eps_H    = get_eps_H(Psi, y, N, Lambda_bare, xi, Lambda0)
    R_over_H2 = 6.0 * (2.0 - eps_H)
    source   = F_eff_prime(Psi, xi, Lambda0) * R_over_H2
    dy_dN    = source - (3.0 - eps_H) * y
    return [y, dy_dN]


# ─────────────────────────────────────────────────────────────────────────────
# INTEGRATION + SELF-CONSISTENT Λ_bare CALIBRATION
# ─────────────────────────────────────────────────────────────────────────────

def integrate_phase3(psi_ini, Lambda_bare, xi, Lambda0, z_ini=1e5, n_points=4000):
    """Integrate Phase 3 background from z_ini to z=0."""
    N_ini  = np.log(1.0 / (1.0 + z_ini))
    N_end  = 0.0
    N_eval = np.linspace(N_ini, N_end, n_points)

    sol = solve_ivp(
        phase3_ode,
        (N_ini, N_end),
        [psi_ini, 0.0],
        args=(Lambda_bare, xi, Lambda0),
        t_eval=N_eval,
        method='RK45',
        rtol=1e-8, atol=1e-12,
        max_step=0.05
    )

    if not sol.success:
        return None

    N_arr   = sol.t
    z_arr   = np.exp(-N_arr) - 1.0
    psi_arr = sol.y[0]
    y_arr   = sol.y[1]

    H_arr, Feff_arr = [], []
    for i in range(len(N_arr)):
        H2  = get_H2(psi_arr[i], y_arr[i], N_arr[i], Lambda_bare, xi, Lambda0)
        Fv  = F_eff(psi_arr[i], xi, Lambda0)
        H_arr.append(H100 * np.sqrt(max(H2, 0.0)))
        Feff_arr.append(Fv)

    H_arr    = np.array(H_arr)
    Feff_arr = np.array(Feff_arr)

    idx = np.argsort(z_arr)
    return {
        'z':    z_arr[idx],
        'N':    N_arr[idx],
        'psi':  psi_arr[idx],
        'y':    y_arr[idx],
        'H':    H_arr[idx],
        'Feff': Feff_arr[idx],
        'psi0': float(psi_arr[-1]),
        'y0':   float(y_arr[-1]),
        'Feff0': float(Feff_arr[-1]),
        'H0':   float(H_arr[-1]),
        'Lambda_bare': Lambda_bare,
        'xi': xi, 'Lambda0': Lambda0,
    }


def calibrate_Lambda_bare(psi0, y0, xi, Lambda0):
    """Compute Λ_bare from today's (Ψ₀, Ψ'₀) to enforce H₀ = h_target."""
    Fv   = F_eff(psi0, xi, Lambda0)
    Fp   = F_eff_prime(psi0, xi, Lambda0)
    coef = 3.0 * Fv + 6.0 * Fp * y0 - 0.5 * y0**2
    return h_target**2 * coef - omh2_m - omh2_r


def integrate_selfconsistent(psi_ini, xi, Lambda0, z_ini=1e5, n_iter=6):
    """
    Self-consistent: iterate Λ_bare until H₀ = h_target.
    Start from Phase 1 Λ_bare as first guess.
    """
    # Initial guess: Phase 1 reference
    Lambda_bare = 3.0 * F0_phase1 * h_target**2 - omh2_m - omh2_r

    bg = None
    for k in range(n_iter):
        bg = integrate_phase3(psi_ini, Lambda_bare, xi, Lambda0, z_ini=z_ini)
        if bg is None:
            return None
        psi0 = bg['psi0']
        y0   = bg['y0']
        Lambda_bare_new = calibrate_Lambda_bare(psi0, y0, xi, Lambda0)
        if abs(Lambda_bare_new - Lambda_bare) < 1e-10:
            break
        Lambda_bare = Lambda_bare_new

    # Final run with calibrated Λ_bare
    return integrate_phase3(psi_ini, Lambda_bare, xi, Lambda0, z_ini=z_ini)


def H_interp(z_query, bg):
    return float(np.interp(z_query, bg['z'], bg['H']))


def theta_star_from_bg(bg):
    """100θ_* = 100 × r_s(z_drag) / D_C(z_*)"""
    H0  = bg['H0']
    h   = H0 / 100.0
    Ogam = 2.469e-5 / h**2

    def rs_integrand(z):
        R   = (3.0 * omh2_b / h**2) / (4.0 * Ogam * (1 + z))
        cs  = c_kms / np.sqrt(3.0 * (1.0 + R))
        Hz  = H_interp(z, bg)
        return cs / Hz if Hz > 0 else 0.0

    rs, _ = quad(rs_integrand, z_drag, 1e4, limit=200, epsrel=1e-5)

    def DC_integrand(z):
        Hz = H_interp(z, bg)
        return c_kms / Hz if Hz > 0 else 0.0

    DC, _ = quad(DC_integrand, 0, z_star, limit=200, epsrel=1e-5)
    return 100.0 * rs / DC if DC > 0 else np.nan


def chi2_DESI(bg):
    H_model = np.array([H_interp(z, bg) for z in DESI_z])
    return float(np.sum(((H_model - DESI_H) / DESI_s)**2))


def w_eff_from_bg(bg):
    """Effective w(z) from Hubble: w = -1 + (2/3)(1+z)/H × dH/dz × 1/(1-Ω_m_eff)"""
    z_arr = bg['z']
    H_arr = bg['H']
    H0    = bg['H0']
    h     = H0 / 100.0
    Omega_m = omh2_m / h**2

    mask = (z_arr > 0.01) & (z_arr < 2.5)
    z_w  = z_arr[mask]
    H_w  = H_arr[mask]

    # numerical dH/dz
    dHdz = np.gradient(H_w, z_w)
    E2   = (H_w / H0)**2
    Om_z = Omega_m * (1 + z_w)**3 / E2

    with np.errstate(divide='ignore', invalid='ignore'):
        numer = (2.0 / 3.0) * (1 + z_w) / H_w * dHdz - 1.0
        denom = 1.0 - Om_z
        w     = np.where(np.abs(denom) > 0.05, numer / denom, np.nan)

    return z_w, w


# ─────────────────────────────────────────────────────────────────────────────
# MAIN
# ─────────────────────────────────────────────────────────────────────────────
print("=" * 72)
print("SIM131 — CMSTG Phase 3: Curvature-Memory Scalar (ξΨR coupling)")
print("=" * 72)
print(f"  ξ_conformal = 1/6 = {XI_CONFORMAL:.6f}  (fixed by symmetry)")
print(f"  Ψ_ini = 0 (sourced from curvature — no tuning)")
print(f"  Baseline: SIM121C DESI tension = {SIM121C_tension:.2f}σ, χ²={SIM121C_chi2:.2f}")
print()


# ─────────────────────────────────────────────────────────────────────────────
# PART A: Reference run — ξ=1/6, Λ₀=0.003, Ψ_ini=0
# ─────────────────────────────────────────────────────────────────────────────
print("─" * 72)
print("Part A: Reference run (ξ=1/6, Λ₀=0.003, Ψ_ini=0)")
print("─" * 72)

xi_ref     = XI_CONFORMAL
Lambda0_ref = 0.003

print(f"  Integrating Phase 3 from z=1e5 to z=0 with zero ICs...")
bg_ref = integrate_selfconsistent(0.0, xi_ref, Lambda0_ref)

if bg_ref is not None:
    psi0_ref   = bg_ref['psi0']
    Feff0_ref  = bg_ref['Feff0']
    H0_ref     = bg_ref['H0']
    Lb_ref     = bg_ref['Lambda_bare']

    psi_CMB  = float(np.interp(z_star, bg_ref['z'], bg_ref['psi']))
    Feff_CMB = float(np.interp(z_star, bg_ref['z'], bg_ref['Feff']))

    print(f"\n  Ψ trajectory:")
    print(f"    Ψ(z_CMB={z_star:.0f}) = {psi_CMB:.6f} M_Pl")
    print(f"    Ψ(z=0)              = {psi0_ref:.6f} M_Pl")
    print(f"    Phase 1 Ψ̄          = {PSI0_phase1:.3f} M_Pl")
    print(f"\n  F_eff trajectory:")
    print(f"    F_eff(z_CMB)        = {Feff_CMB:.6f}")
    print(f"    F_eff(z=0)          = {Feff0_ref:.6f}")
    print(f"    Phase 1 F₀          = {F0_phase1:.6f}")
    print(f"\n  Cosmology:")
    print(f"    H₀                  = {H0_ref:.2f} km/s/Mpc  (target: {h_target*100:.1f})")
    print(f"    Λ_bare              = {Lb_ref:.6f}")

    theta_ref  = theta_star_from_bg(bg_ref)
    chi2_ref   = chi2_DESI(bg_ref)
    tension_ref = np.sqrt(chi2_ref / len(DESI_z))
    chi2_theta_ref = ((theta_ref - theta_obs) / theta_obs_err)**2

    print(f"\n  Observables:")
    print(f"    100θ_*              = {theta_ref:.5f}  (Planck: {theta_obs:.5f}, χ²_θ={chi2_theta_ref:.2f})")
    print(f"    χ²_DESI             = {chi2_ref:.3f}  (SIM121C: {SIM121C_chi2:.3f})")
    print(f"    DESI tension        = {tension_ref:.3f}σ  (SIM121C: {SIM121C_tension:.2f}σ)")
    print(f"    Δχ²_DESI            = {chi2_ref - SIM121C_chi2:+.3f}")

    print(f"\n  H(z) at DESI redshifts:")
    print(f"    {'z':>6}  {'H_obs':>8}  {'H_Phase3':>10}  {'pull':>6}")
    for z_d, Ho, sd in zip(DESI_z, DESI_H, DESI_s):
        Hm = H_interp(z_d, bg_ref)
        print(f"    {z_d:6.3f}  {Ho:8.1f}  {Hm:10.2f}  {(Hm-Ho)/sd:6.2f}")
else:
    print("  FAILED — integration error")
    theta_ref = chi2_ref = tension_ref = chi2_theta_ref = np.nan
    psi0_ref = Feff0_ref = H0_ref = 0.0


# ─────────────────────────────────────────────────────────────────────────────
# PART B: Λ₀ scan (ξ = 1/6 fixed)
# ─────────────────────────────────────────────────────────────────────────────
print()
print("─" * 72)
print("Part B: Λ₀ scan (ξ=1/6 fixed, Ψ_ini=0)")
print("─" * 72)

Lambda0_grid = np.array([0.001, 0.002, 0.003, 0.005, 0.008, 0.010])

B_Lambda0  = []
B_psi0     = []
B_Feff0    = []
B_H0       = []
B_theta    = []
B_chi2     = []
B_tension  = []

print(f"  {'Λ₀':>8}  {'Ψ₀':>8}  {'F_eff₀':>8}  {'H₀':>7}  {'100θ*':>8}  {'tension':>8}")
for L0 in Lambda0_grid:
    bg = integrate_selfconsistent(0.0, XI_CONFORMAL, L0)
    if bg is None:
        print(f"  {L0:8.4f}  FAILED")
        continue
    theta  = theta_star_from_bg(bg)
    if np.isnan(theta):
        continue
    chi2   = chi2_DESI(bg)
    tens   = np.sqrt(chi2 / len(DESI_z))

    B_Lambda0.append(L0)
    B_psi0.append(bg['psi0'])
    B_Feff0.append(bg['Feff0'])
    B_H0.append(bg['H0'])
    B_theta.append(theta)
    B_chi2.append(chi2)
    B_tension.append(tens)

    print(f"  {L0:8.4f}  {bg['psi0']:8.4f}  {bg['Feff0']:8.5f}  "
          f"{bg['H0']:7.2f}  {theta:8.5f}  {tens:8.3f}σ")

B_Lambda0 = np.array(B_Lambda0)
B_psi0    = np.array(B_psi0)
B_Feff0   = np.array(B_Feff0)
B_H0      = np.array(B_H0)
B_theta   = np.array(B_theta)
B_chi2    = np.array(B_chi2)
B_tension = np.array(B_tension)


# ─────────────────────────────────────────────────────────────────────────────
# PART C: ξ scan (Λ₀ = 0.003 fixed)
# ─────────────────────────────────────────────────────────────────────────────
print()
print("─" * 72)
print("Part C: ξ scan (Λ₀=0.003 fixed, Ψ_ini=0)")
print("─" * 72)

xi_grid = np.array([0.05, 1.0/12.0, 1.0/6.0, 1.0/4.0, 1.0/3.0])
xi_labels = ['1/20', '1/12', '1/6 (conf)', '1/4', '1/3']

C_xi      = []
C_psi0    = []
C_Feff0   = []
C_H0      = []
C_theta   = []
C_chi2    = []
C_tension = []

print(f"  {'ξ':>10}  {'Ψ₀':>8}  {'F_eff₀':>8}  {'H₀':>7}  {'100θ*':>8}  {'tension':>8}")
for xi, lbl in zip(xi_grid, xi_labels):
    bg = integrate_selfconsistent(0.0, xi, 0.003)
    if bg is None:
        print(f"  {lbl:>10}  FAILED")
        continue
    theta  = theta_star_from_bg(bg)
    if np.isnan(theta):
        continue
    chi2   = chi2_DESI(bg)
    tens   = np.sqrt(chi2 / len(DESI_z))

    C_xi.append(xi)
    C_psi0.append(bg['psi0'])
    C_Feff0.append(bg['Feff0'])
    C_H0.append(bg['H0'])
    C_theta.append(theta)
    C_chi2.append(chi2)
    C_tension.append(tens)

    print(f"  {lbl:>10}  {bg['psi0']:8.4f}  {bg['Feff0']:8.5f}  "
          f"{bg['H0']:7.2f}  {theta:8.5f}  {tens:8.3f}σ")

C_xi      = np.array(C_xi)
C_psi0    = np.array(C_psi0)
C_Feff0   = np.array(C_Feff0)
C_H0      = np.array(C_H0)
C_theta   = np.array(C_theta)
C_chi2    = np.array(C_chi2)
C_tension = np.array(C_tension)


# ─────────────────────────────────────────────────────────────────────────────
# PART D: Best joint fit and verdict
# ─────────────────────────────────────────────────────────────────────────────
print()
print("─" * 72)
print("Part D: Best joint fit (DESI + CMB)")
print("─" * 72)

# Combine all scan results from B and C
all_xi     = np.concatenate([np.full(len(B_Lambda0), XI_CONFORMAL), C_xi])
all_L0     = np.concatenate([B_Lambda0, np.full(len(C_xi), 0.003)])
all_psi0   = np.concatenate([B_psi0,   C_psi0])
all_Feff0  = np.concatenate([B_Feff0,  C_Feff0])
all_theta  = np.concatenate([B_theta,  C_theta])
all_chi2   = np.concatenate([B_chi2,   C_chi2])
all_tens   = np.concatenate([B_tension, C_tension])

chi2_theta_arr = ((all_theta - theta_obs) / theta_obs_err)**2
chi2_tot_arr   = all_chi2 + chi2_theta_arr

if len(chi2_tot_arr) > 0:
    idx_best  = np.argmin(chi2_tot_arr)
    xi_map    = all_xi[idx_best]
    L0_map    = all_L0[idx_best]
    psi0_map  = all_psi0[idx_best]
    Feff0_map = all_Feff0[idx_best]
    theta_map = all_theta[idx_best]
    chi2_map  = all_chi2[idx_best]
    chi2t_map = chi2_theta_arr[idx_best]
    tens_map  = all_tens[idx_best]
    delta_chi2 = chi2_map - SIM121C_chi2

    print(f"\n  Best joint point:")
    print(f"    ξ = {xi_map:.5f},  Λ₀ = {L0_map:.4f}")
    print(f"    Ψ₀       = {psi0_map:.4f} M_Pl")
    print(f"    F_eff₀   = {Feff0_map:.5f}")
    print(f"    100θ_*   = {theta_map:.5f}  (χ²_θ={chi2t_map:.3f})")
    print(f"    χ²_DESI  = {chi2_map:.3f}  (SIM121C: {SIM121C_chi2:.3f})")
    print(f"    tension  = {tens_map:.3f}σ  (SIM121C: {SIM121C_tension:.2f}σ)")
    print(f"    Δχ²_DESI = {delta_chi2:+.3f}")

    # Re-integrate at MAP for detailed outputs
    bg_map = integrate_selfconsistent(0.0, xi_map, L0_map)
    if bg_map:
        print(f"\n  H(z) at MAP vs DESI:")
        print(f"    {'z':>6}  {'H_obs':>8}  {'H_MAP':>8}  {'pull':>6}")
        for z_d, Ho, sd in zip(DESI_z, DESI_H, DESI_s):
            Hm = H_interp(z_d, bg_map)
            print(f"    {z_d:6.3f}  {Ho:8.1f}  {Hm:8.2f}  {(Hm-Ho)/sd:6.2f}")
    else:
        bg_map = bg_ref
else:
    print("  No valid scan points — using reference run.")
    xi_map = XI_CONFORMAL; L0_map = 0.003; psi0_map = psi0_ref
    Feff0_map = Feff0_ref; theta_map = theta_ref; chi2_map = chi2_ref
    chi2t_map = chi2_theta_ref; tens_map = tension_ref
    delta_chi2 = chi2_ref - SIM121C_chi2
    bg_map = bg_ref

# Verdict
pass_desi  = (tens_map < 2.0) if not np.isnan(tens_map) else False
pass_theta = (chi2t_map < 4.0) if not np.isnan(chi2t_map) else False
beats_ref  = (chi2_map < SIM121C_chi2) if not np.isnan(chi2_map) else False

if pass_desi and pass_theta:
    verdict = "PASS"
elif beats_ref and pass_theta:
    verdict = "PARTIAL"
elif beats_ref:
    verdict = "PARTIAL"
else:
    verdict = "FAIL"


# ─────────────────────────────────────────────────────────────────────────────
# VERDICT
# ─────────────────────────────────────────────────────────────────────────────
print()
print("=" * 72)
print("SIM131 RESULT:")
print()
print(f"  Verdict:       {verdict}")
print(f"  DESI tension:  {tens_map:.3f}σ  (Phase 1 floor: {SIM121C_tension:.2f}σ)")
print(f"  Δχ²_DESI:      {delta_chi2:+.3f}")
print(f"  CMB θ*:        {'PASS' if pass_theta else 'FAIL'}  (χ²_θ={chi2t_map:.2f})")
print(f"  Ψ₀:            {psi0_map:.4f} M_Pl  (grown from Ψ_ini=0)")
print(f"  F_eff₀:        {Feff0_map:.5f}")
print()

if verdict == "PASS":
    print("  CONCLUSION: Phase 3 curvature-memory CMSTG clears the DESI floor.")
    print("  The ξΨR coupling sources Ψ from zero ICs and generates dynamic DE.")
    print("  Phase 3 represents a genuine first-principles improvement over Phase 1.")
elif verdict == "PARTIAL":
    print("  CONCLUSION: Phase 3 improves on the Phase 1 DESI floor.")
    print(f"  Reduction: {SIM121C_tension:.2f}σ → {tens_map:.3f}σ.")
    print("  Further optimisation (SIM132: 2D scan ξ×Λ₀) warranted.")
else:
    print("  CONCLUSION: Curvature-memory scalar does not improve DESI tension.")
    if psi0_ref > 0.1:
        print(f"  Ψ does evolve (Ψ₀={psi0_ref:.3f} M_Pl) but H(z) shape is wrong direction.")
        print("  Consider: negative ξ (anti-conformal), or higher-order F_eff.")
    else:
        print("  Ψ remains frozen — investigate sourcing at early times.")

print("=" * 72)


# ─────────────────────────────────────────────────────────────────────────────
# FIGURES
# ─────────────────────────────────────────────────────────────────────────────
print()
print("Generating figures...")

# Figure 1: Ψ(z) and F_eff(z) for reference run
if bg_ref is not None:
    fig, axes = plt.subplots(1, 2, figsize=(13, 5))
    z_p  = bg_ref['z']
    mask = z_p < 3000

    axes[0].semilogx(z_p[mask]+1, bg_ref['psi'][mask],
                     color='#d73027', lw=2, label='Phase 3 Ψ(z) (ξ=1/6, Λ₀=0.003)')
    axes[0].axhline(PSI0_phase1, color='#2166ac', ls='--', lw=1.5,
                    label=f'Phase 1 Ψ̄={PSI0_phase1} (frozen)')
    axes[0].axhline(0, color='gray', ls=':', lw=1, alpha=0.5)
    axes[0].axvline(z_star+1, color='orange', ls=':', lw=1.5,
                    label=f'$z_*={z_star:.0f}$')
    axes[0].set_xlabel(r'$1+z$')
    axes[0].set_ylabel(r'$\Psi(z)$ [$M_{\rm Pl}$]')
    axes[0].set_title(r'Phase 3: Curvature-memory $\Psi(z)$ from $\Psi_{\rm ini}=0$')
    axes[0].legend(fontsize=9); axes[0].grid(alpha=0.3)

    Feff_arr = bg_ref['Feff']
    axes[1].semilogx(z_p[mask]+1, Feff_arr[mask],
                     color='#d73027', lw=2, label=r'Phase 3 $F_{\rm eff}(z)$')
    axes[1].axhline(F0_phase1, color='#2166ac', ls='--', lw=1.5,
                    label=f'Phase 1 $F_0={F0_phase1:.4f}$ (frozen)')
    axes[1].axhline(0.5, color='gray', ls=':', lw=1.5, label='GR ($F=1/2$)')
    axes[1].axvline(z_star+1, color='orange', ls=':', lw=1.5)
    axes[1].set_xlabel(r'$1+z$')
    axes[1].set_ylabel(r'$F_{\rm eff}(\Psi) = \frac{1+2\Lambda_0\Psi^2}{2} + \xi\Psi$')
    axes[1].set_title(r'Running $F_{\rm eff}(z)$ — Phase 3 vs Phase 1')
    axes[1].legend(fontsize=9); axes[1].grid(alpha=0.3)

    plt.tight_layout()
    for ext in ('pdf', 'png'):
        fig.savefig(os.path.join(OUT, f'sim131_psi_evolution.{ext}'), dpi=150, bbox_inches='tight')
    plt.close(fig)
    print("  Saved sim131_psi_evolution")

# Figure 2: DESI tension landscape — Λ₀ and ξ scans
if len(B_Lambda0) > 1 and len(C_xi) > 1:
    fig, axes = plt.subplots(1, 2, figsize=(13, 5))

    axes[0].plot(B_Lambda0, B_tension, 'o-', color='#d73027', lw=2, ms=6,
                 label='Phase 3 (ξ=1/6 fixed)')
    axes[0].axhline(SIM121C_tension, color='#2166ac', ls='--', lw=1.5,
                    label=f'Phase 1 floor ({SIM121C_tension:.2f}σ)')
    axes[0].axhline(2.0, color='green', ls=':', lw=1.5, label='2σ PASS threshold')
    axes[0].axvline(0.003, color='orange', ls=':', lw=1.5, label='Λ₀=0.003 (Phase 1)')
    axes[0].set_xlabel(r'$\Lambda_0$')
    axes[0].set_ylabel(r'DESI tension [$\sigma$]')
    axes[0].set_title(r'DESI tension vs $\Lambda_0$ (ξ=1/6)')
    axes[0].legend(fontsize=9); axes[0].grid(alpha=0.3)

    axes[1].plot(C_xi, C_tension, 's-', color='#1a9850', lw=2, ms=6,
                 label='Phase 3 (Λ₀=0.003 fixed)')
    axes[1].axhline(SIM121C_tension, color='#2166ac', ls='--', lw=1.5,
                    label=f'Phase 1 floor ({SIM121C_tension:.2f}σ)')
    axes[1].axhline(2.0, color='green', ls=':', lw=1.5, label='2σ PASS threshold')
    axes[1].axvline(XI_CONFORMAL, color='orange', ls=':', lw=1.5,
                    label=f'ξ=1/6 (conformal)')
    axes[1].set_xlabel(r'$\xi$')
    axes[1].set_ylabel(r'DESI tension [$\sigma$]')
    axes[1].set_title(r'DESI tension vs $\xi$ ($\Lambda_0=0.003$)')
    axes[1].legend(fontsize=9); axes[1].grid(alpha=0.3)

    plt.tight_layout()
    for ext in ('pdf', 'png'):
        fig.savefig(os.path.join(OUT, f'sim131_tension_landscape.{ext}'), dpi=150, bbox_inches='tight')
    plt.close(fig)
    print("  Saved sim131_tension_landscape")

# Figure 3: H(z) comparison — Phase 3 vs Phase 1 vs DESI
if bg_ref is not None:
    fig, axes = plt.subplots(1, 2, figsize=(13, 5))
    z_fine = np.linspace(0.01, 2.5, 300)

    H_ref_fine = [H_interp(z, bg_ref) for z in z_fine]

    axes[0].plot(z_fine, H_ref_fine, color='#d73027', lw=2,
                 label=f'Phase 3 (ξ=1/6, Λ₀=0.003)')
    if bg_map is not None and bg_map is not bg_ref:
        H_map_fine = [H_interp(z, bg_map) for z in z_fine]
        axes[0].plot(z_fine, H_map_fine, color='purple', lw=2, ls='--',
                     label=f'Phase 3 MAP (ξ={xi_map:.3f}, Λ₀={L0_map:.3f})')
    axes[0].errorbar(DESI_z, DESI_H, yerr=DESI_s, fmt='ko', ms=7, capsize=4,
                     label='DESI Y1 BAO', zorder=5)
    axes[0].set_xlabel(r'Redshift $z$')
    axes[0].set_ylabel(r'$H(z)$ [km/s/Mpc]')
    axes[0].set_title('Phase 3 H(z) vs DESI')
    axes[0].legend(fontsize=9); axes[0].grid(alpha=0.3)

    # w_eff(z)
    if bg_ref is not None:
        z_w, w_w = w_eff_from_bg(bg_ref)
        axes[1].plot(z_w, w_w, color='#d73027', lw=2, label='Phase 3 (ξ=1/6, Λ₀=0.003)')
        axes[1].axhline(-1.0, color='#2166ac', ls='--', lw=1.5, label='ΛCDM (w=−1)')
        axes[1].set_xlabel(r'Redshift $z$')
        axes[1].set_ylabel(r'$w_{\rm eff}(z)$')
        axes[1].set_title(r'Effective $w(z)$ — Phase 3 curvature memory')
        axes[1].set_ylim(-2.5, 0.5)
        axes[1].legend(fontsize=9); axes[1].grid(alpha=0.3)

    plt.tight_layout()
    for ext in ('pdf', 'png'):
        fig.savefig(os.path.join(OUT, f'sim131_Hz_weff.{ext}'), dpi=150, bbox_inches='tight')
    plt.close(fig)
    print("  Saved sim131_Hz_weff")

# Figure 4: Ψ₀ and F_eff₀ vs Λ₀ (Part B summary)
if len(B_Lambda0) > 1:
    fig, axes = plt.subplots(1, 2, figsize=(13, 5))

    axes[0].plot(B_Lambda0, B_psi0, 'o-', color='#d73027', lw=2, ms=6)
    axes[0].axhline(PSI0_phase1, color='#2166ac', ls='--', lw=1.5,
                    label=f'Phase 1 Ψ̄={PSI0_phase1}')
    axes[0].axvline(0.003, color='orange', ls=':', lw=1.5, label='Phase 1 Λ₀')
    axes[0].set_xlabel(r'$\Lambda_0$')
    axes[0].set_ylabel(r'$\Psi(z=0)$ [$M_{\rm Pl}$]')
    axes[0].set_title(r'Today\'s $\Psi_0$ vs $\Lambda_0$ (ξ=1/6, $\Psi_{\rm ini}=0$)')
    axes[0].legend(fontsize=9); axes[0].grid(alpha=0.3)

    axes[1].plot(B_Lambda0, B_Feff0, 's-', color='#1a9850', lw=2, ms=6)
    axes[1].axhline(F0_phase1, color='#2166ac', ls='--', lw=1.5,
                    label=f'Phase 1 $F_0={F0_phase1:.4f}$')
    axes[1].axhline(0.5, color='gray', ls=':', lw=1.5, label='GR (F=1/2)')
    axes[1].axvline(0.003, color='orange', ls=':', lw=1.5, label='Phase 1 Λ₀')
    axes[1].set_xlabel(r'$\Lambda_0$')
    axes[1].set_ylabel(r'$F_{\rm eff}(z=0)$')
    axes[1].set_title(r'$F_{\rm eff,0}$ vs $\Lambda_0$')
    axes[1].legend(fontsize=9); axes[1].grid(alpha=0.3)

    plt.tight_layout()
    for ext in ('pdf', 'png'):
        fig.savefig(os.path.join(OUT, f'sim131_psi0_Feff0.{ext}'), dpi=150, bbox_inches='tight')
    plt.close(fig)
    print("  Saved sim131_psi0_Feff0")


# ─────────────────────────────────────────────────────────────────────────────
# SAVE JSON
# ─────────────────────────────────────────────────────────────────────────────
results = {
    "sim": "SIM131",
    "phase": "Phase 3",
    "verdict": verdict,
    "description": "Phase 3 curvature-memory scalar: F_eff=(1+2Λ₀Ψ²)/2+ξΨ, ξ=1/6 fixed, Ψ_ini=0",
    "action": "S = ∫d⁴x√g [F_eff(Ψ)·R − ½(∂Ψ)²] + S_SM, F_eff = (1+2Λ₀Ψ²)/2 + ξΨ",
    "field_eom": "Ψ'' + (3−ε_H)Ψ' = (2Λ₀Ψ + ξ)·R/H², R/H²=6(2−ε_H)",
    "xi_conformal": float(XI_CONFORMAL),
    "reference_run": {
        "xi": float(XI_CONFORMAL),
        "Lambda0": 0.003,
        "psi_ini": 0.0,
        "psi0": float(psi0_ref),
        "Feff0": float(Feff0_ref),
        "H0": float(H0_ref),
        "theta_star_100": float(theta_ref),
        "chi2_DESI": float(chi2_ref),
        "DESI_tension_sigma": float(tension_ref),
    },
    "map": {
        "xi": float(xi_map),
        "Lambda0": float(L0_map),
        "psi0": float(psi0_map),
        "Feff0": float(Feff0_map),
        "theta_star_100": float(theta_map),
        "chi2_DESI": float(chi2_map),
        "chi2_theta": float(chi2t_map),
        "DESI_tension_sigma": float(tens_map),
        "delta_chi2_vs_SIM121C": float(delta_chi2),
    },
    "scan_Lambda0": {
        "xi_fixed": float(XI_CONFORMAL),
        "Lambda0": B_Lambda0.tolist(),
        "psi0": B_psi0.tolist(),
        "Feff0": B_Feff0.tolist(),
        "tension_sigma": B_tension.tolist(),
        "chi2_DESI": B_chi2.tolist(),
        "theta_100": B_theta.tolist(),
    },
    "scan_xi": {
        "Lambda0_fixed": 0.003,
        "xi": C_xi.tolist(),
        "psi0": C_psi0.tolist(),
        "Feff0": C_Feff0.tolist(),
        "tension_sigma": C_tension.tolist(),
        "chi2_DESI": C_chi2.tolist(),
        "theta_100": C_theta.tolist(),
    },
    "baseline_SIM121C": {
        "chi2_DESI": SIM121C_chi2,
        "tension_sigma": SIM121C_tension,
        "psi0_frozen": PSI0_phase1,
        "F0_frozen": F0_phase1,
    },
    "pass_desi": bool(pass_desi),
    "pass_theta": bool(pass_theta),
    "beats_SIM121C": bool(beats_ref),
}

with open(os.path.join(OUT, 'sim131_results.json'), 'w') as f:
    json.dump(results, f, indent=2)
print("  Saved sim131_results.json")
print(f"\nAll outputs in: {OUT}")
print("SIM131 complete.")
