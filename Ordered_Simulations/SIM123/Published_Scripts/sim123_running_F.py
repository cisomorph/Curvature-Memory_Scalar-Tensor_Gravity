"""
SIM123 — CMSTG Phase 3: Running F(z) — Full Background ODE Evolution
====================================================================
SIM121C and SIM122 used a frozen-Ψ approximation: F₀ = ½ + Λ₀Ψ₀² = const
at ALL redshifts, including z_CMB = 1090. This forced a structural tension:
matching Planck θ_* required F₀≈0.560, which overshoots DESI H(z) by ~10%.

The frozen-Ψ approximation is WRONG for the CMSTG hilltop potential. In the
SSB potential V_J = λ(Ψ²−v²)², the field starts near the hilltop (Ψ≈0) at
high z and slowly rolls toward Ψ₀ = 2.62 M_Pl today. The effective coupling
F(z) = ½ + Λ₀Ψ²(z) therefore EVOLVES:
  F(z_CMB) ≈ ½ + Λ₀Ψ²_ini  ≈ ½   (near GR if Ψ_ini << Ψ₀)
  F(z=0)   = ½ + Λ₀Ψ₀²     = F₀   (CMSTG modification today)

If Ψ_ini << Ψ₀, the CMB sees approximately GR (F_CMB ≈ ½) → θ_* gives
H₀ ≈ 67 km/s/Mpc (standard value). DESI observes H(z) at low z where F(z)
has evolved from ½ to F₀ — a mild, z-dependent modification. This may allow
simultaneous CMB+DESI consistency.

This simulation integrates the full CMSTG background ODEs:
  3F(Ψ)H² = ρ_m + ρ_r + V(Ψ)        [Friedmann, slow-roll]
  3Hψ̇ ≈ −V′(Ψ) + 24Λ₀Ψ H²          [Klein-Gordon, slow-roll]
  V_J(Ψ) = λ(Ψ²−v²)²

Free parameters:
  Λ₀    ∈ [0.001, 0.010]    coupling (was locked at 0.003)
  Ψ_ini ∈ [0.001, 2.5] M_Pl  initial field value (high-z boundary condition)
  v     = 13.16 M_Pl         (fixed from SIM113)

Fixed: ω_m = 0.1430, ω_b = 0.02237, ω_r = 4.18×10⁻⁵ (Planck 2018)

Outputs:
  • F(z) evolution plot: running vs frozen
  • θ_*(Λ₀, Ψ_ini) landscape
  • DESI tension landscape
  • sim123_results.json
"""

import numpy as np
from scipy.integrate import quad, solve_ivp
from scipy.optimize import brentq
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
# CONSTANTS
# ─────────────────────────────────────────────────────────────────────────────
c_kms    = 2.998e5         # km/s
H100     = 100.0           # reference H [km/s/Mpc]; physical densities in H100² units
omh2_m   = 0.1430          # Ω_m h²
omh2_b   = 0.02237         # Ω_b h²
omh2_r   = 4.18e-5         # Ω_r h²

theta_obs     = 1.04101
theta_obs_err = 0.00029
z_drag = 1059.6
z_star = 1089.8

DESI_z = np.array([0.295, 0.510, 0.706, 0.930, 1.317, 2.330])
DESI_H = np.array([ 81.7,  97.9, 110.7, 128.1, 156.4, 240.8])
DESI_s = np.array([  4.5,   4.4,   6.2,   5.6,   8.6,  11.0])

# SIM113/SIM121C reference values
L0_ref   = 0.003
PSI0_ref = 2.62    # field value today [M_Pl]
V_ref    = 13.16   # SSB potential VEV [M_Pl]
F0_ref   = 0.5 + L0_ref * PSI0_ref**2  # = 0.52059

# SIM121C reference
SIM121C_chi2    = 41.492
SIM121C_tension = 2.630

# ─────────────────────────────────────────────────────────────────────────────
# ODE SYSTEM: CMSTG background in slow-roll
# ─────────────────────────────────────────────────────────────────────────────

def cmstg_ode_slowroll(N, psi_arr, L0, lam_scaled, v):
    """
    First-order slow-roll ODE for CMSTG scalar field.
    State: [psi]  (scalar, not a 2-component system)

    Slow-roll KG equation:
        3H² × dψ/dN ≈ −V′(ψ) + 24Λ₀ψH²
    → dψ/dN = (−V′/H² + 24Λ₀ψ) / 3

    All quantities in H100=1, M_Pl=1 units.
    """
    psi = psi_arr[0]
    a   = np.exp(N)

    rho_m = omh2_m / a**3
    rho_r = omh2_r / a**4
    F     = 0.5 + L0 * psi**2
    V     = lam_scaled * (psi**2 - v**2)**2
    Vp    = 4.0 * lam_scaled * (psi**2 - v**2) * psi

    H_sq = (rho_m + rho_r + V) / (3.0 * F)
    if H_sq <= 0:
        return [0.0]

    dpsi_dN = (-Vp / H_sq + 24.0 * L0 * psi) / 3.0
    return [dpsi_dN]

def integrate_cmstg_background(L0, psi_ini, v, z_ini=1e5, z_end=0.0, n_points=2000):
    """
    Integrate the CMSTG background from z_ini to z_end.
    Returns arrays: z, psi(z), H(z), F(z)
    λ is calibrated to give Ω_DE ≈ 0.685 at z=0 for the given (psi_ini, v).
    """
    N_ini = -np.log(1.0 + z_ini)
    N_end = 0.0

    # Calibrate λ: at z=0 we need 3F(psi0)H0² = rho_m0 + rho_r0 + V(psi0)
    # AND Omega_DE = V(psi0)/(3F0 H0²) ≈ 0.685
    # → V(psi0) = 0.685 × 3F0 × H0²
    # H0² = (rho_m0+rho_r0+V(psi0))/(3F0)
    # V(psi0) = λ(psi0²-v²)²
    # Together: λ = 0.685 × (rho_m0+rho_r0+V) / ((psi0²-v²)²)  [implicit in psi0]
    # We need to find psi0 self-consistently via the ODE, but for a first-order estimate:
    # Use psi0 ≈ PSI0_ref = 2.62 and iterate once.

    psi0_est = PSI0_ref  # estimate of field value today
    rho_m0   = omh2_m    # at z=0 (a=1)
    rho_r0   = omh2_r

    F0_est = 0.5 + L0 * psi0_est**2
    # V(psi0) = 0.685 × 3F0 × H0²
    # H0² = (rho_m0 + rho_r0 + lam*(psi0²-v²)²) / (3F0)
    # → 3F0 H0² = rho_m0 + rho_r0 + lam*(psi0²-v²)²
    # With Omega_DE = lam*(psi0²-v²)² / (3F0 H0²) = 0.685:
    # lam*(psi0²-v²)² = 0.685 × (rho_m0 + rho_r0 + lam*(psi0²-v²)²)
    # lam*(psi0²-v²)²(1-0.685) = 0.685*(rho_m0+rho_r0)
    # lam = 0.685*(rho_m0+rho_r0) / (0.315*(psi0²-v²)²)

    denom_sq = (psi0_est**2 - v**2)**2
    if denom_sq < 1e-10:
        return None  # field at minimum — unphysical for DE

    lam_scaled = 0.685 * (rho_m0 + rho_r0) / (0.315 * denom_sq)

    # Integrate first-order slow-roll ODE: state = [psi]
    N_span = (N_ini, N_end)
    N_eval = np.linspace(N_ini, N_end, n_points)
    sol = solve_ivp(
        cmstg_ode_slowroll,
        N_span,
        [psi_ini],
        args=(L0, lam_scaled, v),
        t_eval=N_eval,
        method='RK45',
        rtol=1e-6, atol=1e-10,
        dense_output=False
    )

    if not sol.success:
        return None

    N_arr   = sol.t
    z_arr   = np.exp(-N_arr) - 1.0
    psi_arr = sol.y[0]

    # Compute H(z) along the solution
    H_arr = []
    F_arr = []
    for i, (N_i, psi_i) in enumerate(zip(N_arr, psi_arr)):
        a_i   = np.exp(N_i)
        rho_m = omh2_m / a_i**3
        rho_r = omh2_r / a_i**4
        V_i   = lam_scaled * (psi_i**2 - v**2)**2
        F_i   = 0.5 + L0 * psi_i**2
        H_sq  = (rho_m + rho_r + V_i) / (3.0 * F_i)
        H_arr.append(H100 * np.sqrt(max(H_sq, 0.0)))
        F_arr.append(F_i)

    H_arr = np.array(H_arr)
    F_arr = np.array(F_arr)

    # Sort by z (increasing)
    idx = np.argsort(z_arr)
    return {
        'z':   z_arr[idx],
        'N':   N_arr[idx],
        'psi': psi_arr[idx],
        'H':   H_arr[idx],
        'F':   F_arr[idx],
        'lam_scaled': lam_scaled,
        'psi0': psi_arr[-1],     # field value today (last point, z=0)
        'F0':   F_arr[-1],       # F today
        'H0':   H_arr[-1],       # H today
    }

def H_interp(z_query, bg):
    """Interpolate H(z) from background solution."""
    return np.interp(z_query, bg['z'], bg['H'])

def theta_star_from_bg(bg):
    """100×θ_* = 100 × r_s(z_drag) / D_C(z_*) from ODE-derived H(z)."""
    H0  = bg['H0']
    omh2_b_loc = omh2_b
    h   = H0 / 100.0
    Ogam = 2.469e-5 / h**2

    def integrand_rs(z):
        R  = (3.0 * omh2_b_loc / h**2) / (4.0 * Ogam * (1+z))
        cs = c_kms / np.sqrt(3.0 * (1.0 + R))
        Hz = H_interp(z, bg)
        if Hz <= 0: return 0.0
        return cs / Hz

    rs, _ = quad(integrand_rs, z_drag, 1e4, limit=200, epsrel=1e-5)

    def integrand_DC(z):
        Hz = H_interp(z, bg)
        if Hz <= 0: return 0.0
        return c_kms / Hz

    DC, _ = quad(integrand_DC, 0, z_star, limit=200, epsrel=1e-5)

    if DC <= 0: return np.nan
    return 100.0 * rs / DC

def chi2_DESI_from_bg(bg):
    """χ²_DESI from ODE-derived H(z)."""
    H_model = np.array([H_interp(z, bg) for z in DESI_z])
    return np.sum(((H_model - DESI_H) / DESI_s)**2)

# ─────────────────────────────────────────────────────────────────────────────
# PART A: Demonstrate F(z) evolution for SIM113 parameters
# ─────────────────────────────────────────────────────────────────────────────
print("=" * 70)
print("SIM123 — CMSTG Phase 3: Running F(z) Full Background Evolution")
print("=" * 70)
print()
print("Part A: F(z) evolution for SIM113 parameters (Λ₀=0.003, Ψ_ini=0.01)...")

bg_run = integrate_cmstg_background(L0=L0_ref, psi_ini=0.01, v=V_ref, z_ini=1e5)

if bg_run is not None:
    # Find F at CMB and today
    z_arr  = bg_run['z']
    F_arr  = bg_run['F']
    psi_arr= bg_run['psi']
    H_arr  = bg_run['H']

    # F at z_CMB
    F_at_CMB   = np.interp(z_star, z_arr, F_arr)
    psi_at_CMB = np.interp(z_star, z_arr, psi_arr)
    F_today    = bg_run['F0']
    psi_today  = bg_run['psi0']
    H_today    = bg_run['H0']

    print(f"\n  Field evolution (Ψ_ini=0.01, Λ₀=0.003, v=13.16):")
    print(f"  Ψ(z_CMB=1090) = {psi_at_CMB:.4f} M_Pl")
    print(f"  Ψ(z=0)        = {psi_today:.4f} M_Pl  (SIM113: {PSI0_ref:.4f})")
    print(f"  F(z_CMB=1090) = {F_at_CMB:.5f}  (GR = 0.50000)")
    print(f"  F(z=0)        = {F_today:.5f}  (frozen-Ψ: {F0_ref:.5f})")
    print(f"  H₀            = {H_today:.2f} km/s/Mpc")

    theta_run = theta_star_from_bg(bg_run)
    chi2_run  = chi2_DESI_from_bg(bg_run)
    tension_run = np.sqrt(chi2_run / len(DESI_z))
    print(f"\n  100θ_* (running F) = {theta_run:.5f}  (Planck: {theta_obs:.5f})")
    print(f"  χ²_DESI = {chi2_run:.3f}  (SIM121C frozen: {SIM121C_chi2:.3f})")
    print(f"  DESI tension = {tension_run:.2f}σ  (SIM121C: {SIM121C_tension:.2f}σ)")
else:
    print("  ODE failed for reference parameters.")
    theta_run = np.nan
    chi2_run  = np.nan
    tension_run = np.nan

# ─────────────────────────────────────────────────────────────────────────────
# PART B: Scan over (Λ₀, Ψ_ini) — tension landscape
# ─────────────────────────────────────────────────────────────────────────────
print()
print("Part B: 2D scan over (Λ₀, Ψ_ini)...")

L0_grid    = np.linspace(0.001, 0.008, 10)
psi_ini_grid = np.array([0.001, 0.01, 0.05, 0.1, 0.3, 0.5, 0.8, 1.0, 1.5, 2.0])

scan_L0      = []
scan_psi_ini = []
scan_theta   = []
scan_chi2    = []
scan_tension = []
scan_F_CMB   = []
scan_psi0    = []
scan_H0      = []

for L0 in L0_grid:
    for psi_ini in psi_ini_grid:
        bg = integrate_cmstg_background(L0=L0, psi_ini=psi_ini, v=V_ref, z_ini=1e5)
        if bg is None:
            continue
        theta = theta_star_from_bg(bg)
        if np.isnan(theta):
            continue
        chi2  = chi2_DESI_from_bg(bg)
        F_CMB = float(np.interp(z_star, bg['z'], bg['F']))
        scan_L0.append(L0)
        scan_psi_ini.append(psi_ini)
        scan_theta.append(theta)
        scan_chi2.append(chi2)
        scan_tension.append(np.sqrt(chi2/len(DESI_z)))
        scan_F_CMB.append(F_CMB)
        scan_psi0.append(bg['psi0'])
        scan_H0.append(bg['H0'])

scan_L0      = np.array(scan_L0)
scan_psi_ini = np.array(scan_psi_ini)
scan_theta   = np.array(scan_theta)
scan_chi2    = np.array(scan_chi2)
scan_tension = np.array(scan_tension)
scan_F_CMB   = np.array(scan_F_CMB)
scan_psi0    = np.array(scan_psi0)
scan_H0      = np.array(scan_H0)

print(f"  Scan complete: {len(scan_chi2)} valid points.")

# Find best point by θ_* proximity first, then DESI
theta_resid  = np.abs(scan_theta - theta_obs)
# Best point with θ_* within 3σ
cmb_mask = theta_resid < 3 * theta_obs_err
if cmb_mask.sum() > 0:
    chi2_sub  = scan_chi2[cmb_mask]
    idx_best  = np.where(cmb_mask)[0][np.argmin(chi2_sub)]
else:
    idx_best  = np.argmin(theta_resid)

L0_best    = scan_L0[idx_best]
psi_best   = scan_psi_ini[idx_best]
theta_best = scan_theta[idx_best]
chi2_best  = scan_chi2[idx_best]
ten_best   = scan_tension[idx_best]
F_CMB_best = scan_F_CMB[idx_best]
psi0_best  = scan_psi0[idx_best]
H0_best    = scan_H0[idx_best]

print(f"\n  Best point (CMB θ_* match + min DESI χ²):")
print(f"  Λ₀={L0_best:.4f}, Ψ_ini={psi_best:.3f}, Ψ₀={psi0_best:.3f}, H₀={H0_best:.2f}")
print(f"  F_CMB={F_CMB_best:.5f},  θ_*={theta_best:.5f} (obs={theta_obs:.5f})")
print(f"  χ²_DESI={chi2_best:.3f}, tension={ten_best:.2f}σ")

# ─────────────────────────────────────────────────────────────────────────────
# PART C: Refine best point — scan finer grid
# ─────────────────────────────────────────────────────────────────────────────
print()
print("Part C: Fine grid around best point...")

L0_fine   = np.linspace(max(0.001, L0_best-0.002), min(0.010, L0_best+0.002), 12)
psi_fine  = np.linspace(max(0.001, psi_best*0.3), min(2.5, psi_best*3.0), 15)

fine_results = []
for L0 in L0_fine:
    for psi_ini in psi_fine:
        bg = integrate_cmstg_background(L0=L0, psi_ini=psi_ini, v=V_ref, z_ini=1e5)
        if bg is None: continue
        theta = theta_star_from_bg(bg)
        if np.isnan(theta): continue
        chi2  = chi2_DESI_from_bg(bg)
        chi2_theta = ((theta - theta_obs) / theta_obs_err)**2
        chi2_tot   = chi2 + chi2_theta
        F_CMB = float(np.interp(z_star, bg['z'], bg['F']))
        fine_results.append((chi2_tot, chi2, chi2_theta, L0, psi_ini,
                             bg['psi0'], bg['H0'], bg['F0'], F_CMB, theta))

fine_results.sort(key=lambda x: x[0])

print(f"  Fine grid: {len(fine_results)} valid points.")
if fine_results:
    print(f"\n  Top 5 joint minima (χ²_DESI + χ²_θ):")
    print(f"  {'chi2_tot':>9} {'chi2_DESI':>9} {'chi2_th':>7} "
          f"{'L0':>7} {'psi_ini':>8} {'psi0':>6} {'H0':>7} {'F_CMB':>7} {'theta':>8}")
    for row in fine_results[:5]:
        print(f"  {row[0]:9.3f} {row[1]:9.3f} {row[2]:7.3f} "
              f"  {row[3]:7.5f} {row[4]:8.3f} {row[5]:6.3f} {row[6]:7.2f} {row[7]+0.0:7.5f} {row[9]:8.5f}")

# ─────────────────────────────────────────────────────────────────────────────
# PART D: Extract best-fit and compare
# ─────────────────────────────────────────────────────────────────────────────
if fine_results:
    best = fine_results[0]
    chi2_tot_map, chi2_DESI_map, chi2_theta_map = best[0], best[1], best[2]
    L0_map, psi_ini_map = best[3], best[4]
    psi0_map, H0_map, F0_map, F_CMB_map, theta_map = best[5], best[6], best[7], best[8], best[9]
    tension_map = np.sqrt(chi2_DESI_map / len(DESI_z))
else:
    chi2_DESI_map = chi2_best
    chi2_theta_map = ((theta_best - theta_obs)/theta_obs_err)**2
    chi2_tot_map = chi2_DESI_map + chi2_theta_map
    L0_map, psi_ini_map = L0_best, psi_best
    psi0_map, H0_map, F0_map, F_CMB_map, theta_map = psi0_best, H0_best, F0_ref, F_CMB_best, theta_best
    tension_map = ten_best

bg_map = integrate_cmstg_background(L0=L0_map, psi_ini=psi_ini_map, v=V_ref, z_ini=1e5)

print()
print("=" * 70)
print("SIM123 RESULT:")
print()
print(f"  MAP parameters:")
print(f"    Λ₀         = {L0_map:.5f}  (ref: {L0_ref:.3f})")
print(f"    Ψ_ini      = {psi_ini_map:.4f} M_Pl  (near hilltop)")
print(f"    Ψ₀ (today) = {psi0_map:.4f} M_Pl  (SIM113: {PSI0_ref:.4f})")
print(f"    H₀         = {H0_map:.2f} km/s/Mpc")
print(f"    F(z_CMB)   = {F_CMB_map:.5f}  (frozen-Ψ approx: {F0_ref:.5f})")
print(f"    F(z=0)     = {F0_map:.5f}")
print()
print(f"  Fit quality:")
print(f"    χ²_DESI / N = {chi2_DESI_map:.3f} / {len(DESI_z)}  ",
      end="")
print("→ PASS" if chi2_DESI_map/len(DESI_z) < 2 else "→ FAIL")
print(f"    χ²_θ        = {chi2_theta_map:.4f}  ",
      end="")
print("→ PASS" if chi2_theta_map < 4 else "→ FAIL")
print(f"    100θ_*      = {theta_map:.5f}  (obs: {theta_obs:.5f})")
print(f"    DESI tension = {tension_map:.2f}σ  (SIM121C frozen: {SIM121C_tension:.2f}σ)")
print(f"    Δχ²_DESI vs SIM121C = {chi2_DESI_map - SIM121C_chi2:+.3f}")

# H(z) comparison
if bg_map is not None:
    print()
    print(f"  H(z) at MAP: CMSTG running F vs DESI ──")
    print(f"  {'z':>6}  {'H_obs':>8}  {'H_MAP':>8}  {'pull':>6}")
    for z_d, Ho, sd in zip(DESI_z, DESI_H, DESI_s):
        Hm = H_interp(z_d, bg_map)
        print(f"  {z_d:6.3f}  {Ho:8.1f}  {Hm:8.2f}  {(Hm-Ho)/sd:6.2f}")

# Verdict
pass_desi   = chi2_DESI_map/len(DESI_z) < 2
pass_theta  = chi2_theta_map < 4
better_121c = chi2_DESI_map < SIM121C_chi2
delta_chi2  = chi2_DESI_map - SIM121C_chi2

if pass_desi and pass_theta:
    verdict = "PASS"
elif better_121c and pass_theta:
    verdict = "PARTIAL"
elif better_121c:
    verdict = "PARTIAL"
else:
    verdict = "FAIL"

print()
print(f"  VERDICT: {verdict}")
if verdict == "PASS":
    print(f"  Running F(z) resolves the frozen-Ψ normalisation problem!")
    print(f"  F_CMB ≈ {F_CMB_map:.4f} (near GR) while F₀ = {F0_map:.4f} today.")
elif verdict == "PARTIAL":
    print(f"  Running F(z) improves over frozen-Ψ (Δχ²={delta_chi2:+.1f}).")
    if not pass_desi:
        print(f"  DESI χ²/N = {chi2_DESI_map/len(DESI_z):.2f} still exceeds 2.")
    if not pass_theta:
        print(f"  CMB θ_* mismatch ({chi2_theta_map:.1f}σ).")
else:
    print(f"  Running F(z) does not resolve the CMB+DESI tension.")
    print(f"  F_CMB = {F_CMB_map:.5f} (frozen approx was F₀ = {F0_ref:.5f}).")
print("=" * 70)

# ─────────────────────────────────────────────────────────────────────────────
# FIGURES
# ─────────────────────────────────────────────────────────────────────────────
print()
print("Generating figures...")

# --- Figure 1: F(z) evolution ---
fig, axes = plt.subplots(1, 2, figsize=(13, 5))

if bg_run is not None:
    z_p    = bg_run['z']
    mask   = z_p < 2000
    ax = axes[0]
    ax.semilogx(z_p[mask]+1, bg_run['F'][mask], color='#d73027', lw=2,
                label=r'CMSTG running $F(z)$, $\Psi_{\rm ini}=0.01$')
    ax.axhline(F0_ref, color='#2166ac', ls='--', lw=1.5,
               label=f'Frozen-$\\Psi$ approx $F_0={F0_ref:.4f}$')
    ax.axhline(0.5, color='gray', ls=':', lw=1.5, label='GR ($F=1/2$)')
    ax.axvline(z_star+1, color='orange', ls=':', lw=1.5, label=f'$z_*={int(z_star)}$')
    ax.set_xlabel(r'$1+z$')
    ax.set_ylabel(r'$F(z) = \frac{1}{2} + \Lambda_0\Psi^2(z)$')
    ax.set_title(r'CMSTG $F(z)$ evolution (running vs frozen)')
    ax.legend(fontsize=9)
    ax.grid(alpha=0.3)

    ax2 = axes[1]
    ax2.semilogx(z_p[mask]+1, bg_run['psi'][mask], color='#d73027', lw=2,
                 label=r'$\Psi(z)$, $\Psi_{\rm ini}=0.01$')
    ax2.axhline(PSI0_ref, color='#2166ac', ls='--', lw=1.5,
                label=f'$\\Psi_0={PSI0_ref}$ (SIM113)')
    ax2.axhline(V_ref, color='green', ls=':', lw=1.5, label=f'$v={V_ref}$ (VEV)')
    ax2.axvline(z_star+1, color='orange', ls=':', lw=1.5, label=f'$z_*$')
    ax2.set_xlabel(r'$1+z$')
    ax2.set_ylabel(r'$\Psi(z)$ [$M_{\rm Pl}$]')
    ax2.set_title(r'CMSTG $\Psi(z)$ field evolution')
    ax2.legend(fontsize=9)
    ax2.grid(alpha=0.3)

plt.tight_layout()
for ext in ('pdf', 'png'):
    fig.savefig(os.path.join(OUT, f'sim123_F_evolution.{ext}'),
                dpi=150, bbox_inches='tight')
plt.close(fig)
print("  Saved sim123_F_evolution")

# --- Figure 2: θ_* and DESI tension landscape ---
if len(scan_chi2) > 5:
    fig, axes = plt.subplots(1, 2, figsize=(13, 5))

    sc1 = axes[0].scatter(scan_psi_ini, scan_theta, c=scan_L0,
                          cmap='RdYlBu', s=60, alpha=0.8)
    axes[0].axhline(theta_obs, color='k', ls='-', lw=1.5, label=f'Planck $\\theta_*={theta_obs}$')
    axes[0].axhline(theta_obs + 2*theta_obs_err, color='k', ls='--', lw=1, alpha=0.5)
    axes[0].axhline(theta_obs - 2*theta_obs_err, color='k', ls='--', lw=1, alpha=0.5)
    axes[0].set_xlabel(r'$\Psi_{\rm ini}$ [$M_{\rm Pl}$]')
    axes[0].set_ylabel(r'$100\theta_*$')
    axes[0].set_title(r'CMB $\theta_*$ vs $\Psi_{\rm ini}$ (colour = $\Lambda_0$)')
    plt.colorbar(sc1, ax=axes[0], label=r'$\Lambda_0$')
    axes[0].legend(fontsize=9)
    axes[0].grid(alpha=0.3)

    sc2 = axes[1].scatter(scan_psi_ini, scan_tension, c=scan_L0,
                          cmap='RdYlBu', s=60, alpha=0.8)
    axes[1].axhline(SIM121C_tension, color='orange', ls='--', lw=1.5,
                    label=f'SIM121C frozen ref ({SIM121C_tension:.2f}σ)')
    axes[1].axhline(2.0, color='green', ls=':', lw=1.5, label='2σ PASS threshold')
    axes[1].set_xlabel(r'$\Psi_{\rm ini}$ [$M_{\rm Pl}$]')
    axes[1].set_ylabel(r'DESI tension [$\sigma$]')
    axes[1].set_title(r'DESI $H(z)$ tension vs $\Psi_{\rm ini}$')
    plt.colorbar(sc2, ax=axes[1], label=r'$\Lambda_0$')
    axes[1].legend(fontsize=9)
    axes[1].grid(alpha=0.3)

    plt.tight_layout()
    for ext in ('pdf', 'png'):
        fig.savefig(os.path.join(OUT, f'sim123_tension_landscape.{ext}'),
                    dpi=150, bbox_inches='tight')
    plt.close(fig)
    print("  Saved sim123_tension_landscape")

# --- Figure 3: H(z) at MAP vs DESI ---
if bg_map is not None:
    fig, ax = plt.subplots(figsize=(8, 5))

    z_fine = np.linspace(0.01, 2.5, 300)
    H_run_fine = [H_interp(z, bg_map) for z in z_fine]

    # Frozen-Ψ SIM121C reference (F₀=0.560, w₀=-0.60, wₐ=+0.49, H₀=76.0)
    def H_frozen(z):
        F0f = 0.5600
        H0f = 76.0
        h   = H0f/100
        Om  = omh2_m/h**2
        Or  = omh2_r/h**2
        ODE = 1-Om-Or
        w0f, waf = -0.6003, 0.4862
        fDE = (1+z)**(3*(1+w0f+waf)) * np.exp(-3*waf*z/(1+z))
        E2  = (0.5/F0f) * (Om*(1+z)**3 + Or*(1+z)**4 + ODE*fDE)
        return H0f * np.sqrt(max(E2, 0.0))

    H_froz_fine = [H_frozen(z) for z in z_fine]

    ax.plot(z_fine, H_run_fine,  color='#d73027', lw=2,
            label=f'SIM123 MAP (running $F$, $H_0={H0_map:.1f}$)')
    ax.plot(z_fine, H_froz_fine, color='#2166ac', lw=2, ls='--',
            label='SIM121C MAP (frozen $F_0=0.560$, $H_0=76.0$)')
    ax.errorbar(DESI_z, DESI_H, yerr=DESI_s, fmt='ko', ms=7, capsize=4,
                label='DESI Y1 BAO')

    ax.set_xlabel(r'Redshift $z$')
    ax.set_ylabel(r'$H(z)$ [km/s/Mpc]')
    ax.set_title('H(z): SIM123 running F vs SIM121C frozen vs DESI')
    ax.legend()
    ax.grid(alpha=0.3)
    plt.tight_layout()
    for ext in ('pdf', 'png'):
        fig.savefig(os.path.join(OUT, f'sim123_Hz.{ext}'), dpi=150, bbox_inches='tight')
    plt.close(fig)
    print("  Saved sim123_Hz")

# --- Figure 4: F_CMB vs Ψ_ini ---
if len(scan_F_CMB) > 0:
    fig, axes = plt.subplots(1, 2, figsize=(13, 5))

    for L0_val in sorted(set(np.round(scan_L0, 4))):
        mask = np.abs(scan_L0 - L0_val) < 1e-5
        if mask.sum() < 2: continue
        psi_s  = scan_psi_ini[mask]
        Fc_s   = scan_F_CMB[mask]
        idx_s  = np.argsort(psi_s)
        axes[0].plot(psi_s[idx_s], Fc_s[idx_s], 'o-', ms=4, lw=1.5,
                     label=f'$\\Lambda_0={L0_val:.3f}$')

    axes[0].axhline(0.5, color='gray', ls=':', lw=1.5, label='GR ($F=1/2$)')
    axes[0].axhline(F0_ref, color='k', ls='--', lw=1.5, label=f'$F_0={F0_ref:.4f}$ (frozen)')
    axes[0].set_xlabel(r'$\Psi_{\rm ini}$ [$M_{\rm Pl}$]')
    axes[0].set_ylabel(r'$F(z_{\rm CMB})$')
    axes[0].set_title(r'$F$ at CMB decoupling vs $\Psi_{\rm ini}$')
    axes[0].legend(fontsize=8, ncol=2)
    axes[0].grid(alpha=0.3)

    axes[1].scatter(scan_F_CMB, scan_theta, c=scan_psi_ini, cmap='plasma', s=50, alpha=0.8)
    axes[1].axhline(theta_obs, color='k', ls='-', lw=1.5, label=f'Planck $\\theta_*$')
    axes[1].axhline(theta_obs + 3*theta_obs_err, color='k', ls='--', lw=1, alpha=0.5)
    axes[1].axhline(theta_obs - 3*theta_obs_err, color='k', ls='--', lw=1, alpha=0.5)
    axes[1].set_xlabel(r'$F(z_{\rm CMB})$')
    axes[1].set_ylabel(r'$100\theta_*$ (running $F$)')
    axes[1].set_title(r'$\theta_*$ vs $F_{\rm CMB}$')
    axes[1].legend(fontsize=9)
    axes[1].grid(alpha=0.3)
    plt.colorbar(axes[1].collections[0], ax=axes[1], label=r'$\Psi_{\rm ini}$')

    plt.tight_layout()
    for ext in ('pdf', 'png'):
        fig.savefig(os.path.join(OUT, f'sim123_F_CMB.{ext}'), dpi=150, bbox_inches='tight')
    plt.close(fig)
    print("  Saved sim123_F_CMB")

# ─────────────────────────────────────────────────────────────────────────────
# SAVE JSON
# ─────────────────────────────────────────────────────────────────────────────
results = {
    "verdict": verdict,
    "map": {
        "L0": float(L0_map),
        "psi_ini": float(psi_ini_map),
        "psi0": float(psi0_map),
        "F0": float(F0_map),
        "F_CMB": float(F_CMB_map),
        "H0": float(H0_map),
        "theta_star": float(theta_map)
    },
    "fit_quality": {
        "chi2_DESI": float(chi2_DESI_map),
        "chi2_DESI_per_N": float(chi2_DESI_map / len(DESI_z)),
        "chi2_theta": float(chi2_theta_map),
        "DESI_tension_sigma": float(tension_map),
        "delta_chi2_vs_SIM121C": float(chi2_DESI_map - SIM121C_chi2)
    },
    "reference_SIM113": {
        "L0": L0_ref, "psi0": PSI0_ref, "F0": F0_ref, "v": V_ref
    },
    "reference_SIM121C": {
        "chi2_DESI": SIM121C_chi2, "tension": SIM121C_tension,
        "F0_frozen": 0.560
    },
    "pass_desi": bool(pass_desi),
    "pass_theta": bool(pass_theta),
    "better_than_SIM121C": bool(better_121c),
    "frozen_psi_approx_valid": bool(abs(F_CMB_map - F0_map) < 0.001)
}

with open(os.path.join(OUT, 'sim123_results.json'), 'w') as f:
    json.dump(results, f, indent=2)
print("  Saved sim123_results.json")

print(f"\nAll outputs in: {OUT}")
print("SIM123 complete.")
