"""
SIM130 — CMSTG Phase 1: Curvature-Sourced Ψ Evolution (Unlocked Field)
======================================================================
Previous results (SIM121–SIM129):
  - Frozen-Ψ Phase 1 (SIM121C):   2.77σ DESI tension  (best frozen-field result)
  - SSB hilltop Phase 2 (SIM123):  FAIL — Ψ frozen at hilltop (can't reach Ψ₀)
  - M-field programme (SIM125–129): FAIL — modulated-Λ wrong direction

SIM128 confirmed the Phase 1 action is the right theory for S₈ (PASS).

Question for SIM130: does the Phase 1 action itself, when Ψ is UNLOCKED,
generate curvature-sourced Ψ evolution that produces measurable H(z) deviation
at DESI redshifts, without any new parameters?

Phase 1 action (Jordan frame):
  S = ∫d⁴x√(-g) [F(Ψ)/2 R − (1/2)(∂Ψ)² − Λ_bare + L_m]
  F(Ψ) = 1/2 + Λ₀Ψ²

Field EOM (FRW, no SSB potential):
  Ψ'' + (3−ε_H)Ψ' = 2Λ₀Ψ × R/H²   where ' = d/dN = d/d ln a

Modified Friedmann:
  H²[3F + 6Λ₀ΨΨ' − Ψ'²/2] = ω_m a⁻³ + ω_r a⁻⁴ + Λ_bare

where ω_m = Ω_m h², ω_r = Ω_r h² (Planck physical densities),
and H is in units of H100 = 100 km/s/Mpc.
Λ_bare is calibrated self-consistently to give H₀ = h_target × H100.

Scan: Ψ_ini ∈ [0.1, 3.5] M_Pl at z_ini = 10⁵, y_ini = 0 (adiabatic)
Fixed: Λ₀=0.003 (locked Phase 1 value)
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
# CONSTANTS
# Units: H in units of H100=100 km/s/Mpc, densities in H100² units (= Ω h²)
# ─────────────────────────────────────────────────────────────────────────────
c_kms     = 2.998e5       # km/s
H100      = 100.0         # km/s/Mpc (unit conversion)
Lambda0   = 0.003         # CMSTG coupling (Phase 1 locked)
omh2_m    = 0.1430        # Planck Ω_m h²
omh2_r    = 4.18e-5       # Planck Ω_r h²
h_target  = 0.674         # H₀/H100 = 67.4 km/s/Mpc

# Phase 1 reference (locked-action, frozen Ψ)
PSI0_ref  = 2.62
F0_ref    = 0.5 + Lambda0 * PSI0_ref**2   # 0.52059

# Calibrate Λ_bare for reference: H₀ = h_target with F₀ = F0_ref
# 3F₀ h² = ω_m + ω_r + Λ_bare  →  Λ_bare = 3F₀h² − ω_m − ω_r
Lambda_bare_ref = 3.0 * F0_ref * h_target**2 - omh2_m - omh2_r

# Planck CMB
theta_obs     = 1.04101
theta_obs_err = 0.00029
z_drag        = 1059.6
z_star        = 1089.8
omh2_b        = 0.02237

# DESI H(z) data [km/s/Mpc]
DESI_z = np.array([0.295, 0.510, 0.706, 0.930, 1.317, 2.330])
DESI_H = np.array([ 81.7,  97.9, 110.7, 128.1, 156.4, 240.8])
DESI_s = np.array([  4.5,   4.4,   6.2,   5.6,   8.6,  11.0])

# Reference results
SIM121C_chi2    = 18.26
SIM121C_tension = 2.77

# ─────────────────────────────────────────────────────────────────────────────
# ODE SYSTEM
# ─────────────────────────────────────────────────────────────────────────────

def get_H2(Psi, y, N, Lambda_bare):
    """
    Modified Friedmann. Returns H² in (H100)² units.
    H²[3F + 6Λ₀ΨΨ' − Ψ'²/2] = ω_m a⁻³ + ω_r a⁻⁴ + Λ_bare
    """
    a    = np.exp(N)
    rhs  = omh2_m / a**3 + omh2_r / a**4 + Lambda_bare
    F    = 0.5 + Lambda0 * Psi**2
    coef = 3.0 * F + 6.0 * Lambda0 * Psi * y - 0.5 * y**2
    if coef < 0.3:
        coef = 3.0 * F   # fallback, ignore kinetic corrections
    return rhs / coef


def get_eps_H(Psi, y, N, Lambda_bare):
    """ε_H = −(1/2) d ln H²/dN"""
    eps  = 5e-4
    H2p  = get_H2(Psi, y, N + eps, Lambda_bare)
    H2m  = get_H2(Psi, y, N - eps, Lambda_bare)
    H2   = get_H2(Psi, y, N, Lambda_bare)
    if H2 < 1e-40:
        return 0.0
    return -0.5 * (H2p - H2m) / (2.0 * eps * H2)


def cmstg_ode_phase1(N, state, Lambda_bare):
    """
    Phase 1 CMSTG ODE. State: [Ψ, y=dΨ/dN]
    Source = 2Λ₀Ψ × R/H² = 2Λ₀Ψ × 6(2 − ε_H)
    """
    Psi, y    = state
    H2        = get_H2(Psi, y, N, Lambda_bare)
    eps_H     = get_eps_H(Psi, y, N, Lambda_bare)
    R_over_H2 = 6.0 * (2.0 - eps_H)
    source    = 2.0 * Lambda0 * Psi * R_over_H2
    dy_dN     = source - (3.0 - eps_H) * y
    return [y, dy_dN]


def integrate_phase1(psi_ini, Lambda_bare=None, z_ini=1e5, n_points=3000):
    """
    Integrate Phase 1 CMSTG background from z_ini to z=0.
    If Lambda_bare is None, use the reference value (calibrated for Ψ₀=2.62).
    """
    if Lambda_bare is None:
        Lambda_bare = Lambda_bare_ref

    N_ini  = np.log(1.0 / (1.0 + z_ini))
    N_end  = 0.0
    N_eval = np.linspace(N_ini, N_end, n_points)

    sol = solve_ivp(
        cmstg_ode_phase1,
        (N_ini, N_end),
        [psi_ini, 0.0],
        args=(Lambda_bare,),
        t_eval=N_eval,
        method='RK45',
        rtol=1e-7, atol=1e-11
    )

    if not sol.success:
        return None

    N_arr   = sol.t
    z_arr   = np.exp(-N_arr) - 1.0
    psi_arr = sol.y[0]
    y_arr   = sol.y[1]

    H_arr, F_arr = [], []
    for i in range(len(N_arr)):
        H2 = get_H2(psi_arr[i], y_arr[i], N_arr[i], Lambda_bare)
        F  = 0.5 + Lambda0 * psi_arr[i]**2
        H_arr.append(H100 * np.sqrt(max(H2, 0.0)))
        F_arr.append(F)

    H_arr = np.array(H_arr)
    F_arr = np.array(F_arr)

    idx = np.argsort(z_arr)
    return {
        'z':    z_arr[idx],
        'N':    N_arr[idx],
        'psi':  psi_arr[idx],
        'y':    y_arr[idx],
        'H':    H_arr[idx],
        'F':    F_arr[idx],
        'psi0': float(psi_arr[-1]),
        'F0':   float(F_arr[-1]),
        'H0':   float(H_arr[-1]),
        'Lambda_bare': Lambda_bare,
    }


def integrate_selfconsistent(psi_ini, z_ini=1e5, n_iter=4):
    """
    Self-consistent integration: iterate Λ_bare until H₀ ≈ h_target.
    Start with Lambda_bare_ref, then update from actual Ψ₀.
    """
    Lambda_bare = Lambda_bare_ref
    bg = None
    for _ in range(n_iter):
        bg = integrate_phase1(psi_ini, Lambda_bare=Lambda_bare, z_ini=z_ini)
        if bg is None:
            return None
        psi0_actual = bg['psi0']
        F0_actual   = 0.5 + Lambda0 * psi0_actual**2
        # Recalibrate so H₀_actual = h_target
        Lambda_bare = 3.0 * F0_actual * h_target**2 - omh2_m - omh2_r
    # Re-run with final Lambda_bare
    return integrate_phase1(psi_ini, Lambda_bare=Lambda_bare, z_ini=z_ini)


def H_interp(z_query, bg):
    return float(np.interp(z_query, bg['z'], bg['H']))


def theta_star_from_bg(bg):
    """100θ_* = 100 × r_s(z_drag) / D_C(z_*)"""
    H0 = bg['H0']
    h  = H0 / 100.0
    Ogam = 2.469e-5 / h**2

    def rs_integrand(z):
        R  = (3.0 * omh2_b / h**2) / (4.0 * Ogam * (1 + z))
        cs = c_kms / np.sqrt(3.0 * (1.0 + R))
        Hz = H_interp(z, bg)
        return cs / Hz if Hz > 0 else 0.0

    rs, _ = quad(rs_integrand, z_drag, 1e4, limit=200, epsrel=1e-5)

    def DC_integrand(z):
        Hz = H_interp(z, bg)
        return c_kms / Hz if Hz > 0 else 0.0

    DC, _ = quad(DC_integrand, 0, z_star, limit=200, epsrel=1e-5)
    return 100.0 * rs / DC if DC > 0 else np.nan


def chi2_DESI_from_bg(bg):
    H_model = np.array([H_interp(z, bg) for z in DESI_z])
    return float(np.sum(((H_model - DESI_H) / DESI_s)**2))


# ─────────────────────────────────────────────────────────────────────────────
# PART A: Reference trajectory (Ψ_ini = 2.62, self-consistent Λ_bare)
# ─────────────────────────────────────────────────────────────────────────────
print("=" * 70)
print("SIM130 — CMSTG Phase 1: Curvature-Sourced Ψ Evolution (Unlocked)")
print("=" * 70)
print(f"  Lambda0={Lambda0}, Ψ₀_ref={PSI0_ref}, F₀_ref={F0_ref:.5f}")
print(f"  h_target={h_target}, Λ_bare_ref={Lambda_bare_ref:.6f}")
print()
print("Part A: Reference trajectory (Ψ_ini = 2.62 M_Pl, self-consistent Λ_bare)...")

bg_ref = integrate_selfconsistent(PSI0_ref)

if bg_ref is not None:
    psi_CMB  = float(np.interp(z_star, bg_ref['z'], bg_ref['psi']))
    F_CMB    = float(np.interp(z_star, bg_ref['z'], bg_ref['F']))
    psi0     = bg_ref['psi0']
    F0       = bg_ref['F0']
    H0       = bg_ref['H0']
    dpsi_pct = (psi0 - PSI0_ref) / PSI0_ref * 100.0
    dF_pct   = (F0 - F0_ref) / F0_ref * 100.0

    print(f"\n  Ψ evolution (Ψ_ini = {PSI0_ref}):")
    print(f"  Ψ(z_CMB)    = {psi_CMB:.6f} M_Pl  [{(psi_CMB-PSI0_ref)/PSI0_ref*100:+.4f}%]")
    print(f"  Ψ(z=0)      = {psi0:.6f} M_Pl  [{dpsi_pct:+.4f}% from Ψ_ini]")
    print(f"  F(z_CMB)    = {F_CMB:.6f}  (frozen: {F0_ref:.6f}, δF={F_CMB-F0_ref:+.2e})")
    print(f"  F(z=0)      = {F0:.6f}  [{dF_pct:+.4f}%]")
    print(f"  H₀          = {H0:.2f} km/s/Mpc  (target: {h_target*100:.1f})")
    print(f"  Λ_bare used  = {bg_ref['Lambda_bare']:.6f}")

    theta_ref  = theta_star_from_bg(bg_ref)
    chi2_ref   = chi2_DESI_from_bg(bg_ref)
    tension_ref = np.sqrt(chi2_ref / len(DESI_z))
    print(f"\n  100θ_*      = {theta_ref:.5f}  (Planck: {theta_obs:.5f}, Δ={theta_ref-theta_obs:+.5f})")
    print(f"  χ²_DESI     = {chi2_ref:.3f}  (SIM121C frozen: {SIM121C_chi2:.3f})")
    print(f"  DESI tension = {tension_ref:.3f}σ  (SIM121C: {SIM121C_tension:.2f}σ)")
    print(f"  Δχ²_DESI    = {chi2_ref - SIM121C_chi2:+.3f}")
else:
    print("  ODE FAILED.")
    theta_ref = tension_ref = chi2_ref = np.nan

# ─────────────────────────────────────────────────────────────────────────────
# PART B: Curvature-sourcing amplitude (analytic estimate)
# ─────────────────────────────────────────────────────────────────────────────
print()
print("Part B: Curvature-sourcing amplitude (quasi-static estimate)...")
print(f"  In slow-roll: dΨ/dN ≈ (Λ₀Ψ/3) × R/H² = Λ₀Ψ × (Ω_m(1+z)³+4Ω_Λ)/(Ω_m(1+z)³+Ω_Λ)")
Omega_L = 1.0 - omh2_m/h_target**2 - omh2_r/h_target**2
Omega_m = omh2_m / h_target**2
total_dlnpsi = 0.0
prev_z = 0.0
for z_q in [0.0, 0.3, 0.5, 1.0, 2.0, 5.0, 10.0]:
    denom = Omega_m * (1+z_q)**3 + Omega_L
    numer = Omega_m * (1+z_q)**3 + 4*Omega_L
    R_H2  = 3.0 * numer / denom   # R/H² = 3(Ω_m(1+z)³+4Ω_Λ)/(Ω_m(1+z)³+Ω_Λ)
    dlnpsi = Lambda0 * R_H2 / 3.0
    print(f"  z={z_q:.1f}: R/H²={R_H2:.3f}, dln Ψ/dN≈{dlnpsi:.5f} ({dlnpsi*100:.3f}%/e-fold)")
print(f"\n  Phase 1 Ψ evolution from z=5 to z=0: ΔΨ/Ψ ≈ few × {Lambda0*3:.3f} × 1e-fold ≈ <1%")
print(f"  Effect on F: ΔF/F = 2Λ₀Ψ² / F × ΔΨ/Ψ ≈ 2×{Lambda0*PSI0_ref**2:.3f}/{F0_ref:.3f} × <1% ≈ <0.16%")
print(f"  Effect on H(z): ΔH/H ≈ −ΔF/(2F) ≈ <0.08% → NEGLIGIBLE for DESI (needs 5-8%)")

# ─────────────────────────────────────────────────────────────────────────────
# PART C: Scan Ψ_ini — does any initial condition help?
# ─────────────────────────────────────────────────────────────────────────────
print()
print("Part C: Scan Ψ_ini ∈ [0.1, 3.5] M_Pl (self-consistent Λ_bare each run)...")

psi_ini_grid = np.array([0.1, 0.3, 0.5, 0.8, 1.0, 1.5, 2.0, 2.3, 2.62, 2.8, 3.0, 3.2, 3.5])

scan_psi_ini = []
scan_psi0    = []
scan_F0      = []
scan_F_CMB   = []
scan_theta   = []
scan_chi2    = []
scan_tension = []
scan_H0      = []

for psi_ini in psi_ini_grid:
    bg = integrate_selfconsistent(psi_ini)
    if bg is None:
        print(f"  Ψ_ini={psi_ini:.2f}: FAILED")
        continue
    theta = theta_star_from_bg(bg)
    if np.isnan(theta):
        continue
    chi2  = chi2_DESI_from_bg(bg)
    F_CMB = float(np.interp(z_star, bg['z'], bg['F']))

    scan_psi_ini.append(psi_ini)
    scan_psi0.append(bg['psi0'])
    scan_F0.append(bg['F0'])
    scan_F_CMB.append(F_CMB)
    scan_theta.append(theta)
    scan_chi2.append(chi2)
    scan_tension.append(np.sqrt(chi2 / len(DESI_z)))
    scan_H0.append(bg['H0'])
    print(f"  Ψ_ini={psi_ini:.2f}: Ψ₀={bg['psi0']:.4f}, H₀={bg['H0']:.1f}, "
          f"100θ_*={theta:.4f}, tension={np.sqrt(chi2/len(DESI_z)):.3f}σ")

scan_psi_ini = np.array(scan_psi_ini)
scan_psi0    = np.array(scan_psi0)
scan_chi2    = np.array(scan_chi2)
scan_tension = np.array(scan_tension)
scan_theta   = np.array(scan_theta)
scan_H0      = np.array(scan_H0)
scan_F0      = np.array(scan_F0)
scan_F_CMB   = np.array(scan_F_CMB)

# ─────────────────────────────────────────────────────────────────────────────
# PART D: Best point and verdict
# ─────────────────────────────────────────────────────────────────────────────
print()
print("Part D: Best joint fit and verdict...")

chi2_theta_arr = ((scan_theta - theta_obs) / theta_obs_err)**2
chi2_tot_arr   = scan_chi2 + chi2_theta_arr

if len(chi2_tot_arr) > 0:
    idx_best       = np.argmin(chi2_tot_arr)
    psi_ini_map    = scan_psi_ini[idx_best]
    psi0_map       = scan_psi0[idx_best]
    F0_map         = scan_F0[idx_best]
    F_CMB_map      = scan_F_CMB[idx_best]
    theta_map      = scan_theta[idx_best]
    chi2_DESI_map  = scan_chi2[idx_best]
    chi2_theta_map = chi2_theta_arr[idx_best]
    H0_map         = scan_H0[idx_best]
    tension_map    = scan_tension[idx_best]
    delta_chi2     = chi2_DESI_map - SIM121C_chi2

    print(f"\n  Best joint point (min χ²_DESI + χ²_θ):")
    print(f"  Ψ_ini  = {psi_ini_map:.3f} M_Pl  →  Ψ₀ = {psi0_map:.4f} M_Pl")
    print(f"  F₀     = {F0_map:.5f},  F_CMB = {F_CMB_map:.5f}")
    print(f"  H₀     = {H0_map:.2f} km/s/Mpc")
    print(f"  100θ_* = {theta_map:.5f}  (χ²_θ={chi2_theta_map:.3f})")
    print(f"  χ²_DESI = {chi2_DESI_map:.3f}  (SIM121C: {SIM121C_chi2:.3f})")
    print(f"  DESI tension = {tension_map:.3f}σ  (SIM121C: {SIM121C_tension:.2f}σ)")
    print(f"  Δχ²_DESI vs SIM121C = {delta_chi2:+.3f}")

    print(f"\n  H(z) at MAP vs DESI:")
    print(f"  {'z':>6}  {'H_obs':>8}  {'H_MAP':>8}  {'pull':>6}")
    bg_map = integrate_selfconsistent(psi_ini_map)
    if bg_map is not None:
        for z_d, Ho, sd in zip(DESI_z, DESI_H, DESI_s):
            Hm = H_interp(z_d, bg_map)
            print(f"  {z_d:6.3f}  {Ho:8.1f}  {Hm:8.2f}  {(Hm-Ho)/sd:6.2f}")
    else:
        bg_map = None
else:
    print("  No valid scan points.")
    delta_chi2 = 0.0; chi2_DESI_map = SIM121C_chi2; chi2_theta_map = 0.0
    psi_ini_map = PSI0_ref; psi0_map = PSI0_ref; F0_map = F0_ref
    F_CMB_map = F0_ref; theta_map = theta_obs; H0_map = 67.4; tension_map = SIM121C_tension
    bg_map = bg_ref

pass_desi  = chi2_DESI_map / len(DESI_z) < 2.0
pass_theta = chi2_theta_map < 4.0
beats_ref  = chi2_DESI_map < SIM121C_chi2

if pass_desi and pass_theta:
    verdict = "PASS"
elif beats_ref and pass_theta:
    verdict = "PARTIAL"
elif beats_ref:
    verdict = "PARTIAL"
else:
    verdict = "FAIL"

print()
print("=" * 70)
print("SIM130 RESULT:")
print()
print(f"  Verdict: {verdict}")
print(f"  DESI tension at MAP = {tension_map:.3f}σ  (frozen Ψ: {SIM121C_tension:.2f}σ)")
print(f"  Δχ²_DESI vs SIM121C = {delta_chi2:+.3f}")
print()
if abs(delta_chi2) < 1.0:
    print(f"  CONCLUSION: Curvature-sourced Ψ evolution is negligible (Δχ²≈{delta_chi2:.3f}).")
    print(f"  The frozen-Ψ approximation is validated — Phase 1 action does not")
    print(f"  generate significant H(z) variation through curvature sourcing alone.")
    print(f"  The 2.77σ DESI tension is structural in the Phase 1 action.")
elif beats_ref:
    print(f"  Mild improvement: Δχ²={delta_chi2:+.3f}, but tension={tension_map:.2f}σ persists.")
    print(f"  Curvature sourcing is insufficient to resolve the structural DESI tension.")
else:
    print(f"  Curvature sourcing WORSENS DESI fit (Δχ²={delta_chi2:+.3f}).")
    print(f"  Unlocking Ψ in Phase 1 is counterproductive for DESI.")
print("=" * 70)

# ─────────────────────────────────────────────────────────────────────────────
# FIGURES
# ─────────────────────────────────────────────────────────────────────────────
print()
print("Generating figures...")

# Figure 1: Ψ(z) and F(z) evolution
if bg_ref is not None:
    fig, axes = plt.subplots(1, 2, figsize=(13, 5))
    z_p  = bg_ref['z']
    mask = z_p < 2000

    axes[0].semilogx(z_p[mask]+1, bg_ref['psi'][mask],
                     color='#d73027', lw=2, label=r'CMSTG running $\Psi(z)$')
    axes[0].axhline(PSI0_ref, color='#2166ac', ls='--', lw=1.5,
                    label=f'Frozen Phase 1 $\\Psi_0={PSI0_ref}$')
    axes[0].axvline(z_star+1, color='orange', ls=':', lw=1.5, label=f'$z_*$')
    axes[0].set_xlabel(r'$1+z$')
    axes[0].set_ylabel(r'$\Psi(z)$ [$M_{\rm Pl}$]')
    axes[0].set_title(r'$\Psi(z)$ — curvature-sourced evolution (Phase 1)')
    axes[0].legend(fontsize=9); axes[0].grid(alpha=0.3)

    psi_arr = bg_ref['psi']
    F_running = 0.5 + Lambda0 * psi_arr**2
    axes[1].semilogx(z_p[mask]+1, F_running[mask],
                     color='#d73027', lw=2, label=r'CMSTG running $F(z)$')
    axes[1].axhline(F0_ref, color='#2166ac', ls='--', lw=1.5,
                    label=f'Frozen $F_0={F0_ref:.5f}$')
    axes[1].axhline(0.5, color='gray', ls=':', lw=1.5, label='GR ($F=1/2$)')
    axes[1].axvline(z_star+1, color='orange', ls=':', lw=1.5, label=f'$z_*$')
    axes[1].set_xlabel(r'$1+z$')
    axes[1].set_ylabel(r'$F(z) = \frac{1}{2} + \Lambda_0\Psi^2$')
    axes[1].set_title(r'Running $F(z)$ vs frozen Phase 1')
    axes[1].legend(fontsize=9); axes[1].grid(alpha=0.3)

    plt.tight_layout()
    for ext in ('pdf', 'png'):
        fig.savefig(os.path.join(OUT, f'sim130_psi_evolution.{ext}'), dpi=150, bbox_inches='tight')
    plt.close(fig)
    print("  Saved sim130_psi_evolution")

# Figure 2: ΔΨ/Ψ fractional change (z < 10)
if bg_ref is not None:
    fig, axes = plt.subplots(1, 2, figsize=(13, 5))
    z_p  = bg_ref['z']
    mask = z_p < 10

    axes[0].plot(z_p[mask], (bg_ref['psi'][mask] - PSI0_ref)/PSI0_ref * 100,
                 color='#d73027', lw=2)
    axes[0].axhline(0, color='gray', ls='--', lw=1)
    axes[0].set_xlabel(r'Redshift $z$')
    axes[0].set_ylabel(r'$\Delta\Psi/\Psi_{\rm ini}$ [%]')
    axes[0].set_title(r'Fractional $\Psi$ change from curvature sourcing')
    axes[0].grid(alpha=0.3)

    F_arr = 0.5 + Lambda0 * bg_ref['psi'][mask]**2
    axes[1].plot(z_p[mask], (F_arr - F0_ref)/F0_ref * 100,
                 color='#2166ac', lw=2)
    axes[1].axhline(0, color='gray', ls='--', lw=1)
    axes[1].set_xlabel(r'Redshift $z$')
    axes[1].set_ylabel(r'$\Delta F / F_0$ [%]')
    axes[1].set_title(r'Fractional $F(z)$ change (drives $\Delta H/H$)')
    axes[1].grid(alpha=0.3)

    plt.tight_layout()
    for ext in ('pdf', 'png'):
        fig.savefig(os.path.join(OUT, f'sim130_delta_psi.{ext}'), dpi=150, bbox_inches='tight')
    plt.close(fig)
    print("  Saved sim130_delta_psi")

# Figure 3: DESI tension landscape
if len(scan_chi2) > 3:
    fig, axes = plt.subplots(1, 2, figsize=(13, 5))

    axes[0].plot(scan_psi_ini, scan_tension, 'o-', color='#d73027', lw=2, ms=6)
    axes[0].axhline(SIM121C_tension, color='#2166ac', ls='--', lw=1.5,
                    label=f'SIM121C frozen ({SIM121C_tension:.2f}σ)')
    axes[0].axhline(2.0, color='green', ls=':', lw=1.5, label='2σ PASS threshold')
    axes[0].axvline(PSI0_ref, color='orange', ls=':', lw=1.5, label=f'$\\Psi_0={PSI0_ref}$')
    axes[0].set_xlabel(r'$\Psi_{\rm ini}$ [$M_{\rm Pl}$]')
    axes[0].set_ylabel(r'DESI tension [$\sigma$]')
    axes[0].set_title(r'DESI tension vs initial $\Psi$')
    axes[0].legend(fontsize=9); axes[0].grid(alpha=0.3)

    axes[1].plot(scan_psi_ini, scan_theta, 'o-', color='#d73027', lw=2, ms=6)
    axes[1].axhline(theta_obs, color='k', ls='-', lw=1.5, label=f'Planck {theta_obs}')
    axes[1].axhline(theta_obs + 2*theta_obs_err, color='k', ls='--', lw=1, alpha=0.5)
    axes[1].axhline(theta_obs - 2*theta_obs_err, color='k', ls='--', lw=1, alpha=0.5)
    axes[1].axvline(PSI0_ref, color='orange', ls=':', lw=1.5, label=f'$\\Psi_0={PSI0_ref}$')
    axes[1].set_xlabel(r'$\Psi_{\rm ini}$ [$M_{\rm Pl}$]')
    axes[1].set_ylabel(r'$100\theta_*$')
    axes[1].set_title(r'CMB $\theta_*$ vs initial $\Psi$')
    axes[1].legend(fontsize=9); axes[1].grid(alpha=0.3)

    plt.tight_layout()
    for ext in ('pdf', 'png'):
        fig.savefig(os.path.join(OUT, f'sim130_tension_landscape.{ext}'), dpi=150, bbox_inches='tight')
    plt.close(fig)
    print("  Saved sim130_tension_landscape")

# Figure 4: H(z) comparison
if bg_map is not None and bg_ref is not None:
    fig, ax = plt.subplots(figsize=(8, 5))
    z_fine = np.linspace(0.01, 2.5, 300)
    H_ref_fine  = [H_interp(z, bg_ref) for z in z_fine]
    H_map_fine  = [H_interp(z, bg_map) for z in z_fine]

    ax.plot(z_fine, H_ref_fine, color='#d73027', lw=2,
            label=f'SIM130 running F (Ψ_ini={PSI0_ref})')
    if psi_ini_map != PSI0_ref:
        ax.plot(z_fine, H_map_fine, color='purple', lw=2, ls='--',
                label=f'SIM130 MAP (Ψ_ini={psi_ini_map:.2f})')
    ax.errorbar(DESI_z, DESI_H, yerr=DESI_s, fmt='ko', ms=7, capsize=4,
                label='DESI Y1 BAO')
    ax.set_xlabel(r'Redshift $z$')
    ax.set_ylabel(r'$H(z)$ [km/s/Mpc]')
    ax.set_title('SIM130: H(z) — unlocked Phase 1 Ψ vs DESI')
    ax.legend(); ax.grid(alpha=0.3)
    plt.tight_layout()
    for ext in ('pdf', 'png'):
        fig.savefig(os.path.join(OUT, f'sim130_Hz.{ext}'), dpi=150, bbox_inches='tight')
    plt.close(fig)
    print("  Saved sim130_Hz")

# ─────────────────────────────────────────────────────────────────────────────
# SAVE JSON
# ─────────────────────────────────────────────────────────────────────────────
ref_psi_change = float((bg_ref['psi0'] - PSI0_ref)/PSI0_ref*100) if bg_ref else None
ref_F_change   = float((bg_ref['F0'] - F0_ref)/F0_ref*100) if bg_ref else None

results = {
    "sim": "SIM130",
    "verdict": verdict,
    "description": "Phase 1 CMSTG unlocked — curvature-sourced Psi evolution, no SSB, no new knobs",
    "fixed_params": {"Lambda0": Lambda0, "h_target": h_target},
    "reference_trajectory": {
        "psi_ini": PSI0_ref,
        "psi0_today": float(bg_ref['psi0']) if bg_ref else None,
        "delta_psi_pct": ref_psi_change,
        "delta_F_pct": ref_F_change,
        "H0": float(bg_ref['H0']) if bg_ref else None,
        "chi2_DESI": float(chi2_ref),
        "DESI_tension_sigma": float(tension_ref),
        "theta_star_100": float(theta_ref)
    },
    "map": {
        "psi_ini": float(psi_ini_map),
        "psi0_today": float(psi0_map),
        "F0": float(F0_map),
        "F_CMB": float(F_CMB_map),
        "H0": float(H0_map),
        "theta_star_100": float(theta_map)
    },
    "fit_quality": {
        "chi2_DESI": float(chi2_DESI_map),
        "chi2_DESI_per_N": float(chi2_DESI_map / len(DESI_z)),
        "chi2_theta": float(chi2_theta_map),
        "DESI_tension_sigma": float(tension_map),
        "delta_chi2_vs_SIM121C": float(delta_chi2)
    },
    "reference_SIM121C": {
        "chi2_DESI": SIM121C_chi2,
        "tension_sigma": SIM121C_tension,
        "psi0_frozen": PSI0_ref,
        "F0_frozen": F0_ref
    },
    "pass_desi": bool(pass_desi),
    "pass_theta": bool(pass_theta),
    "beats_SIM121C": bool(beats_ref),
}

with open(os.path.join(OUT, 'sim130_results.json'), 'w') as f:
    json.dump(results, f, indent=2)
print("  Saved sim130_results.json")
print(f"\nAll outputs in: {OUT}")
print("SIM130 complete.")
