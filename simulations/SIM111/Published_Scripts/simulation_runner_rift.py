#!/usr/bin/env python3
"""
SIM111: Dynamical Ψ — m₀ Scan for DESI w₀-wₐ Tension Resolution
RIFT: Recursive Intelligence-Field Theory

SIM109 found that with m₀ = 0.001 H₀ (slow-roll), RIFT predicts
w₀ = −0.992, wₐ = −0.082: PASS vs Planck (1.3σ) but 3.6σ tension with
DESI Y1 (w₀ = −0.827 ± 0.060, wₐ = −0.75 ± 0.29).

SIM109 concluded: "To match DESI Y1 would require m₀ > H₀ (dynamical Ψ)."

This simulation scans m₀ ∈ [0.001, 10] H₀ and finds:
  1. Whether any m₀ resolves the DESI tension (< 2σ) without breaking Planck (< 2σ).
  2. The optimal m₀ and resulting (w₀, wₐ) CPL parameters.
  3. Physical interpretation — thawing vs frozen quintessence regimes.

Background equations (units 8πG = 1, H₀ = 1):

  Modified Friedmann:
    H²[3(1 + 2Λ₀Ψ²) + 6Λ₀Ψy − ½y²]
      = ρ_m a⁻³ + ρ_r a⁻⁴ + ρ_Λ + ½m₀²Ψ²

  Modified Klein-Gordon (x = ln a):
    d²Ψ/dx² + (3 − ε_H) dΨ/dx + m₀²Ψ/H² = 2Λ₀Ψ R/H²
    where R/H² = 6(2 − ε_H)
"""

import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import curve_fit, minimize_scalar
import json, os, sys

# ─── Fixed cosmological parameters ────────────────────────────────────────────
Lambda0 = 0.003
H0      = 1.0
Omega_m = 0.315
Omega_r = 9.1e-5
Omega_L = 0.685
Psi0_ic = 1.0          # Ψ(z=0); held fixed across scan

rho_m0 = Omega_m * 3.0
rho_r0 = Omega_r * 3.0
rho_L  = Omega_L * 3.0

# DESI Y1 (2024) CPL posteriors
DESI_w0,  DESI_wa  = -0.827, -0.75
DESI_sw0, DESI_swa = 0.060,   0.29

# Planck 2018 w₀CDM
PLANCK_w0,  PLANCK_sw0 = -1.03, 0.03

# ─── m₀ scan grid ─────────────────────────────────────────────────────────────
M0_SCAN = np.array([
    0.001, 0.003, 0.01, 0.03,
    0.1,   0.2,   0.3,  0.5,
    0.7,   1.0,   1.5,  2.0,
    3.0,   5.0,  10.0,
])


# ─── Background solver for a given m₀ ────────────────────────────────────────

def get_H2(Psi, y, x, m0):
    a    = np.exp(x)
    rhs  = rho_m0*a**(-3) + rho_r0*a**(-4) + rho_L + 0.5*m0**2*Psi**2
    coef = 3.0*(1.0 + 2.0*Lambda0*Psi**2) + 6.0*Lambda0*Psi*y - 0.5*y**2
    if coef <= 0.5:
        coef = 3.0
    return rhs / coef


def get_epsilon_H(Psi, y, x, m0):
    eps = 5e-4
    H2p = get_H2(Psi, y, x + eps, m0)
    H2m = get_H2(Psi, y, x - eps, m0)
    H2  = get_H2(Psi, y, x,       m0)
    if H2 < 1e-40:
        return 0.0
    return -0.5 * (H2p - H2m) / (2.0*eps*H2)


def rift_odes(x, state, m0):
    Psi, y = state
    H2     = get_H2(Psi, y, x, m0)
    eps_H  = get_epsilon_H(Psi, y, x, m0)

    R_over_H2 = 6.0 * (2.0 - eps_H)
    source    = 2.0*Lambda0*Psi*R_over_H2 - m0**2*Psi / max(H2, 1e-40)
    dy_dx     = source - (3.0 - eps_H)*y

    return [y, dy_dx]


def cpl(z, w0, wa):
    return w0 + wa * z / (1.0 + z)


def run_m0(m0, verbose=False):
    """
    Integrate RIFT background for given m₀.
    Returns dict with w₀, wₐ, tensions, and diagnostics.
    Returns None if integration fails.
    """
    # Initial condition: generalised slow-roll at z=0
    # 3H·Ψ̇ ≈ (2Λ₀R − m₀²)Ψ  => y₀ = (2Λ₀·12 − m₀²)·Ψ₀ / 3
    # For large m₀ this becomes very negative; cap y₀ to avoid blow-up
    R0_over_H2 = 12.0   # de Sitter-like at z=0
    y0_sr  = (2.0*Lambda0*R0_over_H2 - m0**2) * Psi0_ic / 3.0
    y0     = np.clip(y0_sr, -5.0, 5.0)

    state_ic = [Psi0_ic, y0]
    x_end    = -np.log(6.0)   # z=5

    try:
        sol = solve_ivp(
            lambda x, s: rift_odes(x, s, m0),
            [0.0, x_end],
            state_ic,
            method='DOP853',
            dense_output=True,
            rtol=1e-9, atol=1e-11,
            max_step=0.002,
        )
        if not sol.success:
            return None
    except Exception:
        return None

    # Evaluate on redshift grid
    z_grid = np.array([0.0, 0.1, 0.2, 0.3, 0.5, 0.7, 1.0, 1.5, 2.0, 3.0, 5.0])
    x_grid = -np.log(1.0 + z_grid)
    a_grid = 1.0 / (1.0 + z_grid)

    H_z, H_lcdm_z, Psi_z, y_z = [], [], [], []

    for x in x_grid:
        s    = sol.sol(x)
        P, yy = s[0], s[1]
        H2   = get_H2(P, yy, x, m0)
        H_z.append(np.sqrt(max(H2, 0.0)))
        a    = np.exp(x)
        H_l  = np.sqrt(Omega_m*a**(-3) + Omega_r*a**(-4) + Omega_L)
        H_lcdm_z.append(H_l)
        Psi_z.append(P)
        y_z.append(yy)

    H_z     = np.array(H_z)
    H_lcdm  = np.array(H_lcdm_z)
    Psi_z   = np.array(Psi_z)
    y_z     = np.array(y_z)

    # Effective DE density and w_eff(z)
    rho_m_z = rho_m0 * (1.0 + z_grid)**3
    rho_r_z = rho_r0 * (1.0 + z_grid)**4
    rho_DE  = 3.0*H_z**2 - rho_m_z - rho_r_z

    # Guard: rho_DE must be positive everywhere for a sensible w extraction
    if np.any(rho_DE <= 0):
        return None

    w_eff = np.full_like(z_grid, -1.0)
    for i in range(1, len(z_grid)-1):
        dlna   = np.log(a_grid[i+1]/a_grid[i-1])
        dlnrho = np.log(rho_DE[i+1]/rho_DE[i-1])
        w_eff[i] = -1.0 - dlnrho / (3.0 * dlna)
    w_eff[0]  = w_eff[1]
    w_eff[-1] = w_eff[-2]

    # CPL fit over z ∈ [0, 2]
    mask = z_grid <= 2.0
    try:
        popt, _ = curve_fit(cpl, z_grid[mask], w_eff[mask], p0=[-1.0, 0.0])
        w0_fit, wa_fit = float(popt[0]), float(popt[1])
    except Exception:
        w0_fit, wa_fit = float(w_eff[0]), 0.0

    # Tension calculations
    pull_planck = (w0_fit - PLANCK_w0) / PLANCK_sw0
    pull_desi_w0 = (w0_fit - DESI_w0) / DESI_sw0
    pull_desi_wa = (wa_fit - DESI_wa)  / DESI_swa
    chi2_desi    = pull_desi_w0**2 + pull_desi_wa**2
    desi_tension = np.sqrt(chi2_desi)

    planck_pass = abs(pull_planck) < 2.0
    desi_ok     = desi_tension < 2.0

    delta_Psi = abs(Psi_z[-1] - Psi_z[0]) / abs(Psi_z[0]) * 100.0

    return {
        'm0'              : float(m0),
        'y0_sr'           : float(y0_sr),
        'w0'              : w0_fit,
        'wa'              : wa_fit,
        'pull_planck'     : float(pull_planck),
        'pull_desi_w0'    : float(pull_desi_w0),
        'pull_desi_wa'    : float(pull_desi_wa),
        'desi_tension'    : float(desi_tension),
        'planck_pass'     : bool(planck_pass),
        'desi_ok'         : bool(desi_ok),
        'delta_Psi_pct'   : float(delta_Psi),
        'dPsi_dx_z0'      : float(y_z[0]),
        'Psi_z0'          : float(Psi_z[0]),
        'Psi_z5'          : float(Psi_z[-1]),
        'H_RIFT_z0'       : float(H_z[0]),
        'H_LCDM_z0'       : float(H_lcdm[0]),
        'dH_z0_pct'       : float((H_z[0]-H_lcdm[0])/H_lcdm[0]*100.0),
        'w_eff_grid'      : [float(w) for w in w_eff],
        'z_grid'          : list(z_grid),
        'H_ratio'         : [float(h) for h in H_z/H_lcdm],
    }


# ─── Run scan ─────────────────────────────────────────────────────────────────
print("="*76)
print("SIM111: Dynamical Ψ — m₀ Scan for DESI w₀-wₐ Tension Resolution")
print("="*76)
print(f"Lambda0={Lambda0},  Ψ₀(z=0)={Psi0_ic},  scan {len(M0_SCAN)} m₀ values")
print()
print(f"Targets:")
print(f"  Planck 2018:  w₀ = {PLANCK_w0:.3f} ± {PLANCK_sw0:.3f},  wₐ = 0 (fixed)")
print(f"  DESI Y1 2024: w₀ = {DESI_w0:.3f} ± {DESI_sw0:.3f},  wₐ = {DESI_wa:.2f} ± {DESI_swa:.2f}")
print()

results = []
for m0 in M0_SCAN:
    r = run_m0(m0)
    if r is None:
        print(f"  m₀ = {m0:6.3f} H₀  → integration FAILED")
    else:
        results.append(r)

# ─── Summary table ────────────────────────────────────────────────────────────
print("─── Scan Results ───")
hdr = f"{'m₀/H₀':>8}  {'w₀':>8}  {'wₐ':>7}  {'DESI σ':>8}  {'Planck σ':>9}  {'ΔΨ%':>8}  {'Planck':>7}  {'DESI':>6}"
print(hdr)
print("-"*76)
for r in results:
    p_flag = "PASS" if r['planck_pass'] else "FAIL"
    d_flag = "PASS" if r['desi_ok']    else f"{r['desi_tension']:.1f}σ"
    print(f"{r['m0']:>8.3f}  {r['w0']:>8.4f}  {r['wa']:>7.4f}  "
          f"{r['desi_tension']:>8.2f}  {abs(r['pull_planck']):>9.2f}  "
          f"{r['delta_Psi_pct']:>8.4f}  {p_flag:>7}  {d_flag:>6}")

print()

# ─── Find best-fit m₀ for DESI ────────────────────────────────────────────────
# Minimise DESI tension among Planck-passing results
planck_ok = [r for r in results if r['planck_pass']]
if planck_ok:
    best_desi = min(planck_ok, key=lambda r: r['desi_tension'])
else:
    best_desi = min(results, key=lambda r: r['desi_tension']) if results else None

# Global minimum DESI tension (ignoring Planck)
global_best = min(results, key=lambda r: r['desi_tension']) if results else None

# ─── Fine-grained scan around the best m₀ ────────────────────────────────────
fine_results = []
if best_desi is not None:
    m0_center = best_desi['m0']
    m0_fine   = np.linspace(max(0.001, m0_center*0.3), m0_center*3.0, 30)
    print(f"─── Fine scan around m₀ = {m0_center:.3f} H₀ ───")
    for m0 in m0_fine:
        r = run_m0(m0)
        if r is not None:
            fine_results.append(r)

    if fine_results:
        fine_planck = [r for r in fine_results if r['planck_pass']]
        if fine_planck:
            best_fine = min(fine_planck, key=lambda r: r['desi_tension'])
        else:
            best_fine = min(fine_results, key=lambda r: r['desi_tension'])

        print(f"  Best (Planck-consistent): m₀ = {best_fine['m0']:.4f} H₀  "
              f"w₀ = {best_fine['w0']:+.4f}  wₐ = {best_fine['wa']:+.4f}  "
              f"DESI {best_fine['desi_tension']:.2f}σ  Planck {abs(best_fine['pull_planck']):.2f}σ")
    else:
        best_fine = best_desi
else:
    best_fine = None

print()

# ─── Verdict ──────────────────────────────────────────────────────────────────
print("="*76)
print("VERDICT")
print("="*76)

if best_fine is not None:
    bf = best_fine
    print(f"  Optimal m₀ = {bf['m0']:.4f} H₀")
    print(f"  w₀ (RIFT)  = {bf['w0']:+.4f}")
    print(f"  wₐ (RIFT)  = {bf['wa']:+.4f}")
    print()
    print(f"  vs Planck: pull = {abs(bf['pull_planck']):.2f}σ  → {'PASS' if bf['planck_pass'] else 'FAIL'}")
    print(f"  vs DESI:   {bf['desi_tension']:.2f}σ  → {'PASS (<2σ)' if bf['desi_ok'] else 'TENSION'}")
    print()

    if bf['desi_ok'] and bf['planck_pass']:
        verdict = "PASS — DESI tension resolved within Planck consistency"
    elif bf['desi_ok']:
        verdict = "PARTIAL — DESI resolved but Planck tension introduced"
    else:
        # Report best achievable
        min_tension = min(r['desi_tension'] for r in results if r['planck_pass']) if planck_ok else None
        if min_tension is not None:
            verdict = f"FAIL — minimum DESI tension achievable (Planck-consistent) = {min_tension:.2f}σ"
        else:
            verdict = "FAIL — cannot simultaneously satisfy Planck and reduce DESI tension"

    print(f"  RESULT: {verdict}")
    print()
    print(f"  Physical regime (m₀ = {bf['m0']:.4f} H₀):")
    if bf['m0'] < 0.1:
        regime = "slow-roll / frozen (w ≈ −1, Planck-like)"
    elif bf['m0'] < 1.0:
        regime = "thawing quintessence (intermediate dynamical Ψ)"
    elif bf['m0'] < 3.0:
        regime = "dynamical / fast-roll (m₀ ≳ H₀, Ψ evolves significantly)"
    else:
        regime = "oscillating / matter-like (m₀ >> H₀, Ψ averages to w=0)"
    print(f"  {regime}")
    print(f"  ΔΨ (z=0→5) = {bf['delta_Psi_pct']:.4f}%")
    print(f"  dΨ/dx|z=0  = {bf['dPsi_dx_z0']:.6f}")

    # Tension map summary
    print()
    print("─── Tension landscape (Planck-consistent results) ───")
    print(f"  {'m₀/H₀':>8}  {'DESI σ':>8}  {'w₀':>8}  {'wₐ':>7}")
    print("  " + "-"*38)
    for r in sorted(planck_ok, key=lambda r: r['desi_tension'])[:8]:
        print(f"  {r['m0']:>8.4f}  {r['desi_tension']:>8.2f}  {r['w0']:>8.4f}  {r['wa']:>7.4f}")

# ─── Save diagnostics ─────────────────────────────────────────────────────────
out_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'Outputs')
os.makedirs(out_dir, exist_ok=True)

diag = {
    'sim'         : 'SIM111',
    'description' : 'Dynamical Psi m0 scan for DESI w0-wa tension resolution',
    'Lambda0'     : Lambda0,
    'Psi0_ic'     : Psi0_ic,
    'scan_results': results,
    'fine_results': fine_results,
    'best_m0'     : best_fine,
    'desi_target' : {'w0': DESI_w0, 'wa': DESI_wa, 'sw0': DESI_sw0, 'swa': DESI_swa},
    'planck_target': {'w0': PLANCK_w0, 'sw0': PLANCK_sw0},
}

out_path = os.path.join(out_dir, 'sim111_diagnostics.json')
with open(out_path, 'w') as f:
    json.dump(diag, f, indent=2)

print()
print(f"Diagnostics saved → Outputs/sim111_diagnostics.json")
