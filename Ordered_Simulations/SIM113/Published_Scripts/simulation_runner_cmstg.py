#!/usr/bin/env python3
"""
SIM113: Phase 2 — SSB Hilltop Quintessence V_J = λ(Ψ²−v²)²
CMSTG Phase 2: Curvature-Memory Scalar-Tensor Gravity (Unlocked Action)

SIM112 showed λΨ⁴ gives freezing quintessence (w_a > 0) — wrong sign for DESI.
DESI requires thawing quintessence (w_a < 0): field frozen at high z, rolling now.

Fix derived in Phase 2 Paper I: use SSB (spontaneous symmetry breaking) potential:
  V_J(Ψ) = λ(Ψ² − v²)²

Key properties:
  Hilltop at Ψ=0: V_J(0) = λv⁴ = ρ_DE,0  →  DE energy from hilltop height
  VEV at Ψ=±v:   V_J(±v) = 0              →  true vacuum has zero DE
  G_eff ≈ G_N:   Ψ₀ ≈ 0 today             →  no G_eff modification
  Tachyonic mass: V_E''(0) = −2λv²         →  slow roll AWAY from hilltop (thawing)

Modified equations (units 8πG=1, H₀=1, no bare Λ):

  Friedmann:
    H²[3(1+2Λ₀Ψ²) + 6Λ₀Ψy − ½y²] = ρ_m a⁻³ + ρ_r a⁻⁴ + λ(Ψ²−v²)²

  Klein-Gordon:
    d²Ψ/dx² + (3−ε_H)dΨ/dx + [4λΨ(Ψ²−v²)]/H² = 2Λ₀Ψ R/H²
    Note: dV_J/dΨ = 4λΨ(Ψ²−v²) = 4λΨ³ − 4λv²Ψ

DE scale constraint: λv⁴ = ρ_DE,0 = 0.685 × 3 = 2.055

Strategy:
  1. Fix v, compute λ = ρ_DE,0/v⁴
  2. Scan Ψ₀ ∈ [0, v/2] (field near hilltop today, not yet at VEV)
  3. Self-consistent initial conditions from Friedmann + slow-roll
  4. Integrate backward z=0→5, extract CPL (w₀, w_a)
  5. Find (v, Ψ₀) that satisfies DESI + Planck + G_eff bounds
"""

import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import curve_fit
import json, os

# ─── Fixed parameters ───────────────────��────────────────────���────────────────
Lambda0     = 0.003
Omega_m     = 0.315
Omega_r     = 9.1e-5
rho_DE_target = 0.685 * 3.0    # ρ_DE,0 = λv⁴ in sim units
rho_m0      = Omega_m * 3.0
rho_r0      = Omega_r * 3.0

# Observational targets
PLANCK_w0,  PLANCK_sw0 = -1.03, 0.03
DESI_w0,    DESI_wa    = -0.827, -0.75
DESI_sw0,   DESI_swa   = 0.060, 0.29

# ─── SSB potential ────────────────────────────���───────────────────────────────

def V_J(Psi, lam, v):
    """Jordan-frame SSB potential: λ(Ψ²−v²)²"""
    return lam * (Psi**2 - v**2)**2

def dV_J_dPsi(Psi, lam, v):
    """dV_J/dΨ = 4λΨ(Ψ²−v²)"""
    return 4.0*lam*Psi*(Psi**2 - v**2)

# ─── Background equations ────────────────────────────��────────────────────────

def get_H2(Psi, y, x, lam, v):
    a    = np.exp(x)
    rhs  = rho_m0*a**(-3) + rho_r0*a**(-4) + V_J(Psi, lam, v)
    coef = 3.0*(1+2*Lambda0*Psi**2) + 6*Lambda0*Psi*y - 0.5*y**2
    if coef <= 0.5:
        coef = 3.0
    return rhs / coef


def get_eps_H(Psi, y, x, lam, v):
    eps  = 5e-4
    H2p  = get_H2(Psi, y, x+eps, lam, v)
    H2m  = get_H2(Psi, y, x-eps, lam, v)
    H2   = get_H2(Psi, y, x,     lam, v)
    if H2 < 1e-40:
        return 0.0
    return -0.5*(H2p-H2m)/(2*eps*H2)


def odes(x, state, lam, v):
    Psi, y  = state
    H2      = get_H2(Psi, y, x, lam, v)
    eps_H   = get_eps_H(Psi, y, x, lam, v)
    R_over_H2 = 6*(2 - eps_H)

    # KG: source − restoring − Hubble friction
    # dV_J/dΨ / H² = 4λΨ(Ψ²−v²)/H²
    source  = 2*Lambda0*Psi*R_over_H2 - dV_J_dPsi(Psi, lam, v)/max(H2, 1e-40)
    dy_dx   = source - (3 - eps_H)*y

    return [y, dy_dx]


def cpl(z, w0, wa):
    return w0 + wa*z/(1+z)


# ─── Self-consistent initial conditions ─────────────────────────��─────────────

def find_ic(Psi0, lam, v, n_iter=30):
    """
    Iterate for (H₀=1-satisfying) initial velocity y₀.
    Friedmann: coef × H² = ρ_m + ρ_r + V_J(Ψ₀)
    With H=1: coef = ρ_m + ρ_r + V_J(Ψ₀)
    But coef depends on y₀. Iterate.
    """
    y = 0.0
    for _ in range(n_iter):
        coef    = 3*(1+2*Lambda0*Psi0**2) + 6*Lambda0*Psi0*y - 0.5*y**2
        rhs     = rho_m0 + rho_r0 + V_J(Psi0, lam, v)
        H2_want = rhs / max(coef, 0.5)

        # KG slow-roll at z=0 (ε_H ≈ 0, R/H²=12)
        R0_over_H2 = 12.0
        source = 2*Lambda0*Psi0*R0_over_H2 - dV_J_dPsi(Psi0, lam, v)/max(H2_want, 1e-40)
        y_new  = source / 3.0
        y_new  = np.clip(y_new, -10, 10)

        if abs(y_new - y) < 1e-9:
            break
        y = 0.6*y + 0.4*y_new

    # Check H(z=0) ≈ 1
    H2 = get_H2(Psi0, y, 0.0, lam, v)
    if abs(np.sqrt(max(H2, 0)) - 1.0) > 0.05:
        return None
    return y


# ─── Single-point run ─────────────────────────────────────────────────────────

def run_point(v, Psi0):
    lam  = rho_DE_target / v**4    # DE scale constraint: λv⁴ = ρ_DE,0
    y0   = find_ic(Psi0, lam, v)
    if y0 is None:
        return None

    G_eff_dev = abs(1/(1+2*Lambda0*Psi0**2) - 1.0)

    try:
        sol = solve_ivp(
            lambda x, s: odes(x, s, lam, v),
            [0.0, -np.log(6.0)],
            [Psi0, y0],
            method='DOP853',
            dense_output=True,
            rtol=1e-9, atol=1e-11,
            max_step=0.002,
        )
        if not sol.success:
            return None
    except Exception:
        return None

    z_grid = np.array([0.0, 0.1, 0.2, 0.3, 0.5, 0.7, 1.0, 1.5, 2.0, 3.0, 5.0])
    x_grid = -np.log(1+z_grid)
    a_grid = 1.0/(1+z_grid)

    H_z, H_lcdm_z, Psi_z = [], [], []
    for x in x_grid:
        s = sol.sol(x)
        P, yy = s[0], s[1]
        H2 = get_H2(P, yy, x, lam, v)
        H_z.append(np.sqrt(max(H2, 0.0)))
        a = np.exp(x)
        H_l = np.sqrt(Omega_m*a**(-3) + Omega_r*a**(-4) + 0.685*3.0)
        H_lcdm_z.append(H_l)
        Psi_z.append(P)

    H_z    = np.array(H_z)
    H_lcdm = np.array(H_lcdm_z)
    Psi_z  = np.array(Psi_z)

    rho_m_z = rho_m0*(1+z_grid)**3
    rho_r_z = rho_r0*(1+z_grid)**4
    rho_DE  = 3*H_z**2 - rho_m_z - rho_r_z

    if np.any(rho_DE <= 0):
        return None

    w_eff = np.full_like(z_grid, -1.0)
    for i in range(1, len(z_grid)-1):
        dlna   = np.log(a_grid[i+1]/a_grid[i-1])
        dlnrho = np.log(rho_DE[i+1]/rho_DE[i-1])
        w_eff[i] = -1.0 - dlnrho/(3.0*dlna)
    w_eff[0]  = w_eff[1]
    w_eff[-1] = w_eff[-2]

    mask = z_grid <= 2.0
    try:
        popt, _ = curve_fit(cpl, z_grid[mask], w_eff[mask], p0=[-1.0, -0.5])
        w0_fit, wa_fit = float(popt[0]), float(popt[1])
    except Exception:
        w0_fit, wa_fit = float(w_eff[0]), 0.0

    pull_planck  = (w0_fit - PLANCK_w0) / PLANCK_sw0
    pull_desi_w0 = (w0_fit - DESI_w0)   / DESI_sw0
    pull_desi_wa = (wa_fit - DESI_wa)    / DESI_swa
    desi_tension = float(np.sqrt(pull_desi_w0**2 + pull_desi_wa**2))
    delta_Psi    = abs(Psi_z[-1] - Psi_z[0]) / max(abs(Psi_z[0]), 1e-10) * 100.0

    return {
        'v'            : float(v),
        'lambda'       : float(lam),
        'Psi0'         : float(Psi0),
        'y0'           : float(y0),
        'G_eff_dev_pct': float(G_eff_dev*100),
        'w0'           : w0_fit,
        'wa'           : wa_fit,
        'pull_planck'  : float(pull_planck),
        'pull_desi_w0' : float(pull_desi_w0),
        'pull_desi_wa' : float(pull_desi_wa),
        'desi_tension' : desi_tension,
        'planck_pass'  : bool(abs(pull_planck) < 2.0),
        'desi_ok'      : bool(desi_tension < 2.0),
        'geff_ok'      : bool(G_eff_dev < 0.05),
        'delta_Psi_pct': float(delta_Psi),
        'Psi_z5'       : float(Psi_z[-1]),
        'w_eff_grid'   : [float(w) for w in w_eff],
        'z_grid'       : list(z_grid),
        'H_ratio_z0'   : float(H_z[0]/H_lcdm[0]),
    }


# ─── Scan ─────────────────────────────────────────────────────��───────────────
print("="*80)
print("SIM113: Phase 2 — SSB Hilltop Quintessence V_J = λ(Ψ²−v²)²")
print("="*80)
print(f"DE scale constraint: λv⁴ = ρ_DE,0 = {rho_DE_target:.4f}")
print(f"Targets: Planck w₀={PLANCK_w0}±{PLANCK_sw0}  |  "
      f"DESI w₀={DESI_w0}±{DESI_sw0}, wₐ={DESI_wa}±{DESI_swa}")
print()

# Scan v (VEV scale) and Ψ₀ (field value today)
V_SCAN   = [2.0, 3.0, 4.0, 5.0, 7.0, 10.0, 15.0, 20.0]

all_results = []
best_desi   = None

for v in V_SCAN:
    lam = rho_DE_target / v**4
    tachyon_mass = 2*lam*v**2    # |V_E''(0)|/H₀² should be << 1
    print(f"v={v:5.1f}  λ={lam:.5f}  2λv²={tachyon_mass:.4f}  "
          f"(slow-roll param at hilltop; want << 1)")

    # Scan Ψ₀ ∈ [0, v*0.6] in steps
    psi_scan = np.linspace(0.01, min(v*0.6, 5.0), 30)
    v_results = []
    for Psi0 in psi_scan:
        r = run_point(v, Psi0)
        if r is not None:
            v_results.append(r)
            all_results.append(r)

    if v_results:
        # Show best (Planck-passing or overall best)
        planck_ok = [r for r in v_results if r['planck_pass']]
        if planck_ok:
            br = min(planck_ok, key=lambda r: r['desi_tension'])
        else:
            br = min(v_results, key=lambda r: r['desi_tension'])

        sign_wa = "✓" if br['wa'] < 0 else "✗"
        desi_flag = "PASS" if br['desi_ok'] else f"{br['desi_tension']:.1f}σ"
        print(f"  Best: Ψ₀={br['Psi0']:.3f}  y₀={br['y0']:+.3f}  "
              f"w₀={br['w0']:+.4f}  wₐ={br['wa']:+.4f} {sign_wa}  "
              f"DESI:{desi_flag}  Planck:{abs(br['pull_planck']):.1f}σ  "
              f"G_eff:{br['G_eff_dev_pct']:.2f}%")

        if br['desi_ok'] and br['planck_pass'] and br['geff_ok']:
            if best_desi is None or br['desi_tension'] < best_desi['desi_tension']:
                best_desi = br
    print()

# ─── Fine scan around best ─────────────────────────────���──────────────────────
fine_results = []
if best_desi or all_results:
    ref = best_desi if best_desi else min(all_results, key=lambda r: r['desi_tension'])
    v_c   = ref['v']
    psi_c = ref['Psi0']
    v_fine    = np.linspace(max(1.0, v_c*0.5), v_c*2.0, 20)
    psi_fine  = np.linspace(max(0.01, psi_c*0.3), psi_c*3.0, 20)

    print(f"─── Fine scan around v={v_c:.1f}, Ψ₀={psi_c:.3f} ───")
    for v in v_fine:
        for Psi0 in psi_fine:
            r = run_point(v, Psi0)
            if r is not None:
                fine_results.append(r)

    all_fine_planck = [r for r in fine_results if r['planck_pass'] and r['geff_ok']]
    if all_fine_planck:
        best_fine = min(all_fine_planck, key=lambda r: r['desi_tension'])
        if best_desi is None or best_fine['desi_tension'] < best_desi['desi_tension']:
            best_desi = best_fine
        print(f"  Best fine: v={best_fine['v']:.3f}  Ψ₀={best_fine['Psi0']:.4f}  "
              f"w₀={best_fine['w0']:+.4f}  wₐ={best_fine['wa']:+.4f}  "
              f"DESI {best_fine['desi_tension']:.2f}σ  Planck {abs(best_fine['pull_planck']):.2f}σ")

# ─── Verdict ──────────────────────────────────────────────────────────────────
print()
print("="*80)
print("VERDICT")
print("="*80)

viable = [r for r in all_results + fine_results
          if r['desi_ok'] and r['planck_pass'] and r['geff_ok']]

if viable:
    best = min(viable, key=lambda r: r['desi_tension'])
    verdict = "PASS — SSB hilltop quintessence reproduces DE without bare Λ"
elif best_desi:
    best = best_desi
    verdict = (f"PARTIAL — best: DESI {best['desi_tension']:.2f}σ, "
               f"Planck {abs(best['pull_planck']):.2f}σ")
else:
    pg = [r for r in all_results if r['planck_pass'] and r['geff_ok']]
    best = min(pg, key=lambda r: r['desi_tension']) if pg else (
           min(all_results, key=lambda r: r['desi_tension']) if all_results else None)
    verdict = "FAIL — no viable solution found"

print(f"\n  RESULT: {verdict}")

if best:
    print(f"\n  Optimal parameters:")
    print(f"    v          = {best['v']:.4f} M_Pl  (VEV scale)")
    print(f"    λ          = {best['lambda']:.6f}  (fixed by λv⁴=ρ_DE,0)")
    print(f"    Ψ₀         = {best['Psi0']:.4f} M_Pl  (field today, near hilltop)")
    print(f"    y₀=dΨ/dx   = {best['y0']:+.5f}  (field velocity today)")
    print(f"    2λv²       = {2*best['lambda']*best['v']**2:.4f}  (tachyon mass / H₀²)")
    print(f"    G_eff/G_N  = {1/(1+2*Lambda0*best['Psi0']**2):.5f}  ({best['G_eff_dev_pct']:.3f}% dev)")
    print(f"    w₀ (CMSTG)  = {best['w0']:+.4f}  (DESI: {DESI_w0}±{DESI_sw0})")
    print(f"    wₐ (CMSTG)  = {best['wa']:+.4f}  (DESI: {DESI_wa}±{DESI_swa})")
    print(f"    wₐ sign    = {'CORRECT (negative = thawing)' if best['wa'] < 0 else 'WRONG (positive = freezing)'}")
    print(f"    DESI σ     = {best['desi_tension']:.2f}")
    print(f"    Planck σ   = {abs(best['pull_planck']):.2f}")
    print(f"    ΔΨ(z=0→5)  = {best['delta_Psi_pct']:.2f}%  (Ψ(z=5)={best['Psi_z5']:.4f} M_Pl)")
    print()

    # w_eff table
    print("  w_eff(z) profile:")
    print(f"  {'z':>5}  {'w_eff':>8}")
    for z, w in zip(best['z_grid'], best['w_eff_grid']):
        print(f"  {z:>5.1f}  {w:>8.4f}")

# ─── Save ──────────────────────���───────────────────────────────────��──────────
out_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'Outputs')
os.makedirs(out_dir, exist_ok=True)

diag = {
    'sim'         : 'SIM113',
    'description' : 'Phase 2 SSB hilltop quintessence V_J=lambda(Psi^2-v^2)^2, no bare Lambda',
    'Lambda0'     : Lambda0,
    'rho_DE_target': rho_DE_target,
    'viable'      : viable,
    'best'        : best,
    'all_results' : all_results,
    'fine_results': fine_results,
    'verdict'     : verdict,
}
with open(os.path.join(out_dir, 'sim113_diagnostics.json'), 'w') as f:
    json.dump(diag, f, indent=2)

print(f"\nDiagnostics saved → Outputs/sim113_diagnostics.json")
