#!/usr/bin/env python3
"""
SIM112: Phase 2 Background — Quintessence Dark Energy from λΨ⁴
RIFT Phase 2: Recursive Intelligence-Field Theory (Unlocked Action)

Phase 2 action (Jordan frame, no bare Λ):
  S = ∫d⁴x√(-g) [(M_Pl²+2Λ₀Ψ²)/2 R − ½(∇Ψ)² − λΨ⁴] + S_SM

Phase 1 result (SIM111): locked action cannot resolve DESI tension.
Structural obstruction: ρ_Λ dominates DE at all Planck-consistent m₀.
Resolution: remove bare Λ, derive DE from Einstein-frame plateau of λΨ⁴.

Modified Friedmann (units 8πG=1, H₀=1, no bare Λ):
  H²[3(1+2Λ₀Ψ²) + 6Λ₀Ψy − ½y²] = ρ_m a⁻³ + ρ_r a⁻⁴ + λΨ⁴

Modified Klein-Gordon (Phase 2, no mass term):
  d²Ψ/dx² + (3−ε_H)dΨ/dx + 4λΨ³/H² = 2Λ₀Ψ R/H²

Strategy:
  Scan Ψ₀ ∈ [0.5, 10]. For each Ψ₀, self-consistently solve for (λ, y₀)
  from the Friedmann constraint (H₀=1) and the slow-roll KG condition.
  Integrate backwards z=0→5, extract CPL (w₀, wₐ), report G_eff and tensions.

Key observational bounds:
  G_eff/G_N = 1/(1+2Λ₀Ψ₀²): must be within 5% (Planck/BBN constraints)
  → Ψ₀ ≤ √(0.05/(2Λ₀)) ≈ 2.9 M_Pl  [5% bound]
  → Ψ₀ ≤ √(0.10/(2Λ₀)) ≈ 4.1 M_Pl  [10% bound, conservative]
"""

import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import curve_fit, brentq
import json, os

# ─── Fixed parameters ─────────────────────────────────────────────────────────
Lambda0  = 0.003
Omega_m  = 0.315
Omega_r  = 9.1e-5
# NO Omega_Lambda — dark energy comes from Psi field

rho_m0   = Omega_m * 3.0
rho_r0   = Omega_r * 3.0

# Observational targets
PLANCK_w0,  PLANCK_sw0 = -1.03, 0.03
DESI_w0,    DESI_wa    = -0.827, -0.75
DESI_sw0,   DESI_swa   = 0.060, 0.29

# G_eff tolerance
GEFF_MAX_DEVIATION = 0.10   # 10% tolerance on G_eff/G_N at z=0

# ─── Phase 2 Friedmann solver ─────────────────────────────────────────────────

def get_H2_p2(Psi, y, x, lam):
    """Phase 2 Friedmann: no bare Λ, potential = λΨ⁴."""
    a    = np.exp(x)
    rhs  = rho_m0*a**(-3) + rho_r0*a**(-4) + lam*Psi**4
    coef = 3.0*(1.0 + 2.0*Lambda0*Psi**2) + 6.0*Lambda0*Psi*y - 0.5*y**2
    if coef <= 0.5:
        coef = 3.0
    return rhs / coef


def get_epsilon_H_p2(Psi, y, x, lam):
    eps  = 5e-4
    H2p  = get_H2_p2(Psi, y, x + eps, lam)
    H2m  = get_H2_p2(Psi, y, x - eps, lam)
    H2   = get_H2_p2(Psi, y, x,       lam)
    if H2 < 1e-40:
        return 0.0
    return -0.5*(H2p - H2m)/(2.0*eps*H2)


def rift_odes_p2(x, state, lam):
    """Phase 2 ODE system: Klein-Gordon with 4λΨ³ restoring force."""
    Psi, y = state
    H2     = get_H2_p2(Psi, y, x, lam)
    eps_H  = get_epsilon_H_p2(Psi, y, x, lam)

    R_over_H2 = 6.0*(2.0 - eps_H)
    # Phase 2 KG: no m₀² term; quartic restoring force
    source = 2.0*Lambda0*Psi*R_over_H2 - 4.0*lam*Psi**3/max(H2, 1e-40)
    dy_dx  = source - (3.0 - eps_H)*y

    return [y, dy_dx]


# ─── Self-consistent (λ, y₀) at z=0 ─────────────────────────────────────────

def find_ic(Psi0, n_iter=20):
    """
    Iterate to find self-consistent (lambda, y0) satisfying:
      (1) Friedmann: H(z=0) = 1
      (2) Slow-roll KG at z=0 (ε_H≈0, R/H²=12)

    Returns (lambda, y0) or None if no valid solution.
    """
    R0_over_H2 = 12.0   # de Sitter-like at z=0
    y = 0.0             # initial guess

    for _ in range(n_iter):
        # (1) λ from Friedmann constraint at z=0 with H=1
        coef = 3.0*(1.0 + 2.0*Lambda0*Psi0**2) + 6.0*Lambda0*Psi0*y - 0.5*y**2
        lam  = (coef - rho_m0 - rho_r0) / Psi0**4

        if lam <= 0:
            return None   # negative lambda unphysical

        # (2) slow-roll y from KG: (3−ε_H)y ≈ 2Λ₀Ψ(R/H²) − 4λΨ³/H²
        H2_check = get_H2_p2(Psi0, y, 0.0, lam)
        eps_H    = 0.0   # de Sitter approximation at z=0
        y_new    = (2.0*Lambda0*Psi0*R0_over_H2 - 4.0*lam*Psi0**3/max(H2_check, 1e-40)) / (3.0 - eps_H)
        y_new    = np.clip(y_new, -10.0, 10.0)

        if abs(y_new - y) < 1e-8:
            break
        y = 0.5*y + 0.5*y_new   # damped update for convergence

    # Verify H(z=0) ≈ 1
    H2_final = get_H2_p2(Psi0, y, 0.0, lam)
    if abs(np.sqrt(max(H2_final, 0)) - 1.0) > 0.05:
        return None

    return lam, y


# ─── Background integrator ────────────────────────────────────────────────────

def cpl(z, w0, wa):
    return w0 + wa*z/(1.0 + z)


def run_psi0(Psi0, verbose=False):
    """Run Phase 2 background for given Ψ₀. Returns result dict or None."""
    ic = find_ic(Psi0)
    if ic is None:
        return None
    lam, y0 = ic

    # G_eff at z=0
    G_eff_ratio = 1.0 / (1.0 + 2.0*Lambda0*Psi0**2)
    G_eff_dev   = abs(G_eff_ratio - 1.0)

    try:
        sol = solve_ivp(
            lambda x, s: rift_odes_p2(x, s, lam),
            [0.0, -np.log(6.0)],   # z=0 to z=5
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

    # Evaluate on z-grid
    z_grid = np.array([0.0, 0.1, 0.2, 0.3, 0.5, 0.7, 1.0, 1.5, 2.0, 3.0, 5.0])
    x_grid = -np.log(1.0 + z_grid)
    a_grid = 1.0 / (1.0 + z_grid)

    H_z, H_lcdm_z, Psi_z, y_z = [], [], [], []
    for x in x_grid:
        s      = sol.sol(x)
        P, yy  = s[0], s[1]
        H2     = get_H2_p2(P, yy, x, lam)
        H_z.append(np.sqrt(max(H2, 0.0)))
        a      = np.exp(x)
        H_l    = np.sqrt(Omega_m*a**(-3) + Omega_r*a**(-4) + 0.685*3.0)  # LCDM reference
        H_lcdm_z.append(H_l)
        Psi_z.append(P)
        y_z.append(yy)

    H_z    = np.array(H_z)
    H_lcdm = np.array(H_lcdm_z)
    Psi_z  = np.array(Psi_z)
    y_z    = np.array(y_z)

    # Effective DE density
    rho_m_z = rho_m0*(1.0 + z_grid)**3
    rho_r_z = rho_r0*(1.0 + z_grid)**4
    rho_DE  = 3.0*H_z**2 - rho_m_z - rho_r_z

    if np.any(rho_DE <= 0):
        return None

    # w_eff(z)
    w_eff = np.full_like(z_grid, -1.0)
    for i in range(1, len(z_grid)-1):
        dlna   = np.log(a_grid[i+1]/a_grid[i-1])
        dlnrho = np.log(rho_DE[i+1]/rho_DE[i-1])
        w_eff[i] = -1.0 - dlnrho/(3.0*dlna)
    w_eff[0]  = w_eff[1]
    w_eff[-1] = w_eff[-2]

    # CPL fit
    mask = z_grid <= 2.0
    try:
        popt, _ = curve_fit(cpl, z_grid[mask], w_eff[mask], p0=[-1.0, 0.0])
        w0_fit, wa_fit = float(popt[0]), float(popt[1])
    except Exception:
        w0_fit, wa_fit = float(w_eff[0]), 0.0

    # Tensions
    pull_planck  = (w0_fit - PLANCK_w0) / PLANCK_sw0
    pull_desi_w0 = (w0_fit - DESI_w0)   / DESI_sw0
    pull_desi_wa = (wa_fit - DESI_wa)    / DESI_swa
    desi_tension = float(np.sqrt(pull_desi_w0**2 + pull_desi_wa**2))

    # Einstein-frame plateau height (for reference)
    V_plateau = lam / (4.0*Lambda0**2)   # λ/(4Λ₀²) in M_Pl=1 units
    delta_Psi = abs(Psi_z[-1] - Psi_z[0]) / abs(Psi_z[0]) * 100.0

    return {
        'Psi0'         : float(Psi0),
        'lambda'       : float(lam),
        'y0'           : float(y0),
        'G_eff_ratio'  : float(G_eff_ratio),
        'G_eff_dev_pct': float(G_eff_dev*100.0),
        'w0'           : w0_fit,
        'wa'           : wa_fit,
        'pull_planck'  : float(pull_planck),
        'pull_desi_w0' : float(pull_desi_w0),
        'pull_desi_wa' : float(pull_desi_wa),
        'desi_tension' : desi_tension,
        'planck_pass'  : bool(abs(pull_planck) < 2.0),
        'desi_ok'      : bool(desi_tension < 2.0),
        'geff_ok'      : bool(G_eff_dev < GEFF_MAX_DEVIATION),
        'delta_Psi_pct': float(delta_Psi),
        'dPsi_dx_z0'   : float(y_z[0]),
        'V_plateau'    : float(V_plateau),
        'Psi_plateau'  : float(1.0/np.sqrt(2.0*Lambda0)),
        'H_ratio_z0'   : float(H_z[0]/H_lcdm[0]),
        'w_eff_grid'   : [float(w) for w in w_eff],
        'z_grid'       : list(z_grid),
        'Psi_z'        : [float(p) for p in Psi_z],
    }


# ─── Scan ─────────────────────────────────────────────────────────────────────
PSI0_SCAN = np.concatenate([
    np.linspace(0.5, 3.0, 26),    # fine grid in the G_eff-safe zone
    np.linspace(3.5, 10.0, 14),   # coarser for large Psi0
])

print("="*78)
print("SIM112: Phase 2 Background — Quintessence DE from λΨ⁴ (No Bare Λ)")
print("="*78)
print(f"Scan {len(PSI0_SCAN)} Ψ₀ values ∈ [{PSI0_SCAN[0]:.1f}, {PSI0_SCAN[-1]:.1f}] M_Pl")
print(f"Plateau scale: Ψ* = 1/√(2Λ₀) = {1.0/np.sqrt(2*Lambda0):.2f} M_Pl")
print(f"G_eff 10% bound: Ψ₀ ≤ {np.sqrt(0.10/(2*Lambda0)):.2f} M_Pl")
print()
print(f"Targets:")
print(f"  Planck 2018:  w₀ = {PLANCK_w0:.3f} ± {PLANCK_sw0:.3f}")
print(f"  DESI Y1 2024: w₀ = {DESI_w0:.3f} ± {DESI_sw0:.3f},  wₐ = {DESI_wa:.2f} ± {DESI_swa:.2f}")
print()

results = []
for Psi0 in PSI0_SCAN:
    r = run_psi0(Psi0)
    if r is not None:
        results.append(r)

# ─── Summary table ────────────────────────────────────────────────────────────
print("─── Scan Results ───")
print(f"{'Ψ₀':>7} {'λ':>10} {'y₀':>7} {'G_eff%':>7} {'w₀':>8} {'wₐ':>7} "
      f"{'DESI σ':>8} {'Plan σ':>7} {'G_ok':>5} {'Viable':>7}")
print("-"*80)

viable = []
for r in results:
    g_flag = "OK"  if r['geff_ok']     else "WARN"
    p_flag = "PASS" if r['planck_pass'] else "FAIL"
    d_flag = "PASS" if r['desi_ok']    else f"{r['desi_tension']:.1f}σ"
    all_ok = r['geff_ok'] and r['planck_pass'] and r['desi_ok']
    v_flag = "***" if all_ok else ("P+G" if (r['geff_ok'] and r['planck_pass']) else "")
    if all_ok:
        viable.append(r)

    print(f"{r['Psi0']:>7.3f} {r['lambda']:>10.5f} {r['y0']:>7.3f} "
          f"{r['G_eff_dev_pct']:>7.2f} {r['w0']:>8.4f} {r['wa']:>7.4f} "
          f"{r['desi_tension']:>8.2f} {abs(r['pull_planck']):>7.2f} "
          f"{g_flag:>5} {v_flag:>7}")

print()

# ─── Verdict ──────────────────────────────────────────────────────────────────
print("="*78)
print("VERDICT")
print("="*78)

# Best Planck+G_eff consistent result (minimum DESI tension)
pg_ok = [r for r in results if r['geff_ok'] and r['planck_pass']]

if viable:
    print(f"\n  FULLY VIABLE SOLUTIONS FOUND: {len(viable)}")
    for r in sorted(viable, key=lambda r: r['desi_tension']):
        print(f"  Ψ₀={r['Psi0']:.3f}  λ={r['lambda']:.5f}  "
              f"w₀={r['w0']:+.4f}  wₐ={r['wa']:+.4f}  "
              f"DESI {r['desi_tension']:.2f}σ  G_eff {r['G_eff_dev_pct']:.2f}%")
    best = min(viable, key=lambda r: r['desi_tension'])
    verdict = "PASS — Phase 2 λΨ⁴ action can explain DE without bare Λ"
elif pg_ok:
    best = min(pg_ok, key=lambda r: r['desi_tension'])
    print(f"\n  No fully viable (DESI) solution. Best Planck+G_eff consistent:")
    print(f"  Ψ₀={best['Psi0']:.3f}  λ={best['lambda']:.5f}  "
          f"w₀={best['w0']:+.4f}  wₐ={best['wa']:+.4f}  "
          f"DESI {best['desi_tension']:.2f}σ")
    verdict = f"PARTIAL — Planck+G_eff consistent but DESI tension = {best['desi_tension']:.2f}σ"
else:
    best = min(results, key=lambda r: r['desi_tension']) if results else None
    verdict = "FAIL — No G_eff-consistent, Planck-consistent solution found"

print(f"\n  RESULT: {verdict}")

if best:
    print()
    print(f"  Optimal point:")
    print(f"    Ψ₀          = {best['Psi0']:.4f} M_Pl  (plateau Ψ* = {best['Psi_plateau']:.2f} M_Pl)")
    print(f"    λ           = {best['lambda']:.6f}")
    print(f"    y₀=dΨ/dx    = {best['y0']:+.4f}")
    print(f"    G_eff/G_N   = {best['G_eff_ratio']:.4f}  ({best['G_eff_dev_pct']:.2f}% deviation)")
    print(f"    V_plateau   = {best['V_plateau']:.4f}  (need ρ_DE,0 ≈ 2.055)")
    print(f"    w₀ (RIFT)   = {best['w0']:+.4f}")
    print(f"    wₐ (RIFT)   = {best['wa']:+.4f}")
    print(f"    DESI σ      = {best['desi_tension']:.2f}")
    print(f"    Planck σ    = {abs(best['pull_planck']):.2f}")
    print(f"    ΔΨ (z=0→5)  = {best['delta_Psi_pct']:.2f}%")
    print()
    if best['V_plateau'] < 1.0:
        print(f"  NOTE: Ψ₀ < Ψ*: field is on the SLOPE of the Einstein-frame potential.")
        print(f"  This is the thawing quintessence regime.")
    else:
        print(f"  NOTE: Ψ₀ > Ψ*: field is near or on the PLATEAU.")

# ─── Einstein-frame potential diagnostic ─────────────────────────────────────
if results:
    print()
    print("─── Einstein-frame potential V_E(Ψ) = λΨ⁴/(1+2Λ₀Ψ²)² diagnostic ───")
    Psi_range = np.linspace(0.1, 15.0, 200)
    r_ref     = best if best else results[0]
    lam_ref   = r_ref['lambda']
    V_E       = lam_ref * Psi_range**4 / (1.0 + 2.0*Lambda0*Psi_range**2)**2
    V_plateau = lam_ref / (4.0*Lambda0**2)
    i_half    = np.argmin(abs(V_E - 0.5*V_plateau))
    print(f"  λ = {lam_ref:.6f},  plateau height = {V_plateau:.4f}")
    print(f"  Half-plateau at Ψ ≈ {Psi_range[i_half]:.2f} M_Pl")
    print(f"  Ψ₀ = {r_ref['Psi0']:.3f} M_Pl  →  V_E(Ψ₀)/V_plateau = "
          f"{lam_ref*r_ref['Psi0']**4/(1+2*Lambda0*r_ref['Psi0']**2)**2/V_plateau:.4f}")

# ─── Save ─────────────────────────────────────────────────────────────────────
out_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'Outputs')
os.makedirs(out_dir, exist_ok=True)

diag = {
    'sim'         : 'SIM112',
    'description' : 'Phase 2 background: quintessence DE from lambda*Psi^4, no bare Lambda',
    'Lambda0'     : Lambda0,
    'results'     : results,
    'viable'      : viable,
    'best'        : best,
    'verdict'     : verdict,
    'desi_target' : {'w0': DESI_w0, 'wa': DESI_wa, 'sw0': DESI_sw0, 'swa': DESI_swa},
    'planck_target': {'w0': PLANCK_w0, 'sw0': PLANCK_sw0},
}
out_path = os.path.join(out_dir, 'sim112_diagnostics.json')
with open(out_path, 'w') as f:
    json.dump(diag, f, indent=2)

print(f"\nDiagnostics saved → Outputs/sim112_diagnostics.json")
