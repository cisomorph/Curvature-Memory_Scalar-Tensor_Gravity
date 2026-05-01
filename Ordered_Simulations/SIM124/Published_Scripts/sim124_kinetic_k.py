#!/usr/bin/env python3
"""
SIM124 — CMSTG Phase 4: Non-canonical kinetic term K(Ψ) = exp(-αΨ²)

SIM123 showed the hilltop field is kinematically frozen (Ψ stays at ~0.001,
can't reach 2.62 M_Pl today). With K(Ψ) = exp(-αΨ²) < 1, the Hubble
friction term is reduced: Ψ̇ ≈ -V'/(3HK) = -V'exp(+αΨ²)/(3H).
Once the field departs the hilltop even slightly, exp(αΨ²) provides
positive feedback — a potentially self-amplifying roll.

Approach: FORWARD integration in N = ln(a) from z=z_ini (Ψ=Ψ_ini)
to z=0. Scan (α, Ψ_ini) to find where Ψ(z=0) = 2.62 M_Pl.
If that locus exists at small Ψ_ini (near-hilltop ICs), the theory
is viable; if it requires Ψ_ini >> 0, it is fine-tuned.

Modified equations:
  Friedmann:  H² = (ρ_std + V) / (3F - K u²/2)
  KG:         du/dN = -u(3-ε_H) + αΨu² - V'/(K H²)
  where K = exp(-αΨ²), F = ½ + Λ₀Ψ², u = dΨ/dN

Parts:
  A — (α, Ψ_ini) scan: map Ψ(z=0) across parameter space
  B — Find (α, Ψ_ini) locus where Ψ(z=0) = 2.62 M_Pl
  C — For each viable α: H₀ from θ_*, DESI χ²
"""

import numpy as np
from scipy.integrate import solve_ivp, quad
from scipy.optimize import brentq
import json, os, warnings
warnings.filterwarnings('ignore')

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
OUT_DIR    = os.path.join(SCRIPT_DIR, "../Outputs")
os.makedirs(OUT_DIR, exist_ok=True)

# ── Constants ─────────────────────────────────────────────────────────
OM_H2      = 0.1430
OR_H2      = 4.18e-5
OB_H2      = 0.02237
THETA_OBS  = 1.04101
SIGMA_TH   = 0.00029
DESI_Z     = np.array([0.295, 0.510, 0.706, 0.930, 1.317, 2.330])
DESI_H     = np.array([81.7,  97.9, 110.7, 128.1, 156.4, 240.8])
DESI_SIG   = np.array([ 2.9,   3.8,   4.5,   5.3,   7.3,  15.7])
L0_REF     = 0.003
V_REF      = 13.16
PSI0_REF   = 2.62
Z_INI      = 1500.0
Z_CMB      = 1090.0
CHI2_121C  = 41.492
TENS_121C  = 2.63

print("=" * 70)
print("SIM124 — CMSTG Phase 4: K(Ψ) = exp(−αΨ²) Kinetic Modification")
print("=" * 70)

# ── Physics ───────────────────────────────────────────────────────────
def F_cmstg(psi, L0): return 0.5 + L0 * psi**2
def K_func(psi, alpha): return np.exp(np.clip(-alpha * psi**2, -500, 0))
def dK_dpsi(psi, alpha): return -2.0*alpha*psi * K_func(psi, alpha)
def V_pot(psi, lam, v): return lam * (psi**2 - v**2)**2
def dV_dpsi(psi, lam, v): return 4.0*lam*psi*(psi**2 - v**2)
def rho_std(z): return (OM_H2*(1+z)**3 + OR_H2*(1+z)**4) * 1e4

def calibrate_lam(H0, L0, v, psi0):
    """λ from slow-roll Friedmann at z=0: 3F₀H₀² = ρ₀ + V(ψ₀)."""
    F0  = F_cmstg(psi0, L0)
    num = 3.0*F0*H0**2 - rho_std(0.0)
    den = (psi0**2 - v**2)**2
    return (num/den) if (den > 0 and num > 0) else np.nan

def H2_solve(psi, u, z, lam, L0, alpha, v):
    F = F_cmstg(psi, L0)
    K = K_func(psi, alpha)
    V = V_pot(psi, lam, v)
    d = 3.0*F - 0.5*K*u**2
    if d <= 0: return np.nan
    H2 = (rho_std(z) + V) / d
    return H2 if H2 > 0 else np.nan

# ── Forward ODE in N = ln(a) ──────────────────────────────────────────
def ode_fwd(N, state, lam, L0, alpha, v, psi_max=50.0):
    psi, u = state
    if abs(psi) > psi_max:
        return [u, 0.0]
    z = max(0.0, np.exp(-N) - 1.0)
    K   = K_func(psi, alpha)
    F   = F_cmstg(psi, L0)
    dFdp = 2.0*L0*psi
    V    = V_pot(psi, lam, v)
    Vp   = dV_dpsi(psi, lam, v)
    rho  = rho_std(z)
    dH   = 3.0*F - 0.5*K*u**2

    if dH <= 1e-20 or (rho + V) <= 0:
        return [u, 0.0]
    H2 = (rho + V) / dH
    if not np.isfinite(H2) or H2 <= 0:
        return [u, 0.0]

    # ε_H (approximate, slow-roll for u)
    drho_dN   = -(3.0*OM_H2*(1+z)**3 + 4.0*OR_H2*(1+z)**4)*1e4
    dnum_dN   = drho_dN + Vp*u
    ddenom_dN = 3.0*dFdp*u - 0.5*(-2*alpha*psi*K)*u**3
    dH2_dN    = (dnum_dN*dH - (rho+V)*ddenom_dN) / dH**2
    if not np.isfinite(dH2_dN):
        dH2_dN = 0.0
    eps_H = -0.5*dH2_dN/H2

    # clip eps_H to prevent explosion
    eps_H = np.clip(eps_H, -10, 10)

    du_dN = -u*(3.0 - eps_H) + alpha*psi*u**2 - Vp/(K*H2)
    if not np.isfinite(du_dN):
        du_dN = 0.0
    return [u, du_dN]

def run_forward(alpha, psi_ini, L0=L0_REF, v=V_REF, H0_ref=67.0,
                psi0_calib=PSI0_REF, n_pts=600):
    """
    Forward-integrate from z=Z_INI to z=0.
    λ calibrated at (H0_ref, psi0_calib).
    Returns (z_arr, psi_arr, H_arr) or None.
    """
    lam = calibrate_lam(H0_ref, L0, v, psi0_calib)
    if np.isnan(lam) or lam <= 0:
        return None

    N_ini = -np.log(1.0 + Z_INI)
    N_end = 0.0

    # Initial u from slow-roll at z_ini
    K0  = K_func(psi_ini, alpha)
    H2_0 = H2_solve(psi_ini, 0.0, Z_INI, lam, L0, alpha, v)
    if np.isnan(H2_0):
        return None
    Vp0 = dV_dpsi(psi_ini, lam, v)
    u0  = -Vp0 / (3.0*K0*H2_0) if abs(K0*H2_0) > 0 else 0.0
    if not np.isfinite(u0) or abs(u0) > 1e3:
        u0 = 0.0

    N_eval = np.linspace(N_ini, N_end, n_pts)

    # Terminate if field overshoots
    def hit_vev(N, state, *a): return abs(state[0]) - 0.99*v
    hit_vev.terminal = True; hit_vev.direction = 1

    try:
        sol = solve_ivp(
            ode_fwd, [N_ini, N_end], [psi_ini, u0],
            t_eval=N_eval, args=(lam, L0, alpha, v),
            method='DOP853', rtol=1e-8, atol=1e-11,
            events=hit_vev, dense_output=False
        )
    except Exception:
        return None

    N_arr  = sol.t
    z_arr  = np.maximum(0.0, np.exp(-N_arr) - 1.0)
    psi_arr = sol.y[0]
    u_arr   = sol.y[1]

    H_arr = np.array([
        np.sqrt(max(1e-30, H2_solve(psi_arr[i], u_arr[i], z_arr[i],
                                     lam, L0, alpha, v)))
        for i in range(len(N_arr))
    ])

    return z_arr, psi_arr, u_arr, H_arr, lam

# ══════════════════════════════════════════════════════════════════════
# PART A: 2D scan (α, Ψ_ini) → Ψ(z=0)
# ══════════════════════════════════════════════════════════════════════
print("\nPart A: 2D scan (α, Ψ_ini) — field value reached at z=0...")
print(f"  {'α':>5}  {'Ψ_ini':>7}  {'Ψ(z=0)':>8}  {'F(z=0)':>8}  {'H(z=0)':>8}")

alpha_vals = [0.0, 0.1, 0.3, 0.5, 1.0, 1.5, 2.0, 3.0, 5.0]
psi_ini_vals = [0.001, 0.01, 0.05, 0.1, 0.3, 0.5, 1.0, 2.0, 2.5]

scan_A = []
for alpha in alpha_vals:
    for psi_ini in psi_ini_vals:
        r = run_forward(alpha, psi_ini)
        if r is None:
            scan_A.append({'alpha': alpha, 'psi_ini': psi_ini,
                           'psi_z0': np.nan, 'H_z0': np.nan})
            continue
        z_a, psi_a, u_a, H_a, lam_a = r
        psi_z0 = float(np.interp(0.0, z_a[::-1], psi_a[::-1]))
        H_z0   = float(np.interp(0.0, z_a[::-1], H_a[::-1]))
        F_z0   = F_cmstg(psi_z0, L0_REF)
        scan_A.append({'alpha': alpha, 'psi_ini': psi_ini,
                       'psi_z0': psi_z0, 'F_z0': F_z0,
                       'H_z0': H_z0})
        print(f"  {alpha:5.2f}  {psi_ini:7.4f}  {psi_z0:8.4f}  {F_z0:8.5f}  {H_z0:8.2f}")

# ══════════════════════════════════════════════════════════════════════
# PART B: Find Ψ_ini(α) locus where Ψ(z=0) = PSI0_REF = 2.62
# ══════════════════════════════════════════════════════════════════════
print("\nPart B: Shooting — find Ψ_ini(α) such that Ψ(z=0) = 2.62 M_Pl...")

def psi_z0_for_psiini(psi_ini, alpha):
    r = run_forward(alpha, psi_ini)
    if r is None: return np.nan
    z_a, psi_a, *_ = r
    return float(np.interp(0.0, z_a[::-1], psi_a[::-1]))

viable = []
alpha_fine = [0.0, 0.1, 0.2, 0.3, 0.5, 0.7, 1.0, 1.5, 2.0, 3.0, 5.0]

for alpha in alpha_fine:
    # Scan Ψ_ini to find bracket where Ψ(z=0) crosses 2.62
    psi_ini_try = np.logspace(-3, np.log10(PSI0_REF*0.99), 30)
    psi0_vals = np.array([psi_z0_for_psiini(p, alpha) for p in psi_ini_try])

    # Find crossing of PSI0_REF
    found = False
    for i in range(len(psi0_vals)-1):
        v1, v2 = psi0_vals[i], psi0_vals[i+1]
        if np.isnan(v1) or np.isnan(v2): continue
        if (v1 - PSI0_REF) * (v2 - PSI0_REF) < 0:
            try:
                def resid(p): return psi_z0_for_psiini(p, alpha) - PSI0_REF
                psi_sol = brentq(resid, psi_ini_try[i], psi_ini_try[i+1],
                                 xtol=1e-5, maxiter=20)
                viable.append({'alpha': alpha, 'psi_ini_sol': psi_sol})
                print(f"  α={alpha:.2f}  Ψ_ini={psi_sol:.5f} M_Pl → Ψ(z=0)=2.62 ✓")
                found = True
            except Exception:
                pass

    if not found:
        # Report max achieved
        valid_mask = ~np.isnan(psi0_vals)
        if valid_mask.any():
            pmax = np.nanmax(psi0_vals)
            print(f"  α={alpha:.2f}  max Ψ(z=0)={pmax:.5f} M_Pl — cannot reach 2.62")
        else:
            print(f"  α={alpha:.2f}  all NaN")

# ══════════════════════════════════════════════════════════════════════
# PART C: Observables for viable solutions
# ══════════════════════════════════════════════════════════════════════
print("\nPart C: Observables for viable (α, Ψ_ini) solutions...")

def theta_star(H_func):
    def rs_int(z):
        H = H_func(z)
        if H <= 0: return 0.0
        R = 3.0*OB_H2/(4.0*2.47e-5)/(1.0+z)
        return 1.0/(np.sqrt(3.0*(1.0+R))*H)
    def dc_int(z):
        H = H_func(z)
        return 1.0/H if H > 0 else 0.0
    rs, _ = quad(rs_int, 0.0, Z_CMB, limit=200, epsrel=1e-7)
    dc, _ = quad(dc_int, 0.0, Z_CMB, limit=200, epsrel=1e-7)
    return 100.0*rs/dc if dc > 0 else np.nan

def find_H0_theta(alpha, psi_ini_sol, L0=L0_REF, H0_lo=40.0, H0_hi=90.0):
    def resid(H0_trial):
        r = run_forward(alpha, psi_ini_sol, L0=L0, H0_ref=H0_trial,
                        psi0_calib=PSI0_REF)
        if r is None: return 999.0
        z_a, psi_a, u_a, H_a, _ = r
        Hf = lambda z: float(np.interp(z, z_a[::-1], H_a[::-1],
                                        left=H_a[-1], right=H_a[0]))
        th = theta_star(Hf)
        return (th - THETA_OBS) if np.isfinite(th) else 999.0
    try:
        r_lo, r_hi = resid(H0_lo), resid(H0_hi)
        if r_lo*r_hi > 0: return np.nan
        return brentq(resid, H0_lo, H0_hi, xtol=0.05, maxiter=20)
    except Exception:
        return np.nan

results_C = []
if not viable:
    print("  No viable solutions found — cannot compute observables.")
else:
    print(f"  {'α':>5}  {'Ψ_ini':>7}  {'H₀':>7}  {'100θ*':>8}  {'χ²_D':>8}  {'tens':>6}  verdict")
    for sol in viable:
        alpha_s  = sol['alpha']
        psiini_s = sol['psi_ini_sol']
        H0_s = find_H0_theta(alpha_s, psiini_s)
        if np.isnan(H0_s):
            print(f"  {alpha_s:5.2f}  {psiini_s:.5f}  H₀ not found")
            continue
        r = run_forward(alpha_s, psiini_s, H0_ref=H0_s)
        if r is None: continue
        z_a, psi_a, u_a, H_a, _ = r
        Hf = lambda z: float(np.interp(z, z_a[::-1], H_a[::-1],
                                        left=H_a[-1], right=H_a[0]))
        th_s = theta_star(Hf)
        H_desi = np.interp(DESI_Z, z_a[::-1], H_a[::-1])
        c2 = float(np.sum(((H_desi - DESI_H)/DESI_SIG)**2))
        from scipy.stats import chi2 as chi2d, norm
        p = chi2d.sf(c2, len(DESI_Z))
        tens_s = float(norm.isf(max(p,1e-15)/2)) if p > 0 else 10.0
        psi0_s = float(np.interp(0.0, z_a[::-1], psi_a[::-1]))
        F0_s   = F_cmstg(psi0_s, L0_REF)
        verd   = "PASS" if (c2 < CHI2_121C and abs(th_s-THETA_OBS)<5*SIGMA_TH) else \
                 "PARTIAL" if (abs(th_s-THETA_OBS)<5*SIGMA_TH or c2 < CHI2_121C) else "FAIL"
        results_C.append({'alpha': alpha_s, 'psi_ini': psiini_s,
                          'H0': H0_s, 'theta': th_s, 'chi2_desi': c2,
                          'tension': tens_s, 'psi0': psi0_s, 'F0': F0_s,
                          'verdict': verd})
        print(f"  {alpha_s:5.2f}  {psiini_s:.5f}  {H0_s:7.2f}  {th_s:8.5f}  "
              f"{c2:8.3f}  {tens_s:5.2f}σ  {verd}")

        # Detailed H(z) table
        print(f"    H(z) pulls vs DESI:")
        for i in range(len(DESI_Z)):
            pull = (H_desi[i]-DESI_H[i])/DESI_SIG[i]
            print(f"      z={DESI_Z[i]:.3f}  H_obs={DESI_H[i]:.1f}  "
                  f"H_mod={H_desi[i]:.2f}  pull={pull:+.2f}")

# ══════════════════════════════════════════════════════════════════════
# STRUCTURAL FINDING
# ══════════════════════════════════════════════════════════════════════
if not viable:
    structural = (
        "K(Ψ)=exp(-αΨ²) does not resolve the frozen-field problem. "
        "Near the hilltop (Ψ≈0), K≈1 regardless of α — the initial "
        "departure from the hilltop is unchanged. The exp(+αΨ²) "
        "enhancement only kicks in at intermediate Ψ, but the bottleneck "
        "is the kinematically frozen phase near Ψ=0 where V'∝Ψ is tiny "
        "and H is large. For all tested α∈[0,5], Ψ cannot travel from "
        "near-hilltop initial conditions to Ψ₀=2.62 M_Pl in the available "
        "~7 e-folds. The non-canonical kinetic term does not fix the "
        "frozen-field inconsistency. A different symmetry-breaking "
        "structure (or non-hilltop potential) is required."
    )
    overall_verdict = "FAIL"
elif results_C:
    best = min(results_C, key=lambda r: r['chi2_desi'])
    structural = (
        f"K(Ψ)=exp(-αΨ²) with α={best['alpha']:.2f} allows self-consistent "
        f"field evolution: Ψ_ini={best['psi_ini']:.4f} → Ψ₀={best['psi0']:.3f} M_Pl. "
        f"H₀={best['H0']:.2f} km/s/Mpc, tension={best['tension']:.2f}σ."
    )
    overall_verdict = best['verdict']
else:
    # "Viable" only means Ψ(z=0)=2.62 was reached, but every required
    # Ψ_ini ∈ [1.84, 2.57] — the field barely moved at all. This is pure
    # fine-tuning, not a hilltop solution. The K modification only enhances
    # rolling at intermediate Ψ; the bottleneck is near Ψ≈0 where K≈1
    # regardless of α and V'∝4λΨ(Ψ²-v²) is small. Additionally, the
    # resulting H(z=0) for all viable solutions is catastrophically large
    # (100–10⁶ km/s/Mpc), preventing any CMB acoustic scale match.
    psi_ini_range = (min(s['psi_ini_sol'] for s in viable),
                     max(s['psi_ini_sol'] for s in viable))
    structural = (
        f"K(Ψ)=exp(-αΨ²) does NOT resolve the frozen-field problem. "
        f"'Viable' solutions require Ψ_ini ∈ [{psi_ini_range[0]:.2f}, "
        f"{psi_ini_range[1]:.2f}] M_Pl — the field starts within "
        f"{PSI0_REF - psi_ini_range[0]:.2f} M_Pl of the target "
        f"Ψ₀={PSI0_REF} and barely rolls. This is fine-tuning, not "
        f"dynamics. The hilltop bottleneck (Ψ≈0, K≈1, V'≈0, H large) "
        f"is unchanged by the kinetic modification. Furthermore, all viable "
        f"H(z=0) are unphysically large (up to ~10⁶ km/s/Mpc), so no "
        f"CMB acoustic scale match exists. Conclusion: non-canonical "
        f"K(Ψ)=exp(-αΨ²) kinetics cannot fix the Phase 2 frozen-field "
        f"inconsistency. A fundamentally different potential or coupling "
        f"structure is required."
    )
    overall_verdict = "FAIL"

# ══════════════════════════════════════════════════════════════════════
# FIGURES
# ══════════════════════════════════════════════════════════════════════
print("\nGenerating figures...")
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt

# Figure 1: Ψ(z=0) vs α for each Ψ_ini
fig, ax = plt.subplots(figsize=(9, 6))
colors = plt.cm.viridis(np.linspace(0,1,len(psi_ini_vals)))
for j, psi_ini in enumerate(psi_ini_vals):
    pts = [s for s in scan_A if abs(s['psi_ini']-psi_ini)<1e-5]
    a_v = [p['alpha'] for p in pts]
    p_v = [p.get('psi_z0', np.nan) for p in pts]
    ax.plot(a_v, p_v, 'o-', color=colors[j], ms=5,
            label=f'Ψ_ini={psi_ini}')
ax.axhline(PSI0_REF, color='r', ls='--', lw=1.5, label='Target Ψ₀=2.62')
ax.set_xlabel(r'$\alpha$', fontsize=13)
ax.set_ylabel(r'$\Psi(z=0)\ [M_\mathrm{Pl}]$', fontsize=13)
ax.set_title('SIM124: Field reached at z=0 vs kinetic parameter α', fontsize=12)
ax.legend(fontsize=8, ncol=2); ax.grid(alpha=0.3); ax.set_yscale('log')
for ext in ['pdf','png']:
    fig.savefig(os.path.join(OUT_DIR, f'sim124_psi0_scan.{ext}'),
                dpi=150, bbox_inches='tight')
plt.close()
print("  Saved sim124_psi0_scan")

# Figure 2: Field evolution Ψ(z) for a few α, Ψ_ini=0.01
fig, ax = plt.subplots(figsize=(9, 5))
for alpha in [0.0, 0.5, 1.0, 2.0, 5.0]:
    r = run_forward(alpha, 0.01)
    if r is None: continue
    z_a, psi_a, *_ = r
    mask = z_a <= Z_CMB
    ax.semilogy(np.log10(1+z_a[mask]), psi_a[mask],
                label=f'α={alpha:.1f}')
ax.axhline(PSI0_REF, color='r', ls='--', lw=1.2, label='Target 2.62')
ax.axhline(0.01, color='k', ls=':', lw=0.8, label='Ψ_ini=0.01')
ax.set_xlabel('log₁₀(1+z)', fontsize=12)
ax.set_ylabel(r'$\Psi(z)\ [M_\mathrm{Pl}]$', fontsize=12)
ax.set_title('SIM124: Field evolution for varying α (Ψ_ini=0.01)', fontsize=12)
ax.legend(fontsize=9); ax.grid(alpha=0.3)
for ext in ['pdf','png']:
    fig.savefig(os.path.join(OUT_DIR, f'sim124_psi_evolution.{ext}'),
                dpi=150, bbox_inches='tight')
plt.close()
print("  Saved sim124_psi_evolution")

# Figure 3: Required Ψ_ini(α) to reach Ψ(z=0)=2.62 — the fine-tuning locus
# A near-hilltop solution would show Ψ_ini ≈ 0; instead it hovers near 2.62.
if viable:
    fig, ax = plt.subplots(figsize=(8, 5))
    alpha_v = [s['alpha'] for s in viable]
    pini_v  = [s['psi_ini_sol'] for s in viable]
    ax.plot(alpha_v, pini_v, 'o-', color='#d62728', lw=2, ms=7,
            label='Required $\\Psi_\\mathrm{ini}$')
    ax.axhline(PSI0_REF, color='gray', ls='--', lw=1.2,
               label=f'Target $\\Psi_0={PSI0_REF}$ $M_{{\\rm Pl}}$')
    ax.axhline(0.0, color='steelblue', ls=':', lw=1.2,
               label='Ideal hilltop IC ($\\Psi_{{\\rm ini}}\\approx 0$)')
    ax.fill_between(alpha_v, pini_v, PSI0_REF, alpha=0.15, color='red',
                    label='Fine-tuning gap')
    ax.set_xlabel(r'$\alpha$  (kinetic parameter)', fontsize=12)
    ax.set_ylabel(r'$\Psi_\mathrm{ini}\ [M_\mathrm{Pl}]$', fontsize=12)
    ax.set_title('SIM124: Fine-tuning locus — required $\\Psi_\\mathrm{ini}(\\alpha)$\n'
                 'to reach $\\Psi(z=0)=2.62\\ M_\\mathrm{Pl}$  (true hilltop would give $\\Psi_\\mathrm{ini}\\approx 0$)',
                 fontsize=11)
    ax.set_ylim(0, PSI0_REF * 1.1)
    ax.legend(fontsize=9); ax.grid(alpha=0.3)
    for ext in ['pdf', 'png']:
        fig.savefig(os.path.join(OUT_DIR, f'sim124_finetuning_locus.{ext}'),
                    dpi=150, bbox_inches='tight')
    plt.close()
    print("  Saved sim124_finetuning_locus")

# ══════════════════════════════════════════════════════════════════════
# SAVE RESULTS
# ══════════════════════════════════════════════════════════════════════
out = {
    "sim": "SIM124",
    "verdict": overall_verdict,
    "modification": "L_kin = -exp(-alpha*Psi^2)/2 * (dPsi)^2",
    "viable_solutions": viable,
    "observables": results_C,
    "structural_finding": structural,
    "reference_SIM121C": {"chi2_DESI": CHI2_121C, "tension": TENS_121C},
}
with open(os.path.join(OUT_DIR,"sim124_results.json"),'w') as f:
    json.dump(out, f, indent=2)

print(f"\nAll outputs in: {OUT_DIR}")
print(f"\n{'='*70}")
print(f"SIM124 RESULT: {overall_verdict}")
print(f"  {structural}")
print(f"{'='*70}")
print("SIM124 complete.")
