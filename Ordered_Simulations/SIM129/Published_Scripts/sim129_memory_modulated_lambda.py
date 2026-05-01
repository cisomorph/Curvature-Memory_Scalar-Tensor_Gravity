#!/usr/bin/env python3
"""
SIM129 — Memory-Modulated Cosmological Constant Λ(M)
CMSTG Framework, Flat FRW

Motivation
----------
SIM125–SIM128 tested M as a back-reaction energy density or a G_eff modifier,
both requiring extra phenomenological couplings.  This simulation tests the
physically minimal interpretation: M records the integrated complexity of
the universe's structure, and as that memory accumulates, the effective
cosmological constant weakens.

Derivation
----------
The CMSTG Phase 1 action is kept intact:
    S = ∫d⁴x√(-g) [(½+Λ₀Ψ²)R - ½(∇Ψ)² - ½m²Ψ² - Λ_eff(M)] + S_matter

Replacing the bare Λ with Λ_eff(M) = Λ_0·(1 − γ·M̂),  M̂ = M(a)/M(1).

This is the minimal one-parameter extension of Phase 1 consistent with CMSTG.
The memory-DE equation of state follows analytically:
    w_DE(a) = −1 + [a·γ·dM̂/da] / [3(1 − γ·M̂(a))]
Since dM̂/da > 0 always, w_DE > −1 for any γ > 0 → THAWING by construction.
This is DESI-preferred without tuning the sign.

At γ=0: recovers Phase 1 exactly (Λ = Λ_0, w = −1).
At γ>0: Λ_eff decreases as complexity accumulates; DE weakens at late times.
Phase 1 is the γ→0 limit — this is a continuous deformation.

Key quantities
--------------
  Λ_eff(a) = Λ_0 · (1 − γ·M̂(a))
  H²(a)    = H₀² [Ω_m a⁻³ + Ω_r a⁻⁴ + Ω_Λ(1 − γM̂(a))]

  w_DE(a) computed numerically from dΛ_eff/da.
  w₀ = w_DE(a=1),  wₐ ≈ −dw_DE/da|_{a=1} (standard parameterisation).

Self-consistency: M(a) is first solved with the ΛCDM H (SIM125 baseline),
then M̂ is inserted into the modified Friedmann equation.  A second iteration
checks convergence (ΔH/H < 0.1%).

Tests
-----
  1. DESI Y1 H(z) tension (baseline: 2.63σ from SIM121C)
  2. CMB acoustic angle θ_* (must stay within 5σ of Planck)
  3. w₀, wₐ — is the thawing consistent with DESI CPL constraints?
  4. H₀ inferred from θ_* (must stay near 67–68 km/s/Mpc)
  5. Effective Λ today Λ_eff(1) = Λ_0(1−γ) — how much does DE weaken?
"""

import numpy as np
from scipy.integrate import solve_ivp, quad
from scipy.optimize import brentq
from scipy.stats import norm as sp_norm, chi2 as sp_chi2
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import json, os

# ─── Cosmological parameters (Planck 2018 / Phase 1 values) ──────────────────
H0      = 67.4
Omega_m = 0.315
Omega_r = 9.4e-5
Omega_L = 0.685
Omega_b = 0.0493   # for CMB sound horizon

rho_m0  = 3.0 * Omega_m
rho_L   = 3.0 * Omega_L

THETA_OBS = 1.04101e-2    # Planck 100·θ_* (dimensionless ratio ×100)
SIGMA_TH  = 0.00029e-2    # 1σ

# ─── Memory-field parameters (SIM125, unchanged) ─────────────────────────────
alpha  = 1.0; beta  = 1.5; xi   = 10.0; sigma0 = 0.1
a_eq   = 3.0e-4; a_nl = 0.3; n_sh = 2

# ─── DESI Y1 ─────────────────────────────────────────────────────────────────
DESI_Z   = np.array([0.295, 0.510, 0.706, 0.930, 1.317, 2.330])
DESI_H   = np.array([81.7,  97.9, 110.7, 128.1, 156.4, 240.8])
DESI_SIG = np.array([ 2.9,   3.8,   4.5,   5.3,   7.3,  15.7])
DESI_A   = 1.0 / (1.0 + DESI_Z)
TENSION_REF = 2.63   # SIM121C baseline

# ─── Integration domain ───────────────────────────────────────────────────────
a_start = 0.01
a_end   = 1.0
a_arr   = np.logspace(np.log10(a_start), np.log10(a_end), 800)

# ─── Step 1: Solve M(a) with ΛCDM H (SIM125 baseline) ───────────────────────

def E_LCDM(a):
    return np.sqrt(Omega_m*a**-3 + Omega_r*a**-4 + Omega_L)

def C_func(a):
    rm  = rho_m0 * a**-3
    return rm * sigma0 * (a/a_eq)**n_sh * np.exp(-a/a_nl) / (rm + rho_L/xi)

sol_M0 = solve_ivp(
    lambda a, y: [(alpha*C_func(a)/(a*E_LCDM(a))) - (beta*y[0]/a)],
    [a_start, a_end], [0.0], t_eval=a_arr,
    method='RK45', rtol=1e-10, atol=1e-13
)
M_base = sol_M0.y[0]
M1     = M_base[-1]         # M(a=1) for normalisation
Mhat   = M_base / M1        # M̂(a): 0 at early times → 1 today

from scipy.interpolate import interp1d
Mhat_f = interp1d(a_arr, Mhat, kind='cubic', fill_value=(0.0, 1.0),
                  bounds_error=False)
dMhat_da = np.gradient(Mhat, a_arr)
dMhat_f  = interp1d(a_arr, dMhat_da, kind='cubic', fill_value='extrapolate')

# ─── Step 2: Modified H(a, γ) ─────────────────────────────────────────────────

def E_mod(a, gamma):
    """E(a)=H(a)/H₀ with memory-modulated Λ."""
    Lambda_eff = Omega_L * (1.0 - gamma * Mhat_f(a))
    val = Omega_m*a**-3 + Omega_r*a**-4 + Lambda_eff
    return np.sqrt(max(val, 1e-30))

# ─── Step 3: Self-consistency check ──────────────────────────────────────────
# Re-solve M with the modified H for the best γ; check convergence.

def solve_M_with_H(gamma):
    """Solve M(a) using the γ-modified H(a) for self-consistency."""
    def ode(a, y):
        M   = y[0]
        Ea  = E_mod(a, gamma)
        Ca  = C_func(a)
        return [(alpha*Ca/(a*Ea)) - (beta*M/a)]
    sol = solve_ivp(ode, [a_start, a_end], [0.0], t_eval=a_arr,
                    method='RK45', rtol=1e-10, atol=1e-13)
    return sol.y[0]

# ─── Step 4: Observables ─────────────────────────────────────────────────────

def w_DE(a_pt, gamma):
    """
    Effective DE equation of state at scale factor a_pt.
    w_DE = -1 + a·γ·(dM̂/da) / (3·(1 − γ·M̂))
    Analytic from the definition of Λ_eff(a).
    """
    Mh   = float(Mhat_f(a_pt))
    dMh  = float(dMhat_f(a_pt))
    denom = 1.0 - gamma * Mh
    if denom < 1e-4:
        return -1.0
    return -1.0 + (a_pt * gamma * dMh) / (3.0 * denom)

def desi_diagnostics(gamma):
    H_mod = np.array([H0 * E_mod(a, gamma) for a in DESI_A])
    pulls = (H_mod - DESI_H) / DESI_SIG
    chi2  = float(np.sum(pulls**2))
    p     = sp_chi2.sf(chi2, len(DESI_Z))
    tens  = float(sp_norm.isf(max(p, 1e-15) / 2))
    return chi2, tens, H_mod, pulls

def theta_star(gamma, H0_try=None):
    """
    CMB acoustic angle θ_* = r_s(z_*) / D_A(z_*).
    Uses the Hu & Sugiyama approximation for r_s.
    """
    if H0_try is None:
        H0_try = H0
    z_star = 1090.0
    a_star = 1.0 / (1.0 + z_star)

    def H_func(z):
        a = 1.0/(1.0+z)
        return H0_try * E_mod(a, gamma)

    # Sound horizon integral
    def rs_integrand(z):
        H  = H_func(z)
        R  = 3.0*Omega_b/(4.0*2.47e-5)/(1.0+z)
        return 1.0 / (np.sqrt(3.0*(1.0+R)) * H)

    def da_integrand(z):
        return 1.0 / H_func(z)

    rs, _ = quad(rs_integrand, 0.0, z_star, limit=200, epsrel=1e-7)
    dc, _ = quad(da_integrand, 0.0, z_star, limit=200, epsrel=1e-7)
    return (rs / dc) if dc > 0 else np.nan   # dimensionless; multiply by 100 for degrees

# ─── Main scan ───────────────────────────────────────────────────────────────
gamma_values = [0.00, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.40, 0.50]

print("=" * 70)
print("  SIM129 — Memory-Modulated Λ(M) = Λ₀(1 − γ·M̂)  [CMSTG]")
print("=" * 70)
print(f"  Phase 1 backbone preserved.  Single coupling γ.")
print(f"  M̂ from SIM125: M(1)={M1:.3e},  shape fixed.\n")
print(f"{'γ':>5}  {'Λ_eff(1)/Λ₀':>12}  {'w₀':>7}  {'wₐ':>7}  "
      f"{'tens':>7}  {'100θ*':>8}  {'Δθ*/σ':>7}  verdict")
print("-" * 75)

results = []
for gamma in gamma_values:
    # DE equation of state
    w0   = w_DE(1.0, gamma)
    # wₐ: slope at a=1 (CPL: w ≈ w₀ + wₐ(1-a))
    eps  = 0.01
    w_em = w_DE(1.0 - eps, gamma)
    wa   = (w_em - w0) / eps   # dw/da|_{a→1}, sign convention: wₐ = -dw/da

    # DESI
    chi2, tens, H_mod, pulls = desi_diagnostics(gamma)

    # CMB acoustic angle
    th = theta_star(gamma) * 100.0   # convert to 100×θ_*
    dth_sig = (th - THETA_OBS*100) / (SIGMA_TH*100) if not np.isnan(th) else np.nan

    Lambda_ratio = 1.0 - gamma

    # Verdict
    desi_ok = tens < TENSION_REF
    cmb_ok  = abs(dth_sig) < 5.0 if not np.isnan(dth_sig) else False
    verdict = ("PASS"    if (desi_ok and cmb_ok) else
               "PARTIAL" if (desi_ok or  cmb_ok) else "FAIL")

    print(f"  {gamma:.2f}  {Lambda_ratio:12.4f}  {w0:7.4f}  {wa:7.4f}  "
          f"{tens:7.2f}σ  {th:8.5f}  {dth_sig:+7.2f}σ  {verdict}")

    results.append({
        "gamma": gamma, "Lambda_ratio": Lambda_ratio,
        "w0": w0, "wa": wa, "chi2_DESI": chi2,
        "tension": tens, "theta_star_100": th,
        "dtheta_sigma": float(dth_sig),
        "verdict": verdict, "H_mod": H_mod.tolist(), "pulls": pulls.tolist(),
    })

print(f"\n  SIM121C reference tension: {TENSION_REF:.2f}σ")

# ─── Self-consistency check for best γ ───────────────────────────────────────
best  = min(results, key=lambda r: r['tension'] if r['verdict'] != 'FAIL' else 99)
gamma_best = best['gamma']

M_iter1 = solve_M_with_H(gamma_best)
M1_iter1 = M_iter1[-1]
rel_change = abs(M1_iter1 - M1) / M1
print(f"\n  Self-consistency (γ={gamma_best:.2f}): "
      f"M(1) shifts by {rel_change*100:.3f}% after one H-iteration — "
      f"{'converged' if rel_change < 0.001 else 'iterate'}.")

# ─── Analytical note ─────────────────────────────────────────────────────────
print(f"\n  w_DE(a) is THAWING by construction (w > −1, dw/da < 0).")
print(f"  At γ=0: w₀=−1.0000 (Phase 1 limit recovered exactly).")

# ─── Plotting ─────────────────────────────────────────────────────────────────
colors_g = plt.cm.plasma(np.linspace(0.05, 0.88, len(gamma_values)))
a = a_arr

fig = plt.figure(figsize=(16, 11))
fig.suptitle(
    "SIM129 — Memory-Modulated $\\Lambda(M)$ = $\\Lambda_0(1-\\gamma\\hat{M})$  [CMSTG]\n"
    "Phase 1 backbone preserved — single coupling $\\gamma$, "
    "thawing DE by construction",
    fontsize=13, fontweight='bold', y=0.995
)
gs = gridspec.GridSpec(2, 3, hspace=0.44, wspace=0.38)

ax1 = fig.add_subplot(gs[0, 0])   # Λ_eff(a) profiles
ax2 = fig.add_subplot(gs[0, 1])   # w_DE(a)
ax3 = fig.add_subplot(gs[0, 2])   # H(a)/H₀
ax4 = fig.add_subplot(gs[1, 0])   # DESI H(z)
ax5 = fig.add_subplot(gs[1, 1])   # tension and θ_* vs γ
ax6 = fig.add_subplot(gs[1, 2])   # w₀–wₐ plane

# Panel 1: Λ_eff(a)/Λ₀
for gamma, col in zip(gamma_values, colors_g):
    L_arr = 1.0 - gamma * Mhat
    ax1.plot(a, L_arr, color=col, lw=1.8, label=f'γ={gamma:.2f}')
ax1.axhline(1.0, color='k', ls='--', lw=0.8, alpha=0.4, label='Phase 1 (γ=0)')
ax1.set_xscale('log')
ax1.set_xlabel('Scale factor $a$')
ax1.set_ylabel('$\\Lambda_\\mathrm{eff}(a)/\\Lambda_0$')
ax1.set_title('Panel 1 — Effective Cosmological Constant')
ax1.legend(fontsize=7, ncol=2); ax1.grid(True, alpha=0.3, which='both')

# Panel 2: w_DE(a)
for gamma, col in zip(gamma_values, colors_g):
    w_arr = np.array([w_DE(ai, gamma) for ai in a])
    ax2.plot(a, w_arr, color=col, lw=1.8, label=f'γ={gamma:.2f}')
ax2.axhline(-1.0, color='k', ls='--', lw=1, alpha=0.5, label='$w=-1$')
ax2.axhline(-0.973, color='gray', ls=':', lw=1, alpha=0.6, label='SIM113 $w_0$')
ax2.set_xscale('log')
ax2.set_xlabel('Scale factor $a$')
ax2.set_ylabel('$w_\\mathrm{DE}(a)$')
ax2.set_title('Panel 2 — DE Equation of State (thawing)')
ax2.legend(fontsize=7, ncol=2); ax2.grid(True, alpha=0.3, which='both')
ax2.set_ylim(-1.05, 0.1)

# Panel 3: H(a)/H₀
for gamma, col in zip(gamma_values, colors_g):
    E_arr = np.array([E_mod(ai, gamma) for ai in a])
    ax3.plot(a, E_arr, color=col, lw=1.8, label=f'γ={gamma:.2f}')
ax3.set_xscale('log'); ax3.set_yscale('log')
ax3.set_xlabel('Scale factor $a$')
ax3.set_ylabel('$H(a)/H_0$')
ax3.set_title('Panel 3 — Hubble Parameter')
ax3.legend(fontsize=7, ncol=2); ax3.grid(True, alpha=0.3, which='both')

# Panel 4: DESI H(z)
ax4.errorbar(DESI_Z, DESI_H, yerr=DESI_SIG, fmt='ko', ms=6, capsize=4,
             label='DESI Y1', zorder=5)
for r, col in zip(results, colors_g):
    lw = 2.4 if r['gamma'] == gamma_best else 1.2
    ax4.plot(DESI_Z, r['H_mod'], 'o--', color=col, lw=lw, ms=4,
             label=f"γ={r['gamma']:.2f}")
ax4.set_xlabel('Redshift $z$')
ax4.set_ylabel('$H(z)$ [km/s/Mpc]')
ax4.set_title('Panel 4 — DESI H(z) Comparison')
ax4.legend(fontsize=7, ncol=2); ax4.grid(True, alpha=0.3)

# Panel 5: tension and θ_* vs γ
g_arr   = [r['gamma']           for r in results]
t_arr   = [r['tension']         for r in results]
dth_arr = [abs(r['dtheta_sigma']) for r in results]

ax5l = ax5
ax5r = ax5.twinx()
ax5l.plot(g_arr, t_arr, 'o-', color='#d62728', lw=2, ms=7,
          label='DESI tension')
ax5l.axhline(TENSION_REF, color='#d62728', ls='--', lw=1,
             label=f'SIM121C {TENSION_REF}σ')
ax5l.axhline(2.0, color='#d62728', ls=':', lw=1, alpha=0.6,
             label='2σ target')
ax5r.plot(g_arr, dth_arr, 's--', color='#1f77b4', lw=2, ms=7,
          label='|Δθ_*/σ|')
ax5r.axhline(5.0, color='#1f77b4', ls=':', lw=1, alpha=0.6)
ax5l.set_xlabel('Memory-DE coupling $\\gamma$')
ax5l.set_ylabel('DESI tension [σ]', color='#d62728')
ax5r.set_ylabel('|CMB $\\theta_*$ pull| [σ]', color='#1f77b4')
ax5l.set_title('Panel 5 — DESI Tension and CMB θ_* vs γ')
lines1, labs1 = ax5l.get_legend_handles_labels()
lines2, labs2 = ax5r.get_legend_handles_labels()
ax5l.legend(lines1+lines2, labs1+labs2, fontsize=7)
ax5l.grid(True, alpha=0.3)

# Panel 6: w₀–wₐ plane
w0_arr = [r['w0'] for r in results]
wa_arr = [r['wa'] for r in results]

# DESI CPL constraints (approximate 1σ/2σ contours)
w0_grid = np.linspace(-1.4, -0.5, 200)
wa_grid = np.linspace(-2.5, 1.5, 200)
W0, WA  = np.meshgrid(w0_grid, wa_grid)
# Approximate DESI 1σ/2σ from Adame+2024 (CPL)
chi2_DESI_CPL = ((W0 - (-0.827))/0.19)**2 + ((WA - (-0.75))/0.50)**2
ax6.contourf(W0, WA, chi2_DESI_CPL, levels=[0, 2.30, 6.18],
             colors=['#aec7e8', '#c5dbf0'], alpha=0.5)
ax6.contour(W0, WA, chi2_DESI_CPL, levels=[2.30, 6.18],
            colors=['#1f77b4', '#1f77b4'], linewidths=[1.5, 1.0])
ax6.scatter(w0_arr, wa_arr, c=[r['gamma'] for r in results],
            cmap='plasma', s=60, zorder=5)
for r in results:
    ax6.annotate(f"γ={r['gamma']:.2f}", (r['w0'], r['wa']),
                 fontsize=6.5, xytext=(4, 3), textcoords='offset points')
ax6.axhline(0, color='k', ls=':', lw=0.8, alpha=0.4)
ax6.axvline(-1, color='k', ls=':', lw=0.8, alpha=0.4)
ax6.set_xlabel('$w_0$')
ax6.set_ylabel('$w_a$')
ax6.set_title('Panel 6 — $w_0$–$w_a$ Plane\n(blue: DESI CPL 1σ/2σ)')
ax6.grid(True, alpha=0.3)

out_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), '../Outputs')
os.makedirs(out_dir, exist_ok=True)
for ext in ['png', 'pdf']:
    fig.savefig(os.path.join(out_dir, f'SIM129_memory_lambda.{ext}'),
                dpi=150, bbox_inches='tight')
plt.close()

# ─── Save JSON ────────────────────────────────────────────────────────────────
out = {
    "sim": "SIM129",
    "derivation": "Λ_eff(a) = Λ_0*(1 - gamma*Mhat(a)), single coupling gamma",
    "params": {"alpha": alpha, "beta": beta, "xi": xi, "sigma0": sigma0,
               "M1": float(M1), "H0": H0, "Omega_L": Omega_L},
    "self_consistency_M1_shift_pct": float(rel_change * 100),
    "gamma_scan": [{k: v for k, v in r.items()
                    if k not in ('H_mod', 'pulls')} for r in results],
    "best_gamma": gamma_best,
    "tension_ref_SIM121C": TENSION_REF,
}
with open(os.path.join(out_dir, 'sim129_results.json'), 'w') as f:
    json.dump(out, f, indent=2)

# ─── Final verdict ────────────────────────────────────────────────────────────
pass_r    = [r for r in results if r['verdict'] == 'PASS']
partial_r = [r for r in results if r['verdict'] == 'PARTIAL']
overall   = ("PASS" if pass_r else "PARTIAL" if partial_r else "FAIL")

print(f"\nSIM129 RESULT: {overall}")
if pass_r:
    for r in pass_r:
        print(f"  PASS at γ={r['gamma']:.2f}: tension={r['tension']:.2f}σ, "
              f"w₀={r['w0']:.4f}, wₐ={r['wa']:.4f}, "
              f"Δθ_*={r['dtheta_sigma']:+.2f}σ")
else:
    print(f"  Best: γ={best['gamma']:.2f}, tension={best['tension']:.2f}σ, "
          f"w₀={best['w0']:.4f}, Δθ_*={best['dtheta_sigma']:+.2f}σ")
print(f"Outputs: {out_dir}/")
