#!/usr/bin/env python3
"""
SIM126 — Cosmic Memory Field M(a): Back-reaction on Friedmann Equation
CMSTG Framework, Flat FRW

SIM125 showed M evolves cleanly as a test field. Here M back-reacts on H(a):

    H²(a) = H₀² [Ω_m a⁻³ + Ω_r a⁻⁴ + Ω_Λ + ε·M²(a)]

The ε·M² term treats M as contributing an effective energy density to the
total budget.  ε is the dimensionless back-reaction coupling.  The M ODE
is now self-consistent: H appears in dM/da, and M appears in H².

Key questions:
  1. Does the M back-reaction modify w_DE(a) away from −1?
  2. Does it reduce or worsen the DESI H(z) tension (baseline: 2.63σ from SIM121C)?
  3. What ε range gives cosmologically significant back-reaction?
  4. Is M back-reaction thawing (wₐ < 0, DESI-preferred) or freezing?

Units: natural (8πG=1); densities in H₀²; time in 1/H₀; H → E(a) = H(a)/H₀.
"""

import numpy as np
from scipy.integrate import solve_ivp
from scipy.stats import norm as sp_norm
from scipy.stats import chi2 as sp_chi2
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import json, os

# ─── Cosmological parameters ─────────────────────────────────────────────────
H0      = 67.4
Omega_m = 0.315
Omega_r = 9.4e-5
Omega_L = 0.685

rho_m0  = 3.0 * Omega_m    # 0.945  [H₀²]
rho_L   = 3.0 * Omega_L    # 2.055  [H₀²]

# ─── Memory-field parameters (same as SIM125) ─────────────────────────────────
alpha  = 1.0
beta   = 1.5
xi     = 10.0
sigma0 = 0.1      # fiducial shear amplitude
a_eq   = 3.0e-4
a_nl   = 0.3
n_sh   = 2

# ─── Back-reaction coupling scan ─────────────────────────────────────────────
# Natural scale: ε·M(1)² ≈ Ω_Λ → ε_nat ≈ Ω_Λ / M(1)²
# M(1) ≈ 1.859×10⁴ (from SIM125, σ₀=0.1) → ε_nat ≈ 2.0×10⁻⁹
# We scan 5 decades around ε_nat plus ε=0 (no back-reaction / SIM125 limit)
EPS_NAT   = 2.0e-9
eps_values = [0.0, EPS_NAT * 1e-2, EPS_NAT * 0.1, EPS_NAT, EPS_NAT * 10]
eps_labels = ['ε=0 (SIM125)', 'ε=0.01ε_nat', 'ε=0.1ε_nat', 'ε=ε_nat', 'ε=10ε_nat']

# ─── DESI Y1 data ────────────────────────────────────────────────────────────
DESI_Z   = np.array([0.295, 0.510, 0.706, 0.930, 1.317, 2.330])
DESI_H   = np.array([81.7,  97.9, 110.7, 128.1, 156.4, 240.8])
DESI_SIG = np.array([ 2.9,   3.8,   4.5,   5.3,   7.3,  15.7])
DESI_A   = 1.0 / (1.0 + DESI_Z)

# SIM121C reference tension for comparison
TENSION_REF = 2.63

# ─── Integration domain ───────────────────────────────────────────────────────
a_start = 0.01
a_end   = 1.0
a_eval  = np.logspace(np.log10(a_start), np.log10(a_end), 800)

# ─── Physical functions ───────────────────────────────────────────────────────

def E_BR(a, M, eps):
    """
    Modified E(a) = H(a)/H₀ with M back-reaction.
    The ε·M² term adds M's effective energy density to the Friedmann equation.
    ε=0 recovers SIM125 (test field, no back-reaction).
    """
    val = Omega_m * a**-3 + Omega_r * a**-4 + Omega_L + eps * M**2
    return np.sqrt(max(val, 1e-30))

def shear(a):
    return (a / a_eq)**n_sh * np.exp(-a / a_nl)

def C_func(a):
    """Complexity source — unchanged from SIM125; does not depend on H or M."""
    rm  = rho_m0 * a**-3
    reg = rho_L / xi
    return rm * sigma0 * shear(a) / (rm + reg)

def dM_da(a, M_arr, eps):
    """
    Coupled ODE: M feeds into E_BR, E_BR feeds into dM/da.
    dM/da = α·C/(a·E_BR) − β·M/a
    """
    M  = M_arr[0]
    Ea = E_BR(a, M, eps)
    Ca = C_func(a)
    return [(alpha * Ca / (a * Ea)) - (beta * M / a)]

# ─── Integration ─────────────────────────────────────────────────────────────
print("=" * 60)
print("  SIM126 — Memory Field Back-reaction on H(a)  [CMSTG]")
print("=" * 60)
print(f"  ε_nat = {EPS_NAT:.2e}  (gives ε·M²(1) ≈ Ω_Λ)")
print(f"  σ₀ = {sigma0},  α = {alpha},  β = {beta},  ξ = {xi}")
print()

solutions = {}
for eps, label in zip(eps_values, eps_labels):
    sol = solve_ivp(
        dM_da, [a_start, a_end], [0.0],
        args=(eps,), method='RK45', t_eval=a_eval,
        rtol=1e-9, atol=1e-12,
    )
    solutions[eps] = (sol, label)

# ─── Diagnostics ─────────────────────────────────────────────────────────────

def compute_wDE(a_arr, M_arr, eps):
    """
    Effective DE equation of state from the M back-reaction sector.
    ρ_DE(a) = Ω_Λ + ε·M²   →   w_DE = −1 − (a/3ρ_DE) dρ_DE/da
    At ε=0 this is exactly −1; deviations signal M is doing work on the expansion.
    """
    rho_DE = Omega_L + eps * M_arr**2
    # numerical derivative dρ_DE/da
    drho   = np.gradient(rho_DE, a_arr)
    w_DE   = -1.0 - (a_arr * drho) / (3.0 * rho_DE)
    return rho_DE, w_DE

def desi_tension(a_arr, M_arr, eps):
    """Compute DESI chi² and tension from the modified H(z)."""
    H_mod = np.array([H0 * E_BR(a, float(np.interp(a, a_arr, M_arr)), eps)
                      for a in DESI_A])
    pulls = (H_mod - DESI_H) / DESI_SIG
    chi2  = float(np.sum(pulls**2))
    p     = sp_chi2.sf(chi2, len(DESI_Z))
    tens  = float(sp_norm.isf(max(p, 1e-15) / 2))
    return chi2, tens, H_mod, pulls

# ─── Print results ────────────────────────────────────────────────────────────
print(f"{'Label':<22}  {'M(1)':>10}  {'ε·M²(1)':>10}  {'χ²_DESI':>9}  {'tension':>8}  {'verdict'}")
print("-" * 80)

results_table = []
for eps, label in zip(eps_values, eps_labels):
    sol, lbl = solutions[eps]
    M_final  = sol.y[0, -1]
    eM2      = eps * M_final**2
    chi2, tens, H_mod, pulls = desi_tension(sol.t, sol.y[0], eps)
    _, w_DE  = compute_wDE(sol.t, sol.y[0], eps)
    w0       = float(w_DE[-1])
    w_early  = float(w_DE[0])
    wa_approx = w_early - w0   # sign: positive = thawing-like
    # DESI prefers w0 > -1 (thawing) and wₐ < 0 — check direction
    verdict = ("THAWING" if (w0 > -1.0 and wa_approx < 0) else
               "FREEZING" if (w0 < -1.0 or wa_approx > 0) else "FLAT")
    print(f"  {lbl:<20}  {M_final:>10.3e}  {eM2:>10.4f}  {chi2:>9.3f}  "
          f"{tens:>7.2f}σ  {verdict}  w₀={w0:.4f}")
    results_table.append({
        "label": lbl, "eps": eps, "M_final": float(M_final),
        "eps_M2_final": float(eM2), "chi2_DESI": chi2,
        "tension": tens, "w0_DE": w0, "wa_approx": float(wa_approx),
        "verdict_DE": verdict,
    })

print(f"\n  SIM121C reference tension: {TENSION_REF:.2f}σ")

# ─── DESI pull table for ε=ε_nat ──────────────────────────────────────────────
print(f"\nDESI H(z) pulls for ε=ε_nat:")
print(f"  {'z':>6}  {'H_obs':>7}  {'H_mod':>7}  {'pull':>6}")
sol_nat, _ = solutions[EPS_NAT]
_, _, H_mod_nat, pulls_nat = desi_tension(sol_nat.t, sol_nat.y[0], EPS_NAT)
for i, z in enumerate(DESI_Z):
    print(f"  {z:6.3f}  {DESI_H[i]:7.1f}  {H_mod_nat[i]:7.2f}  {pulls_nat[i]:+6.2f}")

# ─── Plotting ─────────────────────────────────────────────────────────────────
colors_br = ['#888888', '#a6cee3', '#1f78b4', '#e31a1c', '#ff7f00']
a = a_eval

fig = plt.figure(figsize=(15, 11))
fig.suptitle(
    "SIM126 — Memory Field M(a): Back-reaction on Friedmann Equation  [CMSTG]\n"
    r"$H^2(a) = H_0^2[\Omega_m a^{-3} + \Omega_r a^{-4} + \Omega_\Lambda + \varepsilon M^2(a)]$",
    fontsize=13, fontweight='bold', y=0.995
)
gs = gridspec.GridSpec(2, 3, hspace=0.44, wspace=0.38)

ax1 = fig.add_subplot(gs[0, 0])   # H(a)/H₀
ax2 = fig.add_subplot(gs[0, 1])   # M(a)
ax3 = fig.add_subplot(gs[0, 2])   # ε·M²(a)  — back-reaction magnitude
ax4 = fig.add_subplot(gs[1, 0])   # w_DE(a)
ax5 = fig.add_subplot(gs[1, 1])   # DESI H(z) comparison
ax6 = fig.add_subplot(gs[1, 2])   # DESI tension vs ε

for (eps, label), col in zip(zip(eps_values, eps_labels), colors_br):
    sol, lbl = solutions[eps]
    M_arr    = sol.y[0]
    E_arr    = np.array([E_BR(ai, mi, eps) for ai, mi in zip(sol.t, M_arr)])
    eM2_arr  = eps * M_arr**2
    _, w_arr = compute_wDE(sol.t, M_arr, eps)

    lw = 2.2 if eps == EPS_NAT else 1.4
    ls = '-'

    ax1.plot(sol.t, E_arr, color=col, lw=lw, ls=ls, label=lbl)
    ax2.plot(sol.t, M_arr, color=col, lw=lw, ls=ls, label=lbl)
    ax3.plot(sol.t, eM2_arr, color=col, lw=lw, ls=ls, label=lbl)
    ax4.plot(sol.t, w_arr, color=col, lw=lw, ls=ls, label=lbl)

# Panel 1: H(a)/H₀
ax1.set_xscale('log'); ax1.set_yscale('log')
ax1.set_xlabel('Scale factor $a$'); ax1.set_ylabel('$H(a)/H_0$')
ax1.set_title('Panel 1 — Hubble Parameter')
ax1.legend(fontsize=7); ax1.grid(True, alpha=0.3, which='both')

# Panel 2: M(a)
ax2.set_xscale('log')
ax2.set_xlabel('Scale factor $a$'); ax2.set_ylabel('$M(a)$')
ax2.set_title('Panel 2 — Memory Field')
ax2.legend(fontsize=7); ax2.grid(True, alpha=0.3, which='both')

# Panel 3: ε·M²(a)  vs Ω_Λ
ax3.axhline(Omega_L, color='k', ls='--', lw=1, label=f'$\\Omega_\\Lambda={Omega_L}$')
ax3.set_xscale('log')
ax3.set_xlabel('Scale factor $a$'); ax3.set_ylabel(r'$\varepsilon\,M^2(a)$  [$H_0^2$]')
ax3.set_title('Panel 3 — Back-reaction Magnitude')
ax3.legend(fontsize=7); ax3.grid(True, alpha=0.3, which='both')

# Panel 4: w_DE(a)
ax4.axhline(-1.0, color='k', ls=':', lw=1, label='$w=-1$ (ΛCDM)')
ax4.axhline(-0.973, color='gray', ls='--', lw=0.8, label='SIM113 $w_0=-0.973$')
ax4.set_xscale('log')
ax4.set_xlabel('Scale factor $a$')
ax4.set_ylabel('$w_{\\rm DE,eff}(a)$')
ax4.set_title('Panel 4 — Effective DE Equation of State')
ax4.legend(fontsize=7); ax4.grid(True, alpha=0.3, which='both')
ax4.set_ylim(-2.0, 0.5)

# Panel 5: DESI H(z) comparison
ax5.errorbar(DESI_Z, DESI_H, yerr=DESI_SIG, fmt='ko', ms=6, capsize=4,
             label='DESI Y1', zorder=5)
for (eps, label), col in zip(zip(eps_values, eps_labels), colors_br):
    sol, lbl = solutions[eps]
    H_mod = np.array([H0 * E_BR(a_i, float(np.interp(a_i, sol.t, sol.y[0])), eps)
                      for a_i in DESI_A])
    lw = 2.2 if eps == EPS_NAT else 1.2
    ax5.plot(DESI_Z, H_mod, 'o--', color=col, lw=lw, ms=4, label=lbl)
ax5.set_xlabel('Redshift $z$'); ax5.set_ylabel('$H(z)$ [km/s/Mpc]')
ax5.set_title('Panel 5 — DESI H(z) Comparison')
ax5.legend(fontsize=7); ax5.grid(True, alpha=0.3)

# Panel 6: DESI tension vs ε (skip ε=0)
eps_plot   = [r['eps']     for r in results_table if r['eps'] > 0]
tens_plot  = [r['tension'] for r in results_table if r['eps'] > 0]
ax6.semilogx(eps_plot, tens_plot, 'o-', color='#d62728', lw=2, ms=7)
ax6.axhline(TENSION_REF, color='gray', ls='--', lw=1.2,
            label=f'SIM121C baseline {TENSION_REF:.2f}σ')
ax6.axhline(2.0, color='steelblue', ls=':', lw=1,
            label='2σ target')
ax6.set_xlabel(r'Back-reaction coupling $\varepsilon$')
ax6.set_ylabel('DESI tension [σ]')
ax6.set_title('Panel 6 — DESI Tension vs ε')
ax6.legend(fontsize=8); ax6.grid(True, alpha=0.3, which='both')

out_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), '../Outputs')
os.makedirs(out_dir, exist_ok=True)

for ext in ['png', 'pdf']:
    fig.savefig(os.path.join(out_dir, f'SIM126_memory_backreaction.{ext}'),
                dpi=150, bbox_inches='tight')
plt.close()

# ─── Save JSON ────────────────────────────────────────────────────────────────
out = {
    "sim": "SIM126",
    "params": {"alpha": alpha, "beta": beta, "xi": xi, "sigma0": sigma0,
               "eps_nat": EPS_NAT, "H0": H0, "Omega_m": Omega_m,
               "Omega_L": Omega_L},
    "eps_scan": results_table,
    "tension_reference_SIM121C": TENSION_REF,
}
with open(os.path.join(out_dir, 'sim126_results.json'), 'w') as f:
    json.dump(out, f, indent=2)

# ─── Summary verdict ──────────────────────────────────────────────────────────
best = min(results_table, key=lambda r: r['tension'])
nat  = next(r for r in results_table if abs(r['eps'] - EPS_NAT) < 1e-20)
print(f"\nBest tension: {best['tension']:.2f}σ at {best['label']}")
print(f"ε=ε_nat:      {nat['tension']:.2f}σ, w₀={nat['w0_DE']:.4f}, "
      f"DE type={nat['verdict_DE']}")

if best['tension'] < TENSION_REF:
    delta = TENSION_REF - best['tension']
    verdict_str = f"PARTIAL — back-reaction reduces tension by {delta:.2f}σ"
else:
    verdict_str = "FAIL — back-reaction does not reduce DESI tension"

print(f"\nSIM126 RESULT: {verdict_str}")
print(f"Outputs: {out_dir}/")
