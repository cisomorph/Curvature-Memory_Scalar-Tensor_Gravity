#!/usr/bin/env python3
"""
SIM128 — Cosmic Memory Field M(a): Power-Law Growth Suppression Profile
CMSTG Framework, Flat FRW

SIM127 showed G_eff/G = 1−νM̂(a) suppresses S₈ in the right direction
but distorts f·σ₈(z) at intermediate z (χ²_RSD=42.3/14).  The problem:
M̂(a) already grows substantially by a~0.3, so the suppression kicks in
too early, flattening f·σ₈ across all z when RSD data needs ΛCDM-like
shape up to z~0.6.

This simulation tests:
    G_eff(a) / G = 1 − ν · M̂(a)^p

Increasing p concentrates the suppression near a=1 (where M̂→1 rapidly)
while leaving G_eff≈1 at high z where M̂^p → 0.

  p=1  : SIM127 (baseline; too early)
  p>1  : increasingly late-time suppression
  p→∞  : G_eff≈1 everywhere except right at a=1

Strategy:
  2D scan over (p, ν).  For each p, find the ν that minimises a joint
  cost combining the S₈ pull and the RSD χ².  Report the Pareto-optimal
  (p, ν) that simultaneously satisfies the S₈ and RSD constraints.

Physical rationale:
  High p corresponds to M needing to accumulate past a threshold before
  it significantly resists collapse — i.e., structure must be already
  well-formed before memory starts screening gravity.  This is a more
  physically reasonable picture than linear coupling from the start.
"""

import numpy as np
from scipy.integrate import solve_ivp
from scipy.interpolate import interp1d
from scipy.optimize import minimize_scalar
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.colors import LogNorm
import json, os

# ─── Cosmological parameters ─────────────────────────────────────────────────
H0      = 67.4
Omega_m = 0.315
Omega_r = 9.4e-5
Omega_L = 0.685
rho_m0  = 3.0 * Omega_m
rho_L   = 3.0 * Omega_L

SIGMA8_PLANCK = 0.811
S8_PLANCK     = SIGMA8_PLANCK * np.sqrt(Omega_m / 0.3)
S8_LENS       = 0.759     # KiDS-1000
S8_SIGMA      = 0.024

# ─── Memory-field parameters ─────────────────────────────────────────────────
alpha = 1.0; beta = 1.5; xi = 10.0; sigma0 = 0.1
a_eq  = 3.0e-4; a_nl = 0.3; n_sh = 2

# ─── RSD compilation ─────────────────────────────────────────────────────────
RSD_Z   = np.array([0.02, 0.067, 0.10, 0.17, 0.22, 0.25, 0.30,
                    0.37, 0.41,  0.57, 0.60, 0.77, 0.80, 1.40])
RSD_FS8 = np.array([0.428, 0.423, 0.370, 0.510, 0.420, 0.351, 0.407,
                    0.460, 0.450, 0.427, 0.480, 0.490, 0.470, 0.482])
RSD_ERR = np.array([0.047, 0.055, 0.130, 0.060, 0.070, 0.058, 0.055,
                    0.038, 0.040, 0.023, 0.100, 0.080, 0.080, 0.116])
RSD_A   = 1.0 / (1.0 + RSD_Z)
N_RSD   = len(RSD_Z)

# ─── Scan grid ───────────────────────────────────────────────────────────────
p_values  = [1, 2, 3, 5, 8, 12, 20]
nu_values = np.linspace(0.0, 1.5, 61)    # fine ν grid for each p

# ─── Integration domain ───────────────────────────────────────────────────────
a_start = 1e-3
a_end   = 1.0
a_arr   = np.logspace(np.log10(a_start), np.log10(a_end), 800)

# ─── Background ──────────────────────────────────────────────────────────────

def E(a):
    return np.sqrt(Omega_m*a**-3 + Omega_r*a**-4 + Omega_L)

def dE_da(a):
    return (-1.5*Omega_m*a**-4 - 2.0*Omega_r*a**-5) / E(a)

def C_func(a):
    rm  = rho_m0 * a**-3
    return rm * sigma0 * (a/a_eq)**n_sh * np.exp(-a/a_nl) / (rm + rho_L/xi)

# ─── Step 1: Solve M(a) once ─────────────────────────────────────────────────
sol_M  = solve_ivp(lambda a, y: [(alpha*C_func(a)/(a*E(a))) - (beta*y[0]/a)],
                   [a_start, a_end], [0.0], t_eval=a_arr,
                   method='RK45', rtol=1e-10, atol=1e-13)
M_arr      = sol_M.y[0]
M_at_1     = M_arr[-1]
Mhat_arr   = M_arr / M_at_1
Mhat_interp = interp1d(a_arr, Mhat_arr, kind='cubic', fill_value='extrapolate')

# ─── Step 2: Growth solver ───────────────────────────────────────────────────

def solve_growth(nu, p):
    """Solve growth ODE with G_eff/G = 1 − ν·M̂^p."""
    def ode(a, y):
        D, Dp  = y
        Geff   = max(0.01, 1.0 - nu * Mhat_interp(a)**p)
        coeff1 = 3.0/a + dE_da(a)/E(a)
        source = (1.5 * Omega_m / (a**5 * E(a)**2)) * Geff
        return [Dp, source*D - coeff1*Dp]

    sol = solve_ivp(ode, [a_start, a_end], [a_start, 1.0],
                    t_eval=a_arr, method='RK45', rtol=1e-9, atol=1e-12)
    return sol.y[0], sol.y[1]   # D_arr, Dp_arr

# ΛCDM reference
D_lcdm, _ = solve_growth(0.0, 1)
D_lcdm_1  = D_lcdm[-1]

D_interp_lcdm = interp1d(a_arr, D_lcdm, kind='cubic')

def diagnostics(nu, p, D_arr, Dp_arr):
    """Return sigma8, S8, chi2_rsd, fs8_model."""
    D_today  = D_arr[-1]
    D_ratio  = D_today / D_lcdm_1
    sigma8   = SIGMA8_PLANCK * D_ratio
    S8       = sigma8 * np.sqrt(Omega_m / 0.3)

    D_int  = interp1d(a_arr, D_arr,  kind='cubic', fill_value='extrapolate')
    Dp_int = interp1d(a_arr, Dp_arr, kind='cubic', fill_value='extrapolate')
    fs8_model = np.array([(RSD_A[i]/D_int(RSD_A[i]))*Dp_int(RSD_A[i])*sigma8
                          for i in range(N_RSD)])
    chi2 = float(np.sum(((fs8_model - RSD_FS8)/RSD_ERR)**2))
    return sigma8, S8, chi2, fs8_model

# ─── Main scan ───────────────────────────────────────────────────────────────
print("=" * 65)
print("  SIM128 — Memory Growth Profile G_eff/G = 1 − ν·M̂(a)^p")
print("=" * 65)
print(f"  p values: {p_values}")
print(f"  ν grid: 0 → {nu_values[-1]:.1f}  ({len(nu_values)} points)\n")

# Full 2D grid results
grid_S8   = np.full((len(p_values), len(nu_values)), np.nan)
grid_chi2 = np.full((len(p_values), len(nu_values)), np.nan)

for i, p in enumerate(p_values):
    for j, nu in enumerate(nu_values):
        D_arr, Dp_arr = solve_growth(nu, p)
        _, S8, chi2, _ = diagnostics(nu, p, D_arr, Dp_arr)
        grid_S8[i, j]   = S8
        grid_chi2[i, j] = chi2

# ─── Per-p optimal ν (joint minimisation) ────────────────────────────────────
print(f"{'p':>4}  {'ν_opt':>7}  {'σ₈':>7}  {'S₈':>7}  {'ΔS₈/σ':>8}  "
      f"{'χ²_RSD':>9}  {'χ²/N':>6}  verdict")
print("-" * 72)

best_per_p = []
for i, p in enumerate(p_values):
    S8_row   = grid_S8[i]
    chi2_row = grid_chi2[i]

    # Joint cost: penalise S₈ pull and RSD χ²/N equally
    pull_row = (S8_row - S8_LENS) / S8_SIGMA
    cost_row = pull_row**2 + chi2_row / N_RSD

    # Exclude ν where G_eff goes negative (over-suppression)
    valid = (1.0 - nu_values) > 0.05   # G_eff(1) = 1-ν > 0.05
    cost_row[~valid] = np.inf

    j_opt   = np.nanargmin(cost_row)
    nu_opt  = nu_values[j_opt]
    D_opt, Dp_opt = solve_growth(nu_opt, p)
    sigma8, S8, chi2, fs8_mod = diagnostics(nu_opt, p, D_opt, Dp_opt)
    pull  = (S8 - S8_LENS) / S8_SIGMA
    chi2N = chi2 / N_RSD

    s8_ok  = abs(pull) < 2.0
    rsd_ok = chi2N < 2.0
    verdict = ("PASS"    if (s8_ok and rsd_ok) else
               "PARTIAL" if (s8_ok or  rsd_ok) else "FAIL")

    print(f"  {p:2d}  {nu_opt:7.3f}  {sigma8:7.4f}  {S8:7.4f}  {pull:+8.2f}σ  "
          f"{chi2:9.3f}  {chi2N:6.2f}  {verdict}")

    best_per_p.append({
        "p": p, "nu_opt": float(nu_opt), "sigma8": float(sigma8),
        "S8": float(S8), "dS8_sigma": float(pull),
        "chi2_rsd": float(chi2), "chi2N": float(chi2N),
        "verdict": verdict,
        "D_arr": D_opt.tolist(), "Dp_arr": Dp_opt.tolist(),
        "fs8_model": fs8_mod.tolist(),
    })

# ─── Find overall best ───────────────────────────────────────────────────────
pass_list    = [r for r in best_per_p if r['verdict'] == 'PASS']
partial_list = [r for r in best_per_p if r['verdict'] == 'PARTIAL']
best = (min(pass_list,    key=lambda r: r['chi2N']) if pass_list else
        min(partial_list, key=lambda r: abs(r['dS8_sigma']) + r['chi2N']/10)
        if partial_list else min(best_per_p, key=lambda r: r['chi2N']))

print(f"\n  Best overall: p={best['p']}, ν={best['nu_opt']:.3f}  →  "
      f"S₈={best['S8']:.4f} ({best['dS8_sigma']:+.2f}σ), "
      f"χ²/N={best['chi2N']:.2f}  [{best['verdict']}]")

# ─── Plotting ─────────────────────────────────────────────────────────────────
colors_p = plt.cm.viridis(np.linspace(0.05, 0.92, len(p_values)))

fig = plt.figure(figsize=(16, 11))
fig.suptitle(
    "SIM128 — Memory Field: Power-Law Growth Profile  [CMSTG]\n"
    r"$G_\mathrm{eff}(a)/G = 1 - \nu\,\hat{M}(a)^p$"
    "     (p=1 → SIM127 baseline)",
    fontsize=13, fontweight='bold', y=0.995
)
gs = gridspec.GridSpec(2, 3, hspace=0.44, wspace=0.38)

ax1 = fig.add_subplot(gs[0, 0])   # M̂(a)^p profiles
ax2 = fig.add_subplot(gs[0, 1])   # G_eff(a) at best ν for each p
ax3 = fig.add_subplot(gs[0, 2])   # 2D S₈ map (p vs ν)
ax4 = fig.add_subplot(gs[1, 0])   # f·σ₈(z) best solutions
ax5 = fig.add_subplot(gs[1, 1])   # S₈ and χ²/N vs p at optimal ν
ax6 = fig.add_subplot(gs[1, 2])   # 2D χ²/N map (p vs ν)

# Panel 1: M̂^p profiles
for p, col in zip(p_values, colors_p):
    ax1.plot(a_arr, Mhat_arr**p, color=col, lw=1.6, label=f'p={p}')
ax1.axvline(0.5, color='gray', ls=':', lw=1, alpha=0.7, label='a=0.5')
ax1.set_xscale('log')
ax1.set_xlabel('Scale factor $a$')
ax1.set_ylabel(r'$\hat{M}(a)^p$')
ax1.set_title(r'Panel 1 — Profile $\hat{M}^p$')
ax1.legend(fontsize=7, ncol=2); ax1.grid(True, alpha=0.3, which='both')

# Panel 2: G_eff(a) at optimal ν for each p
for r, col in zip(best_per_p, colors_p):
    Geff = 1.0 - r['nu_opt'] * Mhat_arr**r['p']
    ax2.plot(a_arr, Geff, color=col, lw=1.8,
             label=f"p={r['p']}, ν={r['nu_opt']:.2f}")
ax2.axhline(1.0, color='k', ls='--', lw=0.8, alpha=0.4)
ax2.set_xscale('log')
ax2.set_xlabel('Scale factor $a$')
ax2.set_ylabel('$G_\\mathrm{eff}(a)/G$')
ax2.set_title('Panel 2 — Optimal G_eff Profiles')
ax2.legend(fontsize=7, ncol=2); ax2.grid(True, alpha=0.3, which='both')

# Panel 3: 2D S₈ heat map
im3 = ax3.pcolormesh(nu_values, p_values, grid_S8,
                     cmap='RdYlGn', vmin=0.70, vmax=0.85)
ax3.contour(nu_values, p_values, grid_S8,
            levels=[S8_LENS - S8_SIGMA, S8_LENS, S8_LENS + S8_SIGMA],
            colors=['blue', 'cyan', 'blue'], linewidths=[1, 2, 1])
ax3.set_xlabel('Coupling $\\nu$')
ax3.set_ylabel('Power $p$')
ax3.set_title('Panel 3 — $S_8(p, \\nu)$ map\n(cyan = KiDS-1000 target)')
plt.colorbar(im3, ax=ax3, label='$S_8$')

# Panel 4: f·σ₈(z)
ax4.errorbar(RSD_Z, RSD_FS8, yerr=RSD_ERR, fmt='ko', ms=5, capsize=3,
             label='RSD data', zorder=5)
for r, col in zip(best_per_p, colors_p):
    z_fine  = np.linspace(0.01, 1.5, 200)
    a_fine  = 1.0/(1.0+z_fine)
    D_int   = interp1d(a_arr, r['D_arr'],  kind='cubic')
    Dp_int  = interp1d(a_arr, r['Dp_arr'], kind='cubic')
    fs8_arr = np.array([(a_fine[i]/D_int(a_fine[i]))*Dp_int(a_fine[i])*r['sigma8']
                        for i in range(len(a_fine))])
    lw = 2.4 if r['p'] == best['p'] else 1.2
    ax4.plot(z_fine, fs8_arr, color=col, lw=lw,
             label=f"p={r['p']} (ν={r['nu_opt']:.2f})")
ax4.set_xlabel('Redshift $z$')
ax4.set_ylabel('$f\\sigma_8(z)$')
ax4.set_title('Panel 4 — RSD Growth Rate (optimal ν per p)')
ax4.legend(fontsize=7, ncol=2); ax4.grid(True, alpha=0.3)
ax4.set_xlim(0, 1.6)

# Panel 5: S₈ and χ²/N vs p
p_arr_plot = [r['p']      for r in best_per_p]
S8_opt     = [r['S8']     for r in best_per_p]
chi2N_opt  = [r['chi2N']  for r in best_per_p]

ax5l = ax5
ax5r = ax5.twinx()
ax5l.plot(p_arr_plot, S8_opt, 'o-', color='#d62728', lw=2, ms=7, label='$S_8$')
ax5l.axhspan(S8_LENS-S8_SIGMA, S8_LENS+S8_SIGMA,
             alpha=0.15, color='green', label='KiDS-1000 1σ')
ax5l.axhline(S8_LENS, color='green', ls='--', lw=1)
ax5r.plot(p_arr_plot, chi2N_opt, 's--', color='#1f77b4', lw=2, ms=7,
          label='$\\chi^2/N$')
ax5r.axhline(2.0, color='#1f77b4', ls=':', lw=1, alpha=0.7)
ax5l.set_xlabel('Power $p$')
ax5l.set_ylabel('$S_8$', color='#d62728')
ax5r.set_ylabel('$\\chi^2_{\\rm RSD}/N$', color='#1f77b4')
ax5l.set_title('Panel 5 — Optimal S₈ and χ²/N vs p')
lines1, labs1 = ax5l.get_legend_handles_labels()
lines2, labs2 = ax5r.get_legend_handles_labels()
ax5l.legend(lines1+lines2, labs1+labs2, fontsize=7)
ax5l.grid(True, alpha=0.3)

# Panel 6: 2D χ²/N heat map
im6 = ax6.pcolormesh(nu_values, p_values, grid_chi2/N_RSD,
                     cmap='RdYlGn_r', vmin=1.0, vmax=15.0)
ax6.contour(nu_values, p_values, grid_chi2/N_RSD,
            levels=[2.0, 3.0], colors='white', linewidths=1.5)
ax6.set_xlabel('Coupling $\\nu$')
ax6.set_ylabel('Power $p$')
ax6.set_title('Panel 6 — $\\chi^2_{\\rm RSD}/N$ map\n(white = 2 and 3)')
plt.colorbar(im6, ax=ax6, label='$\\chi^2/N$')

out_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), '../Outputs')
os.makedirs(out_dir, exist_ok=True)
for ext in ['png', 'pdf']:
    fig.savefig(os.path.join(out_dir, f'SIM128_memory_growth_profile.{ext}'),
                dpi=150, bbox_inches='tight')
plt.close()

# ─── Save JSON ────────────────────────────────────────────────────────────────
out = {
    "sim": "SIM128",
    "params": {"alpha": alpha, "beta": beta, "xi": xi, "sigma0": sigma0,
               "sigma8_planck": SIGMA8_PLANCK, "S8_kids": S8_LENS},
    "best_overall": {k: v for k, v in best.items()
                     if k not in ('D_arr', 'Dp_arr', 'fs8_model')},
    "best_per_p":   [{k: v for k, v in r.items()
                      if k not in ('D_arr', 'Dp_arr', 'fs8_model')}
                     for r in best_per_p],
}
with open(os.path.join(out_dir, 'sim128_results.json'), 'w') as f:
    json.dump(out, f, indent=2)

overall = ("PASS"    if pass_list else
           "PARTIAL" if partial_list else "FAIL")
print(f"\nSIM128 RESULT: {overall}")
print(f"Outputs: {out_dir}/")
