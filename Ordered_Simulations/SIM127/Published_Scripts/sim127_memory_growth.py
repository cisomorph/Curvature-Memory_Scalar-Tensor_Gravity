#!/usr/bin/env python3
"""
SIM127 — Cosmic Memory Field M(a): Structure Growth Suppression
CMSTG Framework, Flat FRW

SIM126 showed M back-reaction on H(a) is the wrong channel for DESI tension
(gives w₀ < -1, phantom-like).  Here we test whether M acts as a
memory-induced friction on gravitational collapse, suppressing the growth
factor D(a) and reducing the σ₈ / S₈ tension.

Physical motivation:
  M accumulates the integrated complexity of the matter distribution.
  In dense environments where structure has already formed, high M
  represents "structural memory" — a resistance to further collapse,
  analogous to a memory-induced effective screening of gravity.

Implementation:
  Modify the effective gravitational coupling in the growth equation:
    G_eff(a) / G = 1 − ν · M̂(a)
  where M̂(a) = M(a)/M(1) is M normalised to unity today, and ν is the
  memory-gravity coupling.  At early times M̂≈0 so G_eff≈G; at late
  times G_eff is suppressed by the accumulated memory.

  Modified linear growth equation (sub-horizon, scale-independent):
    d²D/da² + [3/a + (1/H)dH/da] dD/da
             − (3Ω_m H₀²)/(2 a⁵ H²) · G_eff(a) · D = 0

Key observables:
  1. Growth factor D(a), normalised to ΛCDM at a_start
  2. Growth rate f(a) = d ln D / d ln a
  3. f·σ₈(z) compared to RSD compilation
  4. σ₈ suppression factor relative to Planck ΛCDM value
  5. S₈ = σ₈√(Ω_m/0.3) — the weak lensing tension diagnostic

Units: natural (8πG=1); H₀=1 time unit.
"""

import numpy as np
from scipy.integrate import solve_ivp
from scipy.interpolate import interp1d
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

rho_m0  = 3.0 * Omega_m
rho_L   = 3.0 * Omega_L

# Planck 2018 normalisation
SIGMA8_PLANCK = 0.811    # CMB-inferred σ₈
S8_PLANCK     = SIGMA8_PLANCK * np.sqrt(Omega_m / 0.3)   # ≈ 0.832

# KiDS-1000 / DES weak lensing target
S8_LENS  = 0.759         # central value (KiDS-1000)
S8_SIGMA = 0.024         # 1σ uncertainty

# ─── Memory field parameters (SIM125 baseline) ────────────────────────────────
alpha  = 1.0
beta   = 1.5
xi     = 10.0
sigma0 = 0.1
a_eq   = 3.0e-4
a_nl   = 0.3
n_sh   = 2

# ─── Memory-gravity coupling scan ────────────────────────────────────────────
# ν=0: standard ΛCDM growth; ν>0: M suppresses G_eff
nu_values = [0.0, 0.05, 0.10, 0.20, 0.40, 0.70, 1.00]

# ─── RSD f·σ₈ compilation (Nesseris et al. 2017 + updates) ──────────────────
RSD_Z    = np.array([0.02, 0.067, 0.10, 0.17, 0.22, 0.25, 0.30,
                     0.37, 0.41,  0.57, 0.60, 0.77, 0.80, 1.40])
RSD_FS8  = np.array([0.428, 0.423, 0.370, 0.510, 0.420, 0.351, 0.407,
                     0.460, 0.450, 0.427, 0.480, 0.490, 0.470, 0.482])
RSD_ERR  = np.array([0.047, 0.055, 0.130, 0.060, 0.070, 0.058, 0.055,
                     0.038, 0.040, 0.023, 0.100, 0.080, 0.080, 0.116])
RSD_A    = 1.0 / (1.0 + RSD_Z)

# ─── Integration domain ───────────────────────────────────────────────────────
a_start = 1e-3     # deep matter era — growth ∝ a there, good IC
a_end   = 1.0
N_pts   = 1000
a_arr   = np.logspace(np.log10(a_start), np.log10(a_end), N_pts)

# ─── Background functions (no M back-reaction — SIM125 test field) ───────────

def E(a):
    """E(a) = H(a)/H₀, standard flat ΛCDM."""
    return np.sqrt(Omega_m * a**-3 + Omega_r * a**-4 + Omega_L)

def dE_da(a):
    """dE/da analytically."""
    E_a = E(a)
    return (-1.5*Omega_m*a**-4 - 2.0*Omega_r*a**-5) / E_a

def shear(a):
    return (a / a_eq)**n_sh * np.exp(-a / a_nl)

def C_func(a):
    rm  = rho_m0 * a**-3
    reg = rho_L / xi
    return rm * sigma0 * shear(a) / (rm + reg)

# ─── Step 1: Solve M(a) once (no back-reaction) ──────────────────────────────

def dM_da_ode(a, y):
    M  = y[0]
    Ea = E(a)
    Ca = C_func(a)
    return [(alpha * Ca / (a * Ea)) - (beta * M / a)]

sol_M = solve_ivp(dM_da_ode, [a_start, a_end], [0.0],
                  t_eval=a_arr, method='RK45', rtol=1e-10, atol=1e-13)
M_arr  = sol_M.y[0]
M_at_1 = M_arr[-1]    # M(a=1) for normalisation

# Normalised memory: M̂(a) = M(a)/M(1)  → 0 at early times, 1 today
Mhat_arr = M_arr / M_at_1
Mhat_interp = interp1d(a_arr, Mhat_arr, kind='cubic', fill_value='extrapolate')

# ─── Step 2: Solve growth equation for each ν ────────────────────────────────

def growth_ode(a, y, nu):
    """
    State: y = [D, dD/da]
    Modified growth equation with G_eff(a)/G = 1 - ν·M̂(a).
    """
    D, Dp = y
    Ea    = E(a)
    dEda  = dE_da(a)
    Geff  = max(0.01, 1.0 - nu * Mhat_interp(a))  # floor at 0.01

    # Coefficient of dD/da from H'/H term
    coeff1 = 3.0/a + dEda/Ea

    # Source term: (3/2) Ω_m H₀² / (a⁵ H²) · G_eff
    source = (1.5 * Omega_m / (a**5 * Ea**2)) * Geff

    dDdp_da = source * D - coeff1 * Dp
    return [Dp, dDdp_da]

# Initial conditions: deep matter era, D ∝ a, dD/da = 1 (normalised)
# Use ΛCDM growth at a_start as reference
D0  = a_start
Dp0 = 1.0

# Solve for each ν and store
growth_solutions = {}
print("=" * 60)
print("  SIM127 — Memory Field: Structure Growth Suppression")
print("=" * 60)
print(f"  M(a=1) = {M_at_1:.4e}  (σ₀={sigma0}, SIM125 baseline)")
print(f"  Planck σ₈ = {SIGMA8_PLANCK},  S₈ = {S8_PLANCK:.3f}")
print(f"  KiDS-1000 S₈ target = {S8_LENS} ± {S8_SIGMA}\n")

for nu in nu_values:
    sol = solve_ivp(growth_ode, [a_start, a_end], [D0, Dp0],
                    t_eval=a_arr, args=(nu,),
                    method='RK45', rtol=1e-10, atol=1e-13)
    growth_solutions[nu] = sol

# ΛCDM reference: D_ΛCDM(1)
D_LCDM_at1 = growth_solutions[0.0].y[0, -1]

# ─── Diagnostics ─────────────────────────────────────────────────────────────

print(f"{'ν':>6}  {'D(1)/D_ΛCDM':>12}  {'σ₈':>7}  {'S₈':>7}  "
      f"{'ΔS₈/σ':>8}  {'χ²_RSD':>9}  verdict")
print("-" * 75)

results = []
for nu in nu_values:
    sol     = growth_solutions[nu]
    D_arr   = sol.y[0]
    Dp_arr  = sol.y[1]
    D_today = D_arr[-1]

    # Growth suppression relative to ΛCDM
    D_ratio = D_today / D_LCDM_at1

    # σ₈ and S₈
    sigma8  = SIGMA8_PLANCK * D_ratio
    S8      = sigma8 * np.sqrt(Omega_m / 0.3)
    dS8_sig = (S8 - S8_LENS) / S8_SIGMA   # pull toward KiDS target

    # f·σ₈(z) at RSD redshifts
    D_interp  = interp1d(a_arr, D_arr,  kind='cubic', fill_value='extrapolate')
    Dp_interp = interp1d(a_arr, Dp_arr, kind='cubic', fill_value='extrapolate')

    fs8_model = np.array([
        (RSD_A[i] / D_interp(RSD_A[i])) * Dp_interp(RSD_A[i]) * sigma8
        for i in range(len(RSD_A))
    ])
    chi2_rsd = float(np.sum(((fs8_model - RSD_FS8) / RSD_ERR)**2))

    # Verdict
    s8_ok  = abs(dS8_sig) < 2.0
    rsd_ok = chi2_rsd < 1.5 * len(RSD_A)
    verdict = ("PASS"    if (s8_ok and rsd_ok) else
               "PARTIAL" if (s8_ok or rsd_ok)  else "FAIL")

    print(f"  {nu:4.2f}  {D_ratio:12.5f}  {sigma8:7.4f}  {S8:7.4f}  "
          f"{dS8_sig:+8.2f}σ  {chi2_rsd:9.3f}  {verdict}")

    results.append({
        "nu": nu, "D_ratio": float(D_ratio),
        "sigma8": float(sigma8), "S8": float(S8),
        "dS8_sigma": float(dS8_sig), "chi2_rsd": float(chi2_rsd),
        "verdict": verdict,
    })

# ─── Find best ν ─────────────────────────────────────────────────────────────
# Best = minimises |S₈ - S₈_lens| while keeping RSD χ² reasonable
best_nu = min(results, key=lambda r: abs(r['dS8_sigma']) + r['chi2_rsd']/100)
pass_nu = [r for r in results if r['verdict'] == 'PASS']

print(f"\n  Best ν = {best_nu['nu']:.2f}: "
      f"σ₈={best_nu['sigma8']:.4f}, S₈={best_nu['S8']:.4f}, "
      f"pull={best_nu['dS8_sigma']:+.2f}σ, χ²_RSD={best_nu['chi2_rsd']:.2f}")
if pass_nu:
    print(f"  PASS solutions: ν ∈ "
          f"[{min(r['nu'] for r in pass_nu):.2f}, "
          f"{max(r['nu'] for r in pass_nu):.2f}]")

# ─── Plotting ─────────────────────────────────────────────────────────────────
colors = plt.cm.plasma(np.linspace(0.05, 0.92, len(nu_values)))
a = a_arr

fig = plt.figure(figsize=(15, 10))
fig.suptitle(
    "SIM127 — Memory Field M(a): Structure Growth Suppression  [CMSTG]\n"
    r"$G_\mathrm{eff}(a)/G = 1 - \nu\,\hat{M}(a)$,  "
    r"$\hat{M} = M(a)/M(1)$",
    fontsize=13, fontweight='bold', y=0.995
)
gs = gridspec.GridSpec(2, 3, hspace=0.44, wspace=0.38)

ax1 = fig.add_subplot(gs[0, 0])   # M̂(a)
ax2 = fig.add_subplot(gs[0, 1])   # G_eff(a) for each ν
ax3 = fig.add_subplot(gs[0, 2])   # D(a)/D_ΛCDM
ax4 = fig.add_subplot(gs[1, 0])   # f·σ₈(z) vs RSD data
ax5 = fig.add_subplot(gs[1, 1])   # σ₈ and S₈ vs ν
ax6 = fig.add_subplot(gs[1, 2])   # RSD χ² vs ν

# Panel 1: M̂(a)
ax1.plot(a, Mhat_arr, 'k-', lw=2)
ax1.axvline(a_nl, color='darkorange', ls=':', lw=1.2, alpha=0.7,
            label=f'$a_{{nl}}={a_nl}$')
ax1.set_xscale('log')
ax1.set_xlabel('Scale factor $a$')
ax1.set_ylabel(r'$\hat{M}(a) = M(a)/M(1)$')
ax1.set_title('Panel 1 — Normalised Memory Field')
ax1.legend(fontsize=8); ax1.grid(True, alpha=0.3, which='both')

# Panel 2: G_eff(a)/G
for nu, col in zip(nu_values, colors):
    Geff = 1.0 - nu * Mhat_arr
    ax2.plot(a, Geff, color=col, lw=1.8, label=f'ν={nu:.2f}')
ax2.axhline(1.0, color='k', ls='--', lw=0.8, alpha=0.5)
ax2.set_xscale('log')
ax2.set_xlabel('Scale factor $a$')
ax2.set_ylabel('$G_\\mathrm{eff}(a)/G$')
ax2.set_title('Panel 2 — Effective Gravitational Coupling')
ax2.legend(fontsize=7.5, ncol=2); ax2.grid(True, alpha=0.3, which='both')

# Panel 3: D(a) / D_ΛCDM(a)
D_lcdm_arr = growth_solutions[0.0].y[0]
for nu, col in zip(nu_values, colors):
    D_nu  = growth_solutions[nu].y[0]
    ratio = D_nu / D_lcdm_arr
    ax3.plot(a, ratio, color=col, lw=1.8, label=f'ν={nu:.2f}')
ax3.axhline(1.0, color='k', ls='--', lw=0.8, alpha=0.5, label='ΛCDM')
ax3.set_xscale('log')
ax3.set_xlabel('Scale factor $a$')
ax3.set_ylabel('$D(a)\\ /\\ D_{\\Lambda\\mathrm{CDM}}(a)$')
ax3.set_title('Panel 3 — Growth Factor Suppression')
ax3.legend(fontsize=7.5, ncol=2); ax3.grid(True, alpha=0.3, which='both')

# Panel 4: f·σ₈(z) vs RSD
ax4.errorbar(RSD_Z, RSD_FS8, yerr=RSD_ERR, fmt='ko', ms=5, capsize=3,
             label='RSD data', zorder=5)
for nu, col, r in zip(nu_values, colors, results):
    sol     = growth_solutions[nu]
    D_int   = interp1d(a_arr, sol.y[0], kind='cubic')
    Dp_int  = interp1d(a_arr, sol.y[1], kind='cubic')
    z_fine  = np.linspace(0.01, 1.5, 200)
    a_fine  = 1.0 / (1.0 + z_fine)
    fs8_arr = np.array([(a_fine[i]/D_int(a_fine[i]))*Dp_int(a_fine[i])*r['sigma8']
                        for i in range(len(a_fine))])
    lw = 2.2 if nu == best_nu['nu'] else 1.2
    ax4.plot(z_fine, fs8_arr, color=col, lw=lw, label=f'ν={nu:.2f}')
ax4.set_xlabel('Redshift $z$')
ax4.set_ylabel('$f\\sigma_8(z)$')
ax4.set_title('Panel 4 — RSD Growth Rate')
ax4.legend(fontsize=7, ncol=2); ax4.grid(True, alpha=0.3)
ax4.set_xlim(0, 1.6)

# Panel 5: σ₈ and S₈ vs ν
nu_arr   = [r['nu']    for r in results]
sig8_arr = [r['sigma8'] for r in results]
S8_arr   = [r['S8']    for r in results]

ax5.plot(nu_arr, sig8_arr, 'o-', color='#1f77b4', lw=2, ms=6, label='$\\sigma_8$')
ax5.plot(nu_arr, S8_arr,   's-', color='#d62728', lw=2, ms=6, label='$S_8$')
ax5.axhline(SIGMA8_PLANCK, color='#1f77b4', ls='--', lw=1,
            label=f'Planck $\\sigma_8={SIGMA8_PLANCK}$')
ax5.axhline(S8_PLANCK, color='#d62728', ls='--', lw=1,
            label=f'Planck $S_8={S8_PLANCK:.3f}$')
ax5.axhspan(S8_LENS - S8_SIGMA, S8_LENS + S8_SIGMA,
            alpha=0.15, color='green', label=f'KiDS-1000 $S_8={S8_LENS}\\pm{S8_SIGMA}$')
ax5.set_xlabel('Memory-gravity coupling $\\nu$')
ax5.set_ylabel('$\\sigma_8$ / $S_8$')
ax5.set_title('Panel 5 — Tension Diagnostics vs ν')
ax5.legend(fontsize=7); ax5.grid(True, alpha=0.3)

# Panel 6: RSD χ² vs ν
chi2_arr = [r['chi2_rsd'] for r in results]
ax6.plot(nu_arr, chi2_arr, 'o-', color='#ff7f0e', lw=2, ms=6)
ax6.axhline(len(RSD_Z), color='gray', ls='--', lw=1,
            label=f'χ²=N={len(RSD_Z)} (reduced=1)')
ax6.axhline(1.5 * len(RSD_Z), color='red', ls=':', lw=1, label='χ²=1.5N (limit)')
ax6.set_xlabel('Memory-gravity coupling $\\nu$')
ax6.set_ylabel('$\\chi^2_{\\rm RSD}$')
ax6.set_title('Panel 6 — RSD Fit Quality vs ν')
ax6.legend(fontsize=8); ax6.grid(True, alpha=0.3)

out_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), '../Outputs')
os.makedirs(out_dir, exist_ok=True)
for ext in ['png', 'pdf']:
    fig.savefig(os.path.join(out_dir, f'SIM127_memory_growth.{ext}'),
                dpi=150, bbox_inches='tight')
plt.close()

# ─── Save JSON ────────────────────────────────────────────────────────────────
out = {
    "sim": "SIM127",
    "params": {"alpha": alpha, "beta": beta, "xi": xi, "sigma0": sigma0,
               "M_at_1": float(M_at_1), "sigma8_planck": SIGMA8_PLANCK,
               "S8_planck": float(S8_PLANCK), "S8_kids": S8_LENS},
    "nu_scan": results,
    "best_nu": best_nu,
}
with open(os.path.join(out_dir, 'sim127_results.json'), 'w') as f:
    json.dump(out, f, indent=2)

# ─── Summary ─────────────────────────────────────────────────────────────────
overall = ("PASS"    if pass_nu else
           "PARTIAL" if best_nu['verdict'] == 'PARTIAL' else "FAIL")

print(f"\nSIM127 RESULT: {overall}")
print(f"  G_eff suppression via M can bring S₈ into KiDS-1000 range "
      f"({'yes' if pass_nu else 'no'}).")
print(f"  Best ν={best_nu['nu']:.2f}: "
      f"S₈={best_nu['S8']:.4f} (KiDS target {S8_LENS}±{S8_SIGMA}), "
      f"χ²_RSD={best_nu['chi2_rsd']:.2f}/{len(RSD_Z)}")
print(f"Outputs: {out_dir}/")
