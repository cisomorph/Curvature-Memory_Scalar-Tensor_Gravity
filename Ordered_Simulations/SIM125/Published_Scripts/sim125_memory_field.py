#!/usr/bin/env python3
"""
SIM125 — Cosmic Memory Field M(a)
CMSTG (Curvature-Memory Scalar-Tensor Gravity) Framework
Flat FRW cosmology: evolve H(a), C(a), and the memory scalar M(a)

Units: natural units with 8πG = 1.
  → ρ_crit,0 = 3 H₀²  (from H² = ρ/3 at a=1)
  → all densities in units of H₀²
  → time in units of 1/H₀, so H → E(a) = H(a)/H₀ (dimensionless)

Physics summary:
  M is a cosmic memory field that accumulates complexity C(a) sourced by
  matter shear, and decays proportionally to the expansion rate H.
  The ODE  dM/dt = α·C − β·M·H  encodes this competition.
  Converting to scale factor:  d/dt = a H d/da  →
    dM/da = α·C/(a·E) − β·M/a   (H₀ cancels throughout)
"""

import numpy as np
from scipy.integrate import solve_ivp
import matplotlib
matplotlib.use('Agg')          # headless-safe backend
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import json, os

# ─── Cosmological parameters (Planck 2018) ──────────────────────────────────
H0      = 67.4          # km/s/Mpc  (sets physical scale; cancels in the ODE)
Omega_m = 0.315
Omega_r = 9.4e-5
Omega_L = 0.685

# ρ_i0 / H₀² = 3 Ω_i  (standard result with 8πG = 1)
rho_m0 = 3.0 * Omega_m   # 0.945  [H₀²]
rho_r0 = 3.0 * Omega_r   # 2.82e-4 [H₀²]
rho_L  = 3.0 * Omega_L   # 2.055  [H₀²]  (constant Λ-term)

# ─── Memory-field parameters ─────────────────────────────────────────────────
alpha = 1.0    # source coupling  (dimensionless)
beta  = 1.5    # Hubble-friction coefficient  (dimensionless)
xi    = 10.0   # regulator that prevents C from diverging in the DE epoch

# ─── Shear model σ(a) ────────────────────────────────────────────────────────
# σ(a) = σ₀ · (a/a_eq)^n · exp(−a/a_nl)
# Physical rationale:
#   (a/a_eq)^n    — shear is suppressed at a << a_eq (radiation domination,
#                   perturbations don't grow); rises through matter era
#   exp(−a/a_nl)  — exponential cutoff once structure goes nonlinear and
#                   virialized halos replace coherent shear flows
a_eq = 3.0e-4   # matter-radiation equality scale factor
a_nl = 0.3      # nonlinear structure formation onset
n_sh = 2        # growth power index

# ─── Integration domain ───────────────────────────────────────────────────────
a_start = 0.01
a_end   = 1.0
a_eval  = np.logspace(np.log10(a_start), np.log10(a_end), 600)

sigma0_values = [0.01, 0.1, 1.0]   # free parameter scan

# ─── Physical functions ───────────────────────────────────────────────────────

def E(a):
    """
    Dimensionless Hubble factor E(a) = H(a)/H₀.
    Standard flat ΛCDM: includes matter, radiation, and Λ.
    """
    return np.sqrt(Omega_m * a**-3 + Omega_r * a**-4 + Omega_L)


def shear(a, sigma0):
    """
    Phenomenological shear proxy σ(a).
    Models the rise of large-scale velocity shear during structure formation
    (grows as a^n in matter era) followed by collapse/virialization cutoff.
    """
    return sigma0 * (a / a_eq)**n_sh * np.exp(-a / a_nl)


def C_func(a, sigma0):
    """
    Complexity source C(a) = ρ_m · σ / (ρ_m + p + ρ_Λ/ξ).
    With p = 0 (pressureless dust) this is ρ_m · σ / (ρ_m + ρ_Λ/ξ).

    The regulator ρ_Λ/ξ:
      - Without it, C → σ when ρ_m dominates (early/matter era, fine).
      - At late times ρ_m → 0, so C → 0 naturally via ρ_m in the numerator.
      - The regulator becomes important when ρ_m ~ ρ_Λ/ξ, providing a smooth
        suppression of C in the DE epoch and preventing numerical instability.
    All densities in units of H₀².
    """
    rm  = rho_m0 * a**-3      # matter density [H₀²]
    reg = rho_L / xi           # Λ-regulator   [H₀²]
    sig = shear(a, sigma0)
    return rm * sig / (rm + reg)


def dM_da(a, M_arr, sigma0):
    """
    ODE RHS for the memory field in scale-factor form.

    From dM/dt = α·C − β·M·H,  substituting d/dt = a·H·d/da:
        dM/da = (α·C − β·M·H) / (a·H)
              = α·C / (a·E)  −  β·M / a        [H₀ = 1 units]

    The two terms have clear physical roles:
      + α·C/(a·E)  : complexity injection, modulated by expansion dilution
      − β·M/a      : Hubble friction; M is erased faster in fast-expansion eras
    """
    M  = M_arr[0]
    Ea = E(a)
    Ca = C_func(a, sigma0)
    return [(alpha * Ca / (a * Ea)) - (beta * M / a)]


# ─── Integration ─────────────────────────────────────────────────────────────
print("=" * 58)
print("  SIM125 — Cosmic Memory Field M(a)  [CMSTG Framework]")
print("=" * 58)
print(f"  α = {alpha},  β = {beta},  ξ = {xi:.0f}")
print(f"  σ(a) model: σ₀·(a/a_eq)^{n_sh}·exp(−a/a_nl),  "
      f"a_eq={a_eq}, a_nl={a_nl}")
print()

M_solutions = {}
for s0 in sigma0_values:
    sol = solve_ivp(
        dM_da,
        [a_start, a_end],
        [0.0],                        # M(a_start) = 0 : field starts at zero
        args=(s0,),
        method='RK45',
        t_eval=a_eval,
        rtol=1e-9,
        atol=1e-12,
    )
    if not sol.success:
        print(f"  WARNING: integrator failed for σ₀={s0}: {sol.message}")
    M_solutions[s0] = sol

# ─── M(a=1) values ───────────────────────────────────────────────────────────
print("M(a=1.0):")
results = {}
for s0 in sigma0_values:
    Mf = M_solutions[s0].y[0, -1]
    print(f"  σ₀ = {s0:5.3f}  →  M(1) = {Mf:.6e}")
    results[f"sigma0_{s0}"] = {"M_final": float(Mf)}

# ─── Monotonicity check ──────────────────────────────────────────────────────
print("\nMonotonicity of M(a):")
for s0 in sigma0_values:
    y = M_solutions[s0].y[0]
    diffs = np.diff(y)
    n_neg = int(np.sum(diffs < 0))
    mono  = n_neg == 0
    tag   = "MONOTONICALLY INCREASING ✓" if mono else f"NOT monotonic  ({n_neg} decreasing steps)"
    print(f"  σ₀ = {s0:5.3f}  →  {tag}")
    results[f"sigma0_{s0}"]["monotonic"] = mono

# ─── Regulator significance check ────────────────────────────────────────────
# Fractional change in C from the regulator = R / (ρ_m + R)
# This exceeds 10% when  ρ_m < 9·R  ⟺  a > (ρ_m0 / (9R))^(1/3)
R_val        = rho_L / xi                      # regulator value [H₀²]
threshold_rm = 9.0 * R_val                     # ρ_m below which effect > 10%
a_reg_cross  = (rho_m0 / threshold_rm)**(1.0 / 3.0)

print("\nRegulator ρ_Λ/ξ significance:")
print(f"  ρ_Λ/ξ = {R_val:.4f} H₀²   (ρ_Λ = {rho_L:.4f}, ξ = {xi:.0f})")
print(f"  Effect > 10% when ρ_m < 9R = {threshold_rm:.4f} H₀²")
if a_reg_cross <= a_end:
    frac = (a_end - a_reg_cross) / (a_end - a_start) * 100.0
    print(f"  Threshold crossed at a = {a_reg_cross:.4f}  "
          f"→ regulator active over last {frac:.1f}% of the a-range")
    print(f"  Physical interpretation: DE-dominated epoch; Λ-floor suppresses C,")
    print(f"  preventing runaway growth as matter dilutes.")
else:
    print(f"  Regulator never significant over a=[{a_start},{a_end}]  (ρ_m always dominates)")

print()

# ─── Plotting ─────────────────────────────────────────────────────────────────
colors = ['#1f77b4', '#ff7f0e', '#2ca02c']   # blue, orange, green
a = a_eval

fig = plt.figure(figsize=(14, 10))
fig.suptitle(
    "SIM125 — Cosmic Memory Field $M(a)$  [CMSTG Framework, Flat FRW]\n"
    r"$dM/da = [\alpha C - \beta M H]\,/\,(aH)$"
    f"     $\\alpha={alpha}$,  $\\beta={beta}$,  $\\xi={int(xi)}$",
    fontsize=13, fontweight='bold', y=0.995
)
gs = gridspec.GridSpec(2, 2, hspace=0.44, wspace=0.36)

# ── Panel 1: H(a)/H₀ ─────────────────────────────────────────────────────────
ax1 = fig.add_subplot(gs[0, 0])
ax1.plot(a, E(a), 'k-', lw=2)
ax1.axvline(a_eq, color='steelblue', ls='--', lw=1.2, alpha=0.8,
            label=f'$a_{{eq}}={a_eq}$  (eq.)')
ax1.axvline(a_nl, color='darkorange', ls=':', lw=1.5, alpha=0.8,
            label=f'$a_{{nl}}={a_nl}$  (nonlinear)')
ax1.set_xscale('log'); ax1.set_yscale('log')
ax1.set_xlabel('Scale factor $a$')
ax1.set_ylabel('$H(a)\,/\,H_0$')
ax1.set_title('Panel 1 — Hubble Parameter')
ax1.legend(fontsize=8); ax1.grid(True, alpha=0.3, which='both')

# ── Panel 2: C(a) ─────────────────────────────────────────────────────────────
ax2 = fig.add_subplot(gs[0, 1])
for s0, col in zip(sigma0_values, colors):
    C_vals = np.array([C_func(ai, s0) for ai in a])
    ax2.plot(a, C_vals, color=col, lw=1.8, label=f'$\\sigma_0={s0}$')
ax2.axvline(a_reg_cross, color='gray', ls='--', lw=1.2, alpha=0.8,
            label=f'reg. threshold  $a={a_reg_cross:.3f}$')
ax2.set_xscale('log'); ax2.set_yscale('log')
ax2.set_xlabel('Scale factor $a$')
ax2.set_ylabel('$C(a)$  [complexity source]')
ax2.set_title('Panel 2 — Complexity Source Term')
ax2.legend(fontsize=8); ax2.grid(True, alpha=0.3, which='both')

# ── Panel 3: M(a) for σ₀=0.1 single curve ───────────────────────────────────
ax3 = fig.add_subplot(gs[1, 0])
sol_mid = M_solutions[0.1]
ax3.plot(sol_mid.t, sol_mid.y[0], color='#ff7f0e', lw=2.2)
ax3.axvline(a_nl, color='darkorange', ls=':', lw=1.5, alpha=0.7,
            label=f'$a_{{nl}}={a_nl}$')
# Annotate the plateau-onset region
M_plateau = sol_mid.y[0, -1]
ax3.axhline(M_plateau, color='gray', ls='--', lw=0.8, alpha=0.6)
ax3.set_xscale('log')
ax3.set_xlabel('Scale factor $a$')
ax3.set_ylabel('$M(a)$  [memory field]')
ax3.set_title(r'Panel 3 — Memory Field ($\sigma_0 = 0.1$)')
ax3.legend(fontsize=8); ax3.grid(True, alpha=0.3, which='both')
ax3.annotate(f"$M(1) = {M_plateau:.4f}$",
             xy=(0.58, 0.12), xycoords='axes fraction', fontsize=9,
             bbox=dict(boxstyle='round,pad=0.3', facecolor='wheat', alpha=0.8))

# ── Panel 4: M(a) comparison across all σ₀ ───────────────────────────────────
ax4 = fig.add_subplot(gs[1, 1])
for s0, col in zip(sigma0_values, colors):
    sol = M_solutions[s0]
    Mf  = sol.y[0, -1]
    ax4.plot(sol.t, sol.y[0], color=col, lw=1.8,
             label=f'$\\sigma_0={s0}$  ($M_{{f}}={Mf:.3e}$)')
ax4.axvline(a_nl, color='gray', ls=':', lw=1, alpha=0.6)
ax4.set_xscale('log')
ax4.set_xlabel('Scale factor $a$')
ax4.set_ylabel('$M(a)$  [memory field]')
ax4.set_title(r'Panel 4 — Memory Field: $\sigma_0$ Comparison')
ax4.legend(fontsize=8.5); ax4.grid(True, alpha=0.3, which='both')

out_dir = os.path.dirname(os.path.abspath(__file__)).replace(
    "Published_Scripts", "Outputs"
)
os.makedirs(out_dir, exist_ok=True)
png_path = os.path.join(out_dir, "SIM125_memory_field.png")
pdf_path = os.path.join(out_dir, "SIM125_memory_field.pdf")

plt.savefig(png_path, dpi=150, bbox_inches='tight')
plt.savefig(pdf_path, bbox_inches='tight')
plt.close()

# ─── Save JSON results ────────────────────────────────────────────────────────
results["params"] = {
    "alpha": alpha, "beta": beta, "xi": xi,
    "a_eq": a_eq, "a_nl": a_nl, "n_shear": n_sh,
    "Omega_m": Omega_m, "Omega_r": Omega_r, "Omega_L": Omega_L,
    "H0_kms_Mpc": H0,
    "regulator_threshold_a": float(a_reg_cross),
}
json_path = os.path.join(out_dir, "sim125_results.json")
with open(json_path, "w") as f:
    json.dump(results, f, indent=2)

print(f"Outputs written to {out_dir}/")
print(f"  → SIM125_memory_field.png")
print(f"  → SIM125_memory_field.pdf")
print(f"  → sim125_results.json")
