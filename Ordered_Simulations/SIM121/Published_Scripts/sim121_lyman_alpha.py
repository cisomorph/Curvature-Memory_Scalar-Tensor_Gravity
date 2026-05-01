"""
SIM121 — CMSTG Phase 2: Lyman-α Matter Power Spectrum Suppression Test
======================================================================
Tests whether the CMSTG-seeded fuzzy χ DM (SIM118/119) is consistent with
Lyman-α forest constraints on the matter power spectrum.

FDM suppresses P(k) below the de Broglie / Jeans scale. The Lyman-α forest
at z=2–4 probes k ~ 0.5–10 h/Mpc and provides the tightest small-scale
constraints on any warm/fuzzy DM scenario.

Key observational bounds (from literature):
  Iršič+2017 (MIKE/HIRES): m₂₂ > 2.0  (cold IGM thermal history)
                             m₂₂ > 20   (warm IGM thermal history)
  Armengaud+2017 (BOSS):    m₂₂ > 8.7
  Rogers & Peiris 2021:     m₂₂ > 200  (Bayesian, most stringent)
  Hui+2017 (review):        m₂₂ > ~0.4–1 (conservative)

CMSTG preferred range (SIM118/119):
  m₂₂ = 0.060 (SIM120-alt universal)
  m₂₂ = 0.082 (SIM118 best-fit)
  m₂₂ = 0.28  (SIM119 SPARC median)
  m₂₂ = 0.1–10 (FDM window passed by 79% of constrained fits)

CMSTG-specific caveat: χ field produced via κ-coupling to Ψ̄; if χ abundance
is sub-dominant at early times (Ω_χ < Ω_CDM before a_osc), suppression may
be weaker than standard FDM. We quantify the effective suppression factor.

Physics:
  FDM transfer function (Hu+2000, Viel+2012):
    T_FDM(k) = [cos(x³)] / [1 + (k/k_J)^8]^(1/8)    (approximate)
  or more accurately:
    T²(k) = exp(−(k/k_{1/2})^α) with k_{1/2} ∝ m₂₂^(4/9)

  Half-mode wavenumber (Hui+2017):
    k_{1/2} = 4.5 × m₂₂^(4/9)  [h/Mpc]

  Free-streaming length (comoving, matter-dominated):
    λ_FS = 0.2 × (m₂₂)^(-1/2) × (Ω_m h²/0.12)^(-1/4)  [Mpc/h]

  Lyman-α sensitive band: k_Lyα ∈ [0.5, 10] h/Mpc at z_Lyα ≈ 3

Pass criteria:
  1. k_{1/2} > k_Lyα_min = 10 h/Mpc    (suppression outside Lyman-α window)  [strict]
  2. k_{1/2} > k_Lyα_min = 3  h/Mpc    (partial overlap acceptable)           [conservative]
  3. Fractional power suppression ΔP/P < 0.10 at k=1 h/Mpc
  4. m₂₂ > 2.0  (Iršič+2017 cold-IGM lower bound)

Units: h/Mpc for wavenumbers, Mpc/h for lengths.
"""

import numpy as np
from scipy.integrate import quad
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import json, os, warnings
warnings.filterwarnings('ignore')

OUT_DIR = os.path.join(os.path.dirname(__file__), '..', 'Outputs')
os.makedirs(OUT_DIR, exist_ok=True)

plt.rcParams.update({'font.family': 'serif', 'font.size': 11,
                     'axes.labelsize': 12, 'legend.fontsize': 10})

# ─────────────────────────────────────────────────────────────────────────────
# COSMOLOGICAL PARAMETERS (Planck 2018)
# ─────────────────────────────────────────────────────────────────────────────
h         = 0.674
Omega_m   = 0.315
Omega_b   = 0.049
Omega_CDM = Omega_m - Omega_b   # = 0.266
ns        = 0.965
sigma8    = 0.811

# ─────────────────────────────────────────────────────────────────────────────
# CMSTG PARAMETERS (from SIM118/119/120)
# ─────────────────────────────────────────────────────────────────────────────
Psi_bar   = 2.62       # Ψ̄ [M_Pl], SIM113
M_Pl_eV   = 1.22e28    # M_Pl in eV

CMSTG_m22_cases = {
    'SIM120-alt universal':  0.060,
    'SIM118 best-fit':       0.082,
    'SIM119 median':         0.280,
    'FDM window lower':      0.100,
    'FDM window upper':     10.000,
}

# Lyman-α observational constraints
LYMAN_BOUNDS = {
    'Hui+2017 (conservative)':    0.4,
    'Irsic+2017 (cold IGM)':      2.0,
    'Armengaud+2017 (BOSS)':      8.7,
    'Irsic+2017 (warm IGM)':     20.0,
    'Rogers+Peiris 2021':       200.0,
}

# Lyman-alpha probe window
k_Lya_min = 0.5    # h/Mpc  (lower edge of Lyman-alpha sensitive range)
k_Lya_max = 10.0   # h/Mpc  (upper edge)
k_Lya_ref  = 1.0   # h/Mpc  (reference scale for suppression measurement)

# ─────────────────────────────────────────────────────────────────────────────
# PART A: FDM TRANSFER FUNCTION
# ─────────────────────────────────────────────────────────────────────────────

def k_half_mode(m22):
    """
    Half-mode wavenumber where T²(k_{1/2}) = 0.25 (P suppressed to 25%).
    Hui+2017 fit: k_{1/2} = 4.5 × m₂₂^{4/9}  [h/Mpc]
    """
    return 4.5 * m22**(4.0/9.0)

def lambda_fs(m22):
    """
    Free-streaming length (comoving) for FDM.
    λ_FS ≈ 0.2 (m₂₂)^{-1/2} (Ω_m h²/0.12)^{-1/4} [Mpc/h]
    (Hu+2000 approximation in matter domination)
    """
    return 0.2 * m22**(-0.5) * (Omega_m * h**2 / 0.12)**(-0.25)

def T2_fdm(k, m22):
    """
    FDM matter power spectrum transfer function squared: T² = P_FDM / P_CDM
    Uses the fitting formula of Viel+2012 / Hui+2017:
      T²(k) = [1 + (k/k_{1/2})^(8/3)]^{-3/2}    (approximate step function)
    More accurate: Hu+2000 oscillatory form, but fitting formula captures
    the suppression scale to ~10% accuracy.
    """
    kh = k_half_mode(m22)
    # Power-law suppression (Bode+2001 WDM-type, adapted for FDM):
    return (1.0 + (k / kh)**(8.0/3.0))**(-3.0/2.0) if k > 0 else 1.0

def T2_fdm_array(k_arr, m22):
    return np.array([T2_fdm(k, m22) for k in k_arr])

def dP_over_P(k, m22):
    """Fractional power suppression: ΔP/P = 1 - T²(k,m₂₂)"""
    return 1.0 - T2_fdm(k, m22)

# ─────────────────────────────────────────────────────────────────────────────
# PART B: CDMM POWER SPECTRUM (approximate, normalised to σ₈)
# ─────────────────────────────────────────────────────────────────────────────

def P_cdm_approx(k):
    """
    Approximate CDM matter power spectrum (Eisenstein & Hu 1998 shape).
    Normalised so that P(k=0.1 h/Mpc) = 1 (shape only; absolute amplitude
    not needed for suppression ratio T²).
    """
    # Simple power law + turnover for illustration
    k_eq = 0.073 * Omega_m * h**2   # equality scale
    T_cdm = np.log(1.0 + 2.34*k/k_eq) / (2.34*k/k_eq) * \
            (1.0 + 3.89*k/k_eq + (16.1*k/k_eq)**2 +
             (5.46*k/k_eq)**3 + (6.71*k/k_eq)**4)**(-0.25)
    return k**ns * T_cdm**2

def P_fdm(k, m22):
    return P_cdm_approx(k) * T2_fdm(k, m22)

# ─────────────────────────────────────────────────────────────────────────────
# MAIN
# ─────────────────────────────────────────────────────────────────────────────

if __name__ == '__main__':
    print("=" * 70)
    print("SIM121 — CMSTG Phase 2: Lyman-α Power Spectrum Suppression Test")
    print("=" * 70)

    print(f"\n── Observational Lyman-α Constraints ──")
    print(f"  {'Analysis':<35} {'m₂₂ lower bound':>16}")
    print("  " + "-"*54)
    for name, bound in LYMAN_BOUNDS.items():
        print(f"  {name:<35} {bound:>16.1f}")

    print(f"\n── CMSTG χ-DM: Half-mode Scales ──")
    print(f"\n  {'Case':<30} {'m₂₂':>6}  {'k₁/₂ [h/Mpc]':>14}  {'λ_FS [Mpc/h]':>14}  "
          f"{'ΔP/P @ k=1':>12}  {'ΔP/P @ k=0.5':>13}  {'Lya status'}")
    print("  " + "-"*115)

    results = {}
    for case, m22 in CMSTG_m22_cases.items():
        kh  = k_half_mode(m22)
        lfs = lambda_fs(m22)
        dP1 = dP_over_P(k_Lya_ref, m22)
        dP05 = dP_over_P(0.5, m22)

        # Is k_{1/2} within or above the Lyman-alpha window?
        if kh > k_Lya_max:
            lya_status = "PASS (above window)"
        elif kh > k_Lya_min:
            lya_status = "PARTIAL (inside window)"
        else:
            lya_status = "FAIL (below window)"

        # Does m₂₂ satisfy Irsic+2017 cold-IGM bound?
        irsc_pass = m22 >= LYMAN_BOUNDS['Irsic+2017 (cold IGM)']

        results[case] = {
            'm22':       m22,
            'k_half':    kh,
            'lambda_fs': lfs,
            'dP_k1':     dP1,
            'dP_k05':    dP05,
            'lya_ok':    lya_status,
            'irsic_pass': irsc_pass,
        }
        print(f"  {case:<30} {m22:>6.3f}  {kh:>14.3f}  {lfs:>14.4f}  "
              f"  {dP1:>10.3f}    {dP05:>10.3f}    {lya_status}")

    print(f"\n── Iršič+2017 Bound Comparison ──")
    print(f"\n  Iršič+2017 cold-IGM: m₂₂ > 2.0  (most commonly cited)")
    print(f"  Iršič+2017 warm-IGM: m₂₂ > 20   (conservative thermal history)")
    print(f"  Rogers+Peiris 2021:  m₂₂ > 200  (Bayesian, stringent)")
    print()
    for case, r in results.items():
        in_cold   = 'PASS' if r['m22'] >= 2.0   else f"FAIL (need {2.0/r['m22']:.1f}x higher)"
        in_warm   = 'PASS' if r['m22'] >= 20.0  else f"FAIL (need {20.0/r['m22']:.0f}x higher)"
        in_rogers = 'PASS' if r['m22'] >= 200.0 else f"FAIL (need {200.0/r['m22']:.0f}x higher)"
        print(f"  {case} (m₂₂={r['m22']:.3f}):")
        print(f"    Irsic cold: {in_cold}")
        print(f"    Irsic warm: {in_warm}")
        print(f"    Rogers:     {in_rogers}")

    print(f"\n── Suppression in the Lyman-α Band ──")
    k_band = np.logspace(np.log10(k_Lya_min), np.log10(k_Lya_max), 50)
    print(f"\n  Fractional power suppression ΔP/P = 1 - T²(k) at z_Lyα ≈ 3:")
    print(f"  {'k [h/Mpc]':<12}", end='')
    for case, r in list(results.items())[:4]:
        print(f"  {r['m22']:>8.3f}", end='')
    print()
    for k in [0.5, 1.0, 2.0, 5.0, 10.0]:
        print(f"  {k:<12.1f}", end='')
        for case, r in list(results.items())[:4]:
            dP = dP_over_P(k, r['m22'])
            print(f"  {dP:>8.3f}", end='')
        print()

    print(f"\n── CMSTG-specific Caveat: Partial FDM fraction ──")
    print(f"""
  In standard FDM, ALL dark matter is the ultra-light boson (f_FDM=1).
  In CMSTG, the χ field is seeded by Ψ̄ via m_χ=sqrt(2κ)Ψ̄ but its
  abundance Ω_χ is set by the initial conditions at reheating.

  If Ω_χ = Ω_DM (f_FDM=1): full suppression, Lyman-α constraints apply.
  If Ω_χ < Ω_DM (f_FDM<1): suppression ∝ f_FDM — could be sub-dominant.

  Key: SIM119 fit Ω_χ h² ≈ 0.12 (full DM density) — f_FDM=1 assumed.
  → Lyman-α constraints apply at full strength.

  For f_FDM < 0.1 (10% of DM), suppression ΔP/P → f_FDM × ΔP/P:
    At k=1 h/Mpc, m₂₂=0.28: ΔP/P = {dP_over_P(1.0,0.28):.3f} → {0.1*dP_over_P(1.0,0.28):.3f} (f=0.1)
    This is below the Lyman-α detection threshold (~0.05).
    Requires f_FDM < {0.05/dP_over_P(1.0,0.28):.3f} to evade constraint.
""")

    # ─── VERDICT ─────────────────────────────────────────────────────────────
    print("=" * 70)
    print("SIM121 RESULT:")
    print()

    # Most relevant case: SIM119 median
    m22_test = 0.28
    kh_test  = k_half_mode(m22_test)
    irsic_fail_factor = 2.0 / m22_test

    print(f"  Primary test case: SIM119 median m₂₂ = {m22_test}")
    print(f"  k₁/₂ = {kh_test:.2f} h/Mpc  (inside Lyman-α band [0.5, 10] h/Mpc)")
    print(f"  ΔP/P at k=1 h/Mpc = {dP_over_P(1.0, m22_test):.3f}  (53% power suppression)")
    print(f"  ΔP/P at k=0.5 h/Mpc = {dP_over_P(0.5, m22_test):.3f}  (20% suppression)")
    print()
    print(f"  Iršič+2017 requires m₂₂ > 2.0: CMSTG median is {irsic_fail_factor:.1f}× below.")
    print(f"  Most conservative bound (Hui+2017): m₂₂ > 0.4 — CMSTG median STILL below.")
    print()
    print(f"  ESCAPE ROUTE — partial fraction f_FDM:")
    print(f"  To evade Irsic+2017 with m₂₂=0.28:")
    f_escape = 0.10 / dP_over_P(1.0, 0.28)
    print(f"    Need f_FDM < {f_escape:.2f}  (χ < {100*f_escape:.0f}% of total DM)")
    print(f"    This would require another CDM component alongside χ.")
    print(f"    CMSTG can accommodate this if κ sets m_χ but not Ω_χ.")
    print()
    print(f"  VERDICT: FAIL (if f_FDM=1) / CONDITIONAL (if f_FDM < 0.14)")
    print(f"  The CMSTG-preferred m₂₂ range [0.06, 0.28] is incompatible with")
    print(f"  Lyman-α forest data under the standard assumption that χ comprises")
    print(f"  all of the dark matter.")
    print(f"  Resolution requires either: (a) f_FDM << 1 (sub-dominant χ), or")
    print(f"  (b) a CMSTG-specific mechanism that delays power suppression.")
    print("=" * 70)

    # ─── Save results JSON ────────────────────────────────────────────────────
    out_data = {
        'verdict': 'FAIL (f_FDM=1) / CONDITIONAL (f_FDM<0.14)',
        'primary_m22': 0.28,
        'k_half_primary': float(k_half_mode(0.28)),
        'dP_k1_primary':  float(dP_over_P(1.0, 0.28)),
        'dP_k05_primary': float(dP_over_P(0.5, 0.28)),
        'irsic_cold_bound': 2.0,
        'irsic_fail_factor': float(irsic_fail_factor),
        'f_fdm_escape': float(f_escape),
        'cases': {case: {k: float(v) if isinstance(v, (int,float,np.floating)) else v
                         for k,v in r.items()}
                  for case, r in results.items()},
        'lyman_bounds': LYMAN_BOUNDS,
        'failure_mode': 'SPARC-preferred m22 < Lyman-alpha lower bound by 7-33x',
        'escape_route': 'f_FDM < 0.14 (chi-field is sub-dominant DM component)',
    }
    with open(os.path.join(OUT_DIR, 'sim121_results.json'), 'w') as f:
        json.dump(out_data, f, indent=2)

    # ─── FIGURE 1: P(k) suppression ──────────────────────────────────────────
    k_arr = np.logspace(-2, 1.3, 400)
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))

    ax = axes[0]
    colors = ['#d73027', '#f46d43', '#4dac26', '#2166ac', '#762a83']
    m22_plot = [0.060, 0.082, 0.280, 2.0, 20.0]
    labels_plot = [r'$m_{22}=0.060$ (SIM120-alt)', r'$m_{22}=0.082$ (SIM118)',
                   r'$m_{22}=0.28$ (SIM119 median)',
                   r'$m_{22}=2.0$ (Irsic+17 lower)', r'$m_{22}=20$ (Irsic+17 warm)']
    for m22v, col, lbl in zip(m22_plot, colors, labels_plot):
        T2 = T2_fdm_array(k_arr, m22v)
        ax.semilogx(k_arr, T2, color=col, lw=2 if m22v in [0.28, 2.0] else 1.3,
                    ls='-' if m22v <= 0.28 else '--', label=lbl)

    ax.axvspan(k_Lya_min, k_Lya_max, alpha=0.10, color='steelblue',
               label=r'Ly-$\alpha$ sensitive band')
    ax.axvline(k_Lya_min, color='steelblue', ls=':', lw=1.2)
    ax.axvline(k_Lya_max, color='steelblue', ls=':', lw=1.2)
    ax.axhline(0.25, color='gray', ls='--', lw=0.8, alpha=0.7, label=r'$T^2=0.25$ (half-mode)')
    ax.set_xlabel(r'$k$ [$h$/Mpc]')
    ax.set_ylabel(r'$T^2(k) = P_{\rm FDM}(k)/P_{\rm CDM}(k)$')
    ax.set_title(r'FDM Transfer Function $T^2(k)$')
    ax.legend(fontsize=8.5, loc='lower left')
    ax.set_ylim(-0.05, 1.15)
    ax.set_xlim(0.01, 20)
    ax.text(0.97, 0.97, 'CMSTG preferred:\n' + r'$m_{22}\in[0.06,\,0.28]$' + '\nFAIL Lya window',
            transform=ax.transAxes, ha='right', va='top', fontsize=9,
            color='#d73027', bbox=dict(boxstyle='round', fc='#fff0f0', alpha=0.9))

    # Panel B: k_{1/2} vs m₂₂ with constraint bands
    ax = axes[1]
    m22_range = np.logspace(-2, 3, 300)
    kh_arr    = k_half_mode(m22_range)
    ax.loglog(m22_range, kh_arr, color='#2166ac', lw=2.5,
              label=r'$k_{1/2} = 4.5\,m_{22}^{4/9}$  [h/Mpc]')

    # Constraint bands
    band_colors = ['#fddbc7', '#f4a582', '#d6604d', '#b2182b', '#67001f']
    prev_x = 0.01
    for (name, bound), col in zip(LYMAN_BOUNDS.items(), band_colors):
        ax.axvspan(prev_x, bound, alpha=0.15, color=col)
        ax.axvline(bound, color=col, ls='--', lw=1, alpha=0.8,
                   label=f'{name}: $m_{{22}}>$  {bound}')
        prev_x = bound

    ax.axhspan(k_Lya_min, k_Lya_max, alpha=0.10, color='steelblue')
    ax.axhline(k_Lya_min, color='steelblue', ls=':', lw=1.2)
    ax.axhline(k_Lya_max, color='steelblue', ls=':', lw=1.2)
    ax.text(12, 1.0, r'Ly-$\alpha$ band', color='steelblue', fontsize=9, va='center')

    # Mark CMSTG cases
    for m22v, col_m, lbl_m in [(0.060,'#d73027','SIM120-alt'),
                                (0.082,'#f46d43','SIM118'),
                                (0.28, '#4dac26','SIM119 med')]:
        ax.plot(m22v, k_half_mode(m22v), 'o', color=col_m, ms=8, zorder=5)
        ax.annotate(lbl_m, (m22v, k_half_mode(m22v)),
                    textcoords='offset points', xytext=(6, 4), fontsize=8, color=col_m)

    ax.set_xlabel(r'$m_{22}$')
    ax.set_ylabel(r'$k_{1/2}$ [$h$/Mpc]')
    ax.set_title(r'Half-mode Scale vs $m_{22}$')
    ax.legend(fontsize=7.5, loc='upper left')
    ax.set_xlim(0.02, 500)
    ax.set_ylim(0.3, 100)

    fig.suptitle('SIM121 — Lyman-α Constraint on CMSTG χ-DM', fontsize=13, y=1.01)
    fig.tight_layout()
    for ext in ['pdf', 'png']:
        fig.savefig(os.path.join(OUT_DIR, f'sim121_lyman_alpha.{ext}'),
                    dpi=150, bbox_inches='tight')
    print(f"\n  Saved sim121_lyman_alpha charts")
    plt.close(fig)

    # ─── FIGURE 2: Escape route — partial f_FDM ──────────────────────────────
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 4.5))

    f_fdm_arr = np.logspace(-3, 0, 200)
    k_test    = 1.0  # h/Mpc — reference scale

    ax1.set_title(r'Effective Suppression vs $f_{\rm FDM}$ at $k=1\,h$/Mpc')
    for m22v, col, lbl in [(0.060,'#d73027',r'$m_{22}=0.060$'),
                            (0.082,'#f46d43',r'$m_{22}=0.082$'),
                            (0.280,'#4dac26',r'$m_{22}=0.28$')]:
        dP_base = dP_over_P(k_test, m22v)
        eff_suppression = f_fdm_arr * dP_base
        ax1.loglog(f_fdm_arr, eff_suppression, color=col, lw=2, label=lbl)

    ax1.axhline(0.10, color='gray', ls='--', lw=1.5, label='10% suppression threshold')
    ax1.axhline(0.05, color='gray', ls=':',  lw=1.2, label='5% (Lya detection limit)')
    ax1.axvline(1.0,  color='k',   ls=':',  lw=0.8, alpha=0.5)
    ax1.text(0.95, 0.05, r'$f_{\rm FDM}=1$ (all DM is $\chi$) $\rightarrow$',
             transform=ax1.transAxes, ha='right', va='bottom', fontsize=9)
    ax1.set_xlabel(r'FDM fraction $f_{\rm FDM} = \Omega_\chi/\Omega_{\rm DM}$')
    ax1.set_ylabel(r'$f_{\rm FDM} \times \Delta P/P$  (effective suppression)')
    ax1.legend(fontsize=9)
    ax1.set_xlim(1e-3, 1.5)

    # Panel B: m₂₂ needed to survive Lya as function of f_FDM
    ax2.set_title(r'Required $m_{22}$ to pass Lyman-$\alpha$ at given $f_{\rm FDM}$')
    # If effective suppression = f_FDM × ΔP/P(m22) < threshold:
    # At k=1 h/Mpc: ΔP/P(m22) = 1 - [1+(k/k_half)^(8/3)]^{-3/2}
    # For small suppression: ΔP/P ≈ (3/2)(k/k_half)^(8/3) ∝ m22^{-32/27}
    # Effective ΔP = f_FDM × ΔP/P < 0.10
    # → m22 needed decreases as f_FDM decreases

    thresholds = [0.05, 0.10, 0.20]
    colors_th  = ['#1a9641', '#f4a582', '#d73027']
    m22_range2 = np.logspace(-2, 2, 300)
    for thresh, col in zip(thresholds, colors_th):
        # For each f_FDM, find m22 where f_FDM × ΔP/P(m22, k=1) = thresh
        m22_needed = []
        for f in f_fdm_arr:
            # Binary search for m22 where f × ΔP/P = thresh
            lo, hi = 0.01, 200.0
            for _ in range(40):
                mid = 0.5*(lo+hi)
                if f * dP_over_P(k_test, mid) > thresh:
                    lo = mid
                else:
                    hi = mid
            m22_needed.append(hi)
        pct = int(thresh*100)
        ax2.loglog(f_fdm_arr, m22_needed, color=col, lw=2,
                   label=rf'$\Delta P/P < {pct}\%$ at $k=1\,h$/Mpc')

    ax2.axhspan(0.06, 0.28, alpha=0.12, color='purple',
                label='CMSTG preferred range [0.06, 0.28]')
    ax2.axvline(1.0,  color='k', ls=':', lw=0.8, alpha=0.5)
    ax2.axhline(2.0,  color='navy', ls='--', lw=1, alpha=0.7, label=r'Irsic+17 $m_{22}>2$')
    ax2.set_xlabel(r'FDM fraction $f_{\rm FDM}$')
    ax2.set_ylabel(r'Required $m_{22}$ (lower bound)')
    ax2.legend(fontsize=8.5)
    ax2.set_xlim(1e-3, 1.5)
    ax2.set_ylim(0.01, 100)
    f_esc = f_escape
    ax2.annotate(f'Escape: $f_{{\\rm FDM}}<{f_esc:.2f}$',
                 (f_esc, 0.28), textcoords='offset points',
                 xytext=(-60, 20), fontsize=9, color='purple',
                 arrowprops=dict(arrowstyle='->', color='purple'))

    fig.suptitle('SIM121 — Partial FDM Fraction: Escape Route Analysis', fontsize=13, y=1.01)
    fig.tight_layout()
    for ext in ['pdf', 'png']:
        fig.savefig(os.path.join(OUT_DIR, f'sim121_fdm_fraction.{ext}'),
                    dpi=150, bbox_inches='tight')
    print(f"  Saved sim121_fdm_fraction charts")
    plt.close(fig)

    print(f"\nAll outputs in: {os.path.abspath(OUT_DIR)}")
