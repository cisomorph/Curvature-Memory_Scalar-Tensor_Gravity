#!/usr/bin/env python3
"""
SIM115: Phase 2 — Gradient Soliton DM (most-RIFT approach)
RIFT Phase 2: Recursive Intelligence-Field Theory (Unlocked Action)

SIM114 found a trilemma for the condensate mechanism:
  - β > β_thresh = 8πΛ₀ = 0.075: condensate forms but G_eff → 0
  - β < β_thresh: Λ₀R curvature term suppresses Ψ → 0; no condensate
  - β_needed for DM (~10⁻⁴) << β_thresh: gap of 3 orders of magnitude

Most-RIFT alternative (this sim): the Λ₀R suppression in dense regions IS the DM mechanism.
In galaxies, curvature drives Ψ BELOW the cosmological background (Ψ_cosmo = 2.62).
This creates a galactic soliton profile: Ψ dips from Ψ_cosmo at large r to Ψ_center at r=0.
The field energy — gradient + potential — constitutes the DM.

This is RIFT-natural because:
  1. The Λ₀ΨR coupling (RIFT's defining term) creates the dip — geometry sources field
  2. The recursion: baryon curvature → Ψ suppression → gradient energy → DM density
  3. G_eff is AUTOMATICALLY safe: Ψ_c ≤ Ψ_cosmo < Ψ_max throughout the galaxy

The field equation (static, spherically symmetric):
  Ψ'' + (2/r)Ψ' = 4λΨ(Ψ²−v²) − 2β·G̃·ρ_b·Ψ + 2Λ₀·Ψ·R
  G̃ ≡ G_gal/c²   [converts M_sun/kpc³ to kpc⁻²]
  R = −8π·G̃·ρ_b  [weak-field Ricci, Ψ contribution negligible for small soliton]

Boundary conditions (gradient soliton regime, β < β_thresh):
  Ψ(r_max) = Ψ_cosmo = 2.618   [cosmological background]
  Ψ'(r_max) = 0                 [field approaches background from outside]
Integrate INWARD → Ψ_center = Ψ(r_min) as output.

DM field energy (excess above cosmological background):
  ρ_DM(r) = (c²/G_gal) × [½(Ψ')² + V_J(Ψ) − V_J(Ψ_cosmo)]
  V_J(Ψ) = λ_kpc × (Ψ² − v²)²

KEY QUESTION: Does the gradient soliton provide enough DM density for flat rotation curves?
  If ρ_DM is sufficient, PASS. If not, quantify the gap.

Expected (from dimensional analysis):
  Source: 2Λ₀·|R|·Ψ_cosmo ~ 1.3×10⁻¹⁰ kpc⁻² at r_d
  δΨ ~ source × r_d² ~ 2×10⁻⁹  (field barely perturbed)
  Ψ' ~ δΨ/r_d ~ 5×10⁻¹⁰ kpc⁻¹
  ρ_grad ~ (c²/G_gal) × (Ψ')²/2 ~ 10⁻³ M_sun/kpc³
  ρ_DM_needed ~ 10⁷ M_sun/kpc³  →  gap ~ 10¹⁰

Root cause: λ_kpc = λ_cosmo × H₀² ~ 3.6×10⁻¹⁹ kpc⁻²  (H₀ suppression squares down)
           (c²/G_gal) × λ_kpc ~ 7×10⁻³  (galactic energy scale / DM energy scale ≪ 1)
"""

import os, json
import numpy as np
from scipy.integrate import solve_ivp, cumulative_trapezoid
from scipy.special import iv, kv
import warnings
warnings.filterwarnings('ignore')

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

BASE    = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
INPUTS  = os.path.join(BASE, 'Inputs')
OUTPUTS = os.path.join(BASE, 'Outputs')
os.makedirs(OUTPUTS, exist_ok=True)

with open(os.path.join(INPUTS, 'sim115_params.json')) as f:
    P = json.load(f)

# ── Physical constants ────────────────────────────────────────────────────
G_gal   = 4.3009e-6    # (km/s)² kpc M_sun⁻¹
c_kms   = 2.998e5      # km/s
G_tilde = G_gal / c_kms**2      # kpc M_sun⁻¹  (ρ unit conversion)
H0_kpc  = 67.4 / (299800.0 * 3086.0)  # = 7.285e-8 kpc⁻¹

# ── Phase 2 parameters ────────────────────────────────────────────────────
P2        = P['phase2_params']
LAMBDA0   = P2['Lambda0']              # 0.003
LAM_CO    = P2['lambda_cosmo']         # 6.856e-5 (cosmological units)
LAM_KPC   = LAM_CO * H0_kpc**2        # 3.64e-19 kpc⁻²
V_SSB     = P2['v_ssb']               # 13.158 (VEV, dimensionless)
PSI_COS   = P2['Psi_cosmo']           # 2.618 (SIM113 best fit)
BETA_THR  = P2['beta_threshold']       # 8π × 0.003 = 0.07540
BETA_SCAN = P2['beta_scan']

# G_eff cosmological constraint (5% → consistent with SIM113 Psi0=2.62, 3.95% dev)
GEF_BOUND  = P['cosmological_constraint']['Geff_max_deviation_pct'] / 100.0  # 0.05
PSI_MAX    = np.sqrt(GEF_BOUND / (2.0 * LAMBDA0))   # = sqrt(0.05/0.006) = 2.887

NUM       = P['numerics']
FLAT_TOL  = P['verdict']['flatness_tolerance']
FLAT_FRAC = P['verdict']['flatness_fraction']

# ── Galaxy ────────────────────────────────────────────────────────────────
GAL       = list(P['galaxies'].values())[0]
GAL_NAME  = list(P['galaxies'].keys())[0]
V_FLAT    = GAL['v_flat_kms']


# ═══════════════════════════════════════════════════════════════════════════
# BARYONIC MASS MODELS
# ═══════════════════════════════════════════════════════════════════════════

def rho_baryon(r, gal):
    M_d, r_d = gal['M_disk_Msun'], gal['r_d_kpc']
    M_g, r_g = gal['M_gas_Msun'],  gal['r_gas_kpc']
    dMdr = (M_d/r_d**2 * r * np.exp(-r/r_d)
          + M_g/r_g**2 * r * np.exp(-r/r_g))
    return dMdr / (4 * np.pi * r**2)


def M_baryon_enclosed(r_arr, gal):
    M_d, r_d = gal['M_disk_Msun'], gal['r_d_kpc']
    M_g, r_g = gal['M_gas_Msun'],  gal['r_gas_kpc']
    enc_d = M_d * (1 - np.exp(-r_arr/r_d) * (1 + r_arr/r_d))
    enc_g = M_g * (1 - np.exp(-r_arr/r_g) * (1 + r_arr/r_g))
    return enc_d + enc_g


def ricci_scalar(r, gal):
    return -8.0 * np.pi * G_tilde * rho_baryon(r, gal)   # < 0 for matter


# ═══════════════════════════════════════════════════════════════════════════
# PHASE 2 FIELD ODE
# ═══════════════════════════════════════════════════════════════════════════

def psi_ode(r, y, beta, gal):
    """
    Ψ'' = -2/r·Ψ' + 4λ_kpc·Ψ(Ψ²−v²) − 2β·G̃·ρ_b·Ψ + 2Λ₀·Ψ·R
    """
    psi, dpsi = y
    rho_b = rho_baryon(r, gal)
    R     = ricci_scalar(r, gal)   # = -8π·G̃·ρ_b < 0

    V_prime = 4.0 * LAM_KPC * psi * (psi**2 - V_SSB**2)  # dV_J/dPsi
    beta_term = -2.0 * beta * G_tilde * rho_b * psi
    curv_term  = 2.0 * LAMBDA0 * psi * R

    d2psi = (-2.0/r * dpsi + V_prime + beta_term + curv_term)
    return [dpsi, d2psi]


def solve_soliton(r_arr, beta, gal):
    """
    Integrate INWARD from r_max with BC: Ψ = Ψ_cosmo, Ψ' = 0.
    Returns (psi, dpsi) arrays on r_arr (increasing order).
    """
    r_max = r_arr[-1]

    sol = solve_ivp(
        psi_ode,
        (r_max, r_arr[0]),    # integrate inward
        [PSI_COS, 0.0],       # Ψ(r_max) = Ψ_cosmo, Ψ'(r_max) = 0
        t_eval=r_arr[::-1],   # decreasing r
        args=(beta, gal),
        method='DOP853',
        rtol=NUM['rtol'], atol=NUM['atol'],
        dense_output=False
    )

    if not sol.success:
        return None, None

    psi  = sol.y[0][::-1]   # restore increasing-r order
    dpsi = sol.y[1][::-1]

    if np.any(~np.isfinite(psi)) or np.any(np.abs(psi) > 1e6):
        return None, None

    return psi, dpsi


# ═══════════════════════════════════════════════════════════════════════════
# DM FIELD ENERGY (GRADIENT + POTENTIAL EXCESS)
# ═══════════════════════════════════════════════════════════════════════════

def rho_dm_field(psi, dpsi):
    """
    DM energy density = excess field energy above cosmological background.
    ρ_DM = (c²/G_gal) × [½(Ψ')² + V_J(Ψ) − V_J(Ψ_cosmo)]
    For Ψ, Ψ_cosmo << v:  V_J(Ψ) ≈ λ_kpc(v²-Ψ²)²  ≈  λ_kpc·v⁴(1 - Ψ²/v²)²
    """
    gradient = 0.5 * dpsi**2
    V_psi    = LAM_KPC * (psi**2 - V_SSB**2)**2
    V_cosmo  = LAM_KPC * (PSI_COS**2 - V_SSB**2)**2
    excess   = gradient + (V_psi - V_cosmo)
    return np.maximum(excess, 0.0) * (c_kms**2 / G_gal)   # M_sun/kpc³


def rho_dm_gradient_only(dpsi):
    """Gradient (kinetic) contribution only."""
    return 0.5 * dpsi**2 * (c_kms**2 / G_gal)


def M_dm_enclosed(r_arr, psi, dpsi):
    rho = rho_dm_field(psi, dpsi)
    return cumulative_trapezoid(4 * np.pi * r_arr**2 * rho, r_arr, initial=0)


def rho_dm_needed(r_arr, gal):
    """Isothermal sphere DM density for flat v_c = v_flat."""
    return gal['v_flat_kms']**2 / (4.0 * np.pi * G_gal * r_arr**2)


# ═══════════════════════════════════════════════════════════════════════════
# G_EFF AND ROTATION CURVE
# ═══════════════════════════════════════════════════════════════════════════

def G_eff_ratio(psi):
    return 1.0 / (1.0 + 2.0 * LAMBDA0 * psi**2)


def rotation_curve(r_arr, psi, dpsi, gal):
    geff  = G_eff_ratio(psi)
    M_bar = M_baryon_enclosed(r_arr, gal)
    M_dm  = M_dm_enclosed(r_arr, psi, dpsi)
    v2    = geff * G_gal * (M_bar + M_dm) / r_arr
    return np.sqrt(np.clip(v2, 0, None)), geff, M_bar, M_dm


def check_flatness(r_arr, vc):
    r_cut = r_arr[0] + FLAT_FRAC * (r_arr[-1] - r_arr[0])
    mask  = r_arr >= r_cut
    if not np.any(mask):
        return False, 1.0
    dev = np.abs(vc[mask] - V_FLAT) / V_FLAT
    return bool(np.all(dev < FLAT_TOL)), float(np.max(dev))


# ═══════════════════════════════════════════════════════════════════════════
# DIMENSIONAL ANALYSIS: the structural gap
# ═══════════════════════════════════════════════════════════════════════════

def structural_gap_analysis(gal):
    """
    Compute the expected DM density from the gradient soliton analytically,
    and the needed DM density, to quantify the gap.
    """
    r_d = gal['r_d_kpc']
    rho_b_rd = rho_baryon(r_d, gal)

    # Source term at r = r_d (max Λ₀R driving)
    source_rd = 2.0 * LAMBDA0 * abs(ricci_scalar(r_d, gal)) * PSI_COS
    # = 2Λ₀ × 8πG̃ × ρ_b × Ψ_cosmo

    # Characteristic field perturbation δΨ ~ source × r_d²
    delta_psi = source_rd * r_d**2

    # Characteristic gradient: Ψ' ~ δΨ / r_d
    psi_prime = delta_psi / r_d

    # Gradient DM density: (c²/G) × (Ψ')²/2
    rho_grad = 0.5 * psi_prime**2 * c_kms**2 / G_gal

    # Needed DM density at r_d
    rho_needed_rd = rho_dm_needed(np.array([r_d]), gal)[0]

    # Gap ratio
    gap = rho_needed_rd / max(rho_grad, 1e-30)

    return {
        'r_d_kpc': r_d,
        'rho_b_at_rd': rho_b_rd,
        'source_term_kpc_inv2': source_rd,
        'delta_Psi_induced': delta_psi,
        'Psi_prime_kpc_inv': psi_prime,
        'rho_DM_gradient': rho_grad,
        'rho_DM_needed': rho_needed_rd,
        'gap_ratio': gap,
        'log10_gap': np.log10(gap),
        'lambda_kpc': LAM_KPC,
        'c2_over_G_times_lambdakpc': LAM_KPC * c_kms**2 / G_gal
    }


# ═══════════════════════════════════════════════════════════════════════════
# MAIN
# ═══════════════════════════════════════════════════════════════════════════

def run():
    print("=" * 72)
    print("SIM115 — Phase 2 RIFT: Gradient Soliton DM (most-RIFT approach)")
    print("=" * 72)
    print(f"  Λ₀       = {LAMBDA0}")
    print(f"  λ_kpc    = {LAM_KPC:.4e} kpc⁻²   [λ_cosmo × H₀²]")
    print(f"  v        = {V_SSB:.4f}  (VEV)")
    print(f"  Ψ_cosmo  = {PSI_COS:.4f}  (SIM113 best fit, BC at r→∞)")
    print(f"  Ψ_max    = {PSI_MAX:.4f}  (G_eff/G within {GEF_BOUND*100:.0f}%)")
    print(f"  β_thresh = {BETA_THR:.5f}  (= 8πΛ₀; soliton regime: β < β_thresh)")
    print(f"  Ψ_cosmo < Ψ_max: G_eff-safe throughout  ✓")
    print()

    diag = {
        'sim': 'SIM115',
        'description': (
            'Phase 2 gradient soliton DM: Psi dips below Psi_cosmo in dense regions. '
            'BC: Psi(r_max)=Psi_cosmo, integrate inward. '
            'DM = excess gradient + potential energy above cosmological background.'
        ),
        'Lambda0': LAMBDA0,
        'lambda_kpc': LAM_KPC,
        'v_ssb': V_SSB,
        'Psi_cosmo': PSI_COS,
        'Psi_max': PSI_MAX,
        'beta_threshold': BETA_THR,
        'profiles': [],
        'structural_gap': None,
        'overall_result': None
    }

    # Set up radial grid
    Nr = GAL['Nr']
    r_arr = np.linspace(GAL['r_min_kpc'], GAL['r_max_kpc'], Nr)
    rho_b = rho_baryon(r_arr, GAL)

    # ── 1. STRUCTURAL GAP ANALYSIS ────────────────────────────────────────
    print("─" * 60)
    print("PART 1: Dimensional analysis — structural gap")
    print("─" * 60)

    gap_info = structural_gap_analysis(GAL)
    diag['structural_gap'] = gap_info

    print(f"\n  At r_d = {gap_info['r_d_kpc']} kpc:")
    print(f"    ρ_b                   = {gap_info['rho_b_at_rd']:.4e} M_sun/kpc³")
    print(f"    Λ₀R source term       = {gap_info['source_term_kpc_inv2']:.4e} kpc⁻²  [drives Ψ dip]")
    print(f"    Induced δΨ            = {gap_info['delta_Psi_induced']:.4e}  [source × r_d²]")
    print(f"    Gradient Ψ'           = {gap_info['Psi_prime_kpc_inv']:.4e} kpc⁻¹")
    print(f"    ρ_DM (gradient)       = {gap_info['rho_DM_gradient']:.4e} M_sun/kpc³")
    print(f"    ρ_DM (needed)         = {gap_info['rho_DM_needed']:.4e} M_sun/kpc³")
    print(f"    GAP RATIO             = {gap_info['gap_ratio']:.4e}  (10^{gap_info['log10_gap']:.1f})")
    print()
    print(f"  Root cause: λ_kpc = λ_cosmo × H₀² = {LAM_CO:.2e} × ({H0_kpc:.2e})² = {LAM_KPC:.2e} kpc⁻²")
    print(f"  Energy scale: (c²/G_gal) × λ_kpc = {LAM_KPC * c_kms**2/G_gal:.4e} M_sun/kpc³ × kpc²")
    print(f"  → Field energy density per (δΨ/kpc)² is {LAM_KPC * c_kms**2/G_gal:.2e}")
    print(f"  → Required (δΨ/kpc)² ~ {gap_info['rho_DM_needed'] * G_gal / c_kms**2:.2e}  but induced ~ {gap_info['Psi_prime_kpc_inv']**2:.2e}")

    # ── 2. NUMERICAL PROFILES ────────────────────────────────────────────
    print("\n" + "─" * 60)
    print("PART 2: Numerical gradient soliton profiles")
    print("─" * 60)

    rho_needed = rho_dm_needed(r_arr, GAL)
    best_profile = None
    best_gap = np.inf

    print(f"\n{'β':>10} | {'Ψ_center':>12} | {'δΨ_max':>12} | {'ρ_DM_peak':>14} | {'ρ_needed_peak':>14} | {'gap':>10} | {'flat':>5}")
    print("-" * 95)

    for beta in BETA_SCAN:
        psi, dpsi = solve_soliton(r_arr, beta, GAL)
        if psi is None:
            print(f"  {beta:10.5f} | INTEGRATION FAILED")
            continue

        rho_dm = rho_dm_field(psi, dpsi)
        rho_dm_grad = rho_dm_gradient_only(dpsi)

        psi_center = float(psi[0])
        delta_psi_max = float(PSI_COS - psi.min())
        rho_dm_peak = float(rho_dm.max())
        rho_needed_at_peak_r = float(rho_needed[np.argmax(rho_dm)])
        gap = float(np.max(rho_needed) / max(rho_dm_peak, 1e-30))

        vc, geff, M_bar, M_dm = rotation_curve(r_arr, psi, dpsi, GAL)
        flat_ok, flat_dev = check_flatness(r_arr, vc)
        geff_max_dev = float(2 * LAMBDA0 * np.max(psi)**2)

        print(f"  {beta:10.5f} | {psi_center:12.6f} | {delta_psi_max:12.4e} | "
              f"{rho_dm_peak:14.4e} | {rho_needed_at_peak_r:14.4e} | "
              f"{gap:10.4e} | {'YES' if flat_ok else 'NO':>5}")

        prof = {
            'beta': beta,
            'Psi_center': psi_center,
            'delta_Psi_max': delta_psi_max,
            'rho_DM_peak': rho_dm_peak,
            'rho_DM_gradient_peak': float(rho_dm_grad.max()),
            'gap_rho': gap,
            'log10_gap': float(np.log10(gap)),
            'flatness_ok': flat_ok,
            'flatness_dev': float(flat_dev),
            'Geff_max_deviation': geff_max_dev,
            'Geff_ok': geff_max_dev < GEF_BOUND
        }
        diag['profiles'].append(prof)

        if gap < best_gap:
            best_gap = gap
            best_profile = (beta, psi.copy(), dpsi.copy(), vc.copy())

    # ── 3. BEST PROFILE DETAILS ───────────────────────────────────────────
    print("\n" + "─" * 60)
    print("PART 3: Best profile analysis")
    print("─" * 60)

    if best_profile is not None:
        beta_best, psi_best, dpsi_best, vc_best = best_profile
        rho_dm_best = rho_dm_field(psi_best, dpsi_best)
        Geff_best   = G_eff_ratio(psi_best)

        print(f"\n  Best β = {beta_best}")
        print(f"  Ψ range: [{psi_best.min():.6f}, {psi_best.max():.6f}]  (Ψ_cosmo = {PSI_COS:.4f})")
        print(f"  G_eff/G range: [{Geff_best.min():.6f}, {Geff_best.max():.6f}]")
        print(f"  Max G_eff deviation: {2*LAMBDA0*np.max(psi_best)**2*100:.3f}%  (limit: {GEF_BOUND*100:.0f}%)")
        print()
        print(f"  DM profile comparison at key radii:")
        print(f"  {'r [kpc]':>10} | {'Ψ_c':>12} | {'ρ_DM_soliton':>16} | {'ρ_DM_needed':>14} | {'ratio':>8}")
        print("  " + "-" * 68)
        r_check = [1.0, 4.0, 10.0, 20.0, 30.0]
        for r_i in r_check:
            idx = np.argmin(np.abs(r_arr - r_i))
            ratio = rho_dm_best[idx] / max(rho_needed[idx], 1e-30)
            print(f"  {r_i:10.1f} | {psi_best[idx]:12.6f} | "
                  f"{rho_dm_best[idx]:16.4e} | {rho_needed[idx]:14.4e} | {ratio:8.4e}")

    # ── 4. FIGURES ────────────────────────────────────────────────────────
    if best_profile is not None:
        _make_figures(r_arr, best_profile, diag, rho_needed)

    # ── 5. STRUCTURAL SUMMARY ─────────────────────────────────────────────
    print("\n" + "─" * 60)
    print("PART 4: Structural summary — Phase 2 DM failure modes")
    print("─" * 60)

    min_gap = min((p['log10_gap'] for p in diag['profiles']), default=99)
    print(f"""
  Phase 2 DM mechanisms tested (SIM114 + SIM115):

  ┌─────────────────────────────────────────────────────────────────┐
  │ Regime          │ G_eff │ ρ_DM/ρ_needed │ Rotation curves      │
  ├─────────────────────────────────────────────────────────────────┤
  │ β > β_thresh    │ FAIL  │ N/A           │ FAIL (G_eff→0 first) │
  │  (condensate)   │ 51%   │               │                      │
  ├─────────────────────────────────────────────────────────────────┤
  │ β < β_thresh    │ PASS  │ ~10^{min_gap:.0f}      │ FAIL (energy too    │
  │  (soliton, SIM115)│       │               │ small by 10^{min_gap:.0f})  │
  └─────────────────────────────────────────────────────────────────┘

  Structural root: λ_kpc = λ_cosmo × H₀² = {LAM_KPC:.2e} kpc⁻²
  The galactic field energy scale is:
    (c²/G_gal) × λ_kpc = {LAM_KPC * c_kms**2 / G_gal:.2e} M_sun/kpc³ per unit Ψ⁴
  Required DM at r~10 kpc: ~ 4×10⁶ M_sun/kpc³
  Available from field: ~ {gap_info['rho_DM_gradient']:.2e} M_sun/kpc³

  BOTH Phase 2 DM mechanisms FAIL.
  The β coupling cannot bridge the gap between:
    (a) cosmological energy density (ρ_DE ~ 1 in H₀=1 units)
    (b) galactic DM energy density (ρ_DM ~ 10⁷ M_sun/kpc³)
  The conversion factor H₀² ~ 5×10⁻¹⁵ kpc⁻² makes the galactic
  manifestation of the cosmological field energy negligible.

  REQUIRED NEW PHYSICS: a DM coupling that operates at galactic scales
  directly, not through the cosmological field energy. Options:
    1. A SEPARATE DM field (not Ψ) seeded by β coupling
    2. A non-minimal kinetic term that amplifies galactic Ψ variations
    3. Reinterpret: DM as a DIFFERENT sector than the Ψ condensate
""")

    # ── 6. VERDICT ────────────────────────────────────────────────────────
    print("=" * 72)
    print("VERDICT")
    print("=" * 72)
    print(f"  Gradient soliton: G_eff OK ✓  but  ρ_DM ~ 10^{min_gap:.0f} × ρ_needed  FAIL ✗")
    print(f"  Best gap: ρ_DM_peak / ρ_needed_peak ~ 10^{min_gap:.1f}")
    print(f"  RESULT: FAIL — DM density {10**(-min_gap):.0e} of required value")
    print(f"  Failure origin: H₀² suppression of λ_kpc renders galactic field")
    print(f"  energy negligible. Not a tuning issue — structural.")

    diag['overall_result'] = 'FAIL'
    diag['failure_mode'] = (
        f'Gradient soliton G_eff-safe but DM density ~10^{min_gap:.0f} of required. '
        f'Root cause: lambda_kpc = lambda_cosmo * H0^2 = {LAM_KPC:.2e} kpc^-2. '
        f'Galactic field energy (c2/G)*lambda_kpc*dPsi^2 negligible vs DM density needed. '
        f'H0^2 suppression is structural — cannot be fixed by tuning beta.'
    )

    out_path = os.path.join(OUTPUTS, 'sim115_diagnostics.json')
    with open(out_path, 'w') as f:
        json.dump(diag, f, indent=2)
    print(f"\n  Diagnostics saved → {out_path}")


def _make_figures(r_arr, best_profile, diag, rho_needed):
    beta_best, psi_best, dpsi_best, vc_best = best_profile
    rho_dm_best = rho_dm_field(psi_best, dpsi_best)
    Geff_best   = G_eff_ratio(psi_best)
    rho_b       = rho_baryon(r_arr, GAL)

    fig, axes = plt.subplots(2, 2, figsize=(12, 10))
    fig.suptitle(f'SIM115 — Phase 2 Gradient Soliton DM ({GAL_NAME}, β_best={beta_best})',
                 fontsize=12)

    colors = plt.cm.viridis(np.linspace(0.1, 0.9, len(diag['profiles'])))

    # Panel (a): Ψ profiles for all β
    ax = axes[0, 0]
    for prof, col in zip(diag['profiles'], colors):
        psi_i, _ = solve_soliton(r_arr, prof['beta'], GAL)
        if psi_i is not None:
            ax.plot(r_arr, psi_i, color=col, lw=1.5, label=f"β={prof['beta']:.4f}")
    ax.axhline(PSI_COS, color='k', ls='--', lw=1.5, label=f'Ψ_cosmo = {PSI_COS:.3f}')
    ax.axhline(PSI_MAX, color='r', ls=':', lw=1, label=f'Ψ_max = {PSI_MAX:.3f}')
    ax.set_xlabel('r [kpc]')
    ax.set_ylabel('Ψ(r)')
    ax.set_title('(a) Soliton profiles Ψ(r)')
    ax.legend(fontsize=7)
    ax.grid(True, alpha=0.3)

    # Panel (b): DM density comparison
    ax = axes[0, 1]
    ax.loglog(r_arr, rho_needed, 'k--', lw=2, label='ρ_DM needed (isothermal)')
    ax.loglog(r_arr, rho_b, 'b-', lw=1.5, alpha=0.5, label='ρ_baryon')
    for prof, col in zip(diag['profiles'], colors):
        psi_i, dpsi_i = solve_soliton(r_arr, prof['beta'], GAL)
        if psi_i is not None:
            rho_i = rho_dm_field(psi_i, dpsi_i)
            rho_i_safe = np.maximum(rho_i, 1e-30)
            ax.loglog(r_arr, rho_i_safe, color=col, lw=1.2,
                      label=f"β={prof['beta']:.4f}  (gap=10^{prof['log10_gap']:.0f})")
    ax.set_xlabel('r [kpc]')
    ax.set_ylabel('ρ [M_sun/kpc³]')
    ax.set_title('(b) DM density: soliton vs needed')
    ax.legend(fontsize=6)
    ax.grid(True, alpha=0.3)

    # Panel (c): G_eff ratio
    ax = axes[1, 0]
    for prof, col in zip(diag['profiles'], colors):
        psi_i, _ = solve_soliton(r_arr, prof['beta'], GAL)
        if psi_i is not None:
            ax.plot(r_arr, G_eff_ratio(psi_i), color=col, lw=1.5,
                    label=f"β={prof['beta']:.4f}")
    ax.axhline(1.0 - GEF_BOUND, color='r', ls='--', lw=1.5,
               label=f'G_eff lower bound ({(1-GEF_BOUND)*100:.0f}%)')
    ax.set_xlabel('r [kpc]')
    ax.set_ylabel('G_eff/G')
    ax.set_title('(c) G_eff/G (all soliton profiles)')
    ax.legend(fontsize=7)
    ax.grid(True, alpha=0.3)
    ax.set_ylim(0.9, 1.01)

    # Panel (d): Rotation curves
    ax = axes[1, 1]
    for prof, col in zip(diag['profiles'], colors):
        psi_i, dpsi_i = solve_soliton(r_arr, prof['beta'], GAL)
        if psi_i is not None:
            vc_i, *_ = rotation_curve(r_arr, psi_i, dpsi_i, GAL)
            ax.plot(r_arr, vc_i, color=col, lw=1.5, label=f"β={prof['beta']:.4f}")
    ax.axhline(V_FLAT, color='k', ls='--', lw=2, label=f'v_flat = {V_FLAT} km/s')
    ax.set_xlabel('r [kpc]')
    ax.set_ylabel('v_c [km/s]')
    ax.set_title('(d) Rotation curves')
    ax.legend(fontsize=7)
    ax.grid(True, alpha=0.3)
    ax.set_ylim(0, None)

    plt.tight_layout()
    fig_path = os.path.join(OUTPUTS, 'sim115_soliton_profiles.pdf')
    plt.savefig(fig_path, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"\n  Figure saved → {fig_path}")


if __name__ == '__main__':
    run()
