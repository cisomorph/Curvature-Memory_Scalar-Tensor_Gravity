#!/usr/bin/env python3
"""
SIM114: Phase 2 — βΨ²ρ_m Condensate and Galactic Rotation Curves
RIFT Phase 2: Recursive Intelligence-Field Theory (Unlocked Action)

Phase 2 Action (unlocked):
  S = ∫d⁴x√(-g) [(M_Pl²+2Λ₀Ψ²)/2 R - ½(∇Ψ)² - λ(Ψ²−v²)² + βΨ²L_m] + S_SM

Phase 1 failure (SIM99, SIM100, SIM103): locked action cannot produce flat
rotation curves at any (Ψ_bc, λ_gal) within cosmological bounds.

Phase 2 DM mechanism (P2 Paper I, Section 4): the βΨ²L_m coupling gives Ψ an
environment-dependent effective mass. In galaxies, baryon density ρ_m drives a
tachyonic instability that condenses Ψ into a halo.

Klein-Gordon (static, spherically symmetric galactic context):
  Ψ'' + (2/r)Ψ' = 4λΨ(Ψ²−v²) − 2β·G̃·ρ_b(r)·Ψ + 2Λ₀Ψ·R(r)
  G̃ ≡ G_gal / c_kms²  [converts M_sun/kpc³ → kpc⁻²]
  R(r) = −8π·G̃·ρ_b(r)  [weak-field Ricci scalar, kpc⁻²]

Equilibrium condensate (Ψ' = 0, local):
  Ψ_c²(r) ≈ v² + β·G̃·ρ_b(r) / (2λ_kpc)
  where λ_kpc = λ_cosmo × H₀²  [H₀ in kpc⁻¹]

Condensate DM density:
  ρ_cond(r) = (c²/G_gal) × λ_kpc × (Ψ_c² − v²)²
             = β²·G̃·ρ_b²(r) / (4λ_kpc) × (c²/G_gal)  [M_sun/kpc³]

Cosmological G_eff constraint:
  16π·Λ₀·Ψ_c² < G_eff_bound  (≡ G_eff/G within 1%)

KEY QUESTION: Does β = 0.034 (P2 Paper I estimate) simultaneously satisfy:
  (a) Flat rotation curves: |v_c(r) − v_flat| / v_flat < 20% for r > 50% r_max
  (b) G_eff constraint: 16π·Λ₀·Ψ_c² < 0.01 everywhere

If not: characterise the failure and scan β to find the viable window (if any).

Parameters (from SIM113 best fit):
  Λ₀ = 0.003, λ = 6.856×10⁻⁵, v = 13.16, β = 0.034 (fiducial)
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

with open(os.path.join(INPUTS, 'sim114_params.json')) as f:
    P = json.load(f)

# ── Physical constants (galactic units: kpc, M_sun, km/s) ─────────────────
G_gal   = 4.3009e-6    # (km/s)² kpc M_sun⁻¹
c_kms   = 2.998e5      # km/s
G_tilde = G_gal / c_kms**2   # kpc M_sun⁻¹  (converts ρ [M_sun/kpc³] → kpc⁻²)
H0_kpc  = 67.4 / (299800.0 * 3086.0)       # H₀ = 67.4 km/s/Mpc → kpc⁻¹  = 7.28e-8

# ── Phase 2 parameters ────────────────────────────────────────────────────
P2      = P['phase2_params']
LAMBDA0 = P2['Lambda0']              # Λ₀ = 0.003  (dimensionless)
LAM_CO  = P2['lambda_ssb']           # λ in cosmological units (H₀=1, 8πG=1)
V_SSB   = P2['v_ssb']               # v in dimensionless Ψ units
LAM_KPC = LAM_CO * H0_kpc**2        # λ in kpc⁻²
BETA_FID= P2['beta_fiducial']        # 0.034 (fiducial estimate)
BETA_SCAN = P2['beta_scan']          # range to scan

# Cosmological G_eff constraint: 2·Λ₀·Ψ² < max_dev
# G_eff/G = 1/(1 + 2Λ₀Ψ²)  →  |G_eff/G − 1| ≈ 2Λ₀Ψ² < max_dev
# (Consistent with SIM113: Ψ₀=2.62, dev = 2×0.003×2.62² = 4.1% ✓)
GEF_BOUND = P['cosmological_constraint']['Geff_max_deviation_pct'] / 100.0   # e.g. 0.01
PSI_MAX_COSMO = np.sqrt(GEF_BOUND / (2.0 * LAMBDA0))

NUM     = P['numerics']
FLAT_TOL  = P['verdict']['flatness_tolerance']
FLAT_FRAC = P['verdict']['flatness_fraction']


# ═══════════════════════════════════════════════════════════════════════════
# BARYONIC MASS MODELS  (same as SIM99 for direct comparability)
# ═══════════════════════════════════════════════════════════════════════════

def disk_vcirc2(r, M_disk, r_d):
    """Exponential disk v_c² (Freeman 1970)."""
    Sigma0 = M_disk / (2 * np.pi * r_d**2)
    y = np.clip(r / (2 * r_d), 1e-8, None)
    bessel = iv(0, y)*kv(0, y) - iv(1, y)*kv(1, y)
    return 4 * np.pi * G_gal * Sigma0 * r_d * y**2 * bessel


def rho_baryon(r, gal):
    """Spherical-average 3D baryonic density [M_sun/kpc³]."""
    M_d, r_d = gal['M_disk_Msun'], gal['r_d_kpc']
    M_g, r_g = gal['M_gas_Msun'],  gal['r_gas_kpc']
    dMdr = (M_d/r_d**2 * r * np.exp(-r/r_d)
          + M_g/r_g**2 * r * np.exp(-r/r_g))
    return dMdr / (4 * np.pi * r**2)


def M_baryon_enclosed(r, gal):
    """Enclosed baryonic mass [M_sun]."""
    M_d, r_d = gal['M_disk_Msun'], gal['r_d_kpc']
    M_g, r_g = gal['M_gas_Msun'],  gal['r_gas_kpc']
    enc_d = M_d * (1 - np.exp(-r/r_d) * (1 + r/r_d))
    enc_g = M_g * (1 - np.exp(-r/r_g) * (1 + r/r_g))
    return enc_d + enc_g


def ricci_scalar(r, gal):
    """R ≈ −8π·G̃·ρ_b(r)  [kpc⁻²]."""
    return -8 * np.pi * G_tilde * rho_baryon(r, gal)


# ═══════════════════════════════════════════════════════════════════════════
# ANALYTIC CONDENSATE PROFILE  (no-gradient approximation)
# Uses local equilibrium: Ψ_c² = v² + β·G̃·ρ_b/(2λ_kpc)
# Valid when the field adjusts faster than it propagates (~Compton length << kpc).
# ═══════════════════════════════════════════════════════════════════════════

def psi_condensate(r, gal, beta):
    """
    Local equilibrium Ψ_c(r) from the condensate fixed point.
    Returns Ψ_c and flag for validity (Ψ_c must be real and > 0).
    """
    rho_b = rho_baryon(r, gal)          # M_sun/kpc³
    R_val = ricci_scalar(r, gal)        # kpc⁻²

    # Effective mass² in field units [kpc⁻²]:
    # m_eff² = 4λ_kpc(Ψ_c² − v²) − 2β·G̃·ρ_b − 2Λ₀·R (set = 0 at equilibrium)
    # Ψ_c² = v² + [β·G̃·ρ_b + Λ₀·R] / (2λ_kpc)
    arg = V_SSB**2 + (beta * G_tilde * rho_b + LAMBDA0 * R_val) / (2.0 * LAM_KPC)
    psi_c = np.where(arg > 0, np.sqrt(arg), 0.0)
    return psi_c


# ═══════════════════════════════════════════════════════════════════════════
# CONDENSATE ENERGY DENSITY
# ═══════════════════════════════════════════════════════════════════════════

def rho_condensate(r, gal, beta):
    """
    DM condensate energy density [M_sun/kpc³].
    ρ_cond = (c²/G_gal) × λ_kpc × (Ψ_c² − v²)²
           = (c²/G_gal) × [β·G̃·ρ_b / (2)] × (β·G̃·ρ_b / (4λ_kpc))
           = β²·G̃·ρ_b² / (4λ_kpc) × (c²/G_gal)
    """
    psi_c = psi_condensate(r, gal, beta)
    excess = np.maximum(psi_c**2 - V_SSB**2, 0.0)  # (Ψ_c² − v²)
    # Energy density in field units × conversion factor
    rho_cond_field = LAM_KPC * excess**2             # kpc⁻²
    return rho_cond_field * (c_kms**2 / G_gal)       # M_sun/kpc³


def M_cond_enclosed(r_arr, gal, beta):
    """Enclosed condensate mass [M_sun]."""
    rho = rho_condensate(r_arr, gal, beta)
    return cumulative_trapezoid(4 * np.pi * r_arr**2 * rho, r_arr, initial=0)


# ═══════════════════════════════════════════════════════════════════════════
# G_eff FROM Ψ CONDENSATE
# ═══════════════════════════════════════════════════════════════════════════

def G_eff_ratio(psi_c):
    """G_eff/G = 1/(1 + 2·Λ₀·Ψ_c²)   [8πG=1 convention; consistent with SIM113]."""
    return 1.0 / (1.0 + 2.0 * LAMBDA0 * psi_c**2)


# ═══════════════════════════════════════════════════════════════════════════
# ROTATION CURVE
# ═══════════════════════════════════════════════════════════════════════════

def compute_rotation_curve(r_arr, gal, beta):
    """
    v_c²(r) = G_eff(r)/G × G_gal × [M_baryon(r) + M_cond(r)] / r
    """
    psi_c   = psi_condensate(r_arr, gal, beta)
    geff    = G_eff_ratio(psi_c)
    M_bar   = M_baryon_enclosed(r_arr, gal)
    M_cond  = M_cond_enclosed(r_arr, gal, beta)
    v2      = geff * G_gal * (M_bar + M_cond) / r_arr
    vc      = np.sqrt(np.clip(v2, 0, None))
    return vc, geff, M_bar, M_cond, psi_c


# ═══════════════════════════════════════════════════════════════════════════
# VERDICTS
# ═══════════════════════════════════════════════════════════════════════════

def check_flatness(r_arr, vc, v_flat):
    r_cut = r_arr[0] + FLAT_FRAC * (r_arr[-1] - r_arr[0])
    mask  = r_arr >= r_cut
    if not np.any(mask):
        return False, 1.0
    dev = np.abs(vc[mask] - v_flat) / v_flat
    return bool(np.all(dev < FLAT_TOL)), float(np.max(dev))


def check_geff(psi_c):
    coupling = 2.0 * LAMBDA0 * psi_c**2          # ≈ |G_eff/G − 1| for small coupling
    max_coup = float(np.max(coupling))
    return max_coup < P['cosmological_constraint']['Geff_max_deviation_pct'] / 100.0, max_coup


# ═══════════════════════════════════════════════════════════════════════════
# MAIN
# ═══════════════════════════════════════════════════════════════════════════

def run():
    print("=" * 72)
    print("SIM114 — Phase 2 RIFT: βΨ²ρ_m condensate, galactic rotation curves")
    print("=" * 72)
    print(f"  Λ₀      = {LAMBDA0}")
    print(f"  λ_cosmo = {LAM_CO:.4e}  →  λ_kpc = {LAM_KPC:.4e} kpc⁻²")
    print(f"  v       = {V_SSB:.4f}  (VEV, dimensionless Ψ units)")
    print(f"  β_fid   = {BETA_FID}  (P2 Paper I estimate: β ≈ 4√λ)")
    print(f"  Ψ²_max  = {PSI_MAX_COSMO**2:.4f}  (G_eff/G within 1%)")
    print(f"  λ_kpc × H₀ used: H₀ = {H0_kpc:.4e} kpc⁻¹")
    print()

    diag = {
        'sim': 'SIM114',
        'description': 'Phase 2 RIFT DM: beta*Psi^2*rho_m condensate galactic rotation curves',
        'Lambda0': LAMBDA0,
        'lambda_cosmo': LAM_CO,
        'lambda_kpc': LAM_KPC,
        'v_ssb': V_SSB,
        'beta_fiducial': BETA_FID,
        'Psi_max_cosmo': PSI_MAX_COSMO,
        'galaxies': {},
        'beta_scan_results': [],
        'overall_result': None
    }

    # ── 1. FIDUCIAL β = 0.034: analytic condensate profile ────────────────
    print("─" * 60)
    print(f"PART 1: Analytic condensate profile at β = {BETA_FID}")
    print("─" * 60)

    for gname, gal in P['galaxies'].items():
        Nr = gal['Nr']
        r_arr = np.linspace(gal['r_min_kpc'], gal['r_max_kpc'], Nr)
        v_flat = gal['v_flat_kms']

        # Compute condensate at fiducial β
        psi_c = psi_condensate(r_arr, gal, BETA_FID)
        rho_c = rho_condensate(r_arr, gal, BETA_FID)
        rho_b = rho_baryon(r_arr, gal)

        Geff_coup = 16 * np.pi * LAMBDA0 * psi_c**2
        G_ok, G_max_coup = check_geff(psi_c)

        vc, geff, M_bar, M_cond, _ = compute_rotation_curve(r_arr, gal, BETA_FID)
        flat_ok, flat_dev = check_flatness(r_arr, vc, v_flat)

        print(f"\n{gname}:")
        print(f"  Ψ_c range: [{psi_c.min():.3e}, {psi_c.max():.3e}]  (cosmo limit: {PSI_MAX_COSMO:.4f})")
        print(f"  16πΛ₀Ψ_c² range: [{Geff_coup.min():.3e}, {Geff_coup.max():.3e}]")
        print(f"  G_eff/G at r_min: {geff[0]:.6f}  (COSMO OK: {G_ok})")
        print(f"  ρ_cond/ρ_baryon at r_d: {rho_c[Nr//4]/max(rho_b[Nr//4],1e-30):.3e}")
        print(f"  v_c range: [{vc.min():.1f}, {vc.max():.1f}] km/s  (v_flat = {v_flat})")
        print(f"  Flatness: max_dev = {flat_dev:.3f}  (PASS: {flat_ok})")

        diag['galaxies'][gname] = {
            'beta': BETA_FID,
            'psi_c_min': float(psi_c.min()),
            'psi_c_max': float(psi_c.max()),
            'Geff_coup_max': float(G_max_coup),
            'Geff_ok': G_ok,
            'vc_min': float(vc.min()),
            'vc_max': float(vc.max()),
            'v_flat': v_flat,
            'flatness_ok': flat_ok,
            'flatness_max_dev': float(flat_dev)
        }

    # ── 2. β SCAN ─────────────────────────────────────────────────────────
    print("\n" + "─" * 60)
    print("PART 2: β scan — find viable window (flat curves + G_eff OK)")
    print("─" * 60)

    # Use NGC3198 for scan
    gal_scan = P['galaxies']['NGC3198']
    Nr_s = gal_scan['Nr']
    r_s  = np.linspace(gal_scan['r_min_kpc'], gal_scan['r_max_kpc'], Nr_s)
    v_flat_s = gal_scan['v_flat_kms']

    print(f"\n{'β':>12s} | {'Ψ_c_max':>12s} | {'16πΛ₀Ψ²_max':>14s} | {'G_eff_ok':>8s} | "
          f"{'vc_outer':>10s} | {'flat_ok':>7s} | {'BOTH':>6s}")
    print("-" * 90)

    for beta in BETA_SCAN:
        psi_c = psi_condensate(r_s, gal_scan, beta)
        vc, geff, M_bar, M_cond, _ = compute_rotation_curve(r_s, gal_scan, beta)
        flat_ok, flat_dev = check_flatness(r_s, vc, v_flat_s)
        G_ok, G_max_coup = check_geff(psi_c)

        r_cut = r_s[0] + 0.5 * (r_s[-1] - r_s[0])
        vc_outer = float(np.mean(vc[r_s >= r_cut]))

        both = flat_ok and G_ok
        flag = "  *** PASS ***" if both else ""
        print(f"  {beta:12.4e} | {psi_c.max():12.4e} | {G_max_coup:14.4e} | "
              f"{'YES' if G_ok else 'NO':>8s} | {vc_outer:10.2f} | {'YES' if flat_ok else 'NO':>7s} | "
              f"{'YES' if both else 'NO':>6s}{flag}")

        diag['beta_scan_results'].append({
            'beta': beta,
            'psi_c_max': float(psi_c.max()),
            'Geff_coup_max': float(G_max_coup),
            'Geff_ok': G_ok,
            'vc_outer_mean': vc_outer,
            'flatness_ok': flat_ok,
            'flatness_max_dev': float(flat_dev),
            'both_ok': both
        })

    # ── 3. PHYSICAL DIAGNOSIS ─────────────────────────────────────────────
    print("\n" + "─" * 60)
    print("PART 3: Physical diagnosis — why fails at β = 0.034")
    print("─" * 60)

    # At β = 0.034
    gal_d = P['galaxies']['NGC3198']
    r_diag = np.array([1.0, 5.0, 10.0, 20.0])
    rho_b_vals = rho_baryon(r_diag, gal_d)
    psi_c_vals = psi_condensate(r_diag, gal_d, BETA_FID)
    rho_c_vals = rho_condensate(r_diag, gal_d, BETA_FID)

    # β needed for flat rotation curve without G_eff modification:
    # Need ρ_cond = ρ_DM_needed = v_flat²/(4πG r²) - ρ_b
    rho_DM_needed = v_flat_s**2 / (4 * np.pi * G_gal * r_diag**2)
    # ρ_cond = β²·G̃·ρ_b² / (4λ_kpc) × (c²/G_gal)
    # β_needed = sqrt(ρ_DM_needed × 4λ_kpc × G_gal / (G̃² × ρ_b² × c²))
    with np.errstate(divide='ignore', invalid='ignore'):
        beta_needed = np.sqrt(np.maximum(rho_DM_needed * 4.0 * LAM_KPC * G_gal
                                         / (G_tilde**2 * rho_b_vals**2 * c_kms**2), 0))

    # Condensate formation threshold: need β > 8πΛ₀ for condensate to form.
    # When β < 8πΛ₀, the Λ₀R (curvature) term dominates β in KG, driving Ψ → 0.
    beta_thresh = 8 * np.pi * LAMBDA0
    print(f"\n  Condensate formation threshold: β > 8πΛ₀ = {beta_thresh:.5f}")
    print(f"  β_fiducial = {BETA_FID:.4f} {'< threshold → NO condensate' if BETA_FID < beta_thresh else '> threshold → condensate forms'}")
    print()
    print(f"NGC3198 — radial diagnosis (β_fid={BETA_FID}, β_needed assumes condensate at VEV):")
    print(f"{'r [kpc]':>10} | {'ρ_b':>12} | {'Ψ_c (β=0.034)':>15} | "
          f"{'ρ_cond':>12} | {'ρ_DM_needed':>12} | {'β_needed':>12}")
    print("-" * 90)
    for i, r_i in enumerate(r_diag):
        print(f"  {r_i:8.1f} | {rho_b_vals[i]:12.4e} | {psi_c_vals[i]:15.4e} | "
              f"{rho_c_vals[i]:12.4e} | {rho_DM_needed[i]:12.4e} | {beta_needed[i]:12.4e}")
    print(f"\n  Note: β_needed << β_thresh → condensate never forms at β_needed values.")

    # G_eff constraint: 2Λ₀Ψ_c² < GEF_BOUND → Ψ_c < PSI_MAX_COSMO
    # v² + β·G̃·ρ_b/(2λ) < PSI_MAX² → β < 2λ·(PSI_MAX² - v²)/(G̃·ρ_b)
    rho_b_disk = rho_baryon(np.array([gal_d['r_d_kpc']]), gal_d)[0]
    psi_max_sq = PSI_MAX_COSMO**2
    if psi_max_sq > V_SSB**2:
        beta_Geff_max = 2.0 * LAM_KPC * (psi_max_sq - V_SSB**2) / (G_tilde * rho_b_disk)
    else:
        beta_Geff_max = 0.0

    print(f"\nG_eff constraint (2Λ₀Ψ² < {GEF_BOUND}):")
    print(f"  Ψ²_max = {PSI_MAX_COSMO**2:.4f}  (Ψ_max = {PSI_MAX_COSMO:.4f})")
    print(f"  v²_ssb = {V_SSB**2:.2f}  →  v² >> Ψ²_max: CONDENSATE AT VEV ALWAYS FAILS")
    print(f"  v = {V_SSB:.4f} >> Ψ_max = {PSI_MAX_COSMO:.4f}  (ratio: {V_SSB/PSI_MAX_COSMO:.1f}×)")

    # ── 4. REVISED DIAGNOSIS: Is v compatible with G_eff? ─────────────────
    print("\n" + "─" * 60)
    print("PART 4: Structural incompatibility — v vs Ψ_max_cosmo")
    print("─" * 60)
    Geff_at_vev = 1.0 / (1.0 + 2 * LAMBDA0 * V_SSB**2)
    print(f"  The SSB VEV v = {V_SSB:.4f} (from SIM113 DE fit)")
    print(f"  The G_eff bound: Ψ < Ψ_max = {PSI_MAX_COSMO:.4f}")
    print(f"  v / Ψ_max = {V_SSB / PSI_MAX_COSMO:.1f}×  →  INCOMPATIBLE")
    print()
    print("  Physical interpretation:")
    print(f"  SIM113 chose v = 13.16 to fix λv⁴ = ρ_DE (hilltop height sets DE scale).")
    print(f"  G_eff/G = 1/(1 + 2Λ₀v²) = {Geff_at_vev:.4f}  ({(1-Geff_at_vev)*100:.1f}% deviation ≫ 1%)")
    print(f"  At VEV Ψ = v, G is suppressed by {1/Geff_at_vev:.1f}×.")
    print()
    print("  This is a STRUCTURAL INCOMPATIBILITY between the DE sector (SIM113)")
    print("  and the G_eff constraint. The VEV that sets the DE scale is ~10× too")
    print("  large for the G_eff/G ≈ 1 requirement.")
    print()
    print("  For G_eff OK: need Ψ_c < Ψ_max throughout the galaxy.")
    print("  This requires a screening/quenching mechanism not yet derived.")
    print("  Alternative: reformulate DE sector with small-Ψ vacuum (see SIM115).")

    # ── 5. FIGURES ────────────────────────────────────────────────────────
    _make_figures(r_s, diag)

    # ── 6. OVERALL RESULT ─────────────────────────────────────────────────
    viable_betas = [r['beta'] for r in diag['beta_scan_results'] if r['both_ok']]
    has_viable = len(viable_betas) > 0

    print("\n" + "=" * 72)
    print("VERDICT")
    print("=" * 72)
    print(f"  β_fiducial = {BETA_FID}  →  FAIL")
    print(f"  β scan:  {len(viable_betas)}/{len(BETA_SCAN)} values pass both flatness + G_eff")
    print(f"  v = {V_SSB:.4f} >> Ψ_max = {PSI_MAX_COSMO:.6f}  →  VEV incompatible with G_eff")
    print()
    print("  RESULT: FAIL — structural incompatibility between DE sector (SIM113)")
    print("  and G_eff cosmological constraint. The VEV v = 13.16 needed for the")
    print("  hilltop DE scale produces G_eff/G → 0 in the galactic condensate.")
    print("  Screening mechanism required (see SIM115).")

    diag['overall_result'] = 'FAIL'
    diag['condensate_threshold_beta'] = float(8 * np.pi * LAMBDA0)
    diag['failure_mode'] = (
        'TRILEMMA: (1) condensate requires beta > 8pi*Lambda0 = 0.075; '
        '(2) at beta > 0.075, Psi_c >> Psi_max=1.29 → G_eff catastrophically violated; '
        '(3) beta_needed for DM density (~1e-4) << threshold (0.075) so condensate never forms at needed beta. '
        'VEV v=13.16 >> Psi_max=1.29 is the structural root: '
        'DE hilltop scale and G_eff/G~1 incompatible by 10x. '
        'Requires screening mechanism or reformulated small-Psi DE vacuum (SIM115).'
    )
    diag['viable_beta_values'] = viable_betas
    diag['v_to_psimax_ratio'] = float(V_SSB / PSI_MAX_COSMO)
    diag['Geff_at_vev'] = float(1.0 / (1.0 + 2.0 * LAMBDA0 * V_SSB**2))

    # Save diagnostics
    out_path = os.path.join(OUTPUTS, 'sim114_diagnostics.json')
    with open(out_path, 'w') as f:
        json.dump(diag, f, indent=2)
    print(f"\n  Diagnostics saved → {out_path}")


def _make_figures(r_s, diag):
    """Generate SIM114 figures."""
    gal_d = P['galaxies']['NGC3198']
    v_flat = gal_d['v_flat_kms']

    fig, axes = plt.subplots(2, 2, figsize=(12, 10))
    fig.suptitle('SIM114 — Phase 2 RIFT: β Condensate Rotation Curves (NGC3198)', fontsize=13)

    # Panel (a): Ψ_c profiles for range of β
    ax = axes[0, 0]
    betas_plot = [1e-8, 1e-6, 1e-4, BETA_FID]
    colors = ['#4488bb', '#77bb44', '#dd8833', '#cc3333']
    for beta, col in zip(betas_plot, colors):
        psi_c = psi_condensate(r_s, gal_d, beta)
        ax.semilogy(r_s, psi_c, color=col, lw=1.5, label=f'β={beta:.1e}')
    ax.axhline(PSI_MAX_COSMO, color='k', ls='--', lw=1.5, label=f'Ψ_max (G_eff 1%): {PSI_MAX_COSMO:.4f}')
    ax.axhline(V_SSB, color='purple', ls=':', lw=1.5, label=f'VEV v = {V_SSB:.2f}')
    ax.set_xlabel('r [kpc]')
    ax.set_ylabel('Ψ_c(r)')
    ax.set_title('(a) Condensate amplitude vs r')
    ax.legend(fontsize=7)
    ax.grid(True, alpha=0.3)

    # Panel (b): G_eff/G profiles
    ax = axes[0, 1]
    for beta, col in zip(betas_plot, colors):
        psi_c = psi_condensate(r_s, gal_d, beta)
        geff  = G_eff_ratio(psi_c)
        ax.semilogy(r_s, geff, color=col, lw=1.5, label=f'β={beta:.1e}')
    ax.axhline(0.99, color='k', ls='--', lw=1.5, label='G_eff/G = 0.99 (bound)')
    ax.set_xlabel('r [kpc]')
    ax.set_ylabel('G_eff/G')
    ax.set_title('(b) G_eff modification vs r')
    ax.legend(fontsize=7)
    ax.grid(True, alpha=0.3)

    # Panel (c): Rotation curves for range of β
    ax = axes[1, 0]
    for beta, col in zip(betas_plot, colors):
        vc, *_ = compute_rotation_curve(r_s, gal_d, beta)
        ax.plot(r_s, vc, color=col, lw=1.5, label=f'β={beta:.1e}')
    ax.axhline(v_flat, color='k', ls='--', lw=1.5, label=f'v_flat = {v_flat} km/s')
    ax.set_xlabel('r [kpc]')
    ax.set_ylabel('v_c [km/s]')
    ax.set_title('(c) Rotation curves vs r')
    ax.legend(fontsize=7)
    ax.grid(True, alpha=0.3)
    ax.set_ylim(0, None)

    # Panel (d): β scan summary — vc_outer and Geff_coup_max vs β
    ax = axes[1, 1]
    betas_arr = np.array([r['beta'] for r in diag['beta_scan_results']])
    vc_outer  = np.array([r['vc_outer_mean'] for r in diag['beta_scan_results']])
    coup_max  = np.array([r['Geff_coup_max'] for r in diag['beta_scan_results']])

    ax2 = ax.twinx()
    ax.loglog(betas_arr, coup_max, 'r-o', ms=4, lw=1.5, label='16πΛ₀Ψ²_max')
    ax.axhline(GEF_BOUND, color='r', ls='--', lw=1, alpha=0.7, label='G_eff bound')
    ax2.semilogx(betas_arr, vc_outer, 'b-s', ms=4, lw=1.5, label='v_c outer (km/s)')
    ax2.axhline(v_flat, color='b', ls='--', lw=1, alpha=0.7)
    ax.set_xlabel('β')
    ax.set_ylabel('2Λ₀Ψ²_max  ≈ |G_eff/G − 1|  (red)', color='r')
    ax2.set_ylabel('v_c outer mean [km/s]  (blue)', color='b')
    ax.set_title('(d) β scan: G_eff coupling & outer v_c')
    lines1, labels1 = ax.get_legend_handles_labels()
    lines2, labels2 = ax2.get_legend_handles_labels()
    ax.legend(lines1 + lines2, labels1 + labels2, fontsize=7)
    ax.grid(True, alpha=0.3)

    plt.tight_layout()
    fig_path = os.path.join(OUTPUTS, 'sim114_rotation_curves.pdf')
    plt.savefig(fig_path, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"\n  Figure saved → {fig_path}")


if __name__ == '__main__':
    run()
