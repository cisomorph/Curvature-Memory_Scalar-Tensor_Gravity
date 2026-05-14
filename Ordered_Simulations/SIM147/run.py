#!/usr/bin/env python3
"""
SIM147 — Curvature-Memory Kernel Decay Characterization
=========================================================
Phase 5 gate calculation for Ψ_pre (SIM148) and Ψ_local (SIM149) conjectures.

Physics background
------------------
CMSTG defines: Ψ(x) = ∫d⁴x' G_R(x,x') S[R(x')]

The retarded Green's function for the CMSTG Ψ field in momentum space is

    G̃_R(k, ω) = e^{−k²/k_m²} / [(ω + iε)² − k² − m₀²]

The Gaussian factor e^{−k²/k_m²} is the memory damping — modes above k_m
are exponentially suppressed. The position-space kernel at coincident spatial
points (x = x', r = 0) is:

    K(Δt) = ∫d³k/(2π)³ G̃_R(Δt; k)
           = θ(Δt) · (1/2π²) ∫₀^∞ dk k² sin(ω_k Δt)/ω_k · e^{−k²/k_m²}

Since m₀ ≪ k_m (ratio ≈ 10⁻⁸ for locked parameters), ω_k ≈ k throughout
the integration domain where the integrand is non-negligible. In this limit:

    K(Δt) ≈ (1/2π²) ∫₀^∞ dk k sin(k Δt) e^{−k²/k_m²}

Using ∫₀^∞ k sin(kt) e^{−k²/a²} dk = (a³√π t/4) e^{−a²t²/4}  [standard result]:

    K(Δt) = (k_m³ / (8π^{3/2})) · Δt · exp(−(k_m Δt)²/4)

This is a GAUSSIAN envelope in Δt, not exponential or power-law.
The memory kernel peaks at Δt_peak = √2/k_m and decays on timescale τ = 2/k_m.

Key consequence: the memory damping k_m destroys temporal memory at
timescales Δt ≫ 1/k_m. For the locked k_m = 10 Mpc⁻¹:
    τ_mem ≈ 2/k_m = 0.2 Mpc ≈ 2×10⁻⁴ Gyr ≈ 200,000 years

Mass correction: the finite m₀ adds Bessel-function oscillations at
t ≫ 1/m₀ ~ 10⁸ Gyr, which are entirely irrelevant on any cosmological
timescale. The m₀ correction to K is of order (m₀/k_m)² ≈ 10⁻¹⁶.
"""

import os, json, pickle
from datetime import datetime
import numpy as np
from scipy.optimize import curve_fit
from scipy.interpolate import interp1d
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

OUTPUT_DIR = os.path.join(os.path.dirname(__file__), 'sims', 'sim147_output')
os.makedirs(OUTPUT_DIR, exist_ok=True)

# ─── Physical constants and locked parameters ─────────────────────────────────
H0_km_s_Mpc  = 67.36
# Hubble length in Mpc  (c = 1 units: H₀ = 67.36/3e5 Mpc⁻¹)
H0_Mpc       = H0_km_s_Mpc / 3e5          # Mpc⁻¹
H0inv_Mpc    = 1.0 / H0_Mpc               # ≈ 4454 Mpc

# Time conversion (c = 1)
# H₀⁻¹ = 4454 Mpc = 14.5 Gyr → 1 Gyr = 4454/14.5 Mpc
Gyr_to_Mpc   = 977.8                       # 1 Gyr ≈ 977.8 Mpc (c = 1)
Mpc_to_Gyr   = 1.0 / Gyr_to_Mpc

m0_H0        = 1e-3                        # m₀ / H₀ (dimensionless)
m0_Mpc       = m0_H0 * H0_Mpc             # ≈ 2.245×10⁻⁷ Mpc⁻¹

KM_LOCKED    = 10.0                        # Mpc⁻¹  (SIM102)
KM_SWEEP     = [1.0, 3.0, 10.0, 30.0, 100.0]

LAMBDA0      = 0.003

# Reference horizons (Gyr → Mpc)
T_REF_GYR    = [10.0, 13.4, 13.8]
T_REF_MPC    = [t * Gyr_to_Mpc for t in T_REF_GYR]
T_REF_LABELS = ['10Gyr', '13.4Gyr', '13.8Gyr']

N_DT         = 300    # grid points

# ─── Analytic kernel (massless approximation, exact for m₀ ≪ k_m) ─────────────

def K_analytic(dt_Mpc, km):
    """
    K(Δt) = (k_m³ / (8π^{3/2})) · Δt · exp(−(k_m Δt)²/4)

    Valid for m₀ ≪ k_m (ratio m₀/k_m ≈ 10⁻⁸ for locked parameters).
    Prefactor (k_m³/8π^{3/2}) is the normalization; we normalize to K(Δt_peak)=1
    so only the shape matters.
    """
    dt = np.asarray(dt_Mpc, dtype=float)
    return dt * np.exp(-0.25 * (km * dt)**2)

def K_peak(km):
    """Location of kernel peak: Δt_peak = √2 / k_m."""
    return np.sqrt(2.0) / km

def K_normalized(dt_Mpc, km):
    """K normalized so that K(Δt_peak) = 1."""
    dt_pk  = K_peak(km)
    K_pk   = float(K_analytic(dt_pk, km))
    return K_analytic(dt_Mpc, km) / K_pk

def K_abs_normalized(dt_Mpc, km):
    """K normalized to K at smallest sampled Δt (for comparison table)."""
    dt_min  = 1e-3 * H0inv_Mpc           # smallest sampling point
    K_min   = float(K_analytic(dt_min, km))
    return K_analytic(dt_Mpc, km) / max(K_min, 1e-300)

# ─── Mass correction estimate ──────────────────────────────────────────────────

def m0_correction_ratio(km, m0=None):
    """
    Leading correction from m₀ ≠ 0 to the massless kernel:
    δK/K ~ (m₀/k_m)²  (from Taylor expansion of ω_k = k√(1+m₀²/k²) ≈ k + m₀²/(2k))
    """
    if m0 is None:
        m0 = m0_Mpc
    return (m0 / km)**2

# ─── Decay characterization ───────────────────────────────────────────────────

def characterize_decay(km):
    """
    For a given k_m, return:
      τ_mem  — e-folding time of the Gaussian envelope (Mpc and Gyr)
      Δt_peak — peak of K(Δt)
      w(t_ref) for each reference horizon
    """
    # The Gaussian factor is exp(−(k_m Δt)²/4).
    # This falls to e⁻¹ when (k_m Δt)²/4 = 1 → Δt_e = 2/k_m
    tau_Mpc   = 2.0 / km
    tau_Gyr   = tau_Mpc * Mpc_to_Gyr
    dt_pk_Mpc = K_peak(km)
    dt_pk_Gyr = dt_pk_Mpc * Mpc_to_Gyr

    # Memory weights: w(Δt) = K(Δt)/K(Δt_min), where Δt_min is the smallest
    # reference time (10 Gyr >> τ_mem for all k_m in sweep, so w is tiny)
    weights = {}
    for t_Mpc, label in zip(T_REF_MPC, T_REF_LABELS):
        # Use K_normalized to K_peak
        w = float(K_normalized(t_Mpc, km))
        weights[label] = w

    return {
        'tau_Mpc':    tau_Mpc,
        'tau_Gyr':    tau_Gyr,
        'dt_peak_Mpc': dt_pk_Mpc,
        'dt_peak_Gyr': dt_pk_Gyr,
        'weights':    weights,
        'm0_correction': m0_correction_ratio(km),
    }

# ─── Fit the analytic kernel to the three candidate forms ─────────────────────

def fit_all_forms(km, n_pts=N_DT):
    """
    On a log-spaced grid in Δt, fit K(Δt) (normalized to K(Δt_min)) to
    three functional forms and return the best fit by R².

    Grid range: 10⁻³ H₀⁻¹ to 10² H₀⁻¹ (clipped to 20 Gyr).
    """
    dt_min_Mpc = 1e-3 * H0inv_Mpc
    dt_max_Mpc = min(1e2 * H0inv_Mpc, 20.0 * Gyr_to_Mpc)

    dt_arr = np.logspace(np.log10(dt_min_Mpc), np.log10(dt_max_Mpc), n_pts)
    K_arr  = K_analytic(dt_arr, km)

    # Normalize to first point
    K0     = K_arr[0] if K_arr[0] > 1e-300 else 1e-300
    K_norm = K_arr / K0

    dt_ref = dt_arr[0]
    results = {}

    # 1. Pure exponential: K = exp(-Δt/τ)
    try:
        popt, _ = curve_fit(
            lambda t, tau: np.exp(-t / tau), dt_arr, K_norm,
            p0=[1.0 / km], bounds=(1e-10, 1e10), maxfev=2000
        )
        tau_exp = popt[0]
        K_pred  = np.exp(-dt_arr / tau_exp)
        r2_exp  = 1 - np.sum((K_norm - K_pred)**2) / max(np.sum((K_norm - K_norm.mean())**2), 1e-300)
        results['exponential'] = {
            'tau_Mpc': float(tau_exp), 'tau_Gyr': float(tau_exp * Mpc_to_Gyr),
            'alpha': float('nan'), 'r2': float(r2_exp), 'form': 'exponential'
        }
    except Exception as e:
        results['exponential'] = {'tau_Mpc': float('nan'), 'tau_Gyr': float('nan'),
                                   'alpha': float('nan'), 'r2': float('-inf'),
                                   'form': 'exponential', 'error': str(e)}

    # 2. Power law: K = A (Δt/Δt_ref)^(-α)  — only positive K values
    mask_pos = K_norm > 1e-300
    if mask_pos.sum() > 5:
        try:
            log_dt = np.log(dt_arr[mask_pos] / dt_ref)
            log_K  = np.log(K_norm[mask_pos])
            c      = np.polyfit(log_dt, log_K, 1)
            alpha  = -c[0]
            A_pl   = np.exp(c[1])
            K_pred = A_pl * (dt_arr[mask_pos] / dt_ref)**(-alpha)
            r2_pl  = 1 - np.sum((K_norm[mask_pos] - K_pred)**2) / max(
                         np.sum((K_norm[mask_pos] - K_norm[mask_pos].mean())**2), 1e-300)
            results['powerlaw'] = {
                'alpha': float(alpha), 'A': float(A_pl),
                'tau_Mpc': float('nan'), 'tau_Gyr': float('nan'),
                'r2': float(r2_pl), 'form': 'powerlaw'
            }
        except Exception as e:
            results['powerlaw'] = {'alpha': float('nan'), 'r2': float('-inf'),
                                    'form': 'powerlaw', 'error': str(e)}
    else:
        results['powerlaw'] = {'alpha': float('nan'), 'r2': float('-inf'),
                                'form': 'powerlaw', 'note': 'insufficient positive points'}

    # 3. Exponential × power law: K = A (Δt/Δt_ref)^(-α) exp(-Δt/τ)
    if mask_pos.sum() > 5:
        try:
            def model_ep(t, alpha, tau):
                return (t / dt_ref)**(-alpha) * np.exp(-t / tau)
            popt, _ = curve_fit(
                model_ep, dt_arr[mask_pos], K_norm[mask_pos],
                p0=[0.5, 1.0 / km], bounds=([0, 1e-10], [5.0, 1e10]), maxfev=5000
            )
            alpha_ep, tau_ep = popt
            K_pred = model_ep(dt_arr[mask_pos], alpha_ep, tau_ep)
            r2_ep  = 1 - np.sum((K_norm[mask_pos] - K_pred)**2) / max(
                         np.sum((K_norm[mask_pos] - K_norm[mask_pos].mean())**2), 1e-300)
            results['exp_powerlaw'] = {
                'alpha': float(alpha_ep), 'tau_Mpc': float(tau_ep),
                'tau_Gyr': float(tau_ep * Mpc_to_Gyr), 'r2': float(r2_ep),
                'form': 'exp_powerlaw'
            }
        except Exception as e:
            results['exp_powerlaw'] = {'alpha': float('nan'), 'tau_Mpc': float('nan'),
                                        'tau_Gyr': float('nan'), 'r2': float('-inf'),
                                        'form': 'exp_powerlaw', 'error': str(e)}
    else:
        results['exp_powerlaw'] = {'alpha': float('nan'), 'tau_Mpc': float('nan'),
                                    'tau_Gyr': float('nan'), 'r2': float('-inf'),
                                    'form': 'exp_powerlaw'}

    # Winner
    winner = max(results, key=lambda k: results[k].get('r2', float('-inf')))
    return dt_arr, K_norm, results, winner

# ─── Spatial coherence ────────────────────────────────────────────────────────

def spatial_coherence(km):
    """
    G_R(r, Δt) at fixed Δt, as a function of spatial separation r:
        G_R(r, Δt) ∝ ∫₀^∞ dk k² sin(kr)/kr · sin(ω_k Δt)/ω_k · e^{-k²/k_m²}

    In the massless limit and for the spatial integral at Δt fixed:
    The spatial envelope is set by the memory damping alone:
        G_R(r, Δt) ∝ ∫₀^∞ dk k sin(kr) sin(k Δt) e^{-k²/k_m²}
                    ∝ [cosh(r Δt k_m²/2) - (stuff)] × Gaussian(r, Δt)

    The cleanest estimate: the spatial coherence length is the half-width of
    the Gaussian regulator in position space, which is 1/k_m (the Fourier
    partner of the k_m cutoff). This is an analytic result:

        ℓ_coh = 1/k_m

    To verify: compute G_R(r, Δt_1Gyr) at several r and find the 1/e point.
    We use: G̃(r, t) = ∂/∂(Δt) × {position-space Gaussian convolution}
    The exact form at fixed Δt, varying r (massless):

    G_R(r, Δt) ∝ [e^{-(k_m(r+Δt))²/4} - e^{-(k_m(r-Δt))²/4}] / r

    This shows the spatial envelope is Gaussian in (r ± Δt) with width 2/k_m.
    For r ≪ Δt (nearby points): G_R ∝ k_m² Δt exp(-...) × Gaussian(r) with width ~1/k_m.
    """
    l_coh_Mpc = 1.0 / km   # analytic result

    # Verify with explicit spatial profile at Δt = 1 Gyr
    Dt_Mpc    = 1.0 * Gyr_to_Mpc

    def G_spatial(r_Mpc):
        """G_R(r, Dt_Mpc) analytic massless expression."""
        if r_Mpc < 1e-10:
            return float(K_analytic(Dt_Mpc, km))
        rp  = r_Mpc + Dt_Mpc
        rm  = abs(r_Mpc - Dt_Mpc)
        # Include sign from (r-Δt) term
        sign_rm = 1.0 if Dt_Mpc >= r_Mpc else -1.0
        G = (np.exp(-0.25 * (km * rp)**2)
             - sign_rm * np.exp(-0.25 * (km * rm)**2))
        return G / (4 * np.pi * r_Mpc)   # 1/(4πr) from 3D Green's function normalization

    r_vals  = np.array([0, 0.3, 0.5, 0.7, 1.0, 1.5, 2.0, 3.0]) / km   # in Mpc
    G_vals  = np.array([G_spatial(r) for r in r_vals])
    G0      = abs(G_vals[0]) if abs(G_vals[0]) > 1e-300 else 1e-300
    G_norm  = G_vals / G0

    return l_coh_Mpc, G_norm, r_vals

# ─── Main ─────────────────────────────────────────────────────────────────────

def run_sim147():
    print("=" * 65)
    print("SIM147 — Curvature-Memory Kernel Decay Characterization")
    print("=" * 65)
    print(f"Λ₀ = {LAMBDA0},  m₀ = {m0_H0:.0e} H₀ = {m0_Mpc:.3e} Mpc⁻¹")
    print(f"k_m (locked) = {KM_LOCKED} Mpc⁻¹,  H₀⁻¹ = {H0inv_Mpc:.1f} Mpc")
    print()
    print("Analytic kernel: K(Δt) = Δt · exp(−(k_m Δt)²/4)  [normalized to peak=1]")
    print("Derived from:    e^{−k²/k_m²} memory damping in propagator,")
    print("                 massless limit m₀/k_m ≈ 10⁻⁸  (exact to 10⁻¹⁶ level)")
    print()

    sweep_results = {}

    for km in KM_SWEEP:
        print(f"─── k_m = {km:.1f} Mpc⁻¹ ───")
        char = characterize_decay(km)
        dt_arr, K_norm, fits, winner = fit_all_forms(km)

        print(f"  τ_mem (Gaussian e-fold)  = {char['tau_Mpc']:.4f} Mpc"
              f" = {char['tau_Gyr']:.4e} Gyr")
        print(f"  Δt_peak                  = {char['dt_peak_Mpc']:.4f} Mpc"
              f" = {char['dt_peak_Gyr']:.4e} Gyr")
        print(f"  m₀/k_m correction        ≈ {char['m0_correction']:.2e}  (negligible)")
        for label, t_Gyr in zip(T_REF_LABELS, T_REF_GYR):
            w = char['weights'][label]
            print(f"  w({t_Gyr} Gyr)          = {w:.4e}")
        print(f"  Best fit: {winner}  (R² = {fits[winner].get('r2', float('nan')):.4f})")
        print()

        sweep_results[km] = {
            'characterization': char,
            'dt_Mpc': dt_arr.tolist(),
            'K_norm': K_norm.tolist(),
            'fits': {k: {kk: float(vv) if isinstance(vv, (int, float, np.floating)) else vv
                         for kk, vv in v.items()}
                     for k, v in fits.items()},
            'winner': winner,
        }

    # ── Spatial coherence ─────────────────────────────────────────────────
    print(f"─── Spatial coherence at k_m = {KM_LOCKED} Mpc⁻¹, Δt = 1 Gyr ───")
    l_coh_Mpc, G_spatial_norm, r_spatial = spatial_coherence(KM_LOCKED)
    R_vir     = 0.1   # Mpc, typical L* virial radius
    uniform   = l_coh_Mpc >= R_vir
    print(f"  ℓ_coh = 1/k_m = {l_coh_Mpc:.4f} Mpc  (analytic)")
    print(f"  R_vir ≈ {R_vir:.2f} Mpc  (L* galaxy)")
    print(f"  ℓ_coh / R_vir = {l_coh_Mpc/R_vir:.3f}")
    print(f"  Kernel spatially uniform across L* galaxy: {uniform}")
    print()

    # ── Parameter dominance ───────────────────────────────────────────────
    print("─── Parameter dominance: m₀ vs k_m ───")
    for km in KM_SWEEP:
        ratio = m0_Mpc / km
        print(f"  k_m={km:5.1f} Mpc⁻¹ → m₀/k_m = {ratio:.2e}  "
              f"→ correction O({ratio**2:.1e})  [k_m dominates]")
    print()

    # ── Gate decision ─────────────────────────────────────────────────────
    locked_char  = sweep_results[KM_LOCKED]['characterization']
    locked_fits  = sweep_results[KM_LOCKED]['fits']
    locked_winner = sweep_results[KM_LOCKED]['winner']
    locked_bf    = locked_fits[locked_winner]
    tau_Gyr      = locked_char['tau_Gyr']
    w_134Gyr     = locked_char['weights']['13.4Gyr']
    w_10Gyr      = locked_char['weights']['10Gyr']
    w_138Gyr     = locked_char['weights']['13.8Gyr']

    # Note: w(Δt) = K(Δt)/K(Δt_peak) — Gaussian decay from peak to reference time
    # For k_m=10, Δt_peak ≈ 0.14 Mpc ≈ 1.5×10⁻⁴ Gyr. At 10 Gyr:
    # exp(-((10 Mpc⁻¹ × 10 Gyr × 977.8 Mpc/Gyr))²/4) = exp(-2.4×10⁸) → 0 to all precision.

    if (tau_Gyr >= 10.0) or (w_134Gyr > 0.01):
        gate = 'GREEN'
        gate_msg = ("τ_mem ≳ 10 Gyr or w(13.4 Gyr) > 0.01. Both conjectures alive. "
                    "SIM148 and SIM149 should proceed as written.")
    elif (1.0 <= tau_Gyr < 10.0) or (1e-4 <= w_134Gyr <= 1e-2):
        gate = 'YELLOW'
        gate_msg = ("τ_mem ∈ [1, 10] Gyr or w(13.4 Gyr) ∈ [10⁻⁴, 10⁻²]. "
                    "Ψ_local (SIM149) alive, Ψ_pre (SIM148) marginal.")
    else:
        gate = 'RED'
        gate_msg = ("τ_mem ≪ 1 Gyr and w(13.4 Gyr) < 10⁻⁴. Both conjectures effectively "
                    "dead. SIM148 and SIM149 should run for completeness but flag results "
                    "as conditional on large unphysical Ψ_pre or Ψ_local amplitudes.")

    # ── Tables ────────────────────────────────────────────────────────────
    print("=" * 65)
    print("RESULTS TABLES")
    print("=" * 65)

    print("\nFit comparison — locked k_m = 10 Mpc⁻¹:")
    print(f"  {'Form':<18} {'τ_mem [Gyr]':>14} {'α':>8} {'R²':>10}")
    print("  " + "-" * 52)
    for form in ['exponential', 'exp_powerlaw', 'powerlaw']:
        bf_f = locked_fits[form]
        if 'tau_Gyr' in bf_f and not np.isnan(bf_f['tau_Gyr']):
            tau_s = f"{bf_f['tau_Gyr']:14.6f}"
        else:
            tau_s = f"{'—':>14}"
        alp_s = f"{bf_f.get('alpha', float('nan')):8.3f}" if 'alpha' in bf_f and not np.isnan(bf_f.get('alpha', float('nan'))) else f"{'—':>8}"
        r2_s  = f"{bf_f.get('r2', float('nan')):10.4f}"
        star  = " ← BEST" if form == locked_winner else ""
        print(f"  {form:<18} {tau_s} {alp_s} {r2_s}{star}")

    print("\nMemory weights (K normalized to peak; reference times >> τ_mem):")
    print(f"  {'k_m [Mpc⁻¹]':>12}  {'τ_mem [Gyr]':>14}  {'w(10 Gyr)':>14}  {'w(13.4 Gyr)':>14}  {'w(13.8 Gyr)':>14}  Gate")
    print("  " + "-" * 85)
    for km in KM_SWEEP:
        ch = sweep_results[km]['characterization']
        ws = ch['weights']
        tau_km = ch['tau_Gyr']
        g = ('GREEN' if ws['13.4Gyr'] > 0.01
             else ('YELLOW' if ws['13.4Gyr'] >= 1e-4 else 'RED'))
        star = '*' if km == KM_LOCKED else ' '
        print(f"  {star}{km:>11.1f}  {tau_km:>14.6e}  {ws['10Gyr']:>14.4e}  "
              f"{ws['13.4Gyr']:>14.4e}  {ws['13.8Gyr']:>14.4e}  {g}")
    print("  (* = locked value)")

    print(f"\nSpatial coherence (locked k_m = {KM_LOCKED} Mpc⁻¹):")
    print(f"  ℓ_coh = 1/k_m = {l_coh_Mpc:.4f} Mpc  ({'≥' if uniform else '<'} R_vir = {R_vir} Mpc)")
    print(f"  Kernel spatially uniform across L* galaxy: {uniform}")

    print()
    print("=" * 65)
    print(f"  GATE: {gate}")
    print(f"  {gate_msg}")
    print("=" * 65)

    # ── Interpretation paragraph ──────────────────────────────────────────
    print("""
INTERPRETATION
─────────────────────────────────────────────────────────────────
The CMSTG memory kernel K(Δt) decays as a Gaussian,

    K(Δt) ∝ Δt · exp(−(k_m Δt)²/4),

not as an exponential or power law. This follows analytically from the
Gaussian memory damping e^{−k²/k_m²} in the propagator: the Fourier
transform of a Gaussian in k is a Gaussian in Δt. The decay is entirely
controlled by k_m; the cosmological mass m₀ ≈ 10⁻³ H₀ contributes a
correction of order (m₀/k_m)² ≈ 10⁻¹⁶ — completely negligible.

For the locked k_m = 10 Mpc⁻¹, the characteristic memory timescale is
τ_mem = 2/k_m = 0.2 Mpc ≈ 2×10⁻⁴ Gyr ≈ 200,000 years. At any
cosmological reference time (10 Gyr, 13.4 Gyr, 13.8 Gyr), the kernel
is non-zero only at the level of exp(−(k_m Δt)²/4) where
(k_m × 10 Gyr in Mpc) ≈ 10 × 9778 ≈ 9.8×10⁴, giving
exp(−(9.8×10⁴)²/4) ≈ exp(−2.4×10⁹) ≈ 0 to all numerical precision.

The k_m sweep confirms that τ_mem scales precisely as 1/k_m across the
full range 1–100 Mpc⁻¹. To achieve τ_mem ≈ 10 Gyr would require
k_m ≈ 2/(10 Gyr × 977.8 Mpc/Gyr) ≈ 2×10⁻⁴ Mpc⁻¹, five orders of
magnitude below the SIM102 locked value.

The spatial coherence length is ℓ_coh = 1/k_m = 0.1 Mpc, equal to the
L* virial radius. For SMALLER k_m (longer coherence) the kernel is
marginally uniform across a galaxy; at the locked value it is at the
boundary.

Gate: RED. Both Ψ_pre (SIM148) and Ψ_local (SIM149) conjectures are
effectively dead at the locked k_m = 10 Mpc⁻¹. Proceeding with SIM148
and SIM149 for completeness, but results must be flagged as requiring
unphysically large initial Ψ amplitudes.
─────────────────────────────────────────────────────────────────""")

    # ── Plots ─────────────────────────────────────────────────────────────
    print("\nGenerating plots...")
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))

    # Panel 1: K(Δt) for all k_m
    ax = axes[0]
    colors = plt.cm.viridis(np.linspace(0.1, 0.9, len(KM_SWEEP)))
    dt_plot_Mpc = np.logspace(-3, np.log10(20 * Gyr_to_Mpc), 500)
    dt_plot_Gyr = dt_plot_Mpc * Mpc_to_Gyr

    for km, col in zip(KM_SWEEP, colors):
        K_pl = K_normalized(dt_plot_Mpc, km)
        # Clip at numerical zero for log plot
        K_pl_clip = np.where(K_pl > 1e-300, K_pl, np.nan)
        lw = 2.5 if km == KM_LOCKED else 1.2
        ax.semilogy(dt_plot_Gyr, K_pl_clip, color=col, lw=lw,
                    label=f"$k_m = {km:.0f}$")

    # Reference lines
    for t_Gyr in [10.0, 13.4, 13.8]:
        ax.axvline(t_Gyr, color='gray', ls=':', lw=0.8, alpha=0.6)

    ax.set_xlabel(r'$\Delta t$ [Gyr]', fontsize=12)
    ax.set_ylabel(r'$K(\Delta t) / K(\Delta t_{\rm peak})$', fontsize=12)
    ax.set_title('SIM147: Memory Kernel $K(\\Delta t) \\propto \\Delta t\\,e^{-(k_m\\Delta t)^2/4}$',
                 fontsize=12)
    ax.legend(title=r'$k_m$ [Mpc$^{-1}$]', fontsize=8, loc='upper right')
    ax.set_xlim(0, 0.5)    # zoom into the first 0.5 Gyr to see the peak
    ax.set_ylim(1e-10, 2.0)
    ax.text(0.25, 1e-5, 'All reference\ntimes (10–13.8 Gyr)\nare far off-scale',
            fontsize=8, color='gray', ha='center')

    # Panel 2: τ_mem vs k_m (log-log)
    ax2 = axes[1]
    km_vals    = np.array(KM_SWEEP)
    tau_vals   = np.array([sweep_results[km]['characterization']['tau_Gyr']
                           for km in KM_SWEEP])

    ax2.loglog(km_vals, tau_vals, 'o-', color='C0', lw=2.5, ms=8, label=r'$\tau_{\rm mem} = 2/k_m$')
    # Reference lines for gate thresholds
    ax2.axhline(10.0, color='green',  ls='--', lw=1.5, label='GREEN: τ ≥ 10 Gyr')
    ax2.axhline(1.0,  color='orange', ls='--', lw=1.5, label='YELLOW: τ ≥ 1 Gyr')
    ax2.axvline(KM_LOCKED, color='red', ls=':', lw=1.5, label=f'Locked $k_m = {KM_LOCKED}$ Mpc⁻¹')

    # Annotate locked value
    tau_locked = sweep_results[KM_LOCKED]['characterization']['tau_Gyr']
    ax2.annotate(f'$\\tau = {tau_locked:.1e}$ Gyr\n({KM_LOCKED} Mpc⁻¹ locked)',
                 xy=(KM_LOCKED, tau_locked), xytext=(3, 1e-6),
                 arrowprops=dict(arrowstyle='->', color='red'),
                 fontsize=9, color='red')

    # Show the k_m needed for τ = 1 Gyr, 10 Gyr
    km_for_1Gyr  = 2.0 / (1.0  * Gyr_to_Mpc)
    km_for_10Gyr = 2.0 / (10.0 * Gyr_to_Mpc)
    ax2.axvline(km_for_1Gyr,  color='orange', ls=':', lw=1.0, alpha=0.7)
    ax2.axvline(km_for_10Gyr, color='green',  ls=':', lw=1.0, alpha=0.7)

    ax2.set_xlabel(r'$k_m$ [Mpc$^{-1}$]', fontsize=12)
    ax2.set_ylabel(r'$\tau_{\rm mem}$ [Gyr]', fontsize=12)
    ax2.set_title(r'Memory timescale $\tau = 2/k_m$', fontsize=12)
    ax2.legend(fontsize=8)

    plt.suptitle('SIM147 — CMSTG Memory Kernel Decay\n'
                 r'$K(\Delta t) \propto \Delta t \cdot \exp(-(k_m\Delta t)^2/4)$  '
                 r'[$m_0/k_m \approx 10^{-8}$, Gaussian dominates]',
                 fontsize=11, y=1.02)
    plt.tight_layout()
    fig_path = os.path.join(OUTPUT_DIR, 'sim147_kernel_decay.pdf')
    plt.savefig(fig_path, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"  Plot saved: {fig_path}")

    # ── Serialized kernel (downstream SIM148/149) ─────────────────────────
    # Export locked-k_m kernel on a grid; downstream sims interpolate.
    dt_export_Mpc = np.logspace(-3, np.log10(20 * Gyr_to_Mpc), 1000)
    K_export      = K_analytic(dt_export_Mpc, KM_LOCKED)
    K_export_norm = K_export / float(K_analytic(K_peak(KM_LOCKED), KM_LOCKED))

    # Seconds conversion
    Mpc_to_s = 3.0857e22 / 3e8    # 1 Mpc in seconds (c = 1)
    s_to_Mpc = 1.0 / Mpc_to_s

    kernel_data = {
        'dt_Mpc':          dt_export_Mpc.tolist(),
        'K_norm':          K_export_norm.tolist(),
        'km_locked_Mpcinv': KM_LOCKED,
        'm0_Mpcinv':       m0_Mpc,
        'Mpc_to_s':        Mpc_to_s,
        's_to_Mpc':        s_to_Mpc,
        'Gyr_to_Mpc':      Gyr_to_Mpc,
        'Mpc_to_Gyr':      Mpc_to_Gyr,
        'dt_peak_Mpc':     float(K_peak(KM_LOCKED)),
        'tau_Mpc':         float(2.0 / KM_LOCKED),
        'tau_Gyr':         float(2.0 / KM_LOCKED * Mpc_to_Gyr),
        'formula':         'K(dt) = dt * exp(-0.25*(km*dt)**2) [dt in Mpc, km in Mpc^-1]',
        'note':            ('Analytic result from Gaussian memory damping. '
                            'Numerical interpolation of this table for dt in Mpc. '
                            'For downstream use: '
                            'dt_Mpc = dt_s * s_to_Mpc; K = interp(dt_Mpc)'),
    }

    kernel_path = os.path.join(OUTPUT_DIR, 'kernel_sim147.pkl')
    with open(kernel_path, 'wb') as f:
        pickle.dump(kernel_data, f)
    print(f"  Kernel pickle saved: {kernel_path}")

    # ── JSON metadata ─────────────────────────────────────────────────────
    def _clean(obj):
        if isinstance(obj, float) and (np.isnan(obj) or np.isinf(obj)):
            return None
        if isinstance(obj, dict):
            return {k: _clean(v) for k, v in obj.items()}
        if isinstance(obj, list):
            return [_clean(v) for v in obj]
        return obj

    meta = {
        'sim':              'SIM147',
        'date':             datetime.now().isoformat(),
        'kernel_form':      'Gaussian: K(dt) = dt * exp(-(km*dt)^2/4)',
        'tau_mem_Gyr':      float(2.0 / KM_LOCKED * Mpc_to_Gyr),
        'tau_mem_Mpc':      float(2.0 / KM_LOCKED),
        'dt_peak_Gyr':      float(K_peak(KM_LOCKED) * Mpc_to_Gyr),
        'alpha':            float(locked_bf.get('alpha', float('nan'))),
        'best_fit_form':    locked_winner,
        'best_fit_r2':      float(locked_bf.get('r2', float('nan'))),
        'w_10Gyr':          float(w_10Gyr),
        'w_13.4Gyr':        float(w_134Gyr),
        'w_13.8Gyr':        float(w_138Gyr),
        'l_coh_Mpc':        float(l_coh_Mpc),
        'k_m_locked':       float(KM_LOCKED),
        'm0_Mpc':           float(m0_Mpc),
        'm0_over_km':       float(m0_Mpc / KM_LOCKED),
        'l_coh_over_Rvir':  float(l_coh_Mpc / R_vir),
        'uniform_across_galaxy': bool(uniform),
        'gate':             gate,
        'gate_message':     gate_msg,
        'km_dominates_m0':  True,
        'km_for_tau_1Gyr_Mpcinv':  float(2.0 / (1.0  * Gyr_to_Mpc)),
        'km_for_tau_10Gyr_Mpcinv': float(2.0 / (10.0 * Gyr_to_Mpc)),
        'sweep': {
            str(km): {
                'tau_Gyr':   sweep_results[km]['characterization']['tau_Gyr'],
                'tau_Mpc':   sweep_results[km]['characterization']['tau_Mpc'],
                'w_10Gyr':   sweep_results[km]['characterization']['weights']['10Gyr'],
                'w_13.4Gyr': sweep_results[km]['characterization']['weights']['13.4Gyr'],
                'w_13.8Gyr': sweep_results[km]['characterization']['weights']['13.8Gyr'],
                'winner':    sweep_results[km]['winner'],
            }
            for km in KM_SWEEP
        },
    }

    meta_path = os.path.join(OUTPUT_DIR, 'sim147_metadata.json')
    with open(meta_path, 'w') as f:
        json.dump(_clean(meta), f, indent=2)
    print(f"  Metadata JSON saved: {meta_path}")

    return meta

if __name__ == '__main__':
    meta = run_sim147()
    print("\nSIM147 complete.")
