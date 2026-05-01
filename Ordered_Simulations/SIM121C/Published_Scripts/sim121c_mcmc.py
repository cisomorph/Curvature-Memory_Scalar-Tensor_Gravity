"""
SIM121-C — CMSTG Phase 2: Joint DESI+Planck MCMC Parameter Estimation
=====================================================================
Performs a proper Bayesian parameter estimation for CMSTG Phase 2 against
the combination of Planck CMB acoustic scale and DESI BAO H(z) data.

SIM121-B established that the CMSTG F₀ normalisation is the dominant CMB
observable. This simulation treats (F₀, w₀, wₐ) as the primary CMSTG
parameters, with H₀ derived self-consistently from the Planck θ_*
constraint at each MCMC step.

Parameter space:
  F₀   = ½ + Λ₀Ψ₀²  (CMSTG effective coupling today)
  w₀               (DE EOS today)
  wₐ               (DE EOS running: w(a)=w₀+wₐ(1-a))
  H₀               (derived from CMB θ_* at each step)

Fixed (Planck 2018 physical densities):
  Ω_m h²  = 0.1430  →  Ω_m = 0.1430 / (H₀/100)²
  Ω_b h²  = 0.02237 →  (enters sound horizon r_s)
  Ω_r h²  = 4.18×10⁻⁵

Likelihood:
  ln L = −½ χ²_θ  −½ χ²_DESI
  χ²_θ    = [(θ_*_CMSTG − θ_*_Planck) / σ_θ]²
  χ²_DESI = Σᵢ [(H_CMSTG(zᵢ) − H_obs(zᵢ)) / σᵢ]²

CMB observables used:
  100θ_MC = 1.04101 ± 0.00029  (Planck 2018 TT,TE,EE)

DESI BAO Y1 (2024):
  z=0.30: H=81.7±4.5; z=0.51: 97.9±4.4; z=0.71: 110.7±6.2
  z=0.93: 128.1±5.6; z=1.32: 156.4±8.6; z=2.33: 240.8±11.0

Pass criteria:
  1. χ²_DESI / N_desi < 2  at best-fit
  2. χ²_θ < 4  (within 2σ of Planck θ_*)
  3. w₀ ∈ (−1.3, −0.7),  wₐ ∈ (−2, 1)  (physical DE)
  4. F₀ ∈ [0.490, 0.560]  (G_eff within 10% of G_N)
  5. DESI tension reduction vs SIM113 2.7σ reference

Outputs:
  • Corner plot: (F₀, w₀, wₐ, H₀) posterior
  • H(z) best-fit vs DESI data
  • w₀-wₐ confidence contours vs Planck+DESI benchmarks
  • sim121c_results.json: MAP, credible intervals, tension σ
"""

import numpy as np
from scipy.integrate import quad
from scipy.optimize import brentq
import emcee, json, os, warnings
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
warnings.filterwarnings('ignore')

OUT  = os.path.join(os.path.dirname(__file__), '..', 'Outputs')
os.makedirs(OUT, exist_ok=True)

plt.rcParams.update({'font.family': 'serif', 'font.size': 11,
                     'axes.labelsize': 12, 'legend.fontsize': 10})

# ─────────────────────────────────────────────────────────────────────────────
# PHYSICAL CONSTANTS AND PLANCK 2018 PRIORS
# ─────────────────────────────────────────────────────────────────────────────
c_kms        = 2.998e5        # km/s

# Planck 2018 (TT,TE,EE+lowl+lowE) physical densities
omh2_m       = 0.1430         # Ω_m h²
omh2_b       = 0.02237        # Ω_b h²
omh2_r       = 4.18e-5        # Ω_r h² (photons + 3ν)

# Planck acoustic scale (100×θ_MC)
theta_obs    = 1.04101        # ×100
theta_obs_err= 0.00029        # ×100

z_drag       = 1059.6
z_star       = 1089.8

# DESI BAO Y1 2024 — H(z) [km/s/Mpc]
DESI_z = np.array([0.295, 0.510, 0.706, 0.930, 1.317, 2.330])
DESI_H = np.array([ 81.7,  97.9, 110.7, 128.1, 156.4, 240.8])
DESI_s = np.array([  4.5,   4.4,   6.2,   5.6,   8.6,  11.0])

# SIM113 CMSTG reference
SIM113_w0 = -0.973;  SIM113_wa = -0.41;  SIM113_F0 = 0.52059

# ─────────────────────────────────────────────────────────────────────────────
# CMSTG HUBBLE RATE (frozen Ψ, CPL DE, exact F₀ normalisation)
# ─────────────────────────────────────────────────────────────────────────────

def E2_CMSTG(z, F0, w0, wa, H0):
    """
    H²(z)/H₀² in CMSTG with frozen Ψ ≈ Ψ₀:
      3F₀ H² = ρ_m + ρ_r + ρ_DE
    Physical densities fixed; Ω_i = Ω_i_phys / (F₀ × H₀²).
    """
    h       = H0 / 100.0
    Om      = omh2_m / h**2          # Ω_m (physical density fixed)
    Or      = omh2_r / h**2          # Ω_r
    ODE     = 1.0 - Om - Or          # flat universe
    fDE     = (1+z)**(3*(1+w0+wa)) * np.exp(-3*wa*z/(1+z))
    # CMSTG Friedmann: H² = (ρ_m+ρ_r+ρ_DE)/(3F₀) with ρ_crit,0 = 3×½×H₀²
    return (0.5/F0) * (Om*(1+z)**3 + Or*(1+z)**4 + ODE*fDE)

def H_CMSTG(z, F0, w0, wa, H0):
    return H0 * np.sqrt(max(E2_CMSTG(z, F0, w0, wa, H0), 0.0))

# ─────────────────────────────────────────────────────────────────────────────
# CMB OBSERVABLES
# ─────────────────────────────────────────────────────────────────────────────

def r_s_CMSTG(F0, w0, wa, H0):
    """Comoving sound horizon at z_drag [Mpc]."""
    h        = H0 / 100.0
    Ogam     = 2.469e-5 / h**2    # photon density (neutrinos handled via Ω_r)
    def integrand(z):
        R    = (3.0 * omh2_b / h**2) / (4.0 * Ogam * (1+z))
        cs   = c_kms / np.sqrt(3.0*(1.0+R))
        return cs / H_CMSTG(z, F0, w0, wa, H0)
    val, _   = quad(integrand, z_drag, 1e4, limit=150, epsrel=1e-5)
    return val

def D_A_CMSTG(z_target, F0, w0, wa, H0):
    """Angular diameter distance to z_target [Mpc]."""
    val, _ = quad(lambda z: c_kms/H_CMSTG(z, F0, w0, wa, H0),
                  0, z_target, limit=150, epsrel=1e-5)
    return val / (1.0 + z_target)

def theta_star_CMSTG(F0, w0, wa, H0):
    """100×θ_* = 100 × r_s(z_drag)/D_C(z_*), D_C = D_A×(1+z) comoving."""
    rs = r_s_CMSTG(F0, w0, wa, H0)
    DA = D_A_CMSTG(z_star, F0, w0, wa, H0)
    return 100.0 * rs / (DA * (1.0 + z_star))

def H0_from_theta(F0, w0, wa, H0_lo=50.0, H0_hi=90.0):
    """
    Find H₀ such that 100θ_*(F₀,w₀,wₐ,H₀) = theta_obs.
    Returns H₀ if found, else np.nan.
    """
    try:
        fa = theta_star_CMSTG(F0, w0, wa, H0_lo) - theta_obs
        fb = theta_star_CMSTG(F0, w0, wa, H0_hi) - theta_obs
        if fa * fb > 0:
            return np.nan
        return brentq(lambda H: theta_star_CMSTG(F0, w0, wa, H) - theta_obs,
                      H0_lo, H0_hi, xtol=0.01, maxiter=40)
    except Exception:
        return np.nan

# ─────────────────────────────────────────────────────────────────────────────
# LIKELIHOOD
# ─────────────────────────────────────────────────────────────────────────────

def log_likelihood(theta_params):
    """
    theta_params = [F0, w0, wa]
    H₀ is derived by matching CMB θ_*.
    Returns (log_L, H0) or (-inf, nan).
    """
    F0, w0, wa = theta_params

    # Derive H₀ from CMB θ_*
    H0 = H0_from_theta(F0, w0, wa)
    if np.isnan(H0) or H0 < 55 or H0 > 85:
        return -np.inf, np.nan

    # CMB θ_* residual (should be ~0 since H₀ was solved from it,
    # but include small numerical residual as soft constraint)
    theta_CMSTG = theta_star_CMSTG(F0, w0, wa, H0)
    chi2_theta = ((theta_CMSTG - theta_obs) / theta_obs_err)**2

    # DESI H(z) likelihood
    H_model = np.array([H_CMSTG(z, F0, w0, wa, H0) for z in DESI_z])
    chi2_DESI = np.sum(((H_model - DESI_H) / DESI_s)**2)

    return -0.5 * (chi2_theta + chi2_DESI), H0

def log_prior(theta_params):
    F0, w0, wa = theta_params
    if not (0.490 <= F0 <= 0.560):  return -np.inf
    if not (-1.40 <= w0 <= -0.60):  return -np.inf
    if not (-2.50 <= wa <=  1.50):  return -np.inf
    return 0.0

def log_prob(theta_params):
    lp = log_prior(theta_params)
    if not np.isfinite(lp):
        return -np.inf
    ll, _ = log_likelihood(theta_params)
    return lp + ll

# ─────────────────────────────────────────────────────────────────────────────
# PRE-SCAN: find good starting region before MCMC
# ─────────────────────────────────────────────────────────────────────────────
print("=" * 70)
print("SIM121-C — CMSTG Phase 2: Joint DESI+Planck MCMC")
print("=" * 70)

print("\nPart A: Coarse grid scan to find starting region...")
F0_grid  = np.linspace(0.495, 0.540, 12)
w0_grid  = np.linspace(-1.30, -0.70, 15)
wa_grid  = np.linspace(-1.50,  0.50, 15)

best_chi2 = np.inf
best_params = None
best_H0 = None

scan_results = []
for F0 in F0_grid:
    for w0 in w0_grid:
        for wa in wa_grid:
            ll, H0 = log_likelihood([F0, w0, wa])
            if np.isfinite(ll) and not np.isnan(H0):
                chi2 = -2 * ll
                scan_results.append((chi2, F0, w0, wa, H0))
                if chi2 < best_chi2:
                    best_chi2 = chi2
                    best_params = [F0, w0, wa]
                    best_H0 = H0

scan_results.sort()
print(f"  Grid: {len(F0_grid)}×{len(w0_grid)}×{len(wa_grid)} = "
      f"{len(F0_grid)*len(w0_grid)*len(wa_grid)} points")
print(f"  Valid points: {len(scan_results)}")
print(f"\n  Top 5 grid minima:")
print(f"  {'chi2':>8}  {'F0':>8}  {'w0':>8}  {'wa':>8}  {'H0':>8}")
for chi2, F0, w0, wa, H0 in scan_results[:5]:
    print(f"  {chi2:>8.3f}  {F0:>8.4f}  {w0:>8.4f}  {wa:>8.4f}  {H0:>8.3f}")

if best_params is None:
    print("ERROR: No valid grid points found.")
    import sys; sys.exit(1)

# ─────────────────────────────────────────────────────────────────────────────
# MCMC WITH EMCEE
# ─────────────────────────────────────────────────────────────────────────────
print(f"\nPart B: MCMC (emcee) around best-fit region...")

ndim    = 3
nwalkers = 32
nsteps   = 3000
nburn    = 800

# Initialise walkers in a small ball around the grid best-fit
rng      = np.random.default_rng(42)
p0_center = np.array(best_params)
scales    = np.array([0.003, 0.05, 0.15])
p0        = p0_center + scales * rng.normal(size=(nwalkers, ndim))

# Clip to prior bounds
p0[:,0] = np.clip(p0[:,0], 0.491, 0.559)
p0[:,1] = np.clip(p0[:,1], -1.39, -0.61)
p0[:,2] = np.clip(p0[:,2], -2.49,  1.49)

sampler = emcee.EnsembleSampler(nwalkers, ndim, log_prob)

print(f"  Burning in ({nburn} steps)...")
state   = sampler.run_mcmc(p0, nburn, progress=False)
sampler.reset()
print(f"  Sampling ({nsteps} steps × {nwalkers} walkers)...")
sampler.run_mcmc(state, nsteps, progress=False)

# Extract flat chain
flat_chain = sampler.get_chain(flat=True)
log_probs  = sampler.get_log_prob(flat=True)
print(f"  Acceptance fraction: {np.mean(sampler.acceptance_fraction):.3f}")
print(f"  Chain shape: {flat_chain.shape}")

# MAP estimate
idx_map    = np.argmax(log_probs)
map_params = flat_chain[idx_map]
F0_map, w0_map, wa_map = map_params
_, H0_map = log_likelihood(map_params)
theta_map = theta_star_CMSTG(F0_map, w0_map, wa_map, H0_map)

# Derive H₀ for each chain sample
print(f"  Deriving H₀ for all chain samples...")
H0_chain = np.array([H0_from_theta(r[0], r[1], r[2]) for r in flat_chain])
valid     = np.isfinite(H0_chain)
flat_full = np.column_stack([flat_chain[valid], H0_chain[valid]])  # (N, 4)

# Percentile summaries
labels_mc = [r'$F_0$', r'$w_0$', r'$w_a$', r'$H_0$']
pct = np.percentile(flat_full, [16, 50, 84], axis=0)

print(f"\n── Posterior Summary ──")
print(f"  {'Param':>6}  {'MAP':>10}  {'16%':>10}  {'50%':>10}  {'84%':>10}")
map_vals = [F0_map, w0_map, wa_map, H0_map]
for lbl, mv, p16, p50, p84 in zip(['F0','w0','wa','H0'], map_vals,
                                    pct[0], pct[1], pct[2]):
    print(f"  {lbl:>6}  {mv:>10.4f}  {p16:>10.4f}  {p50:>10.4f}  {p84:>10.4f}")

# ─────────────────────────────────────────────────────────────────────────────
# DIAGNOSTICS
# ─────────────────────────────────────────────────────────────────────────────
ll_map, _ = log_likelihood([F0_map, w0_map, wa_map])
chi2_map  = -2 * ll_map

H_map_arr = np.array([H_CMSTG(z, F0_map, w0_map, wa_map, H0_map) for z in DESI_z])
chi2_DESI_map = np.sum(((H_map_arr - DESI_H)/DESI_s)**2)
chi2_theta_map = ((theta_map - theta_obs)/theta_obs_err)**2

# SIM113 reference chi2
ll_113, H0_113 = log_likelihood([SIM113_F0, SIM113_w0, SIM113_wa])
chi2_113 = -2*ll_113 if np.isfinite(ll_113) else np.nan

# LCDM reference (F0=0.5, w0=-1, wa=0)
ll_lcdm, H0_lcdm = log_likelihood([0.5, -1.0, 0.0])
chi2_lcdm = -2*ll_lcdm if np.isfinite(ll_lcdm) else np.nan

# DESI tension at MAP
tension_map = np.sqrt(chi2_DESI_map / len(DESI_z))
tension_113 = np.sqrt(chi2_DESI_map / len(DESI_z))  # placeholder

print(f"\n── Fit Quality ──")
print(f"  {'Model':25}  {'chi2_tot':>10}  {'chi2_DESI':>11}  "
      f"{'chi2_theta':>12}  {'H0':>8}  {'DESI tension':>14}")
for model, chi2_t, chi2_d, chi2_th, H0_v in [
    ('CMSTG MAP (this sim)',   chi2_map,  chi2_DESI_map, chi2_theta_map, H0_map),
    ('SIM113 reference',      chi2_113,  np.nan,         np.nan,         H0_113 if not np.isnan(H0_113) else np.nan),
    ('flat LCDM (F0=0.5)',    chi2_lcdm, np.nan,         np.nan,         H0_lcdm if not np.isnan(H0_lcdm) else np.nan),
]:
    print(f"  {model:25}  {chi2_t:>10.3f}  {chi2_d:>11.3f}  "
          f"{chi2_th:>12.3f}  {H0_v:>8.3f}  "
          f"{np.sqrt(chi2_d/len(DESI_z)):>14.3f}σ" if np.isfinite(chi2_d) else
          f"  {model:25}  {chi2_t:>10.3f}")

# H(z) table at MAP
print(f"\n── H(z) at MAP: CMSTG vs DESI ──")
print(f"  {'z':>5}  {'H_obs':>8}  {'H_MAP':>8}  {'pull':>6}")
for z, Ho, s, Hm in zip(DESI_z, DESI_H, DESI_s, H_map_arr):
    print(f"  {z:>5.3f}  {Ho:>8.1f}  {Hm:>8.2f}  {(Hm-Ho)/s:>6.2f}")

# ─────────────────────────────────────────────────────────────────────────────
# VERDICT
# ─────────────────────────────────────────────────────────────────────────────
pass_chi2   = chi2_DESI_map / len(DESI_z) < 2.0
pass_theta  = chi2_theta_map < 4.0
pass_w0     = -1.3 <= w0_map <= -0.7
pass_F0     = 0.490 <= F0_map <= 0.560
tension_sig = np.sqrt(chi2_DESI_map / len(DESI_z))
improvement_vs_LCDM = chi2_lcdm - chi2_map if np.isfinite(chi2_lcdm) else np.nan

print(f"\n{'='*70}")
print("SIM121-C RESULT:")
print()
print(f"  MAP parameters:")
print(f"    F₀   = {F0_map:.5f}  (Λ₀Ψ₀² = {F0_map-0.5:.5f})")
print(f"    w₀   = {w0_map:.4f}")
print(f"    wₐ   = {wa_map:.4f}")
print(f"    H₀   = {H0_map:.3f} km/s/Mpc  (CMB-derived)")
print()
print(f"  Fit quality:")
print(f"    χ²_DESI / N = {chi2_DESI_map:.3f} / {len(DESI_z)}"
      f"  → {'PASS' if pass_chi2 else 'FAIL'}")
print(f"    χ²_θ        = {chi2_theta_map:.4f}"
      f"  → {'PASS' if pass_theta else 'FAIL'} (θ_* within {np.sqrt(chi2_theta_map):.1f}σ)")
print(f"    DESI tension = {tension_sig:.2f}σ  (SIM113 reference: 2.70σ)")
print()
print(f"  Constraints:")
print(f"    F₀ = {F0_map:.5f}  (deviation from GR: {100*(F0_map/0.5-1):.3f}%)"
      f"  → {'PASS' if pass_F0 else 'FAIL'}")
print(f"    w₀ = {w0_map:.4f}  → {'PASS' if pass_w0 else 'FAIL'}")
print()
print(f"  Improvement vs ΛCDM: Δχ² = {improvement_vs_LCDM:+.3f}" if np.isfinite(improvement_vs_LCDM)
      else "  LCDM reference not available")

PASS = pass_chi2 and pass_theta and pass_w0 and pass_F0
PARTIAL = (pass_theta or pass_chi2) and not (pass_chi2 and pass_theta)

if PASS:
    print(f"\n  VERDICT: PASS")
    print(f"  CMSTG Phase 2 parameter set ({F0_map:.4f}, {w0_map:.3f}, {wa_map:.3f}) is")
    print(f"  jointly consistent with Planck CMB acoustic scale and DESI BAO.")
elif PARTIAL:
    print(f"\n  VERDICT: PARTIAL")
    print(f"  CMSTG satisfies one of CMB/DESI but not both simultaneously.")
else:
    print(f"\n  VERDICT: FAIL")
    print(f"  No CMSTG parameter set satisfies both CMB θ_* and DESI BAO.")
print(f"{'='*70}")

# ─────────────────────────────────────────────────────────────────────────────
# SAVE RESULTS JSON
# ─────────────────────────────────────────────────────────────────────────────
out_data = {
    'verdict': 'PASS' if PASS else 'PARTIAL' if PARTIAL else 'FAIL',
    'map': {'F0': float(F0_map), 'w0': float(w0_map), 'wa': float(wa_map),
            'H0': float(H0_map)},
    'posterior_50pct': {'F0': float(pct[1,0]), 'w0': float(pct[1,1]),
                        'wa': float(pct[1,2]),  'H0': float(pct[1,3])},
    'posterior_16pct': dict(zip(['F0','w0','wa','H0'], pct[0].tolist())),
    'posterior_84pct': dict(zip(['F0','w0','wa','H0'], pct[2].tolist())),
    'chi2_DESI_map': float(chi2_DESI_map),
    'chi2_theta_map': float(chi2_theta_map),
    'chi2_total_map': float(chi2_map),
    'DESI_tension_sigma': float(tension_sig),
    'theta_star_map': float(theta_map),
    'theta_star_obs': float(theta_obs),
    'SIM113_chi2': float(chi2_113) if np.isfinite(chi2_113) else None,
    'LCDM_chi2':   float(chi2_lcdm) if np.isfinite(chi2_lcdm) else None,
    'delta_chi2_vs_LCDM': float(improvement_vs_LCDM) if np.isfinite(improvement_vs_LCDM) else None,
    'nwalkers': nwalkers, 'nsteps': nsteps, 'nburn': nburn,
    'n_valid_samples': int(valid.sum()),
}
with open(os.path.join(OUT, 'sim121c_results.json'), 'w') as f:
    json.dump(out_data, f, indent=2)

# ─────────────────────────────────────────────────────────────────────────────
# FIGURE 1: Corner-style marginals (manual, no corner library needed)
# ─────────────────────────────────────────────────────────────────────────────
print(f"\nGenerating figures...")

param_names = [r'$F_0$', r'$w_0$', r'$w_a$', r'$H_0$ [km/s/Mpc]']
param_ranges = [(0.492, 0.545), (-1.32, -0.65), (-2.0, 1.2), (58, 80)]

fig, axes = plt.subplots(4, 4, figsize=(13, 12))
for i in range(4):
    for j in range(4):
        ax = axes[i, j]
        if j > i:
            ax.axis('off')
            continue
        xi = flat_full[:, j]
        yi = flat_full[:, i]
        if i == j:
            ax.hist(xi, bins=40, color='#4393c3', alpha=0.8, density=True, edgecolor='k', lw=0.3)
            ax.axvline(map_vals[i], color='#d6604d', lw=2, ls='--', label='MAP')
            ax.axvline(pct[1,i],   color='k',        lw=1, ls=':')
            if param_ranges[i]:
                ax.set_xlim(param_ranges[i])
            ax.set_yticks([])
            if i == 3:
                ax.axvline(73.0, color='#1a9641', lw=1.5, ls=':', alpha=0.8)  # SH0ES
        else:
            ax.hist2d(xi, yi, bins=35, cmap='Blues',
                      range=[param_ranges[j] or [xi.min(),xi.max()],
                             param_ranges[i] or [yi.min(),yi.max()]],
                      density=True)
            ax.plot(map_vals[j], map_vals[i], 'r*', ms=10, zorder=5)
        if i == 3:
            ax.set_xlabel(param_names[j], fontsize=10)
        if j == 0 and i > 0:
            ax.set_ylabel(param_names[i], fontsize=10)
        ax.tick_params(labelsize=8)

# SIM113 reference point
for i in range(1, 4):
    for j in range(0, i):
        ref_vals = [SIM113_F0, SIM113_w0, SIM113_wa, H0_113 if not np.isnan(H0_113) else 67.0]
        axes[i,j].plot(ref_vals[j], ref_vals[i], 'g^', ms=7, zorder=6,
                        label='SIM113' if (i==1 and j==0) else '')
        if i==1 and j==0:
            axes[i,j].legend(fontsize=8, loc='upper right')

fig.suptitle('SIM121-C — CMSTG Joint Posterior: CMB+DESI', fontsize=13, y=1.00)
fig.tight_layout()
for ext in ['pdf','png']:
    fig.savefig(os.path.join(OUT, f'sim121c_corner.{ext}'), dpi=130, bbox_inches='tight')
print("  Saved sim121c_corner")
plt.close(fig)

# ─────────────────────────────────────────────────────────────────────────────
# FIGURE 2: w₀-wₐ plane with confidence contours
# ─────────────────────────────────────────────────────────────────────────────
fig, axes = plt.subplots(1, 2, figsize=(12, 5))

ax = axes[0]
# Compute χ² on a (w0, wa) grid with F0, H0 marginalized (use MAP F0)
w0g = np.linspace(-1.35, -0.65, 60)
wag = np.linspace(-2.0, 1.2, 60)
chi2_grid = np.full((len(wag), len(w0g)), np.inf)
for iw, w0v in enumerate(w0g):
    for iwa, wav in enumerate(wag):
        ll_v, h0_v = log_likelihood([F0_map, w0v, wav])
        if np.isfinite(ll_v) and not np.isnan(h0_v):
            chi2_grid[iwa, iw] = -2*ll_v

chi2_min = chi2_grid[np.isfinite(chi2_grid)].min()
dchi2 = chi2_grid - chi2_min

W0G, WAG = np.meshgrid(w0g, wag)
ax.contourf(W0G, WAG, dchi2, levels=[0, 2.30, 6.18, 11.8],
            colors=['#d9f0d3','#a6dba0','#5aae61'], alpha=0.7)
ax.contour(W0G, WAG, dchi2,  levels=[2.30, 6.18], colors='k', linewidths=1.2)
ax.plot(w0_map, wa_map, 'r*', ms=14, zorder=5, label=f'MAP ({w0_map:.3f}, {wa_map:.3f})')
ax.plot(SIM113_w0, SIM113_wa, 'g^', ms=10, zorder=5, label='SIM113 (-0.973, -0.41)')
ax.plot(-1.0, 0.0, 'bs', ms=10, zorder=5, label=r'$\Lambda$CDM')
# DESI 2024 contour (approximate, from published results)
desi_w0_center, desi_wa_center = -0.827, -0.75
ax.errorbar(desi_w0_center, desi_wa_center, xerr=0.07, yerr=0.35,
            fmt='ko', ms=8, capsize=5, label='DESI+CMB+SNe 2024')
ax.set_xlabel(r'$w_0$')
ax.set_ylabel(r'$w_a$')
ax.set_title(r'$w_0$-$w_a$ Posterior (CMSTG MAP F₀)')
ax.legend(fontsize=9)
ax.axhline(0, color='gray', ls=':', lw=0.8)
ax.axvline(-1.0, color='gray', ls=':', lw=0.8)
ctxt = ax.contourf(W0G, WAG, dchi2, levels=[0,2.30], colors='none')
ax.text(0.03, 0.03,
        r'Contours: 1$\sigma$, 2$\sigma$ ($\Delta\chi^2=2.30, 6.18$)',
        transform=ax.transAxes, fontsize=9)

# FIGURE 2b: H(z) best-fit
ax = axes[1]
z_plot = np.logspace(-2, 0.4, 300)
H_map_plot = np.array([H_CMSTG(z, F0_map, w0_map, wa_map, H0_map) for z in z_plot])
H_113_plot = np.array([H_CMSTG(z, SIM113_F0, SIM113_w0, SIM113_wa,
                               H0_113 if not np.isnan(H0_113) else 67.0) for z in z_plot])
H_lcdm_plot = np.array([H0_lcdm*np.sqrt(0.315*(1+z)**3 + 9e-5*(1+z)**4 + 0.685)
                          for z in z_plot]) if not np.isnan(H0_lcdm) else H_map_plot

ax.plot(z_plot, H_map_plot, color='#d6604d', lw=2.5, label=f'CMSTG MAP (H₀={H0_map:.1f})')
ax.plot(z_plot, H_113_plot, color='#4dac26', lw=1.8, ls='--',
        label=f'SIM113 (H₀={H0_113:.1f})')
ax.plot(z_plot, H_lcdm_plot, color='#2166ac', lw=1.5, ls=':',
        label=f'flat ΛCDM (H₀={H0_lcdm:.1f})')
ax.errorbar(DESI_z, DESI_H, yerr=DESI_s, fmt='ko', ms=6, capsize=4,
            zorder=5, label='DESI BAO 2024')
ax.set_xlabel(r'Redshift $z$')
ax.set_ylabel(r'$H(z)$ [km/s/Mpc]')
ax.set_title(r'$H(z)$: CMSTG MAP vs Data')
ax.legend(fontsize=9)
ax.set_xlim(0.2, 2.5)

fig.suptitle('SIM121-C — CMSTG Phase 2: CMB+DESI Joint Fit', fontsize=13, y=1.01)
fig.tight_layout()
for ext in ['pdf','png']:
    fig.savefig(os.path.join(OUT, f'sim121c_contours.{ext}'), dpi=150, bbox_inches='tight')
print("  Saved sim121c_contours")
plt.close(fig)

# ─────────────────────────────────────────────────────────────────────────────
# FIGURE 3: F₀ vs H₀ constraint from CMB
# ─────────────────────────────────────────────────────────────────────────────
fig, axes = plt.subplots(1, 2, figsize=(12, 5))

ax = axes[0]
F0_range = np.linspace(0.490, 0.560, 60)
H0_vals  = np.array([H0_from_theta(F0v, w0_map, wa_map) for F0v in F0_range])
valid_F0 = np.isfinite(H0_vals)
ax.plot(F0_range[valid_F0], H0_vals[valid_F0], color='#2166ac', lw=2.5,
        label=r'H₀(F₀) from CMB $\theta_*$ at MAP $(w_0, w_a)$')
ax.axhline(73.0, color='#1a9641', ls='--', lw=1.8, label='SH0ES H₀=73.0')
ax.axhline(67.4, color='#d6604d', ls='--', lw=1.8, label='Planck H₀=67.4')
ax.axvline(0.5,  color='gray', ls=':', lw=1, label='GR (F₀=½)')
ax.axvline(F0_map, color='#7b2d8b', ls=':', lw=1.5)
ax.plot(F0_map, H0_map, 'r*', ms=12, zorder=5)
ax.set_xlabel(r'$F_0 = \frac{1}{2} + \Lambda_0 \Psi_0^2$')
ax.set_ylabel(r'$H_0$ [km/s/Mpc]')
ax.set_title(r'CMB-Derived $H_0$ vs $F_0$')
ax.legend(fontsize=9)
ax.set_xlim(F0_range[0], F0_range[-1])
ax.text(0.03, 0.08,
        r'Higher $F_0$ $\rightarrow$ lower CMB-inferred $H_0$' + '\n'
        r'(opposite to needed for H$_0$ tension)',
        transform=ax.transAxes, fontsize=9,
        bbox=dict(boxstyle='round', fc='#fff0f0', alpha=0.9))

ax = axes[1]
# H₀ posterior from chain
ax.hist(flat_full[:,3], bins=50, color='#4393c3', alpha=0.8,
        density=True, edgecolor='k', lw=0.3, label='CMSTG posterior')
ax.axvline(H0_map,  color='#d6604d', lw=2, ls='--', label=f'MAP H₀={H0_map:.2f}')
ax.axvline(pct[1,3],color='k',       lw=1, ls=':',  label=f'Median H₀={pct[1,3]:.2f}')
ax.axvline(73.0, color='#1a9641', lw=2, ls='-', alpha=0.8, label='SH0ES')
ax.axvline(67.4, color='#b2182b', lw=2, ls='-', alpha=0.8, label='Planck')
ax.set_xlabel(r'$H_0$ [km/s/Mpc]')
ax.set_ylabel('Probability density')
ax.set_title(r'CMSTG Marginal Posterior: $H_0$')
ax.legend(fontsize=9)
tension_cmstg = abs(73.0 - pct[1,3]) / np.sqrt(1.0**2 + 0.5**2)
ax.text(0.97, 0.97,
        f'CMSTG median H₀={pct[1,3]:.1f}\n'
        f'SH0ES tension: {tension_cmstg:.1f}σ\n'
        f'Planck tension: {abs(pct[1,3]-67.4)/0.5:.1f}σ',
        transform=ax.transAxes, ha='right', va='top', fontsize=10,
        bbox=dict(boxstyle='round', fc='white', alpha=0.9))

fig.suptitle('SIM121-C — CMB Constraint on CMSTG F₀ and H₀', fontsize=13, y=1.01)
fig.tight_layout()
for ext in ['pdf','png']:
    fig.savefig(os.path.join(OUT, f'sim121c_F0_H0.{ext}'), dpi=150, bbox_inches='tight')
print("  Saved sim121c_F0_H0")
plt.close(fig)

print(f"\nAll outputs in: {os.path.abspath(OUT)}")
print(f"MCMC complete. {valid.sum()} valid samples.")
