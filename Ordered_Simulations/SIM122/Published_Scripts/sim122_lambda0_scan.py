"""
SIM122 — CMSTG Phase 3: Unlocked Λ₀ Scan
=========================================
SIM121C established that with Λ₀=0.003 locked, the CMB acoustic scale
and DESI H(z) cannot be simultaneously satisfied: the F₀ normalisation
forces H(z) ~10% too high for DESI.

This simulation unlocks Λ₀ and asks: is there a Λ₀ > 0 (non-trivial CMSTG)
that simultaneously satisfies Planck θ_* and DESI H(z)?

Physics:
  F₀ = ½ + Λ₀ × Ψ₀²,  Ψ₀ = 2.62 M_Pl  (SIM113 best-fit, held fixed)
  For Λ₀=0 → F₀=0.5 → ΛCDM (no CMSTG modification).
  For Λ₀=0.003 → F₀=0.521 (Phase 1 lock / SIM113).
  For Λ₀=0.0073 → F₀=0.560 (SIM121C CMB-derived MAP).

Free parameters:
  Λ₀  ∈ [0.0001, 0.010]   (was fixed at 0.003)
  w₀  ∈ [−1.40, −0.60]
  wₐ  ∈ [−2.50,  1.50]
  H₀  derived from Planck θ_* at each step

Fixed (Planck 2018 physical densities):
  Ω_m h² = 0.1430,  Ω_b h² = 0.02237,  Ω_r h² = 4.18×10⁻⁵

Likelihood:
  ln L = −½ χ²_θ  −½ χ²_DESI

Pass criteria:
  1. χ²_DESI / N < 2   at MAP
  2. χ²_θ < 4          (within 2σ of Planck θ_*)
  3. wₐ < 0            (thawing, DESI preference)
  4. F₀ ∈ [0.500, 0.560]
  5. Λ₀ > 0            (non-trivial CMSTG)
  6. DESI tension < 2.63σ  (better than SIM121C)

Outputs:
  • Corner plot: (Λ₀, w₀, wₐ, H₀, F₀)
  • Λ₀ tension profile: DESI tension vs Λ₀ (marginalised over w₀, wₐ)
  • H(z) best-fit vs DESI at MAP
  • sim122_results.json
"""

import numpy as np
from scipy.integrate import quad
from scipy.optimize import brentq, minimize_scalar
import emcee, json, os, warnings
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
warnings.filterwarnings('ignore')

OUT = os.path.join(os.path.dirname(__file__), '..', 'Outputs')
os.makedirs(OUT, exist_ok=True)

plt.rcParams.update({'font.family': 'serif', 'font.size': 11,
                     'axes.labelsize': 12, 'legend.fontsize': 10})

# ─────────────────────────────────────────────────────────────────────────────
# CONSTANTS AND DATA
# ─────────────────────────────────────────────────────────────────────────────
c_kms       = 2.998e5
omh2_m      = 0.1430
omh2_b      = 0.02237
omh2_r      = 4.18e-5

theta_obs     = 1.04101
theta_obs_err = 0.00029

z_drag = 1059.6
z_star = 1089.8

# DESI BAO Y1 2024 H(z) [km/s/Mpc]
DESI_z = np.array([0.295, 0.510, 0.706, 0.930, 1.317, 2.330])
DESI_H = np.array([ 81.7,  97.9, 110.7, 128.1, 156.4, 240.8])
DESI_s = np.array([  4.5,   4.4,   6.2,   5.6,   8.6,  11.0])

# CMSTG Phase 2 reference values
PSI0      = 2.62          # Ψ₀ in M_Pl units (SIM113)
L0_locked = 0.003         # Phase 1/2 locked value
F0_locked = 0.5 + L0_locked * PSI0**2   # = 0.52059

# SIM121C reference
SIM121C_tension = 2.6297
SIM121C_chi2    = 41.492
SIM121C_L0      = L0_locked

# ─────────────────────────────────────────────────────────────────────────────
# CMSTG HUBBLE RATE
# ─────────────────────────────────────────────────────────────────────────────

def F0_from_L0(L0):
    return 0.5 + L0 * PSI0**2

def E2_CMSTG(z, L0, w0, wa, H0):
    """H²(z)/H₀² in CMSTG: 3F₀H² = ρ with F₀ = ½ + Λ₀Ψ₀²."""
    F0  = F0_from_L0(L0)
    h   = H0 / 100.0
    Om  = omh2_m / h**2
    Or  = omh2_r / h**2
    ODE = 1.0 - Om - Or
    fDE = (1+z)**(3*(1+w0+wa)) * np.exp(-3*wa*z/(1+z))
    return (0.5/F0) * (Om*(1+z)**3 + Or*(1+z)**4 + ODE*fDE)

def H_CMSTG(z, L0, w0, wa, H0):
    return H0 * np.sqrt(max(E2_CMSTG(z, L0, w0, wa, H0), 0.0))

# ─────────────────────────────────────────────────────────────────────────────
# CMB OBSERVABLES
# ─────────────────────────────────────────────────────────────────────────────

def r_s_CMSTG(L0, w0, wa, H0):
    """Comoving sound horizon at z_drag [Mpc]."""
    h    = H0 / 100.0
    Ogam = 2.469e-5 / h**2
    def integrand(z):
        R  = (3.0 * omh2_b / h**2) / (4.0 * Ogam * (1+z))
        cs = c_kms / np.sqrt(3.0*(1.0+R))
        return cs / H_CMSTG(z, L0, w0, wa, H0)
    val, _ = quad(integrand, z_drag, 1e4, limit=150, epsrel=1e-5)
    return val

def D_A_CMSTG(z_target, L0, w0, wa, H0):
    """Proper angular diameter distance [Mpc]."""
    val, _ = quad(lambda z: c_kms / H_CMSTG(z, L0, w0, wa, H0),
                  0, z_target, limit=150, epsrel=1e-5)
    return val / (1.0 + z_target)

def theta_star_CMSTG(L0, w0, wa, H0):
    """100×θ_* = 100 × r_s(z_drag) / D_C(z_*) comoving."""
    rs = r_s_CMSTG(L0, w0, wa, H0)
    DA = D_A_CMSTG(z_star, L0, w0, wa, H0)
    return 100.0 * rs / (DA * (1.0 + z_star))

def H0_from_theta(L0, w0, wa, H0_lo=50.0, H0_hi=95.0):
    """Derive H₀ by matching 100θ_* = theta_obs."""
    try:
        fa = theta_star_CMSTG(L0, w0, wa, H0_lo) - theta_obs
        fb = theta_star_CMSTG(L0, w0, wa, H0_hi) - theta_obs
        if fa * fb > 0:
            return np.nan
        return brentq(lambda H: theta_star_CMSTG(L0, w0, wa, H) - theta_obs,
                      H0_lo, H0_hi, xtol=0.01, maxiter=50)
    except Exception:
        return np.nan

# ─────────────────────────────────────────────────────────────────────────────
# LIKELIHOOD
# ─────────────────────────────────────────────────────────────────────────────

def log_likelihood(params):
    """params = [L0, w0, wa]. Returns (log_L, H0)."""
    L0, w0, wa = params
    H0 = H0_from_theta(L0, w0, wa)
    if np.isnan(H0) or H0 < 50 or H0 > 90:
        return -np.inf, np.nan

    theta_CMSTG = theta_star_CMSTG(L0, w0, wa, H0)
    chi2_theta = ((theta_CMSTG - theta_obs) / theta_obs_err)**2

    H_model  = np.array([H_CMSTG(z, L0, w0, wa, H0) for z in DESI_z])
    chi2_DESI = np.sum(((H_model - DESI_H) / DESI_s)**2)

    return -0.5 * (chi2_theta + chi2_DESI), H0

def log_prior(params):
    L0, w0, wa = params
    if not (0.0001 <= L0 <= 0.010): return -np.inf
    if not (-1.40 <= w0 <= -0.60):  return -np.inf
    if not (-2.50 <= wa <=  1.50):  return -np.inf
    return 0.0

def log_prob(params):
    lp = log_prior(params)
    if not np.isfinite(lp):
        return -np.inf
    ll, _ = log_likelihood(params)
    return lp + ll

# ─────────────────────────────────────────────────────────────────────────────
# PART A: Λ₀ TENSION PROFILE (1D scan, marginalise w₀/wₐ)
# ─────────────────────────────────────────────────────────────────────────────
print("=" * 70)
print("SIM122 — CMSTG Phase 3: Unlocked Λ₀ Scan")
print("=" * 70)
print()
print("Part A: Λ₀ tension profile (1D scan, best w₀/wₐ at each Λ₀)...")

L0_grid   = np.linspace(0.0001, 0.010, 30)
w0_grid   = np.linspace(-1.35, -0.65, 20)
wa_grid   = np.linspace(-2.0,   1.0,  20)

profile_L0      = []
profile_chi2    = []
profile_tension = []
profile_F0      = []
profile_w0best  = []
profile_wabest  = []
profile_H0best  = []

for L0 in L0_grid:
    best_chi2 = np.inf
    best_row  = None
    for w0 in w0_grid:
        for wa in wa_grid:
            ll, H0 = log_likelihood([L0, w0, wa])
            if not np.isfinite(ll) or np.isnan(H0):
                continue
            H_model   = np.array([H_CMSTG(z, L0, w0, wa, H0) for z in DESI_z])
            chi2_DESI = np.sum(((H_model - DESI_H) / DESI_s)**2)
            if chi2_DESI < best_chi2:
                best_chi2 = chi2_DESI
                best_row  = (w0, wa, H0)
    if best_row is not None:
        w0b, wab, H0b = best_row
        tension = np.sqrt(best_chi2 / len(DESI_z))
        profile_L0.append(L0)
        profile_chi2.append(best_chi2)
        profile_tension.append(tension)
        profile_F0.append(F0_from_L0(L0))
        profile_w0best.append(w0b)
        profile_wabest.append(wab)
        profile_H0best.append(H0b)

profile_L0      = np.array(profile_L0)
profile_chi2    = np.array(profile_chi2)
profile_tension = np.array(profile_tension)
profile_F0      = np.array(profile_F0)

best_idx  = np.argmin(profile_chi2)
L0_opt    = profile_L0[best_idx]
chi2_opt  = profile_chi2[best_idx]
F0_opt    = profile_F0[best_idx]
w0_opt    = profile_w0best[best_idx]
wa_opt    = profile_wabest[best_idx]
H0_opt    = profile_H0best[best_idx]
ten_opt   = profile_tension[best_idx]

print(f"\n  Λ₀ scan complete.")
print(f"  Minimum χ²_DESI at Λ₀ = {L0_opt:.5f}  (F₀ = {F0_opt:.4f})")
print(f"  Best-fit: w₀={w0_opt:.3f}, wₐ={wa_opt:.3f}, H₀={H0_opt:.2f}")
print(f"  χ²_DESI = {chi2_opt:.3f} / {len(DESI_z)}  →  tension ≈ {ten_opt:.2f}σ")
print(f"  SIM121C reference: χ²={SIM121C_chi2:.3f}, tension={SIM121C_tension:.2f}σ  (Λ₀=0.003)")

# ─────────────────────────────────────────────────────────────────────────────
# PART B: COARSE GRID SCAN around optimum
# ─────────────────────────────────────────────────────────────────────────────
print()
print("Part B: Coarse grid scan (Λ₀, w₀, wₐ) around optimum region...")

L0_lo = max(0.0001, L0_opt - 0.003)
L0_hi = min(0.010,  L0_opt + 0.003)

L0_scan = np.linspace(L0_lo, L0_hi, 10)
w0_scan = np.linspace(-1.35, -0.65, 15)
wa_scan = np.linspace(-2.0,   1.0,  15)

grid_results = []
for L0 in L0_scan:
    for w0 in w0_scan:
        for wa in wa_scan:
            ll, H0 = log_likelihood([L0, w0, wa])
            if not np.isfinite(ll) or np.isnan(H0):
                continue
            H_model   = np.array([H_CMSTG(z, L0, w0, wa, H0) for z in DESI_z])
            chi2_DESI = np.sum(((H_model - DESI_H) / DESI_s)**2)
            F0        = F0_from_L0(L0)
            grid_results.append((chi2_DESI, L0, w0, wa, H0, F0))

grid_results.sort(key=lambda x: x[0])
valid = grid_results
print(f"  Valid grid points: {len(valid)}")
print(f"\n  Top 5 grid minima:")
print(f"  {'chi2':>8} {'L0':>8} {'F0':>7} {'w0':>7} {'wa':>7} {'H0':>8}")
for row in valid[:5]:
    print(f"  {row[0]:8.3f} {row[1]:8.5f} {row[5]:7.4f} {row[2]:7.3f} {row[3]:7.3f} {row[4]:8.3f}")

if len(valid) == 0:
    print("  ERROR: no valid grid points found. Exiting.")
    import sys; sys.exit(1)

# ─────────────────────────────────────────────────────────────────────────────
# PART C: MCMC
# ─────────────────────────────────────────────────────────────────────────────
print()
print("Part C: MCMC (emcee) over (Λ₀, w₀, wₐ)...")

best_grid = valid[0]
p0_center = np.array([best_grid[1], best_grid[2], best_grid[3]])
ndim, nwalkers = 3, 32
sigma_init = np.array([0.0005, 0.05, 0.10])
p0 = p0_center + sigma_init * np.random.randn(nwalkers, ndim)

# Clip to prior
p0[:, 0] = np.clip(p0[:, 0], 0.0001, 0.010)
p0[:, 1] = np.clip(p0[:, 1], -1.40, -0.60)
p0[:, 2] = np.clip(p0[:, 2], -2.50,  1.50)

sampler = emcee.EnsembleSampler(nwalkers, ndim, log_prob)

n_burn   = 600
n_sample = 2400
print(f"  Burning in ({n_burn} steps)...")
burn_state = sampler.run_mcmc(p0, n_burn, progress=False)
sampler.reset()
print(f"  Sampling ({n_sample} steps × {nwalkers} walkers)...")
sampler.run_mcmc(burn_state, n_sample, progress=False)

acc = np.mean(sampler.acceptance_fraction)
print(f"  Acceptance fraction: {acc:.3f}")

flat_chain = sampler.get_chain(flat=True)
print(f"  Chain shape: {flat_chain.shape}")

# Derive H₀ and F₀ for each sample
print("  Deriving H₀ and F₀ for chain samples...")
n_derive  = min(len(flat_chain), 5000)
idx_d     = np.random.choice(len(flat_chain), n_derive, replace=False)
H0_chain  = []
F0_chain  = []
valid_idx = []
for i in idx_d:
    L0, w0, wa = flat_chain[i]
    H0 = H0_from_theta(L0, w0, wa)
    if not np.isnan(H0) and 50 < H0 < 90:
        H0_chain.append(H0)
        F0_chain.append(F0_from_L0(L0))
        valid_idx.append(i)

chain_valid = flat_chain[valid_idx]
H0_chain    = np.array(H0_chain)
F0_chain    = np.array(F0_chain)
full_chain  = np.column_stack([chain_valid, H0_chain, F0_chain])
# columns: L0, w0, wa, H0, F0

pct = np.percentile(full_chain, [16, 50, 84], axis=0)
labels_all = ['L0', 'w0', 'wa', 'H0', 'F0']

print(f"\n── Posterior Summary ──")
print(f"  {'Param':>8}  {'MAP':>9}  {'16%':>9}  {'50%':>9}  {'84%':>9}")

# Compute MAP from chain (max log_prob)
log_probs = np.array([log_prob(flat_chain[i]) for i in valid_idx])
map_idx   = np.argmax(log_probs)
map_params = full_chain[map_idx]

for k, lab in enumerate(labels_all):
    print(f"  {lab:>8}  {map_params[k]:9.4f}  {pct[0,k]:9.4f}  {pct[1,k]:9.4f}  {pct[2,k]:9.4f}")

# MAP values
L0_map = map_params[0]
w0_map = map_params[1]
wa_map = map_params[2]
H0_map = map_params[3]
F0_map = map_params[4]

# ─────────────────────────────────────────────────────────────────────────────
# PART D: EVALUATE MAP
# ─────────────────────────────────────────────────────────────────────────────
H_model_map  = np.array([H_CMSTG(z, L0_map, w0_map, wa_map, H0_map) for z in DESI_z])
chi2_DESI_map = np.sum(((H_model_map - DESI_H) / DESI_s)**2)
theta_map    = theta_star_CMSTG(L0_map, w0_map, wa_map, H0_map)
chi2_theta_map = ((theta_map - theta_obs) / theta_obs_err)**2
tension_map  = np.sqrt(chi2_DESI_map / len(DESI_z))

print(f"\n── H(z) at MAP: CMSTG vs DESI ──")
print(f"  {'z':>6}  {'H_obs':>8}  {'H_MAP':>8}  {'pull':>6}")
for z, Ho, Hm, s in zip(DESI_z, DESI_H, H_model_map, DESI_s):
    print(f"  {z:6.3f}  {Ho:8.1f}  {Hm:8.2f}  {(Hm-Ho)/s:6.2f}")

# LCDM reference (Λ₀=0, same w0/wa)
L0_lcdm = 1e-9
H_lcdm  = np.array([H_CMSTG(z, L0_lcdm, -1.0, 0.0,
                            H0_from_theta(L0_lcdm, -1.0, 0.0)) for z in DESI_z])
H0_lcdm = H0_from_theta(L0_lcdm, -1.0, 0.0)
chi2_lcdm = np.sum(((H_lcdm - DESI_H) / DESI_s)**2) if not np.isnan(H0_lcdm) else np.nan
tension_lcdm = np.sqrt(chi2_lcdm / len(DESI_z)) if not np.isnan(chi2_lcdm) else np.nan

# SIM121C reference (locked Λ₀=0.003, MAP from that run)
L0_121c = 0.003
w0_121c, wa_121c = -0.6003, 0.4862
H0_121c = H0_from_theta(L0_121c, w0_121c, wa_121c)
H_121c  = np.array([H_CMSTG(z, L0_121c, w0_121c, wa_121c, H0_121c) for z in DESI_z]) if not np.isnan(H0_121c) else np.full(6, np.nan)
chi2_121c = np.sum(((H_121c - DESI_H) / DESI_s)**2) if not np.isnan(H0_121c) else np.nan

# ─────────────────────────────────────────────────────────────────────────────
# PART E: VERDICT
# ─────────────────────────────────────────────────────────────────────────────
print()
print("=" * 70)
print("SIM122 RESULT:")
print()
print(f"  MAP parameters:")
print(f"    Λ₀   = {L0_map:.6f}  (locked was: {L0_locked:.3f})")
print(f"    F₀   = {F0_map:.5f}  (SIM113: {F0_locked:.5f})")
print(f"    w₀   = {w0_map:.4f}")
print(f"    wₐ   = {wa_map:.4f}")
print(f"    H₀   = {H0_map:.3f} km/s/Mpc  (CMB-derived)")
print()
print(f"  Fit quality:")
print(f"    χ²_DESI / N = {chi2_DESI_map:.3f} / {len(DESI_z)}  ",
      end="")
print("→ PASS" if chi2_DESI_map/len(DESI_z) < 2 else "→ FAIL")
print(f"    χ²_θ        = {chi2_theta_map:.4f}  ",
      end="")
print("→ PASS" if chi2_theta_map < 4 else "→ FAIL")
print(f"    DESI tension = {tension_map:.2f}σ  (SIM121C: {SIM121C_tension:.2f}σ)")

delta_chi2 = chi2_DESI_map - SIM121C_chi2
print(f"    Δχ²_DESI vs SIM121C = {delta_chi2:+.3f}")

print()
print(f"  Constraints:")
print(f"    F₀ = {F0_map:.5f}  (deviation from GR: {abs(F0_map-0.5)/0.5*100:.2f}%)",
      end="  ")
print("→ PASS" if 0.500 <= F0_map <= 0.560 else "→ FAIL")
print(f"    wₐ = {wa_map:.4f}  (thawing requires wₐ < 0)",
      end="  ")
print("→ PASS" if wa_map < 0 else "→ FAIL")
print(f"    Λ₀ = {L0_map:.6f}  (non-trivial CMSTG: Λ₀ > 0)",
      end="  ")
print("→ PASS" if L0_map > 0.0005 else "→ MARGINAL")

if not np.isnan(chi2_lcdm):
    print(f"\n  ΛCDM reference (Λ₀→0): χ²_DESI={chi2_lcdm:.3f}, tension={tension_lcdm:.2f}σ")

# Determine verdict
pass_desi   = chi2_DESI_map/len(DESI_z) < 2
pass_theta  = chi2_theta_map < 4
pass_thaw   = wa_map < 0
pass_F0     = 0.500 <= F0_map <= 0.560
pass_L0     = L0_map > 0.0005
better_121c = chi2_DESI_map < SIM121C_chi2

n_pass = sum([pass_desi, pass_theta, pass_thaw, pass_F0, pass_L0])
if n_pass == 5 and pass_desi:
    verdict = "PASS"
elif n_pass >= 3 and better_121c:
    verdict = "PARTIAL"
else:
    verdict = "FAIL"

print(f"\n  VERDICT: {verdict}")
if verdict == "PASS":
    print(f"  CMSTG with Λ₀={L0_map:.5f} satisfies all criteria.")
    print(f"  Unlocking Λ₀ resolves the Phase 2 F₀ normalisation problem.")
elif verdict == "PARTIAL":
    print(f"  Partial improvement over SIM121C (Δχ²={delta_chi2:+.1f}).")
    print(f"  Some criteria fail — see analysis above.")
else:
    print(f"  Unlocking Λ₀ does not resolve the CMB+DESI tension.")
    print(f"  Phase 3 must consider an extended kinetic sector.")
print("=" * 70)

# ─────────────────────────────────────────────────────────────────────────────
# FIGURES
# ─────────────────────────────────────────────────────────────────────────────
print()
print("Generating figures...")

# --- Figure 1: Λ₀ tension profile ---
fig, axes = plt.subplots(1, 2, figsize=(12, 5))

ax = axes[0]
ax.plot(profile_L0, profile_chi2, 'o-', color='#2166ac', lw=2, ms=5)
ax.axvline(L0_locked, color='gray', ls='--', lw=1.5, label=r'$\Lambda_0=0.003$ (locked)')
ax.axvline(L0_opt,    color='#d73027', ls='-',  lw=1.5, label=r'$\Lambda_0^{\rm opt}$='+f'{L0_opt:.4f}')
ax.axhline(SIM121C_chi2, color='orange', ls=':', lw=1.5, label=f'SIM121C ref ({SIM121C_chi2:.1f})')
ax.axhline(2*len(DESI_z), color='green', ls=':', lw=1.5, label=f'PASS threshold ({2*len(DESI_z)})')
ax.set_xlabel(r'$\Lambda_0$')
ax.set_ylabel(r'$\chi^2_{\rm DESI}$ (best $w_0, w_a$)')
ax.set_title(r'DESI $\chi^2$ vs $\Lambda_0$')
ax.legend(fontsize=9)
ax.grid(alpha=0.3)

ax2 = axes[1]
ax2.plot(profile_F0, profile_chi2, 's-', color='#4dac26', lw=2, ms=5)
ax2.axvline(F0_locked, color='gray', ls='--', lw=1.5, label=f'SIM113 $F_0$={F0_locked:.4f}')
ax2.axvline(F0_opt,    color='#d73027', ls='-', lw=1.5, label=f'Optimal $F_0$={F0_opt:.4f}')
ax2.axhline(SIM121C_chi2, color='orange', ls=':', lw=1.5, label='SIM121C ref')
ax2.axhline(2*len(DESI_z), color='green', ls=':', lw=1.5, label='PASS threshold')
ax2.set_xlabel(r'$F_0 = \frac{1}{2} + \Lambda_0\bar\Psi^2$')
ax2.set_ylabel(r'$\chi^2_{\rm DESI}$ (best $w_0, w_a$)')
ax2.set_title(r'DESI $\chi^2$ vs $F_0$')
ax2.legend(fontsize=9)
ax2.grid(alpha=0.3)

plt.tight_layout()
for ext in ('pdf', 'png'):
    fig.savefig(os.path.join(OUT, f'sim122_tension_profile.{ext}'),
                dpi=150, bbox_inches='tight')
plt.close(fig)
print("  Saved sim122_tension_profile")

# --- Figure 2: H(z) comparison at MAP ---
fig, ax = plt.subplots(figsize=(8, 5))

z_fine = np.linspace(0.01, 2.5, 200)

H_map_fine   = [H_CMSTG(z, L0_map, w0_map, wa_map, H0_map) for z in z_fine]
H_lock_fine  = [H_CMSTG(z, L0_locked, w0_121c, wa_121c,
                        H0_from_theta(L0_locked, w0_121c, wa_121c)) for z in z_fine]
H0_lcdm2     = H0_from_theta(1e-9, -1.0, 0.0)
if not np.isnan(H0_lcdm2):
    H_lcdm_fine = [H_CMSTG(z, 1e-9, -1.0, 0.0, H0_lcdm2) for z in z_fine]
    ax.plot(z_fine, H_lcdm_fine, 'k--', lw=1.5, label=r'$\Lambda$CDM', alpha=0.6)

ax.plot(z_fine, H_map_fine,  color='#d73027', lw=2,
        label=fr'SIM122 MAP ($\Lambda_0$={L0_map:.5f}, $H_0$={H0_map:.1f})')
ax.plot(z_fine, H_lock_fine, color='#2166ac', lw=2, ls=':',
        label=fr'SIM121C MAP ($\Lambda_0$=0.003, $H_0$={H0_121c:.1f})')
ax.errorbar(DESI_z, DESI_H, yerr=DESI_s, fmt='ko', ms=7, capsize=4,
            label='DESI Y1 BAO', zorder=5)

ax.set_xlabel(r'Redshift $z$')
ax.set_ylabel(r'$H(z)$ [km/s/Mpc]')
ax.set_title(r'$H(z)$ at SIM122 MAP vs SIM121C vs DESI')
ax.legend()
ax.grid(alpha=0.3)
plt.tight_layout()
for ext in ('pdf', 'png'):
    fig.savefig(os.path.join(OUT, f'sim122_Hz.{ext}'), dpi=150, bbox_inches='tight')
plt.close(fig)
print("  Saved sim122_Hz")

# --- Figure 3: 2D posterior (Λ₀ vs w₀, Λ₀ vs wₐ) ---
fig, axes = plt.subplots(1, 2, figsize=(12, 5))

L0_samp = full_chain[:, 0]
w0_samp = full_chain[:, 1]
wa_samp = full_chain[:, 2]

axes[0].hexbin(L0_samp, w0_samp, gridsize=40, cmap='Blues', mincnt=1)
axes[0].axvline(L0_locked, color='gray', ls='--', lw=1.5, label=r'$\Lambda_0$ locked')
axes[0].axvline(L0_map,    color='red',  ls='-',  lw=1.5, label='MAP')
axes[0].set_xlabel(r'$\Lambda_0$')
axes[0].set_ylabel(r'$w_0$')
axes[0].set_title(r'Posterior: $\Lambda_0$ vs $w_0$')
axes[0].legend(fontsize=9)

axes[1].hexbin(L0_samp, wa_samp, gridsize=40, cmap='Greens', mincnt=1)
axes[1].axvline(L0_locked, color='gray', ls='--', lw=1.5, label=r'$\Lambda_0$ locked')
axes[1].axvline(L0_map,    color='red',  ls='-',  lw=1.5, label='MAP')
axes[1].axhline(0, color='black', ls=':', lw=1, label=r'$w_a=0$ (freezing boundary)')
axes[1].set_xlabel(r'$\Lambda_0$')
axes[1].set_ylabel(r'$w_a$')
axes[1].set_title(r'Posterior: $\Lambda_0$ vs $w_a$')
axes[1].legend(fontsize=9)

plt.tight_layout()
for ext in ('pdf', 'png'):
    fig.savefig(os.path.join(OUT, f'sim122_posterior_L0.{ext}'),
                dpi=150, bbox_inches='tight')
plt.close(fig)
print("  Saved sim122_posterior_L0")

# --- Figure 4: w₀-wₐ plane with DESI+Planck benchmarks ---
fig, ax = plt.subplots(figsize=(8, 6))

from matplotlib.patches import Ellipse
ax.hexbin(w0_samp, wa_samp, gridsize=40, cmap='Purples', mincnt=1, alpha=0.8)
ax.scatter(w0_map, wa_map, color='red', s=120, zorder=10, label=f'SIM122 MAP')
ax.scatter(-0.6003, 0.4862, color='blue', s=80, marker='s', zorder=10,
           label='SIM121C MAP')
ax.scatter(-0.973, -0.41, color='green', s=80, marker='^', zorder=10,
           label='SIM113 MAP')
# DESI central values
ax.scatter(-0.76, -0.79, color='orange', s=80, marker='*', zorder=10,
           label='DESI best-fit')
ax.axvline(-1.0, color='k', ls='--', lw=1, alpha=0.5)
ax.axhline(0,    color='k', ls=':', lw=1, alpha=0.5)
ax.set_xlabel(r'$w_0$')
ax.set_ylabel(r'$w_a$')
ax.set_title(r'$w_0$--$w_a$ posterior: SIM122 vs references')
ax.legend(fontsize=9)
ax.grid(alpha=0.3)
plt.tight_layout()
for ext in ('pdf', 'png'):
    fig.savefig(os.path.join(OUT, f'sim122_w0wa.{ext}'), dpi=150, bbox_inches='tight')
plt.close(fig)
print("  Saved sim122_w0wa")

# ─────────────────────────────────────────────────────────────────────────────
# SAVE JSON
# ─────────────────────────────────────────────────────────────────────────────
results = {
    "verdict": verdict,
    "map": {
        "L0": float(L0_map),
        "F0": float(F0_map),
        "w0": float(w0_map),
        "wa": float(wa_map),
        "H0": float(H0_map)
    },
    "posterior_50pct": {k: float(v) for k, v in zip(labels_all, pct[1])},
    "posterior_16pct": {k: float(v) for k, v in zip(labels_all, pct[0])},
    "posterior_84pct": {k: float(v) for k, v in zip(labels_all, pct[2])},
    "chi2_DESI_map": float(chi2_DESI_map),
    "chi2_theta_map": float(chi2_theta_map),
    "DESI_tension_sigma": float(tension_map),
    "SIM121C_chi2": float(SIM121C_chi2),
    "SIM121C_tension": float(SIM121C_tension),
    "delta_chi2_vs_SIM121C": float(delta_chi2),
    "L0_optimal_profile": float(L0_opt),
    "L0_locked": float(L0_locked),
    "n_valid_samples": len(valid_idx),
    "acceptance_fraction": float(acc),
    "pass_desi": bool(pass_desi),
    "pass_thawing": bool(pass_thaw),
    "pass_F0": bool(pass_F0)
}
with open(os.path.join(OUT, 'sim122_results.json'), 'w') as f:
    json.dump(results, f, indent=2)
print("  Saved sim122_results.json")

print(f"\nAll outputs in: {OUT}")
print("SIM122 complete.")
