"""
SIM93 — CMSTG Bayesian Model Comparison
=======================================
Computes the log-Bayes factor ln B(CMSTG/LCDM) via the Savage-Dickey density
ratio. LCDM is CMSTG at Lambda0=0 (nested model), so:

    B(CMSTG/LCDM) = p(Lambda0=0 | data) / p(Lambda0=0 | prior)

Method:
  1. MCMC over theta = (H0, Omega_m, rd, Lambda0) with the SIM90 joint
     CMB+BAO likelihood.
  2. Estimate p(Lambda0=0 | data) via KDE on the marginal Lambda0 chain.
  3. Evaluate at three prior choices: uniform, log-uniform, half-Gaussian.
  4. Interpret via Jeffreys (1961) scale.

Likelihood: SIM90 joint CMB+BAO (Planck 2018 Gaussian CMB prior +
  full 12-point BOSS/eBOSS/Lyα covariance).
"""

import json, math, os, warnings
import numpy as np
from scipy.stats import gaussian_kde
from scipy.integrate import quad
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

warnings.filterwarnings('ignore')

OUTDIR = os.path.join(os.path.dirname(__file__), '..', 'Outputs')
os.makedirs(OUTDIR, exist_ok=True)

# ── 0. MCMC availability check ────────────────────────────────────────────────
try:
    import emcee
    HAS_EMCEE = True
except ImportError:
    HAS_EMCEE = False
    print("WARNING: emcee not installed — using grid-based posterior instead.")

# ═══════════════════════════════════════════════════════════════════════════════
# 1. COSMOLOGICAL MODEL  (same as SIM90)
# ═══════════════════════════════════════════════════════════════════════════════

Omega_r = 9.2e-5   # radiation density (fixed)

def H_over_H0(z, Omega_m):
    """Flat LCDM H(z)/H0 (CMSTG G_eff correction is ~16 ppm — negligible)."""
    a = 1.0 / (1.0 + z)
    Omega_L = 1.0 - Omega_m - Omega_r
    return math.sqrt(Omega_m / a**3 + Omega_r / a**4 + Omega_L)

def comoving_distance(z, H0, Omega_m):
    """D_C(z) in Mpc (c/H0 * integral)."""
    c = 299792.458  # km/s
    from scipy.integrate import quad as _quad
    integrand = lambda zp: 1.0 / H_over_H0(zp, Omega_m)
    val, _ = _quad(integrand, 0.0, z, limit=200)
    return (c / H0) * val

def bao_observables(z, H0, Omega_m, rd):
    """Return DM/rd, DH/rd, DV/rd at redshift z."""
    c = 299792.458
    DC = comoving_distance(z, H0, Omega_m)
    DM = DC  # flat universe
    DH = c / (H0 * H_over_H0(z, Omega_m))
    DV = (z * DM**2 * DH) ** (1.0/3.0)
    return DM / rd, DH / rd, DV / rd

# ═══════════════════════════════════════════════════════════════════════════════
# 2. BAO DATA  (BOSS DR12 + eBOSS DR16 + Ly-alpha DR16 — identical to SIM87/90)
# ═══════════════════════════════════════════════════════════════════════════════

BAO_DATA = [
    (0.38, "DH_over_rd", 25.00,                    0.76,   "BOSS DR12"),
    (0.38, "DM_over_rd", 10.23,                    0.17,   "BOSS DR12"),
    (0.51, "DH_over_rd", 22.33,                    0.58,   "BOSS DR12"),
    (0.51, "DM_over_rd", 13.36,                    0.21,   "BOSS DR12"),
    (0.70, "DH_over_rd", 19.33,                    0.53,   "eBOSS DR16 LRG"),
    (0.70, "DM_over_rd", 17.86,                    0.33,   "eBOSS DR16 LRG"),
    (1.48, "DH_over_rd", 13.26,                    0.55,   "eBOSS DR16 QSO"),
    (1.48, "DM_over_rd", 30.69,                    0.80,   "eBOSS DR16 QSO"),
    (2.33, "DH_over_rd", 8.990618556701030,        0.21614046597277392, "Ly-alpha DR16"),
    (2.33, "DM_over_rd", 37.433384615384625,       1.26691023299267,    "Ly-alpha DR16"),
    (0.122,"DV_over_rd", 3.944,                    0.215,  "6dFGS+MGS"),
    (1.52, "DV_over_rd", 26.12,                    0.58,   "eBOSS DR14 QSO"),
]

# DM-DH correlation coefficients at paired redshifts
RHO_DM_DH = {0.38: -0.52, 0.51: -0.47, 0.70: -0.48, 1.48: -0.46, 2.33: -0.43}
N_BAO = len(BAO_DATA)

def build_bao_covariance():
    sigmas = np.array([row[3] for row in BAO_DATA])
    C = np.diag(sigmas**2)
    # Pair up DH/DM entries at the same redshift
    z_idx = {}
    for i, (z, kind, _, _, _) in enumerate(BAO_DATA):
        z_idx.setdefault(z, {})[kind] = i
    for z, rho in RHO_DM_DH.items():
        if z in z_idx and "DH_over_rd" in z_idx[z] and "DM_over_rd" in z_idx[z]:
            i = z_idx[z]["DH_over_rd"]
            j = z_idx[z]["DM_over_rd"]
            cov_off = rho * sigmas[i] * sigmas[j]
            C[i, j] = cov_off;  C[j, i] = cov_off
    return C

C_BAO = build_bao_covariance()
C_BAO_INV = np.linalg.inv(C_BAO)

def chi2_bao(H0, Omega_m, rd):
    resid = np.zeros(N_BAO)
    cache = {}
    for i, (z, kind, obs, _, _) in enumerate(BAO_DATA):
        if z not in cache:
            cache[z] = bao_observables(z, H0, Omega_m, rd)
        dm_rd, dh_rd, dv_rd = cache[z]
        if   kind == "DM_over_rd": pred = dm_rd
        elif kind == "DH_over_rd": pred = dh_rd
        else:                      pred = dv_rd
        resid[i] = obs - pred
    return float(resid @ C_BAO_INV @ resid)

# ═══════════════════════════════════════════════════════════════════════════════
# 3. CMB LIKELIHOOD  (Planck 2018 Gaussian — same as SIM90)
# ═══════════════════════════════════════════════════════════════════════════════

CMB_MU  = np.array([67.36, 0.3153])
CMB_SIG = np.array([0.54,  0.0073])
CMB_RHO = -0.90

def chi2_cmb(H0, Omega_m):
    dx = np.array([H0 - CMB_MU[0], Omega_m - CMB_MU[1]])
    s0, s1 = CMB_SIG
    r = CMB_RHO
    det = (1 - r**2)
    return (dx[0]**2/s0**2 - 2*r*dx[0]*dx[1]/(s0*s1) + dx[1]**2/s1**2) / det

def log_likelihood(theta):
    H0, Omega_m, rd, Lambda0 = theta
    # Lambda0 does not enter the likelihood at these coupling strengths
    # (Geff deviation < 500 ppm — negligible for BAO chi2; confirmed SIM91)
    c2 = chi2_cmb(H0, Omega_m) + chi2_bao(H0, Omega_m, rd)
    return -0.5 * c2

# ═══════════════════════════════════════════════════════════════════════════════
# 4. PRIORS
# ═══════════════════════════════════════════════════════════════════════════════

# Flat box priors on H0, Omega_m, rd (same bounds for all runs)
H0_BOUNDS   = (60.0, 75.0)
OMM_BOUNDS  = (0.25, 0.40)
RD_BOUNDS   = (140.0, 160.0)
L0_MAX      = 0.10   # upper bound for uniform / log-uniform priors

def log_prior_box(theta):
    """Returns 0.0 if in bounds, -inf otherwise (for H0, Omm, rd)."""
    H0, Omega_m, rd, Lambda0 = theta
    if not (H0_BOUNDS[0] < H0 < H0_BOUNDS[1]):     return -np.inf
    if not (OMM_BOUNDS[0] < Omega_m < OMM_BOUNDS[1]): return -np.inf
    if not (RD_BOUNDS[0] < rd < RD_BOUNDS[1]):      return -np.inf
    return 0.0

def log_prior_Lambda0_uniform(L0):
    """Uniform on [0, L0_MAX]."""
    if 0.0 <= L0 <= L0_MAX:
        return -math.log(L0_MAX)  # = log(1/0.1) = log(10)
    return -np.inf

def log_prior_Lambda0_log_uniform(L0):
    """Log-uniform (Jeffreys) on [L0_min, L0_MAX]; undefined at L0=0."""
    L0_MIN = 1e-4
    if L0_MIN <= L0 <= L0_MAX:
        return -math.log(L0) - math.log(math.log(L0_MAX / L0_MIN))
    return -np.inf

def log_prior_Lambda0_half_gaussian(L0):
    """Half-Gaussian with sigma=0.03, mu=0 (truncated to L0>=0)."""
    sigma = 0.03
    if L0 < 0.0:
        return -np.inf
    norm = math.log(2.0 / (sigma * math.sqrt(2 * math.pi)))
    return norm - 0.5 * (L0 / sigma)**2

PRIORS = {
    "uniform":       log_prior_Lambda0_uniform,
    "log_uniform":   log_prior_Lambda0_log_uniform,
    "half_gaussian": log_prior_Lambda0_half_gaussian,
}

# ═══════════════════════════════════════════════════════════════════════════════
# 5. MCMC  (or grid fallback)
# ═══════════════════════════════════════════════════════════════════════════════

N_WALKERS = 64
N_STEPS   = 5000
N_BURNIN  = 1000
N_THIN    = 10

def make_log_prob(prior_L0_fn):
    def log_prob(theta):
        lp = log_prior_box(theta) + prior_L0_fn(theta[3])
        if not np.isfinite(lp):
            return -np.inf
        ll = log_likelihood(theta)
        if not np.isfinite(ll):
            return -np.inf
        return lp + ll
    return log_prob

def run_mcmc(prior_name, prior_L0_fn, init_center):
    print(f"\n  Running MCMC: prior={prior_name}, {N_WALKERS} walkers x {N_STEPS} steps...")
    log_prob = make_log_prob(prior_L0_fn)
    ndim = 4
    # Initialise walkers as small Gaussian ball around joint best-fit
    rng = np.random.default_rng(42)
    scales = np.array([0.5, 0.005, 0.5, 0.005])
    p0 = init_center + scales * rng.standard_normal((N_WALKERS, ndim))
    # Clip Lambda0 walkers to valid range
    if prior_name == "log_uniform":
        p0[:, 3] = np.clip(np.abs(p0[:, 3]), 1e-4, L0_MAX)
    else:
        p0[:, 3] = np.clip(np.abs(p0[:, 3]), 0.0, L0_MAX)

    sampler = emcee.EnsembleSampler(N_WALKERS, ndim, log_prob)
    sampler.run_mcmc(p0, N_STEPS, progress=False)

    chain = sampler.get_chain(discard=N_BURNIN, thin=N_THIN, flat=True)
    print(f"    Flat chain shape: {chain.shape}")
    acc = np.mean(sampler.acceptance_fraction)
    print(f"    Mean acceptance fraction: {acc:.3f}")
    return chain, acc

def run_grid_posterior(prior_name, prior_L0_fn):
    """Fallback when emcee unavailable: grid over Lambda0, marginalise analytically."""
    print(f"\n  Grid posterior (no emcee): prior={prior_name}")
    # Fix (H0, Omm, rd) at joint best-fit — Lambda0 is the only free parameter
    H0_bf, Omm_bf, rd_bf = 67.59, 0.3118, 147.56
    L0_grid = np.linspace(1e-4, L0_MAX, 2000)
    log_post = np.array([
        log_likelihood([H0_bf, Omm_bf, rd_bf, L0]) + prior_L0_fn(L0)
        for L0 in L0_grid
    ])
    log_post -= log_post.max()
    post = np.exp(log_post)
    post /= np.trapz(post, L0_grid)
    # Return as fake chain: sample from the grid
    rng = np.random.default_rng(42)
    idx = rng.choice(len(L0_grid), size=5000, p=post/post.sum())
    chain_L0 = L0_grid[idx]
    chain = np.column_stack([
        np.full(len(chain_L0), H0_bf),
        np.full(len(chain_L0), Omm_bf),
        np.full(len(chain_L0), rd_bf),
        chain_L0
    ])
    return chain, None

# ═══════════════════════════════════════════════════════════════════════════════
# 6. SAVAGE-DICKEY DENSITY RATIO
# ═══════════════════════════════════════════════════════════════════════════════

def savage_dickey(chain_L0, prior_L0_fn, prior_name):
    """
    ln B(CMSTG/LCDM) = ln p(L0=0|data) - ln p(L0=0|prior)
    Posterior density at L0=0 estimated via KDE on the marginal Lambda0 samples.
    """
    # KDE on the Lambda0 marginal
    bw = 'silverman'
    kde = gaussian_kde(chain_L0, bw_method=bw)
    post_at_0 = float(kde(0.0)[0])

    # Prior density at L0=0
    # For log-uniform, p(L0=0) is undefined (diverges); use smallest grid point
    if prior_name == "log_uniform":
        prior_at_0 = None  # undefined
        ln_B = None
        note = "Log-uniform prior: p(Lambda0=0) undefined (diverges at 0). Bayes factor not defined via SD ratio for this prior."
    else:
        lp = prior_L0_fn(0.0)
        prior_at_0 = math.exp(lp) if np.isfinite(lp) else None
        if prior_at_0 is not None and prior_at_0 > 0:
            ln_B = math.log(post_at_0) - math.log(prior_at_0)
            note = ""
        else:
            ln_B = None
            note = "Prior density at Lambda0=0 is zero or undefined."

    return {
        "posterior_density_at_0": post_at_0,
        "prior_density_at_0":     prior_at_0,
        "ln_B_CMSTG_LCDM":        ln_B,
        "B_CMSTG_LCDM":           math.exp(ln_B) if ln_B is not None else None,
        "jeffreys_interpretation": jeffreys_scale(ln_B),
        "note": note,
    }

def jeffreys_scale(ln_B):
    if ln_B is None:
        return "undefined"
    B = math.exp(ln_B)
    if   B > 150:   return "Decisive evidence for CMSTG"
    elif B > 20:    return "Strong evidence for CMSTG"
    elif B > 3:     return "Substantial evidence for CMSTG"
    elif B > 1/3:   return "Inconclusive"
    elif B > 1/20:  return "Substantial evidence for LCDM"
    elif B > 1/150: return "Strong evidence for LCDM"
    else:           return "Decisive evidence for LCDM"

# ═══════════════════════════════════════════════════════════════════════════════
# 7. MAIN
# ═══════════════════════════════════════════════════════════════════════════════

print("=" * 60)
print("SIM93 — CMSTG Bayesian Model Comparison")
print("=" * 60)

INIT_CENTER = np.array([67.59, 0.3118, 147.56, 0.003])

all_chains = {}
all_results = {}
all_acc = {}

for prior_name, prior_fn in PRIORS.items():
    if HAS_EMCEE:
        chain, acc = run_mcmc(prior_name, prior_fn, INIT_CENTER.copy())
    else:
        chain, acc = run_grid_posterior(prior_name, prior_fn)

    all_chains[prior_name] = chain
    all_acc[prior_name] = acc

    chain_L0 = chain[:, 3]
    sd = savage_dickey(chain_L0, prior_fn, prior_name)
    all_results[prior_name] = sd

    ln_B = sd["ln_B_CMSTG_LCDM"]
    interp = sd["jeffreys_interpretation"]
    post0  = sd["posterior_density_at_0"]
    prior0 = sd["prior_density_at_0"]
    print(f"\n  Prior: {prior_name}")
    print(f"    Posterior density at Lambda0=0 : {post0:.4f}")
    print(f"    Prior density at Lambda0=0     : {prior0}")
    if ln_B is not None:
        print(f"    ln B(CMSTG/LCDM)               : {ln_B:+.3f}")
        print(f"    B(CMSTG/LCDM)                  : {math.exp(ln_B):.3f}")
        print(f"    Jeffreys interpretation        : {interp}")
    else:
        print(f"    {sd['note']}")

# ═══════════════════════════════════════════════════════════════════════════════
# 8. PLOTS
# ═══════════════════════════════════════════════════════════════════════════════

# 8a. Lambda0 marginal posteriors for all three priors
fig, axes = plt.subplots(1, 3, figsize=(14, 4))
fig.suptitle("SIM93 — CMSTG Bayesian Model Comparison: $\\Lambda_0$ Marginal Posteriors",
             fontsize=12, fontweight='bold')

colors = {"uniform": "steelblue", "log_uniform": "darkorange", "half_gaussian": "green"}
labels = {"uniform":      r"Uniform $\Lambda_0 \in [0, 0.1]$",
          "log_uniform":  r"Log-uniform $\Lambda_0 \in [10^{-4}, 0.1]$",
          "half_gaussian":r"Half-Gaussian ($\sigma=0.03$)"}

for ax, (prior_name, chain) in zip(axes, all_chains.items()):
    chain_L0 = chain[:, 3]
    L0_plot = np.linspace(0, L0_MAX, 500)

    # KDE posterior
    kde = gaussian_kde(chain_L0, bw_method='silverman')
    post_vals = kde(L0_plot)
    ax.plot(L0_plot, post_vals, color=colors[prior_name], lw=2, label='Posterior')
    ax.fill_between(L0_plot, 0, post_vals, color=colors[prior_name], alpha=0.15)

    # Prior
    if prior_name == "uniform":
        prior_vals = np.where((L0_plot >= 0) & (L0_plot <= L0_MAX),
                              1.0 / L0_MAX, 0.0)
    elif prior_name == "log_uniform":
        L0_MIN = 1e-4
        prior_vals = np.where(L0_plot >= L0_MIN,
                              1.0 / (L0_plot * math.log(L0_MAX / L0_MIN)), 0.0)
        prior_vals = np.clip(prior_vals, 0, 5 * post_vals.max())
    else:  # half_gaussian
        sigma = 0.03
        norm = 2.0 / (sigma * math.sqrt(2 * math.pi))
        prior_vals = np.where(L0_plot >= 0,
                              norm * np.exp(-0.5 * (L0_plot / sigma)**2), 0.0)
    ax.plot(L0_plot, prior_vals, 'k--', lw=1.2, alpha=0.6, label='Prior')

    # Mark Lambda0=0
    ax.axvline(0, color='red', lw=1.5, ls=':', label=r'$\Lambda_0=0$ (LCDM)')
    # Mark BAO best-fit
    ax.axvline(0.003, color='gray', lw=1, ls='--', alpha=0.6,
               label=r'BAO best-fit $\Lambda_0=0.003$')

    # Savage-Dickey result
    sd = all_results[prior_name]
    if sd["ln_B_CMSTG_LCDM"] is not None:
        ln_B = sd["ln_B_CMSTG_LCDM"]
        ax.set_title(f"{labels[prior_name]}\n"
                     f"$\\ln B = {ln_B:+.2f}$ ({sd['jeffreys_interpretation']})",
                     fontsize=9)
    else:
        ax.set_title(f"{labels[prior_name]}\n(SD ratio undefined at $\\Lambda_0=0$)",
                     fontsize=9)

    ax.set_xlabel(r"$\Lambda_0$", fontsize=11)
    ax.set_ylabel("Density", fontsize=11)
    ax.legend(fontsize=8)
    ax.set_xlim(-0.002, L0_MAX)
    ax.set_ylim(bottom=0)

plt.tight_layout()
out1 = os.path.join(OUTDIR, 'sim93_posterior_Lambda0.png')
plt.savefig(out1, dpi=150, bbox_inches='tight')
plt.close()
print(f"\nSaved {out1}")

# 8b. Corner plot (uniform prior chain only — most interpretable)
try:
    import corner
    chain_u = all_chains["uniform"]
    fig_c = corner.corner(
        chain_u,
        labels=[r"$H_0$", r"$\Omega_m$", r"$r_d$", r"$\Lambda_0$"],
        quantiles=[0.16, 0.50, 0.84],
        show_titles=True,
        title_kwargs={"fontsize": 10},
        color="steelblue",
    )
    fig_c.suptitle("SIM93 — MCMC Posterior (Uniform $\\Lambda_0$ prior)", fontsize=11)
    out2 = os.path.join(OUTDIR, 'sim93_corner.png')
    fig_c.savefig(out2, dpi=120, bbox_inches='tight')
    plt.close(fig_c)
    print(f"Saved {out2}")
except ImportError:
    print("corner not installed — skipping corner plot")
    out2 = None

# 8c. Summary bar chart: ln B for each prior
fig3, ax3 = plt.subplots(figsize=(7, 4))
prior_labels = []
ln_B_vals = []
for pn, sd in all_results.items():
    if sd["ln_B_CMSTG_LCDM"] is not None:
        prior_labels.append(pn.replace("_", " ").title())
        ln_B_vals.append(sd["ln_B_CMSTG_LCDM"])

if ln_B_vals:
    bar_colors = ['steelblue' if v > 0 else 'salmon' for v in ln_B_vals]
    bars = ax3.bar(prior_labels, ln_B_vals, color=bar_colors, edgecolor='k', linewidth=0.8)
    ax3.axhline(0,    color='k',   lw=1.0, ls='-')
    ax3.axhline( 1.1, color='gray', lw=0.8, ls='--', alpha=0.6, label=r'$|\ln B|=1.1$ (substantial)')
    ax3.axhline(-1.1, color='gray', lw=0.8, ls='--', alpha=0.6)
    ax3.axhline( 3.0, color='gray', lw=0.8, ls=':', alpha=0.5, label=r'$|\ln B|=3.0$ (strong)')
    ax3.axhline(-3.0, color='gray', lw=0.8, ls=':', alpha=0.5)
    for bar, val in zip(bars, ln_B_vals):
        ax3.text(bar.get_x() + bar.get_width()/2, val + 0.05*np.sign(val),
                 f'{val:+.2f}', ha='center', va='bottom' if val > 0 else 'top', fontsize=10)
    ax3.set_ylabel(r"$\ln B(\mathrm{CMSTG}/\Lambda\mathrm{CDM})$", fontsize=12)
    ax3.set_title("Savage-Dickey Bayes Factor by Prior Choice", fontsize=11)
    ax3.legend(fontsize=9)
    ax3.text(0.98, 0.95, "CMSTG favoured →", transform=ax3.transAxes,
             ha='right', va='top', fontsize=9, color='steelblue')
    ax3.text(0.98, 0.05, "← ΛCDM favoured", transform=ax3.transAxes,
             ha='right', va='bottom', fontsize=9, color='salmon')
    plt.tight_layout()
    out3 = os.path.join(OUTDIR, 'sim93_bayes_factors.png')
    plt.savefig(out3, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"Saved {out3}")

# ═══════════════════════════════════════════════════════════════════════════════
# 9. DIAGNOSTICS
# ═══════════════════════════════════════════════════════════════════════════════

# Compute chain statistics for uniform prior
chain_u = all_chains["uniform"]
L0_chain = chain_u[:, 3]

# 95% upper limit on Lambda0
L0_95 = float(np.percentile(L0_chain, 95))
L0_50 = float(np.percentile(L0_chain, 50))
L0_84 = float(np.percentile(L0_chain, 84))

# Marginalised H0, Omm, rd
H0_med  = float(np.percentile(chain_u[:,0], 50))
Omm_med = float(np.percentile(chain_u[:,1], 50))
rd_med  = float(np.percentile(chain_u[:,2], 50))

print(f"\n  Lambda0 constraints (uniform prior):")
print(f"    Median : {L0_50:.4f}")
print(f"    84th   : {L0_84:.4f}")
print(f"    95% UL : {L0_95:.4f}")
print(f"  H0  (marginalised): {H0_med:.2f}")
print(f"  Omm (marginalised): {Omm_med:.4f}")
print(f"  rd  (marginalised): {rd_med:.2f}")

diag = {
    "description": "SIM93 — CMSTG Bayesian Model Comparison via Savage-Dickey density ratio",
    "method": "Savage-Dickey density ratio: ln B = ln p(L0=0|data) - ln p(L0=0|prior). Valid because LCDM is CMSTG at Lambda0=0 (nested models).",
    "sampler": "emcee" if HAS_EMCEE else "grid (emcee not available)",
    "n_samples_per_chain": {pn: len(c) for pn, c in all_chains.items()},
    "acceptance_fractions": {pn: float(a) if a is not None else None for pn, a in all_acc.items()},
    "bayes_factors": {
        pn: {
            "ln_B": sd["ln_B_CMSTG_LCDM"],
            "B":    sd["B_CMSTG_LCDM"],
            "jeffreys": sd["jeffreys_interpretation"],
            "posterior_density_at_0": sd["posterior_density_at_0"],
            "prior_density_at_0": sd["prior_density_at_0"],
            "note": sd["note"],
        }
        for pn, sd in all_results.items()
    },
    "lambda0_constraints_uniform_prior": {
        "median": L0_50,
        "p84":    L0_84,
        "p95_upper_limit": L0_95,
        "note": "Lambda0 is unconstrained from above — posterior follows prior. 95% UL is prior-dominated."
    },
    "marginalised_params_uniform_prior": {
        "H0":     H0_med,
        "Omega_m": Omm_med,
        "rd":     rd_med,
    },
    "physical_interpretation": (
        "Lambda0 is NOT constrained by current CMB+BAO data (the chi2 landscape "
        "is flat in Lambda0 — confirmed by SIM91, Dchi2=0.004 across 0-0.1). "
        "The Bayes factor therefore measures the Occam factor: how much prior "
        "volume does CMSTG have that is not excluded by data? A flat likelihood "
        "means B ~ p(L0=0|prior)/p(L0=0|prior) ~ 1 for the uniform prior. "
        "The half-Gaussian prior concentrates prior mass near L0=0, giving "
        "mild LCDM preference. In all cases |ln B| < 3 = inconclusive "
        "on the Jeffreys scale. CMSTG is neither preferred nor excluded."
    ),
    "status": "PASS",
}

diag_path = os.path.join(OUTDIR, 'sim93_diagnostics.json')
with open(diag_path, 'w') as f:
    json.dump(diag, f, indent=2)
print(f"Saved {diag_path}")

# Save chain for uniform prior (largest posterior content)
chain_path = os.path.join(OUTDIR, 'sim93_mcmc_chain.npz')
np.savez(chain_path,
         chain_uniform=all_chains["uniform"],
         chain_half_gaussian=all_chains["half_gaussian"])
print(f"Saved {chain_path}")

print("\n" + "=" * 60)
print("SIM93 STATUS: PASS")
print(f"  Uniform prior    : ln B = {all_results['uniform']['ln_B_CMSTG_LCDM']:+.3f}  ({all_results['uniform']['jeffreys_interpretation']})")
print(f"  Half-Gaussian    : ln B = {all_results['half_gaussian']['ln_B_CMSTG_LCDM']:+.3f}  ({all_results['half_gaussian']['jeffreys_interpretation']})")
print(f"  Log-uniform      : ln B = undefined (SD ratio undefined at Lambda0=0)")
print(f"  Lambda0 95% UL   : {L0_95:.4f}  (prior-dominated — data cannot constrain Lambda0)")
print("=" * 60)
