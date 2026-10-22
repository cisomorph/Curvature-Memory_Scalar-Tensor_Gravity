"""
SIM121-C — CMSTG Phase 2: Joint DESI+Planck MCMC Parameter Estimation
=====================================================================
Performs a proper Bayesian parameter estimation for CMSTG Phase 2 against
the combination of Planck CMB acoustic scale and DESI BAO data.

SIM121-B established that the CMSTG F₀ normalisation is the dominant CMB
observable. This simulation treats (F₀, w₀, wₐ) as the primary CMSTG
parameters, with H₀ derived self-consistently from the Planck θ_*
constraint at each MCMC step.

Parameter space:
  F₀   = ½ + Λ₀Ψ₀²  (CMSTG effective coupling, Ψ frozen)
  w₀               (DE EOS today)           — phenomenological (CPL), not derived
  wₐ               (DE EOS running)         — phenomenological (CPL), not derived
  H₀               (derived from CMB θ_* at each step)

Fixed (Planck 2018 physical densities):
  Ω_m h²  = 0.1430,  Ω_b h² = 0.02237,  Ω_r h² = 4.18×10⁻⁵

Likelihood:
  ln L = −½ χ²_θ  −½ χ²_DESI
  χ²_θ    = [(100θ_*,CMSTG − 1.04110) / 0.00031]²   (≈0: H₀ is solved from it)
  χ²_DESI = official DESI Gaussian BAO likelihood (D_M/r_d, D_H/r_d, D_V/r_d,
            full covariance; DR1 by default, DR2 with --release DR2)

AUDIT 2026-10-22 (see AUDIT_NOTES.md). Changes from the published version:
  1. DESI data: the hardcoded "H(z)" table (81.7, 97.9, 110.7, 128.1, 156.4,
     240.8) is not DESI data. Replaced with the official DESI BAO likelihood.
  2. r_s was integrated only to z = 1e4 (r_d = 122.6 instead of 147.1 Mpc).
     Now integrated to z = ∞.
  3. θ_* used r_s(z_drag); now r_s(z_*). θ_obs was Planck 100θ_MC (+BAO);
     now 100θ_* = 1.04110 ± 0.00031 (TT,TE,EE+lowE+lensing).
  4. Flatness closure: Ω_DE = 1 − Ω_m − Ω_r gave H(z=0) = H₀·√(0.5/F₀) ≠ H₀.
     Now Ω_DE = 2F₀ − Ω_m − Ω_r so the H₀ parameter is the Hubble constant.
  5. "DESI tension" = √(χ²/N) is not a significance. Now reported: χ², dof,
     p-value, Gaussian-equivalent σ, and Δχ² against ΛCDM and w₀wₐCDM run
     through the same pipeline.
  6. Prior w₀ + wₐ < 0 added (DE must not dominate at early times).
  7. H₀ is stored as an emcee blob instead of being re-derived afterwards.

Pass criteria (unchanged in intent):
  1. χ²_DESI / N_DESI < 2  at best fit
  2. χ²_θ < 4
  3. w₀ ∈ (−1.3, −0.7)
  4. F₀ ∈ [0.490, 0.560]
"""

import numpy as np
import emcee, json, os, sys, argparse, warnings
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.optimize import minimize
warnings.filterwarnings('ignore')

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', '..', 'common'))
import cmstg_obs as obs

ap = argparse.ArgumentParser()
ap.add_argument('--release', default='DR1', choices=['DR1', 'DR2'])
ap.add_argument('--nsteps', type=int, default=3000)
ap.add_argument('--nburn', type=int, default=800)
args, _ = ap.parse_known_args()

OUT = os.path.join(HERE, '..', 'Outputs')
if args.release != 'DR1':
    OUT = os.path.join(OUT, args.release)
os.makedirs(OUT, exist_ok=True)

plt.rcParams.update({'font.family': 'serif', 'font.size': 11,
                     'axes.labelsize': 12, 'legend.fontsize': 10})

# ─────────────────────────────────────────────────────────────────────────────
# PHYSICAL CONSTANTS AND PLANCK 2018 PRIORS
# ─────────────────────────────────────────────────────────────────────────────
omh2_m        = obs.OMEGA_M
omh2_b        = obs.OMEGA_B
omh2_r        = obs.OMEGA_R
theta_obs     = obs.THETA_STAR_100
theta_obs_err = obs.THETA_STAR_100_ERR

DESI = obs.DesiBAO(args.release)

# SIM113 CMSTG reference
SIM113_w0 = -0.973;  SIM113_wa = -0.41;  SIM113_F0 = 0.52059

# Published (pre-audit) values, for the before/after record
PUBLISHED = dict(map=dict(F0=0.5599845, w0=-0.6002629, wa=0.4862087, H0_label=75.957),
                 chi2_DESI=41.4917, DESI_tension_sigma_old_metric=2.6297, verdict='PARTIAL')

# ─────────────────────────────────────────────────────────────────────────────
# CMSTG HUBBLE RATE (frozen Ψ, CPL DE, exact F₀ normalisation)
# ─────────────────────────────────────────────────────────────────────────────

def make_H(F0, w0, wa, H0):
    """
    H(z) [km/s/Mpc] with frozen Ψ ≈ Ψ₀:  6F₀H² = ρ_m + ρ_r + ρ_DE  (GR: 3H² = ρ),
    i.e. H² = (0.5/F₀)·H100²·(ω_m a⁻³ + ω_r a⁻⁴) + DE, with the DE amount fixed by
    H(z=0) = H₀ (flatness).
    """
    h = H0 / 100.0
    Om, Or = omh2_m / h**2, omh2_r / h**2
    ODE = 2.0 * F0 - Om - Or               # GR-normalised; gives E(0) = 1 exactly
    g = 0.5 / F0

    def H(z):
        z = np.asarray(z, dtype=float)
        a = 1.0 / (1.0 + z)
        fDE = a**(-3.0 * (1.0 + w0 + wa)) * np.exp(-3.0 * wa * (1.0 - a))
        return H0 * np.sqrt(g * (Om / a**3 + Or / a**4 + ODE * fDE))
    return H


def H0_from_theta(F0, w0, wa, H0_lo=50.0, H0_hi=90.0):
    """H₀ such that 100θ_*(F₀,w₀,wₐ,H₀) = θ_obs, or nan."""
    h = obs.h_from_theta(lambda hh: make_H(F0, w0, wa, 100.0 * hh),
                         h_lo=H0_lo / 100.0, h_hi=H0_hi / 100.0)
    return 100.0 * h if np.isfinite(h) else np.nan


def evaluate(F0, w0, wa):
    """Full evaluation at one parameter point (H₀ from θ_*)."""
    H0 = H0_from_theta(F0, w0, wa)
    if np.isnan(H0) or H0 < 55 or H0 > 85:
        return None
    H = make_H(F0, w0, wa, H0)
    rd = obs.r_drag(H)
    th = obs.theta_star_100(H)
    return dict(H0=H0, H=H, rd=rd, theta=th, chi2_theta=obs.chi2_theta(th),
                chi2_DESI=DESI.chi2(H, rd))

# ─────────────────────────────────────────────────────────────────────────────
# LIKELIHOOD
# ─────────────────────────────────────────────────────────────────────────────

def log_likelihood(theta_params):
    F0, w0, wa = theta_params
    e = evaluate(F0, w0, wa)
    if e is None:
        return -np.inf, np.nan
    return -0.5 * (e['chi2_theta'] + e['chi2_DESI']), e['H0']

def log_prior(theta_params):
    F0, w0, wa = theta_params
    if not (0.490 <= F0 <= 0.560):  return -np.inf
    if not (-1.40 <= w0 <= -0.60):  return -np.inf
    if not (-2.50 <= wa <=  1.50):  return -np.inf
    if w0 + wa >= 0.0:              return -np.inf
    return 0.0

def log_prob(theta_params):
    lp = log_prior(theta_params)
    if not np.isfinite(lp):
        return -np.inf, np.nan
    ll, H0 = log_likelihood(theta_params)
    return lp + ll, H0


def best_fit(free, x0, fixed):
    """Minimise χ²_tot over the named free parameters (others fixed)."""
    names = ['F0', 'w0', 'wa']
    def f(x):
        p = dict(fixed); p.update(zip(free, x))
        v = [p[n] for n in names]
        if not np.isfinite(log_prior(v)):
            return 1e10
        ll, _ = log_likelihood(v)
        return -2 * ll if np.isfinite(ll) else 1e10
    if not free:
        return f([]), dict(fixed)
    r = minimize(f, x0, method='Nelder-Mead',
                 options=dict(xatol=1e-6, fatol=1e-6, maxiter=4000))
    p = dict(fixed); p.update(zip(free, r.x))
    return r.fun, p


if __name__ == '__main__':
    # ─────────────────────────────────────────────────────────────────────────
    # PRE-SCAN: find good starting region before MCMC
    # ─────────────────────────────────────────────────────────────────────────
    print("=" * 70)
    print(f"SIM121-C — CMSTG Phase 2: Joint DESI+Planck MCMC   [DESI {args.release}, audited]")
    print("=" * 70)

    print("\nPart A: Coarse grid scan to find starting region...")
    F0_grid  = np.linspace(0.495, 0.555, 13)
    w0_grid  = np.linspace(-1.30, -0.62, 15)
    wa_grid  = np.linspace(-2.40,  0.50, 15)

    scan_results = []
    for F0 in F0_grid:
        for w0 in w0_grid:
            for wa in wa_grid:
                if not np.isfinite(log_prior([F0, w0, wa])):
                    continue
                ll, H0 = log_likelihood([F0, w0, wa])
                if np.isfinite(ll):
                    scan_results.append((-2 * ll, F0, w0, wa, H0))

    scan_results.sort()
    print(f"  Valid points: {len(scan_results)}")
    print(f"\n  Top 5 grid minima:")
    print(f"  {'chi2':>8}  {'F0':>8}  {'w0':>8}  {'wa':>8}  {'H0':>8}")
    for chi2, F0, w0, wa, H0 in scan_results[:5]:
        print(f"  {chi2:>8.3f}  {F0:>8.4f}  {w0:>8.4f}  {wa:>8.4f}  {H0:>8.3f}")
    if not scan_results:
        print("ERROR: No valid grid points found.")
        sys.exit(1)
    best_params = list(scan_results[0][1:4])

    # ─────────────────────────────────────────────────────────────────────────
    # MCMC WITH EMCEE
    # ─────────────────────────────────────────────────────────────────────────
    print(f"\nPart B: MCMC (emcee) around best-fit region...")
    ndim, nwalkers = 3, 32
    nsteps, nburn  = args.nsteps, args.nburn

    rng       = np.random.default_rng(42)
    p0_center = np.array(best_params)
    scales    = np.array([0.003, 0.05, 0.15])
    p0        = p0_center + scales * rng.normal(size=(nwalkers, ndim))
    p0[:,0] = np.clip(p0[:,0], 0.491, 0.559)
    p0[:,1] = np.clip(p0[:,1], -1.39, -0.61)
    p0[:,2] = np.clip(p0[:,2], -2.49,  1.49)
    for k in range(nwalkers):                       # keep every walker inside the prior
        while not np.isfinite(log_prior(p0[k])):
            p0[k] = p0_center + scales * rng.normal(size=ndim)

    sampler = emcee.EnsembleSampler(nwalkers, ndim, log_prob,
                                    blobs_dtype=[('H0', float)])
    print(f"  Burning in ({nburn} steps)...")
    state = sampler.run_mcmc(p0, nburn, progress=False)
    sampler.reset()
    print(f"  Sampling ({nsteps} steps × {nwalkers} walkers)...")
    sampler.run_mcmc(state, nsteps, progress=False)

    flat_chain = sampler.get_chain(flat=True)
    log_probs  = sampler.get_log_prob(flat=True)
    H0_chain   = sampler.get_blobs(flat=True)['H0']
    print(f"  Acceptance fraction: {np.mean(sampler.acceptance_fraction):.3f}")
    try:
        tau = sampler.get_autocorr_time(quiet=True)
        print(f"  Autocorrelation times: {np.round(tau, 1)}  (chain length {nsteps})")
    except Exception:
        tau = np.full(ndim, np.nan)

    idx_map = np.argmax(log_probs)
    F0_map, w0_map, wa_map = flat_chain[idx_map]
    e_map = evaluate(F0_map, w0_map, wa_map)
    H0_map, theta_map = e_map['H0'], e_map['theta']

    valid     = np.isfinite(H0_chain)
    flat_full = np.column_stack([flat_chain[valid], H0_chain[valid]])
    pct = np.percentile(flat_full, [16, 50, 84], axis=0)

    print(f"\n── Posterior Summary ──")
    print(f"  {'Param':>6}  {'MAP':>10}  {'16%':>10}  {'50%':>10}  {'84%':>10}")
    map_vals = [F0_map, w0_map, wa_map, H0_map]
    for lbl, mv, p16, p50, p84 in zip(['F0','w0','wa','H0'], map_vals, pct[0], pct[1], pct[2]):
        print(f"  {lbl:>6}  {mv:>10.4f}  {p16:>10.4f}  {p50:>10.4f}  {p84:>10.4f}")

    # ─────────────────────────────────────────────────────────────────────────
    # DIAGNOSTICS — every model through the same pipeline
    # ─────────────────────────────────────────────────────────────────────────
    fits = {}
    for key, free, x0, fixed in [
            ('LCDM',             [],                 [],                       dict(F0=0.5, w0=-1.0, wa=0.0)),
            ('Phase1_frozen_F0', [],                 [],                       dict(F0=0.5 + 0.003 * 2.62**2, w0=-1.0, wa=0.0)),
            ('SIM113_point',     [],                 [],                       dict(F0=SIM113_F0, w0=SIM113_w0, wa=SIM113_wa)),
            ('F0_only',          ['F0'],             [0.51],                   dict(w0=-1.0, wa=0.0)),
            ('w0waCDM',          ['w0', 'wa'],       [w0_map, wa_map],         dict(F0=0.5)),
            ('CMSTG_full',       ['F0', 'w0', 'wa'], [F0_map, w0_map, wa_map], dict())]:
        c2, p = best_fit(free, x0, fixed)
        e = evaluate(p['F0'], p['w0'], p['wa'])
        if e is None:
            fits[key] = dict(params=p, valid=False)
            continue
        s = obs.summarize_fit(e['chi2_DESI'], DESI.n, len(free), e['chi2_theta'], label=key)
        s.update(params={k: float(v) for k, v in p.items()}, H0=float(e['H0']),
                 r_drag=float(e['rd']), theta_star_100=float(e['theta']), valid=True)
        fits[key] = s

    ref = fits['LCDM']
    for key, s in fits.items():
        if s.get('valid') and key != 'LCDM':
            k = s['n_params']
            d = ref['chi2_DESI'] + ref['chi2_theta'] - s['chi2_DESI'] - s['chi2_theta']
            s['delta_chi2_vs_LCDM'] = float(d)          # > 0: better than ΛCDM
            s['preference_over_LCDM_sigma'] = float(obs.delta_chi2_to_sigma(d, k)) if k > 0 else None
    w0wa = fits['w0waCDM']; full = fits['CMSTG_full']
    d_F0 = w0wa['chi2_DESI'] + w0wa['chi2_theta'] - full['chi2_DESI'] - full['chi2_theta']
    full['delta_chi2_vs_w0waCDM'] = float(d_F0)
    full['F0_preference_sigma'] = float(obs.delta_chi2_to_sigma(d_F0, 1))

    print(f"\n── Fit quality (DESI {args.release} + θ_*) ──")
    print(f"  {'Model':20} {'k':>2} {'χ²_DESI':>8} {'dof':>4} {'p':>9} {'GOF σ':>6} "
          f"{'Δχ² vs ΛCDM':>12} {'H0':>7}  params")
    for key, s in fits.items():
        if not s.get('valid'):
            print(f"  {key:20} invalid (no θ* solution)")
            continue
        d = s.get('delta_chi2_vs_LCDM', 0.0)
        print(f"  {key:20} {s['n_params']:>2} {s['chi2_DESI']:>8.2f} {s['dof']:>4} {s['p_value']:>9.2e} "
              f"{s['gof_sigma']:>6.2f} {d:>+12.2f} {s['H0']:>7.2f}  "
              + ", ".join(f"{k}={v:.4f}" for k, v in s['params'].items()))
    print(f"\n  F₀ adds Δχ² = {d_F0:+.2f} beyond w₀wₐCDM → {full['F0_preference_sigma']:.2f}σ for F₀ ≠ ½")

    print(f"\n── DESI points at MAP ──")
    print(f"  {'z':>6} {'quantity':>11} {'data':>8} {'σ':>6} {'MAP':>8} {'pull':>6}")
    for row in DESI.table(e_map['H'], e_map['rd']):
        print(f"  {row['z']:>6.3f} {row['quantity']:>11} {row['data']:>8.3f} {row['sigma']:>6.3f} "
              f"{row['model']:>8.3f} {row['pull']:>+6.2f}")

    # ─────────────────────────────────────────────────────────────────────────
    # VERDICT
    # ─────────────────────────────────────────────────────────────────────────
    chi2_DESI_map  = e_map['chi2_DESI']
    chi2_theta_map = e_map['chi2_theta']
    pass_chi2   = chi2_DESI_map / DESI.n < 2.0
    pass_theta  = chi2_theta_map < 4.0
    pass_w0     = -1.3 <= w0_map <= -0.7
    pass_F0     = 0.490 <= F0_map <= 0.560
    gof_map     = obs.summarize_fit(chi2_DESI_map, DESI.n, 3, chi2_theta_map)

    PASS = pass_chi2 and pass_theta and pass_w0 and pass_F0
    PARTIAL = (pass_theta or pass_chi2) and not PASS
    verdict = 'PASS' if PASS else 'PARTIAL' if PARTIAL else 'FAIL'

    print(f"\n{'='*70}")
    print("SIM121-C RESULT:")
    print(f"  MAP: F₀ = {F0_map:.5f}, w₀ = {w0_map:.4f}, wₐ = {wa_map:.4f}, H₀ = {H0_map:.3f} (CMB-derived)")
    print(f"  χ²_DESI = {chi2_DESI_map:.3f} for N = {DESI.n}, dof = {gof_map['dof']}: "
          f"p = {gof_map['p_value']:.3g} → {gof_map['gof_sigma']:.2f}σ   "
          f"[old metric √(χ²/N) = {gof_map['rms_pull_old_metric']:.2f}]")
    print(f"  criteria: χ²/N<2 {'PASS' if pass_chi2 else 'FAIL'} | θ_* {'PASS' if pass_theta else 'FAIL'} | "
          f"w₀ {'PASS' if pass_w0 else 'FAIL'} | F₀ {'PASS' if pass_F0 else 'FAIL'}")
    print(f"  VERDICT: {verdict}")
    print(f"  Published (pre-audit): F₀={PUBLISHED['map']['F0']:.4f}, 'H₀'={PUBLISHED['map']['H0_label']:.2f}, "
          f"χ²={PUBLISHED['chi2_DESI']:.2f} on the non-DESI table, "
          f"'tension'={PUBLISHED['DESI_tension_sigma_old_metric']:.2f}")
    print(f"{'='*70}")

    # ─────────────────────────────────────────────────────────────────────────
    # SAVE RESULTS JSON
    # ─────────────────────────────────────────────────────────────────────────
    out_data = {
        'verdict': verdict,
        'desi_release': args.release,
        'map': {'F0': float(F0_map), 'w0': float(w0_map), 'wa': float(wa_map), 'H0': float(H0_map)},
        'posterior_50pct': dict(zip(['F0','w0','wa','H0'], pct[1].tolist())),
        'posterior_16pct': dict(zip(['F0','w0','wa','H0'], pct[0].tolist())),
        'posterior_84pct': dict(zip(['F0','w0','wa','H0'], pct[2].tolist())),
        'chi2_DESI_map': float(chi2_DESI_map),
        'N_DESI': DESI.n,
        'gof_map': gof_map,
        'chi2_theta_map': float(chi2_theta_map),
        'theta_star_map': float(theta_map),
        'theta_star_obs': float(theta_obs),
        'r_drag_map': float(e_map['rd']),
        'fits_same_pipeline': fits,
        'desi_points_map': DESI.table(e_map['H'], e_map['rd']),
        'nwalkers': nwalkers, 'nsteps': nsteps, 'nburn': nburn,
        'acceptance_fraction': float(np.mean(sampler.acceptance_fraction)),
        'autocorr_time': [float(t) for t in tau],
        'n_valid_samples': int(valid.sum()),
        'published_pre_audit': PUBLISHED,
        'audit_note': 'Audited 2026-10-22: official DESI likelihood, r_s to infinity, '
                      'theta*=r_s(z*)/D_M(z*), H0 closure fix, proper significance. See AUDIT_NOTES.md.',
        'phenomenological_inputs': ['w0', 'wa (CPL dark energy, not derived from the action)'],
    }
    with open(os.path.join(OUT, 'sim121c_results.json'), 'w') as f:
        json.dump(out_data, f, indent=2)

    # ─────────────────────────────────────────────────────────────────────────
    # FIGURE 1: Corner-style marginals
    # ─────────────────────────────────────────────────────────────────────────
    print(f"\nGenerating figures...")
    param_names = [r'$F_0$', r'$w_0$', r'$w_a$', r'$H_0$ [km/s/Mpc]']
    param_ranges = [(0.490, 0.560), (-1.40, -0.60), (-2.5, 1.0),
                    (min(55, np.percentile(flat_full[:, 3], 0.5)),
                     max(80, np.percentile(flat_full[:, 3], 99.5)))]

    fig, axes = plt.subplots(4, 4, figsize=(13, 12))
    for i in range(4):
        for j in range(4):
            ax = axes[i, j]
            if j > i:
                ax.axis('off'); continue
            xi, yi = flat_full[:, j], flat_full[:, i]
            if i == j:
                ax.hist(xi, bins=40, range=param_ranges[i], color='#4393c3', alpha=0.8,
                        density=True, edgecolor='k', lw=0.3)
                ax.axvline(map_vals[i], color='#d6604d', lw=2, ls='--')
                ax.axvline(pct[1, i], color='k', lw=1, ls=':')
                ax.set_xlim(param_ranges[i]); ax.set_yticks([])
            else:
                ax.hist2d(xi, yi, bins=35, cmap='Blues',
                          range=[param_ranges[j], param_ranges[i]], density=True)
                ax.plot(map_vals[j], map_vals[i], 'r*', ms=10, zorder=5)
            if i == 3: ax.set_xlabel(param_names[j], fontsize=10)
            if j == 0 and i > 0: ax.set_ylabel(param_names[i], fontsize=10)
            ax.tick_params(labelsize=8)
    fig.suptitle(f'SIM121-C (audited) — CMSTG joint posterior: θ* + DESI {args.release}',
                 fontsize=13, y=1.00)
    fig.tight_layout()
    for ext in ['pdf', 'png']:
        fig.savefig(os.path.join(OUT, f'sim121c_corner.{ext}'), dpi=130, bbox_inches='tight')
    plt.close(fig)
    print("  Saved sim121c_corner")

    # ─────────────────────────────────────────────────────────────────────────
    # FIGURE 2: w₀-wₐ plane + DESI distances relative to ΛCDM
    # ─────────────────────────────────────────────────────────────────────────
    fig, axes = plt.subplots(1, 2, figsize=(13, 5))
    ax = axes[0]
    w0g = np.linspace(-1.35, -0.62, 45)
    wag = np.linspace(-2.45, 0.6, 45)
    chi2_grid = np.full((len(wag), len(w0g)), np.inf)
    for iw, w0v in enumerate(w0g):
        for iwa, wav in enumerate(wag):
            if w0v + wav >= 0:
                continue
            ll_v, _ = log_likelihood([F0_map, w0v, wav])
            if np.isfinite(ll_v):
                chi2_grid[iwa, iw] = -2 * ll_v
    dchi2 = chi2_grid - chi2_grid[np.isfinite(chi2_grid)].min()
    W0G, WAG = np.meshgrid(w0g, wag)
    ax.contourf(W0G, WAG, dchi2, levels=[0, 2.30, 6.18, 11.8],
                colors=['#d9f0d3', '#a6dba0', '#5aae61'], alpha=0.7)
    ax.contour(W0G, WAG, dchi2, levels=[2.30, 6.18], colors='k', linewidths=1.2)
    ax.plot(w0_map, wa_map, 'r*', ms=14, zorder=5, label=f'MAP ({w0_map:.3f}, {wa_map:.3f})')
    ax.plot(SIM113_w0, SIM113_wa, 'g^', ms=10, zorder=5, label='SIM113 (−0.973, −0.41)')
    ax.plot(-1.0, 0.0, 'bs', ms=10, zorder=5, label=r'$\Lambda$CDM')
    ax.set_xlabel(r'$w_0$'); ax.set_ylabel(r'$w_a$')
    ax.set_title(rf'$w_0$–$w_a$ ($F_0$ = MAP {F0_map:.4f})')
    ax.legend(fontsize=9)
    ax.axhline(0, color='gray', ls=':', lw=0.8); ax.axvline(-1.0, color='gray', ls=':', lw=0.8)

    ax = axes[1]
    H_l = make_H(0.5, -1.0, 0.0, fits['LCDM']['H0']); rd_l = obs.r_drag(H_l)
    m_l = DESI.model(H_l, rd_l)
    zz = np.linspace(0.05, 2.5, 120)

    def curves(Hf, rd):
        DMg = np.array([obs.comoving_distance(Hf, z, n=400) for z in zz])
        DH = obs.C_KMS / Hf(zz)
        return {'DM_over_rs': DMg / rd, 'DH_over_rs': DH / rd,
                'DV_over_rs': (zz * DMg**2 * DH)**(1 / 3) / rd}
    c_map, c_l = curves(e_map['H'], e_map['rd']), curves(H_l, rd_l)
    for col, kind, mk in [('#d6604d', 'DM_over_rs', 'o'), ('#2166ac', 'DH_over_rs', 's'),
                          ('#4dac26', 'DV_over_rs', '^')]:
        sel = [i for i, k in enumerate(DESI.kind) if k == kind]
        if sel:
            ax.errorbar(DESI.z[sel], DESI.value[sel] / m_l[sel], yerr=DESI.sigma[sel] / m_l[sel],
                        fmt=mk, color=col, ms=6, capsize=3,
                        label=f'DESI {args.release} {kind.split("_")[0]}')
        ax.plot(zz, c_map[kind] / c_l[kind], color=col, lw=1.8)
    ax.axhline(1.0, color='k', lw=1, ls=':')
    ax.set_xlabel('Redshift $z$'); ax.set_ylabel(r'$(D/r_d)\,/\,(D/r_d)_{\Lambda{\rm CDM}}$')
    ax.set_title('DESI BAO relative to ΛCDM (lines: CMSTG MAP)')
    ax.legend(fontsize=8)
    fig.suptitle(f'SIM121-C (audited) — θ* + DESI {args.release}', fontsize=13, y=1.01)
    fig.tight_layout()
    for ext in ['pdf', 'png']:
        fig.savefig(os.path.join(OUT, f'sim121c_contours.{ext}'), dpi=150, bbox_inches='tight')
    plt.close(fig)
    print("  Saved sim121c_contours")

    # ─────────────────────────────────────────────────────────────────────────
    # FIGURE 3: F₀ vs H₀ constraint from CMB
    # ─────────────────────────────────────────────────────────────────────────
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    ax = axes[0]
    F0_range = np.linspace(0.490, 0.560, 40)
    H0_vals = np.array([H0_from_theta(F0v, w0_map, wa_map) for F0v in F0_range])
    ok = np.isfinite(H0_vals)
    ax.plot(F0_range[ok], H0_vals[ok], color='#2166ac', lw=2.5,
            label=r'$H_0(F_0)$ from $\theta_*$ at MAP $(w_0, w_a)$')
    ax.axhline(73.0, color='#1a9641', ls='--', lw=1.8, label='SH0ES 73.0')
    ax.axhline(67.4, color='#d6604d', ls='--', lw=1.8, label='Planck 67.4')
    ax.axvline(0.5, color='gray', ls=':', lw=1, label='GR ($F_0$=½)')
    ax.plot(F0_map, H0_map, 'r*', ms=12, zorder=5)
    ax.set_xlabel(r'$F_0 = \frac{1}{2} + \Lambda_0 \Psi_0^2$'); ax.set_ylabel(r'$H_0$ [km/s/Mpc]')
    ax.set_title(r'CMB-derived $H_0$ vs $F_0$'); ax.legend(fontsize=9)
    ax = axes[1]
    ax.hist(flat_full[:, 3], bins=50, color='#4393c3', alpha=0.8, density=True, edgecolor='k', lw=0.3)
    ax.axvline(H0_map, color='#d6604d', lw=2, ls='--', label=f'MAP {H0_map:.2f}')
    ax.axvline(pct[1, 3], color='k', lw=1, ls=':', label=f'median {pct[1, 3]:.2f}')
    ax.axvline(73.0, color='#1a9641', lw=2, alpha=0.8, label='SH0ES')
    ax.axvline(67.4, color='#b2182b', lw=2, alpha=0.8, label='Planck')
    ax.set_xlabel(r'$H_0$ [km/s/Mpc]'); ax.set_ylabel('Probability density')
    ax.set_title(r'Marginal posterior: $H_0$'); ax.legend(fontsize=9)
    fig.suptitle(f'SIM121-C (audited) — CMB constraint on F₀ and H₀ (DESI {args.release})',
                 fontsize=13, y=1.01)
    fig.tight_layout()
    for ext in ['pdf', 'png']:
        fig.savefig(os.path.join(OUT, f'sim121c_F0_H0.{ext}'), dpi=150, bbox_inches='tight')
    plt.close(fig)
    print("  Saved sim121c_F0_H0")

    print(f"\nAll outputs in: {os.path.abspath(OUT)}")
    print(f"MCMC complete. {valid.sum()} valid samples.")
