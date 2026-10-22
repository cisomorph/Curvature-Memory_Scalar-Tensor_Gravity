"""
SIM122 — CMSTG Phase 3: Unlocked Λ₀ Scan
=========================================
Question (as published): is there a Λ₀ > 0 that simultaneously satisfies Planck θ_*
and DESI, with F₀ = ½ + Λ₀Ψ₀² (Ψ₀ = 2.62 frozen) and CPL dark energy (w₀, wₐ)?

Free parameters: Λ₀ ∈ [0.0001, 0.010], w₀ ∈ [−1.40, −0.60], wₐ ∈ [−2.50, 1.50]
(+ w₀ + wₐ < 0); H₀ derived from θ_*.  Fixed: Planck ω_m, ω_b, ω_r.

AUDIT 2026-10-22 (see AUDIT_NOTES.md). Same corrections as SIM121C (official DESI
likelihood, r_s to ∞, θ* = r_s(z*)/D_M(z*), Planck 100θ* = 1.04110 ± 0.00031, closure so
that H₀ is the Hubble constant, proper statistics). Model code is shared with SIM121C.

Key result of the audit: with Ψ frozen, F₀ (hence Λ₀) is exactly degenerate with H₀ for
θ* + BAO — a constant G rescaling with fixed ω_b/ω_γ rescales H(z), r_s and every distance
by the same factor once H₀ is re-solved from θ*. Λ₀ is therefore unconstrained, and any
improvement over ΛCDM comes entirely from the phenomenological (w₀, wₐ).

Pass criteria (as published, with criterion 6 restated):
  1. χ²_DESI/N < 2   2. χ²_θ < 4   3. wₐ < 0   4. F₀ ∈ [0.500, 0.560]   5. Λ₀ > 0
  6. better than ΛCDM in the same pipeline (published: "tension < 2.63σ")
  + audit flag: is Λ₀ actually constrained by the data?
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
sys.path.insert(0, os.path.join(HERE, '..', '..', 'SIM121C', 'Published_Scripts'))
import cmstg_obs as obs
import sim121c_mcmc as s121          # shared corrected model (make_H, evaluate)

ap = argparse.ArgumentParser()
ap.add_argument('--release', default='DR1', choices=['DR1', 'DR2'])
args, _ = ap.parse_known_args()
s121.DESI = obs.DesiBAO(args.release)
OUT = os.path.join(HERE, '..', 'Outputs') + ('' if args.release == 'DR1' else '/' + args.release)
os.makedirs(OUT, exist_ok=True)

PSI0 = 2.62
F0_of = lambda L0: 0.5 + L0 * PSI0**2

def log_prior(p):
    L0, w0, wa = p
    if not (0.0001 <= L0 <= 0.010): return -np.inf
    if not (-1.40 <= w0 <= -0.60): return -np.inf
    if not (-2.50 <= wa <= 1.50): return -np.inf
    if w0 + wa >= 0: return -np.inf
    return 0.0

def log_prob(p):
    lp = log_prior(p)
    if not np.isfinite(lp): return -np.inf, np.nan
    e = s121.evaluate(F0_of(p[0]), p[1], p[2])
    if e is None: return -np.inf, np.nan
    return -0.5 * (e['chi2_theta'] + e['chi2_DESI']), e['H0']

print("=" * 70)
print(f"SIM122 — Unlocked Λ₀ scan   [DESI {args.release}, audited]")
print("=" * 70)
ref = obs.lcdm_reference(args.release)
print(f"  ΛCDM reference: χ² = {ref['chi2_DESI']:.2f}/{ref['N_DESI']} (GOF {ref['gof_sigma']:.2f}σ), H₀ = {ref['H0']:.2f}")

# Part A: χ² profile in Λ₀ (minimised over w₀, wₐ)
print("\nPart A: profile over Λ₀")
profile = []
x0 = [-0.7, -1.0]
for L0 in np.linspace(0.0001, 0.010, 12):
    f = lambda q: (-2 * log_prob([L0, q[0], q[1]])[0]) if np.isfinite(log_prob([L0, q[0], q[1]])[0]) else 1e10
    r = minimize(f, x0, method='Nelder-Mead', options=dict(xatol=1e-5, fatol=1e-6, maxiter=3000))
    e = s121.evaluate(F0_of(L0), r.x[0], r.x[1])
    profile.append(dict(Lambda0=float(L0), F0=F0_of(L0), w0=float(r.x[0]), wa=float(r.x[1]),
                        chi2_tot=float(r.fun), H0=float(e['H0']) if e else None))
    print(f"  Λ₀={L0:.4f}  F₀={F0_of(L0):.4f}  min χ²_tot={r.fun:8.4f}  (w₀={r.x[0]:+.3f}, wₐ={r.x[1]:+.3f}, H₀={e['H0'] if e else np.nan:.2f})")
spread = max(p['chi2_tot'] for p in profile) - min(p['chi2_tot'] for p in profile)
print(f"  spread of the profile over Λ₀: Δχ² = {spread:.2e}  → Λ₀ is {'UNCONSTRAINED' if spread < 0.1 else 'constrained'}")

# Part B: MCMC
print("\nPart B: MCMC")
rng = np.random.default_rng(7)
best = min(profile, key=lambda p: p['chi2_tot'])
center = np.array([0.005, best['w0'], best['wa']])
p0 = center + np.array([0.002, 0.03, 0.08]) * rng.normal(size=(32, 3))
for k in range(32):
    while not np.isfinite(log_prior(p0[k])):
        p0[k] = center + np.array([0.002, 0.03, 0.08]) * rng.normal(size=3)
sampler = emcee.EnsembleSampler(32, 3, log_prob, blobs_dtype=[('H0', float)])
st = sampler.run_mcmc(p0, 600, progress=False); sampler.reset()
sampler.run_mcmc(st, 2400, progress=False)
chain = sampler.get_chain(flat=True); lp = sampler.get_log_prob(flat=True); H0c = sampler.get_blobs(flat=True)['H0']
imap = np.argmax(lp); L0m, w0m, wam = chain[imap]
em = s121.evaluate(F0_of(L0m), w0m, wam)
pct = np.percentile(np.column_stack([chain, H0c]), [16, 50, 84], axis=0)
prior_sd = 0.0099 / np.sqrt(12)
post_sd = np.std(chain[:, 0])
print(f"  acceptance {np.mean(sampler.acceptance_fraction):.3f}")
print(f"  Λ₀ posterior: median {pct[1,0]:.4f} (16–84%: {pct[0,0]:.4f}–{pct[2,0]:.4f}); sd {post_sd:.4f} vs flat-prior sd {prior_sd:.4f}")
gof = obs.summarize_fit(em['chi2_DESI'], s121.DESI.n, 3, em['chi2_theta'])
dchi2 = ref['chi2_DESI'] - em['chi2_DESI'] - em['chi2_theta']
crit = dict(chi2N=em['chi2_DESI'] / s121.DESI.n < 2, theta=em['chi2_theta'] < 4, wa_neg=wam < 0,
            F0_range=0.5 <= F0_of(L0m) <= 0.56, L0_pos=L0m > 0, beats_LCDM=dchi2 > 0)
constrained = post_sd < 0.7 * prior_sd
verdict = ('PASS' if all(crit.values()) else 'PARTIAL' if (crit['theta'] or crit['chi2N']) else 'FAIL')
verdict_note = '' if constrained else ' — but Λ₀ is UNCONSTRAINED (exact F₀–H₀ degeneracy): no evidence for CMSTG; gain is from (w₀, wₐ)'
print("=" * 70)
print(f"SIM122 RESULT (audited): {verdict}{verdict_note}")
print(f"  MAP Λ₀={L0m:.4f}, w₀={w0m:.3f}, wₐ={wam:.3f}, H₀={em['H0']:.2f}; χ²_DESI={em['chi2_DESI']:.2f}, "
      f"GOF {gof['gof_sigma']:.2f}σ, Δχ² vs ΛCDM = {dchi2:+.2f}   (published: FAIL)")
print("=" * 70)

fig, ax = plt.subplots(1, 2, figsize=(12, 4.5))
ax[0].plot([p['Lambda0'] for p in profile], [p['chi2_tot'] for p in profile], 'o-')
ax[0].set_xlabel('Λ₀'); ax[0].set_ylabel('min χ²_tot over (w₀, wₐ)'); ax[0].set_title('Profile: flat ⇒ Λ₀ unconstrained')
ax[1].hist(chain[:, 0], bins=40, density=True, alpha=0.7); ax[1].axhline(1 / 0.0099, color='k', ls=':', label='flat prior')
ax[1].set_xlabel('Λ₀'); ax[1].legend(); ax[1].set_title('Λ₀ posterior')
fig.suptitle(f'SIM122 (audited) — DESI {args.release} + θ*'); fig.tight_layout()
for ext in ('pdf', 'png'):
    fig.savefig(os.path.join(OUT, f'sim122_lambda0_profile.{ext}'), dpi=150, bbox_inches='tight')
plt.close(fig)

with open(os.path.join(OUT, 'sim122_results.json'), 'w') as f:
    json.dump(dict(verdict=verdict, verdict_note=verdict_note.strip(' —'), desi_release=args.release,
                   map=dict(Lambda0=float(L0m), w0=float(w0m), wa=float(wam), H0=float(em['H0']), F0=F0_of(L0m)),
                   posterior_16_50_84=dict(zip(['Lambda0', 'w0', 'wa', 'H0'], pct.T.tolist())),
                   Lambda0_constrained=bool(constrained), Lambda0_profile=profile, profile_spread=spread,
                   gof_map=gof, delta_chi2_vs_LCDM=dchi2, criteria=crit, lcdm_reference=ref,
                   published_pre_audit=dict(verdict='FAIL'),
                   audit_note='Audited 2026-10-22; see AUDIT_NOTES.md.'), f, indent=2, default=float)
print("SIM122 complete.")
