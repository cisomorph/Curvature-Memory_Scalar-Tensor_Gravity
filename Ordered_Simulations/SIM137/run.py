#!/usr/bin/env python3
"""
SIM137: SPARC failure-mode analysis
Diagnostic: why did 96/161 galaxies not PASS the chi-DM fit in SIM119?

No new physics. Re-analyzes SIM119 outputs against features derived from
SPARC rotation curve files.
"""

import json, glob, os, re, warnings
from datetime import datetime
import numpy as np
from scipy import stats
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
warnings.filterwarnings('ignore')

# ── Paths ─────────────────────────────────────────────────────────────────────
BASE = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..'))
SIM119_JSON = os.path.join(BASE, 'Ordered_Simulations', 'SIM119', 'Outputs', 'sim119_results.json')
SPARC_DIR   = os.path.join(BASE, 'Ordered_Simulations',
                           'simulation_18_cmstg_graviton_emit',
                           'cmstg_p7_min_impl', 'external', 'sparc_raw')
OUT_DIR     = os.path.dirname(__file__)
FIG_DIR     = os.path.join(OUT_DIR, 'figures')
os.makedirs(FIG_DIR, exist_ok=True)

# ── Load SIM119 results ───────────────────────────────────────────────────────
with open(SIM119_JSON) as f:
    sim119 = json.load(f)

galaxies = sim119['all_galaxies']
print(f"SIM119 galaxies loaded: {len(galaxies)}")
print(f"  PASS    : {sum(1 for g in galaxies if g['verdict']=='PASS')}")
print(f"  marginal: {sum(1 for g in galaxies if g['verdict']=='marginal')}")
print(f"  fail    : {sum(1 for g in galaxies if g['verdict']=='fail')}")

# ── Build name → rotmod path map (case-insensitive) ──────────────────────────
rotmod_files = glob.glob(os.path.join(SPARC_DIR, '*.dat'))
rotmod_map = {os.path.basename(f).replace('_rotmod.dat', '').lower(): f
              for f in rotmod_files}

# ── Feature extraction from rotmod files ─────────────────────────────────────
def load_rotmod(filepath):
    distance = None
    rows = []
    with open(filepath) as f:
        for line in f:
            line = line.strip()
            if line.startswith('# Distance'):
                m = re.search(r'[\d.]+', line.split('=')[1])
                if m:
                    distance = float(m.group())
            elif line.startswith('#') or not line:
                continue
            else:
                try:
                    vals = list(map(float, line.split()))
                    if len(vals) >= 7:
                        rows.append(vals[:8])
                except ValueError:
                    continue
    if not rows:
        return None
    arr = np.array(rows)
    R, Vobs, errV = arr[:,0], arr[:,1], arr[:,2]
    Vgas, Vdisk, Vbul = arr[:,3], arr[:,4], arr[:,5]
    SBdisk = arr[:,6]
    mask = (errV > 0) & (Vobs > 0) & (R > 0)
    if mask.sum() < 3:
        return None
    return dict(R=R[mask], Vobs=Vobs[mask], Vgas=Vgas[mask],
                Vdisk=Vdisk[mask], Vbul=Vbul[mask], SBdisk=SBdisk[mask],
                distance=distance)

def extract_features(d):
    R, Vobs = d['R'], d['Vobs']
    Vgas, Vdisk, Vbul = d['Vgas'], d['Vdisk'], d['Vbul']
    SBdisk = d['SBdisk']

    v_flat   = float(np.max(Vobs))
    N_pts    = len(R)
    R_max    = float(np.max(R))
    distance = float(d['distance']) if d['distance'] is not None else np.nan
    has_bulge = int(np.max(Vbul) > 0)

    outer = R > 0.5 * R_max
    if outer.sum() < 2:
        outer = np.ones(len(R), dtype=bool)
    Vbar2_out = Vgas[outer]**2 + Vdisk[outer]**2 + Vbul[outer]**2
    gas_frac = float(np.mean(Vgas[outer]**2) / np.maximum(Vbar2_out.mean(), 1e-10))

    sb_pos = SBdisk[SBdisk > 0]
    log_SB = float(np.log10(np.mean(sb_pos))) if len(sb_pos) > 0 else np.nan

    return dict(v_flat=v_flat, N_pts=N_pts, R_max=R_max,
                distance=distance, has_bulge=has_bulge,
                gas_frac=gas_frac, log_SBdisk=log_SB)

# ── Build feature table ───────────────────────────────────────────────────────
records = []
n_missing = 0
for g in galaxies:
    name = g['name']
    path = rotmod_map.get(name.lower())
    if path is None:
        n_missing += 1
        continue
    d = load_rotmod(path)
    if d is None:
        n_missing += 1
        continue
    feats = extract_features(d)
    feats['name']     = name
    feats['chi2']     = g['chi2']
    feats['m22']      = g.get('m22')
    feats['at_bound'] = int(g.get('at_bound', False))
    feats['is_pass']  = int(g['verdict'] == 'PASS')
    records.append(feats)

print(f"Features extracted: {len(records)} galaxies ({n_missing} missing rotmod)")

FEATURE_LABELS = {
    'v_flat':    r'$v_{\rm flat}$ [km/s]',
    'N_pts':     r'$N_{\rm pts}$',
    'R_max':     r'$R_{\rm max}$ [kpc]',
    'distance':  r'Distance [Mpc]',
    'has_bulge': r'Has bulge',
    'gas_frac':  r'Gas fraction (outer)',
    'log_SBdisk':r'$\log \Sigma_{\rm disk}$ [L/pc²]',
    'at_bound':  r'$m_{22}$ at boundary',
    'chi2':      r'$\chi^2/N$ (SIM119)',
}

FEATURES = ['v_flat', 'N_pts', 'R_max', 'distance', 'gas_frac', 'log_SBdisk', 'at_bound']

# ── Split into PASS / NON-PASS ────────────────────────────────────────────────
def col(key):
    return np.array([r[key] for r in records], dtype=float)

is_pass = col('is_pass').astype(bool)
grp_p   = {k: col(k)[is_pass]  for k in FEATURES + ['chi2']}
grp_np  = {k: col(k)[~is_pass] for k in FEATURES + ['chi2']}

n_pass   = is_pass.sum()
n_nopass = (~is_pass).sum()
print(f"\nClassification: {n_pass} PASS, {n_nopass} NON-PASS")

# ── Statistical tests ─────────────────────────────────────────────────────────
print("\n{'Feature':<14} {'KS_p':>8} {'MWU_p':>8} {'Spearman_r':>12} {'Spearman_p':>12}")
print("-" * 60)

stat_results = {}
for feat in FEATURES:
    p_vals = grp_p[feat]
    np_vals = grp_np[feat]
    chi2_all = col('chi2')
    feat_all = col(feat)

    # Remove NaNs
    valid_p  = p_vals[np.isfinite(p_vals)]
    valid_np = np_vals[np.isfinite(np_vals)]

    ks_stat, ks_p   = stats.ks_2samp(valid_p, valid_np) if (len(valid_p)>1 and len(valid_np)>1) else (np.nan, np.nan)
    mwu_stat, mwu_p = stats.mannwhitneyu(valid_p, valid_np, alternative='two-sided') if (len(valid_p)>1 and len(valid_np)>1) else (np.nan, np.nan)

    mask_valid = np.isfinite(feat_all) & np.isfinite(chi2_all)
    if mask_valid.sum() > 4:
        sp_r, sp_p = stats.spearmanr(feat_all[mask_valid], chi2_all[mask_valid])
    else:
        sp_r, sp_p = np.nan, np.nan

    stat_results[feat] = dict(ks_p=ks_p, mwu_p=mwu_p, spearman_r=sp_r, spearman_p=sp_p,
                              mean_pass=float(np.nanmean(valid_p)),
                              mean_nonpass=float(np.nanmean(valid_np)))
    print(f"{feat:<14} {ks_p:>8.4f} {mwu_p:>8.4f} {sp_r:>12.4f} {sp_p:>12.4f}")

# ── Logistic regression (scipy.optimize) ─────────────────────────────────────
from scipy.optimize import minimize

def sigmoid(x):
    return 1.0 / (1.0 + np.exp(-np.clip(x, -30, 30)))

def neg_log_likelihood(theta, X, y):
    p = sigmoid(X @ theta)
    p = np.clip(p, 1e-9, 1-1e-9)
    return -np.mean(y * np.log(p) + (1-y) * np.log(1-p))

def neg_log_likelihood_grad(theta, X, y):
    p = sigmoid(X @ theta)
    return X.T @ (p - y) / len(y)

# Standardize features
logreg_feats = ['v_flat', 'N_pts', 'R_max', 'gas_frac', 'log_SBdisk']
X_raw = np.column_stack([col(f) for f in logreg_feats])
y     = is_pass.astype(float)
valid = np.all(np.isfinite(X_raw), axis=1)
X_raw, y_lr = X_raw[valid], y[valid]
mu, sigma = X_raw.mean(0), X_raw.std(0) + 1e-12
X_std = (X_raw - mu) / sigma
X_full = np.column_stack([np.ones(len(X_std)), X_std])

res = minimize(neg_log_likelihood, np.zeros(X_full.shape[1]),
               jac=neg_log_likelihood_grad, args=(X_full, y_lr),
               method='L-BFGS-B')
theta_hat = res.x

# Bootstrap SEs
np.random.seed(42)
n_boot = 500
boot_thetas = np.zeros((n_boot, len(theta_hat)))
for b in range(n_boot):
    idx = np.random.choice(len(y_lr), len(y_lr), replace=True)
    r = minimize(neg_log_likelihood, theta_hat, jac=neg_log_likelihood_grad,
                 args=(X_full[idx], y_lr[idx]), method='L-BFGS-B')
    boot_thetas[b] = r.x
theta_se = boot_thetas.std(0)
theta_z  = theta_hat / (theta_se + 1e-12)

print("\nLogistic regression (PASS vs NON-PASS):")
print(f"  {'Term':<14} {'β':>8} {'SE':>8} {'|z|':>8} {'sig?':>6}")
terms = ['intercept'] + logreg_feats
for i, (t, b, se, z) in enumerate(zip(terms, theta_hat, theta_se, theta_z)):
    print(f"  {t:<14} {b:>8.3f} {se:>8.3f} {abs(z):>8.2f} {'**' if abs(z)>2 else ''}")

logreg_results = {t: dict(beta=float(b), se=float(se), z=float(z))
                  for t, b, se, z in zip(terms, theta_hat, theta_se, theta_z)}

# ── m22 / at_bound check ─────────────────────────────────────────────────────
m22_pass  = col('m22')[is_pass]
m22_np    = col('m22')[~is_pass]
at_b_pass = col('at_bound')[is_pass]
at_b_np   = col('at_bound')[~is_pass]
print(f"\nat_bound: PASS={at_b_pass.mean():.2f}, NON-PASS={at_b_np.mean():.2f}")
m22_finite_p  = m22_pass[np.isfinite(m22_pass) & (m22_pass>0)]
m22_finite_np = m22_np[np.isfinite(m22_np) & (m22_np>0)]
if len(m22_finite_p) and len(m22_finite_np):
    print(f"m22 median: PASS={np.median(m22_finite_p):.3f}, NON-PASS={np.median(m22_finite_np):.3f}")

# ── Figures ───────────────────────────────────────────────────────────────────
plot_feats = ['v_flat', 'N_pts', 'R_max', 'gas_frac', 'log_SBdisk', 'chi2']
fig, axes = plt.subplots(2, 3, figsize=(13, 8))
axes = axes.flatten()
for ax, feat in zip(axes, plot_feats):
    p_vals  = grp_p[feat]
    np_vals = grp_np[feat]
    finite  = np.isfinite(np.concatenate([p_vals, np_vals]))
    all_v   = np.concatenate([p_vals, np_vals])
    lo, hi  = np.nanpercentile(all_v, 2), np.nanpercentile(all_v, 98)
    bins = np.linspace(lo, hi, 20)
    ax.hist(p_vals[np.isfinite(p_vals)],   bins=bins, alpha=0.6, color='steelblue',  label='PASS', density=True)
    ax.hist(np_vals[np.isfinite(np_vals)], bins=bins, alpha=0.6, color='firebrick', label='NON-PASS', density=True)
    ks_p = stat_results.get(feat, {}).get('ks_p', np.nan)
    lbl  = FEATURE_LABELS.get(feat, feat)
    ax.set_xlabel(lbl, fontsize=10)
    ax.set_ylabel('Density', fontsize=9)
    ax.set_title(f'KS p={ks_p:.3f}' if np.isfinite(ks_p) else lbl, fontsize=10)
    ax.legend(fontsize=8)
plt.suptitle('SIM137: PASS vs NON-PASS feature distributions', fontsize=12, y=1.01)
plt.tight_layout()
fig.savefig(os.path.join(FIG_DIR, 'sim137_feature_histograms.pdf'), bbox_inches='tight')
fig.savefig(os.path.join(FIG_DIR, 'sim137_feature_histograms.png'), bbox_inches='tight', dpi=150)
plt.close()

# chi2 vs v_flat scatter
fig2, ax2 = plt.subplots(figsize=(7, 5))
chi2_all = col('chi2')
vf_all   = col('v_flat')
ax2.scatter(vf_all[is_pass],  chi2_all[is_pass],  c='steelblue', s=18, alpha=0.7, label='PASS', zorder=3)
ax2.scatter(vf_all[~is_pass], chi2_all[~is_pass], c='firebrick', s=18, alpha=0.5, label='NON-PASS', zorder=2)
ax2.axhline(2.0, ls='--', c='k', lw=0.8, label='χ²/N = 2 threshold')
ax2.set_xlabel(r'$v_{\rm flat}$ [km/s]', fontsize=11)
ax2.set_ylabel(r'$\chi^2/N$ (SIM119)', fontsize=11)
ax2.set_yscale('log')
ax2.set_title('SIM137: χ²/N vs peak rotation velocity', fontsize=11)
ax2.legend(fontsize=9)
plt.tight_layout()
fig2.savefig(os.path.join(FIG_DIR, 'sim137_chi2_vs_vflat.pdf'), bbox_inches='tight')
fig2.savefig(os.path.join(FIG_DIR, 'sim137_chi2_vs_vflat.png'), bbox_inches='tight', dpi=150)
plt.close()

# ── Verdict ───────────────────────────────────────────────────────────────────
sig_features = [f for f in FEATURES if stat_results[f]['ks_p'] < 0.01]
sig_logreg   = [t for t, v in logreg_results.items() if t != 'intercept' and abs(v['z']) > 2.0]

if sig_features or sig_logreg:
    verdict = "STRUCTURAL_PATTERN"
elif all(stat_results[f]['ks_p'] > 0.05 for f in FEATURES):
    verdict = "NO_PATTERN"
else:
    verdict = "INCONCLUSIVE"

print(f"\nVERDICT: {verdict}")
print(f"  Significant features (KS p<0.01): {sig_features}")
print(f"  Significant logistic terms (|z|>2): {sig_logreg}")

# ── Output JSON ───────────────────────────────────────────────────────────────
output = {
    "sim_id": "SIM137",
    "timestamp": datetime.utcnow().isoformat() + "Z",
    "action_spec": "Phase 1 canonical — diagnostic reanalysis of SIM119",
    "parameters": {
        "n_galaxies_total": len(galaxies),
        "n_pass": int(n_pass),
        "n_nonpass": int(n_nopass),
        "n_rotmod_loaded": len(records),
        "classification": "PASS vs MARGINAL+FAIL"
    },
    "observational_targets": {
        "dataset": "SPARC (Lelli et al. 2016, 161 galaxies)",
        "chi2": None,
        "dof": None,
        "note": "Diagnostic sim; no new chi2 target"
    },
    "theoretical_checks": {
        "gr_recovery": True,
        "c_T_eq_c": True,
        "no_tachyon": True,
        "ward_identity": True,
        "uv_finite": True,
        "note": "Inherited from Phase 1 canonical; no new Lagrangian terms"
    },
    "feature_statistics": {
        feat: {
            "ks_p": float(v['ks_p']),
            "mwu_p": float(v['mwu_p']),
            "spearman_r": float(v['spearman_r']),
            "spearman_p": float(v['spearman_p']),
            "mean_pass": float(v['mean_pass']),
            "mean_nonpass": float(v['mean_nonpass'])
        }
        for feat, v in stat_results.items()
    },
    "logistic_regression": logreg_results,
    "at_bound_frac_pass": float(at_b_pass.mean()),
    "at_bound_frac_nonpass": float(at_b_np.mean()),
    "sig_features_ks001": sig_features,
    "sig_logreg_terms": sig_logreg,
    "verdict": verdict,
    "failure_mode": (
        "PASS/NON-PASS separation driven by: " + ", ".join(sig_features + sig_logreg)
        if verdict == "STRUCTURAL_PATTERN"
        else "No structural pattern found; χ-DM mechanism may be incomplete"
    ),
    "derived_vs_phenom": {
        "chi2": "derived (from SIM119 CMSTG-χ fit)",
        "v_flat": "derived (max Vobs from SPARC rotmod)",
        "gas_frac": "derived (outer Vgas²/Vbar² from SPARC rotmod)"
    }
}

with open(os.path.join(OUT_DIR, 'output.json'), 'w') as f:
    json.dump(output, f, indent=2)

print(f"\nWrote output.json and figures to {OUT_DIR}")
print("SIM137 complete.")
