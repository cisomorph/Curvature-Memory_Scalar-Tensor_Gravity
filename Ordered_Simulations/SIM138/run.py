#!/usr/bin/env python3
"""
SIM138: DESI per-bin sensitivity decomposition (audited)
Diagnostic: which DESI tracer bins drive the misfit of each model?

AUDIT 2026-10-22 (see AUDIT_NOTES.md). The published version compared models to a
hardcoded "DESI H(z)" table (81.7, 97.9, 110.7, 128.1, 156.4, 240.8) described as
"converted from D_M/r_d, D_H/r_d (Adame et al. 2024)". The conversion is wrong
(z=0.706 and z=0.93 off by 3σ and 6σ; z=0.295 has no H measurement), which is why it
concluded "DESI prefers H₀ ≈ 74". This version uses the official DESI DR1 and DR2
Gaussian likelihoods. DESI's covariance is block-diagonal per tracer bin, so the χ²
decomposes exactly into per-bin contributions Δ_bᵀC_b⁻¹Δ_b.

Models (all with Planck ω_m, ω_b; H₀ from θ*):
  ΛCDM; Phase 1 frozen F₀ (identical to ΛCDM for BAO+θ*, SIM121C);
  Phase 1 unlocked Ψ_ini = 2.62 (SIM130); w₀wₐCDM best fit (SIM121C).
"""

import json, os, sys, warnings
from datetime import datetime
import numpy as np
from scipy.optimize import minimize
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
warnings.filterwarnings('ignore')

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', 'common'))
import cmstg_obs as obs
from st_solver import ScalarTensorBG

FIG_DIR = os.path.join(HERE, 'figures')
os.makedirs(FIG_DIR, exist_ok=True)

def blocks(lik):
    """Group data indices by redshift (one block per tracer bin)."""
    out = {}
    for i, z in enumerate(lik.z):
        out.setdefault(round(float(z), 3), []).append(i)
    return out

def decompose(lik, Hf):
    rd = obs.r_drag(Hf)
    d = lik.model(Hf, rd) - lik.value
    res = {}
    for z, idx in blocks(lik).items():
        C = lik.cov[np.ix_(idx, idx)]
        res[z] = dict(chi2=float(d[idx] @ np.linalg.solve(C, d[idx])),
                      quantities=[lik.kind[i] for i in idx],
                      pulls=[float(d[i] / lik.sigma[i]) for i in idx])
    total = float(d @ lik.icov @ d)
    return res, total, rd

models = {}
h_l = obs.h_from_theta(lambda h: obs.H_lcdm(h))
models['ΛCDM'] = obs.H_lcdm(h_l)
m130 = ScalarTensorBG(F=lambda p, N: 0.5 + 0.003 * p * p, F_psi=lambda p, N: 0.006 * p)
models['Phase 1 unlocked (Ψ_ini=2.62)'] = m130.H_func(m130.theta_matched(2.62, 0.0))
def c2w(p, lik):
    w0, wa = p
    if w0 + wa >= 0: return 1e10
    h = obs.h_from_theta(lambda x: obs.H_lcdm(x, w0=w0, wa=wa))
    if not np.isfinite(h): return 1e10
    H = obs.H_lcdm(h, w0=w0, wa=wa); return lik.chi2(H, obs.r_drag(H))

out = dict(sim_id='SIM138', timestamp=datetime.utcnow().isoformat() + 'Z', releases={})
print("=" * 72)
print("SIM138: DESI per-bin decomposition   [audited]")
print("=" * 72)
for rel in ('DR1', 'DR2'):
    lik = obs.DesiBAO(rel)
    bw = minimize(c2w, [-0.6, -1.3], args=(lik,), method='Nelder-Mead', options=dict(xatol=1e-4, fatol=1e-6))
    hw = obs.h_from_theta(lambda x: obs.H_lcdm(x, w0=bw.x[0], wa=bw.x[1]))
    ms = dict(models); ms[f'w0waCDM (w0={bw.x[0]:.2f}, wa={bw.x[1]:.2f})'] = obs.H_lcdm(hw, w0=bw.x[0], wa=bw.x[1])
    print(f"\nDESI {rel}:   per-bin χ² (block covariance)")
    zs = sorted(blocks(lik).keys())
    print(f"  {'model':42s} " + " ".join(f"z={z:<5}" for z in zs) + "   total")
    rel_out = {}
    for name, Hf in ms.items():
        res, total, rd = decompose(lik, Hf)
        rel_out[name] = dict(per_bin=res, total=total, r_drag=rd, H0=float(Hf(np.array([0.0]))[0]))
        print(f"  {name:42s} " + " ".join(f"{res[z]['chi2']:7.2f}" for z in zs) + f"   {total:6.2f}")
    out['releases'][rel] = rel_out
    # figure: per-bin χ²
    fig, ax = plt.subplots(figsize=(9, 4.5))
    wdt = 0.8 / len(ms)
    for j, (name, r) in enumerate(rel_out.items()):
        ax.bar(np.arange(len(zs)) + j * wdt, [r['per_bin'][z]['chi2'] for z in zs], wdt, label=f"{name} (Σ={r['total']:.1f})")
    ax.set_xticks(np.arange(len(zs)) + 0.4 - wdt / 2); ax.set_xticklabels([f'z={z}' for z in zs])
    ax.set_ylabel('χ² contribution'); ax.set_title(f'SIM138 (audited) — DESI {rel} per tracer bin'); ax.legend(fontsize=7)
    fig.tight_layout()
    for ext in ('pdf', 'png'):
        fig.savefig(os.path.join(FIG_DIR, f'sim138_per_bin_{rel}.{ext}'), dpi=150, bbox_inches='tight')
    plt.close(fig)

out['verdict'] = 'DIAGNOSTIC'
out['published_pre_audit'] = dict(verdict='DISTRIBUTED', chi2_phase1=13.15,
                                  note='non-DESI H(z) table; "DESI prefers H0≈74" was an artifact of that table')
out['audit_note'] = 'Audited 2026-10-22; see AUDIT_NOTES.md.'
with open(os.path.join(HERE, 'output.json'), 'w') as f:
    json.dump(out, f, indent=2, default=float)
print("\nSIM138 complete.")
