"""
SIM132 — CMSTG Phase 3: Subtractive Quadratic F (Local Deser-Woodard Analog)
============================================================================
SIM131 showed that additive F_eff = ½(1+2Λ₀Ψ²)+ξΨ grows as Ψ grows.
Here the quadratic backbone is flipped:
  F_eff = ½(1 − 2Λ₀Ψ²) + ξΨ,   attractor at Ψ* = ξ/(2Λ₀),  F_eff(Ψ*) = ½ + ξ²/(4Λ₀)
Gravity flips (F_eff = 0) at Ψ_flip = [ξ + √(ξ²+4Λ₀)]/(4Λ₀) — integration stops there.

Action (8πG = 1, GR is F = ½):
  S = ∫d⁴x√(-g) [ F_eff(Ψ)·R − ½(∂Ψ)² ] + S_SM + Λ_bare
Corrected equations (audit):
  Friedmann: H²(2F + 2F'Ψ' − Ψ'²/6) = ω_m a⁻³ + ω_r a⁻⁴ + Λ_bare
  Field:     Ψ'' + (3−ε)Ψ' = F'·6(2−ε),  ε constraint-consistent
Scan: Λ₀ ∈ {0.20, 0.33, 0.50, 1.00}, ξ = 1/6, Ψ_ini = 0, H₀ = 67.4 (faithful) or from θ*.

AUDIT 2026-10-22: same corrections as SIM131 (Friedmann with a GR limit, ε including
field evolution, official DESI likelihood, r_s to ∞, θ* at z*, proper statistics, ΛCDM
baseline in the same pipeline, faithful + θ*-matched modes). See AUDIT_NOTES.md.
"""

import numpy as np
import os, sys, warnings
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
warnings.filterwarnings('ignore')

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', '..', 'common'))
import cmstg_obs as obs
from st_solver import ScalarTensorBG, verdict_from
import audit_driver as drv

OUT = os.path.join(HERE, '..', 'Outputs')
os.makedirs(OUT, exist_ok=True)

h_target, XI = 0.674, 1.0 / 6.0
PUBLISHED = dict(best=dict(L0=0.2, tension=3.83, chi2=88.12, theta=0.93175), verdict='FAIL')

def model(L0):
    return ScalarTensorBG(F=lambda p, N: 0.5 * (1 - 2 * L0 * p * p) + XI * p,
                          F_psi=lambda p, N: -2.0 * L0 * p + XI)

print("=" * 72)
print("SIM132 — Phase 3: Subtractive quadratic F (DW analog)   [audited]")
print("=" * 72)
cases = []
for L0 in [0.20, 0.33, 0.50, 1.00]:
    print(f"  Λ₀={L0}: Ψ* = {XI/(2*L0):.4f}, F(Ψ*) = {0.5 + XI**2/(4*L0):.5f}, "
          f"Ψ_flip = {(XI + np.sqrt(XI**2 + 4*L0))/(4*L0):.4f}")
    c = drv.run_case(model(L0), 0.0, 0.0, h_target, 1, f'Λ₀={L0}')
    c['Lambda0'] = L0
    cases.append(c); print(drv.line(c))

bf, bm = drv.best(cases, 'faithful'), drv.best(cases, 'theta_matched')
v = verdict_from(bf['faithful'] if bf else None, bm['theta_matched'] if bm else None, 1)
if bm:
    for rel in ('DR1', 'DR2'):
        s = bm['theta_matched'][rel]
        print(f"  θ*-matched best ({bm['label']}), DESI {rel}: χ² = {s['chi2_DESI']:.2f}, "
              f"GOF {s['gof_sigma']:.2f}σ, Δχ² vs ΛCDM = {s['delta_chi2_tot_vs_LCDM']:+.2f}")
print("=" * 72)
print(f"SIM132 RESULT (audited): faithful = {v['faithful']},  θ*-matched = {v['theta_matched']}   "
      f"(published: {PUBLISHED['verdict']}, 'tension' {PUBLISHED['best']['tension']}σ)")
print("=" * 72)

fig, ax = plt.subplots(1, 2, figsize=(13, 5))
for c in cases:
    b = c['bg_faithful']
    if b is not None:
        m = b['z'] < 3000
        ax[0].semilogx(b['z'][m] + 1, b['psi'][m], label=c['label'])
        ax[1].semilogx(b['z'][m] + 1, b['F'][m], label=c['label'])
ax[0].set_xlabel('1+z'); ax[0].set_ylabel('Ψ'); ax[0].legend(); ax[0].grid(alpha=.3)
ax[1].axhline(0.5, color='gray', ls=':'); ax[1].set_xlabel('1+z'); ax[1].set_ylabel('F_eff'); ax[1].legend(); ax[1].grid(alpha=.3)
fig.suptitle('SIM132 (audited) — subtractive F, Ψ_ini = 0 (faithful H₀)')
fig.tight_layout()
for ext in ('pdf', 'png'):
    fig.savefig(os.path.join(OUT, f'sim132_evolution.{ext}'), dpi=150, bbox_inches='tight')
plt.close(fig)

drv.save_json(os.path.join(OUT, 'sim132_results.json'), dict(
    sim='SIM132', phase='Phase 3', verdict=v['theta_matched'], verdict_modes=v,
    action='S=∫d⁴x√g[F_eff(Ψ)R - ½(∂Ψ)²]+S_SM+Λ_bare, F_eff=½(1-2Λ₀Ψ²)+ξΨ', xi=XI,
    lcdm_reference={rel: obs.lcdm_reference(rel) for rel in ('DR1', 'DR2')},
    scan=[drv.strip(c) for c in cases],
    best_faithful=drv.strip(bf) if bf else None, best_theta_matched=drv.strip(bm) if bm else None,
    published_pre_audit=PUBLISHED, audit_note='Audited 2026-10-22; see AUDIT_NOTES.md.'))
print("SIM132 complete.")
