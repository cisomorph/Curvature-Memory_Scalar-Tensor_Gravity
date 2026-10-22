"""
SIM134 — CMSTG Phase 3: Exponential Coupling (Dilaton-type)
===========================================================
  F_eff = ½·exp(−2Λ₀Ψ²) + ξΨ  — no gravity flip; attractor where F_eff'(Ψ*) = 0:
      −2Λ₀Ψ*·exp(−2Λ₀Ψ*²) + ξ = 0

Action (8πG = 1, GR is F = ½):
  S = ∫d⁴x√(-g) [ F_eff(Ψ)·R − ½(∂Ψ)² ] + S_SM + Λ_bare
Corrected equations (audit):
  Friedmann: H²(2F + 2F'Ψ' − Ψ'²/6) = ω_m a⁻³ + ω_r a⁻⁴ + Λ_bare
  Field:     Ψ'' + (3−ε)Ψ' = F'·6(2−ε),  ε constraint-consistent
Scan: Λ₀ ∈ {0.20, 0.33, 0.50, 1.00}, ξ = 1/6, Ψ_ini = 0.

AUDIT 2026-10-22: same corrections as SIM131 (see AUDIT_NOTES.md).
"""

import numpy as np
import os, sys, warnings
from scipy.optimize import brentq
warnings.filterwarnings('ignore')

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', '..', 'common'))
import cmstg_obs as obs
from st_solver import ScalarTensorBG, verdict_from
import audit_driver as drv

OUT = os.path.join(HERE, '..', 'Outputs')
os.makedirs(OUT, exist_ok=True)

h_target, XI = 0.674, 1.0 / 6.0
PUBLISHED = dict(best=dict(L0=0.2, tension=3.84, chi2=88.65, theta=0.93117), verdict='FAIL')

def F(p, L0):  return 0.5 * np.exp(-2 * L0 * p * p) + XI * p
def Fp(p, L0): return -2 * L0 * p * np.exp(-2 * L0 * p * p) + XI

def model(L0):
    return ScalarTensorBG(F=lambda p, N: F(p, L0), F_psi=lambda p, N: Fp(p, L0))

print("=" * 72)
print("SIM134 — Phase 3: Exponential (dilaton-type) coupling   [audited]")
print("=" * 72)
cases = []
for L0 in [0.20, 0.33, 0.50, 1.00]:
    try:   # stable attractor: first root of F' below the maximum of 2Λ₀x·e^{−2Λ₀x²} at x = 1/(2√Λ₀)
        ps = brentq(lambda x: Fp(x, L0), 1e-6, 0.5 / np.sqrt(L0))
        print(f"  Λ₀={L0}: Ψ* = {ps:.4f}, F(Ψ*) = {F(ps, L0):.5f}")
    except ValueError:
        print(f"  Λ₀={L0}: no attractor")
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
print(f"SIM134 RESULT (audited): faithful = {v['faithful']},  θ*-matched = {v['theta_matched']}   "
      f"(published: {PUBLISHED['verdict']}, 'tension' {PUBLISHED['best']['tension']}σ)")
print("=" * 72)

drv.save_json(os.path.join(OUT, 'sim134_results.json'), dict(
    sim='SIM134', phase='Phase 3', verdict=v['theta_matched'], verdict_modes=v,
    action='S=∫d⁴x√g[F_eff(Ψ)R - ½(∂Ψ)²]+S_SM+Λ_bare, F_eff=½exp(-2Λ₀Ψ²)+ξΨ', xi=XI,
    lcdm_reference={rel: obs.lcdm_reference(rel) for rel in ('DR1', 'DR2')},
    scan=[drv.strip(c) for c in cases],
    best_faithful=drv.strip(bf) if bf else None, best_theta_matched=drv.strip(bm) if bm else None,
    published_pre_audit=PUBLISHED, audit_note='Audited 2026-10-22; see AUDIT_NOTES.md.'))
print("SIM134 complete.")
