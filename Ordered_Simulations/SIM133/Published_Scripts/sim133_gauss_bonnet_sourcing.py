"""
SIM133 — CMSTG Phase 3: Gauss-Bonnet Sourcing (Direction B)
===========================================================
Source Ψ from the Gauss–Bonnet invariant 𝒢 instead of R. In FRW
𝒢 = 24H²(Ḣ + H²) = 24H⁴(1−ε): negative in matter domination, positive under Λ.

Action (8πG = 1, GR is F = ½):
  S = ∫d⁴x√(-g) [ ½(1+2Λ₀Ψ²)R − ½(∂Ψ)² + α_GB·Ψ·𝒢 ] + S_SM + Λ_bare
  α̃ = α_GB·H₀² (dimensionless), α_GB = α̃/h_target² in H100⁻² units
Corrected equations (audit):
  Friedmann: H²(2F + 2F'Ψ' − Ψ'²/6) + 8α_GB Ψ' H⁴ = ω_m a⁻³ + ω_r a⁻⁴ + Λ_bare
             (the GB term follows from the lapse variation of √-g α_GB Ψ 𝒢
              = −8 α_GB Ψ̇ ȧ³/N³ after integration by parts; the published sim
              omitted it: "Modified Friedmann: same Phase 1 form")
  Field:     Ψ'' + (3−ε)Ψ' = F'·6(2−ε) + α_GB·24H²(1−ε),  ε constraint-consistent
Scan: α̃ ∈ {1e-8, 1e-7, 1e-6, 1e-5, 1e-4}, Ψ_ini = 0, Λ₀ = 0.003.

AUDIT 2026-10-22: GB Friedmann term added; otherwise the same corrections as SIM131
(see AUDIT_NOTES.md).
"""

import numpy as np
import os, sys, warnings
warnings.filterwarnings('ignore')

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', '..', 'common'))
import cmstg_obs as obs
from st_solver import ScalarTensorBG, verdict_from
import audit_driver as drv

OUT = os.path.join(HERE, '..', 'Outputs')
os.makedirs(OUT, exist_ok=True)

h_target, Lambda0 = 0.674, 0.003
PUBLISHED = dict(best=dict(alpha=1e-8, tension=9.68, chi2=562.5, theta=63.03), verdict='FAIL')

def model(alpha_tilde):
    a_gb = alpha_tilde / h_target**2
    return ScalarTensorBG(F=lambda p, N: 0.5 + Lambda0 * p * p, F_psi=lambda p, N: 2.0 * Lambda0 * p,
                          f_gb_psi=lambda p: a_gb)

print("=" * 72)
print("SIM133 — Phase 3: Gauss–Bonnet sourcing   [audited]")
print("=" * 72)
cases = []
for at in [1e-8, 1e-7, 1e-6, 1e-5, 1e-4]:
    c = drv.run_case(model(at), 0.0, 0.0, h_target, 1, f'α̃={at:.0e}')
    c['alpha_tilde'] = at
    b = c['bg_faithful']
    if b is not None:
        print(f"    α̃={at:.0e}: Ψ(z*) = {np.interp(1089.92, b['z'], b['psi']):.4g}, Ψ₀ = {b['psi0']:.4g}, F₀ = {b['F0']:.4g}")
    cases.append(c); print(drv.line(c))

bf, bm = drv.best(cases, 'faithful'), drv.best(cases, 'theta_matched')
v = verdict_from(bf['faithful'] if bf else None, bm['theta_matched'] if bm else None, 1)
if bm:
    for rel in ('DR1', 'DR2'):
        s = bm['theta_matched'][rel]
        print(f"  θ*-matched best ({bm['label']}), DESI {rel}: χ² = {s['chi2_DESI']:.2f}, "
              f"GOF {s['gof_sigma']:.2f}σ, Δχ² vs ΛCDM = {s['delta_chi2_tot_vs_LCDM']:+.2f}")
print("=" * 72)
print(f"SIM133 RESULT (audited): faithful = {v['faithful']},  θ*-matched = {v['theta_matched']}   "
      f"(published: {PUBLISHED['verdict']}, 100θ* {PUBLISHED['best']['theta']})")
print("=" * 72)

drv.save_json(os.path.join(OUT, 'sim133_results.json'), dict(
    sim='SIM133', phase='Phase 3', verdict=v['theta_matched'], verdict_modes=v,
    action='S=∫d⁴x√g[½(1+2Λ₀Ψ²)R - ½(∂Ψ)² + α_GB·Ψ·𝒢] + S_SM + Λ_bare', Lambda0=Lambda0,
    lcdm_reference={rel: obs.lcdm_reference(rel) for rel in ('DR1', 'DR2')},
    scan=[drv.strip(c) for c in cases],
    best_faithful=drv.strip(bf) if bf else None, best_theta_matched=drv.strip(bm) if bm else None,
    published_pre_audit=PUBLISHED, audit_note='Audited 2026-10-22 (GB Friedmann term added); see AUDIT_NOTES.md.'))
print("SIM133 complete.")
