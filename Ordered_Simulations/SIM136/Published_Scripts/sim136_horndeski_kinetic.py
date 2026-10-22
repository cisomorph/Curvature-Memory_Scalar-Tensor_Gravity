"""
SIM136 — CMSTG Phase 3: Horndeski non-minimal kinetic coupling
=============================================================
Action (8πG = 1, GR is F = ½):
  S = ∫d⁴x√g [ (1+2Λ₀Ψ²)/2 · R − ½(∂Ψ)² + α·G^{μν}∂_μΨ∂_νΨ ] + S_SM + Λ_bare
In Horndeski form this is G₄ = F(Ψ) + αX (X = ½Ψ̇² in FRW).

Corrected equations (audit; derived by lapse variation and checked against the
general Horndeski background equations; ρ in GR units 3H² = ρ, densities ω):
  Friedmann: 6FH² + 6HḞ = ρ + ½Ψ̇²(1 + 18αH²)
             → H²(2F + 2F'Ψ' − Ψ'²/6) − 3αH⁴Ψ'² = ω_m a⁻³ + ω_r a⁻⁴ + Λ_bare
  Field:     (1 + 6αH²)Ψ̈ + [3H + 6α(3H³ + 2HḢ)]Ψ̇ = F'R
             → Ψ'' = [12F' − (3 + 18αH²)Ψ' − ε(6F' − 12αH²Ψ')]/(1 + 6αH²) + εΨ'
  Tensor speed: c_T² = (F + αX)/(F − αX)  ⇒  c_T² − 1 ≈ αH²Ψ'²/F (GW170817: < ~6e-15)
The published sim used ½(1+6αH²)Ψ̇² (energy 3× too small: the term scales as N⁻³ under
the lapse) and left the Ψ equation "unchanged at leading order"; c_T was not checked.
Tests: A — Ψ_ini = 0; B — Ψ_ini = 2.62; α ∈ {1e-8 … 1e-3} [H100⁻²].

AUDIT 2026-10-22: corrected Friedmann and field equation, constraint-consistent ε,
official DESI likelihood, r_s to ∞, θ* at z*, proper statistics, faithful + θ*-matched,
c_T check (RESEARCH_RULES 2.3.2).
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
h_target, Lambda0, Psi_bar = 0.674, 0.003, 2.62
PUBLISHED = dict(best=dict(alpha=1e-8, tension=3.53, chi2=74.65, theta=0.92066), verdict='FAIL')


class HorndeskiBG(ScalarTensorBG):
    def __init__(self, alpha, **kw):
        super().__init__(F=lambda p, N: 0.5 + Lambda0 * p * p, F_psi=lambda p, N: 2 * Lambda0 * p, **kw)
        self.alpha = alpha

    def H2(self, p, y, N, Lb):
        c = self.coef(p, y, N)
        r = self.rho_hat(N, Lb, p)
        b = 3.0 * self.alpha * y * y               # b X² − c X + r = 0
        if b < 1e-300:
            return r / c if c > self.coef_floor else np.nan
        disc = c * c - 4.0 * b * r
        if disc < 0 or c <= 0:
            return np.nan
        X = 2.0 * r / (c + np.sqrt(disc))           # root → r/c as α → 0
        return X if X > 0 else np.nan

    def AB(self, p, y, N, H2):
        Fp = self.F_psi(p, N); k = 1.0 + 6.0 * self.alpha * H2
        if k <= 0:
            raise ValueError('ghost')
        A = (12.0 * Fp - (3.0 + 18.0 * self.alpha * H2) * y) / k
        B = (-6.0 * Fp + 12.0 * self.alpha * H2 * y) / k + y
        return A, B

    def cT2_minus_1_today(self, bg):
        X = 0.5 * (bg['H0'] / 100.0)**2 * bg['y0']**2
        F = bg['F0']
        return 2 * self.alpha * X / (F - self.alpha * X)


cases = []
print("=" * 72)
print("SIM136 — Phase 3: Horndeski G^{μν}∂Ψ∂Ψ coupling   [audited]")
print("=" * 72)
for test, psi_ini in (('A', 0.0), ('B', Psi_bar)):
    for alpha in [1e-8, 1e-7, 1e-6, 1e-5, 1e-4, 1e-3]:
        m = HorndeskiBG(alpha)
        try:
            c = drv.run_case(m, psi_ini, 0.0, h_target, 1, f'{test}: α={alpha:.0e}')
        except ValueError:
            c = dict(label=f'{test}: α={alpha:.0e}', faithful=None, theta_matched=None, bg_faithful=None, bg_matched=None)
        c['alpha'], c['test'] = alpha, test
        cts = []
        for key in ('bg_faithful', 'bg_matched'):
            if c.get(key) is not None:
                cts.append(m.cT2_minus_1_today(c[key]))
        c['cT2_minus_1_today'] = cts
        c['gw170817_ok'] = bool(cts) and all(abs(x) < 6e-15 for x in cts)
        cases.append(c)
        print(drv.line(c), f"| c_T²−1 = {cts[0]:.1e}" if cts else "")

bf, bm = drv.best(cases, 'faithful'), drv.best(cases, 'theta_matched')
v = verdict_from(bf['faithful'] if bf else None, bm['theta_matched'] if bm else None, 1)
# c_T = c is a hard requirement (RESEARCH_RULES 2.3.2)
def gw_gate(verdict_str, case):
    return verdict_str if (case is None or case['gw170817_ok']) else f'FAIL (c_T: {verdict_str} on cosmology)'
v = dict(faithful=gw_gate(v['faithful'], bf), theta_matched=gw_gate(v['theta_matched'], bm))
print("=" * 72)
print(f"SIM136 RESULT (audited): faithful = {v['faithful']}, θ*-matched = {v['theta_matched']}   (published: FAIL)")
print("=" * 72)
drv.save_json(os.path.join(OUT, 'sim136_results.json'), dict(
    sim='SIM136', phase='Phase 3', verdict=v['theta_matched'], verdict_modes=v, Lambda0=Lambda0,
    action='S=∫d⁴x√g[(1+2Λ₀Ψ²)/2·R − ½(∂Ψ)² + α·G^{μν}∂_μΨ∂_νΨ]+S_SM+Λ_bare',
    lcdm_reference={rel: obs.lcdm_reference(rel) for rel in ('DR1', 'DR2')},
    scan=[drv.strip(c) for c in cases],
    best_faithful=drv.strip(bf) if bf else None, best_theta_matched=drv.strip(bm) if bm else None,
    published_pre_audit=PUBLISHED, audit_note='Audited 2026-10-22; see AUDIT_NOTES.md.'))
print("SIM136 complete.")
