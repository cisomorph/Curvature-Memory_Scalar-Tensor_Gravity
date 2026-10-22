"""
SIM135 — CMSTG Phase 3: Bi-scalar (Ψ frozen at Phase 1 canonical, φ sourced by curvature)
===========================================================================================
Ψ is locked at Ψ̄ = 2.62 (F₀ = ½ + Λ₀Ψ̄² = 0.5206, treated as a constant); a second
scalar φ couples linearly to R and starts from φ_ini = 0:
  S = ∫d⁴x√(-g) [ (F₀ + ξ_φ φ)·R − ½(∂φ)² ] + S_SM + Λ_bare
Corrected equations (audit; GR-normalised densities ω):
  Friedmann: H²(2F + 2ξφ' − φ'²/6) = ω_m a⁻³ + ω_r a⁻⁴ + Λ_bare,   F = F₀ + ξφ
  Field:     φ'' + (3−ε)φ' = ξ·6(2−ε),  ε constraint-consistent
(published: H²(3F + 6ξφ' − ½φ'²) = ω… — no GR limit, F'φ' and kinetic terms doubled.)
Scan: ξ ∈ {0.01, 0.03, 0.05, 0.08, 0.10, 0.15}, H₀ = 67.4 (faithful) or from θ*.

AUDIT 2026-10-22: corrected Friedmann/ε, official DESI likelihood, r_s to ∞, θ* at z*,
proper statistics, ΛCDM baseline, faithful + θ*-matched modes, local-gravity checks.
"""

import numpy as np
import os, sys, warnings
warnings.filterwarnings('ignore')

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', '..', 'common'))
import cmstg_obs as obs
from mf_solver import MultiFieldBG
import audit_driver as drv

OUT = os.path.join(HERE, '..', 'Outputs')
os.makedirs(OUT, exist_ok=True)
h_target = 0.674
F0 = 0.5 + 0.003 * 2.62**2
PUBLISHED = dict(best=dict(xi=0.01, tension=3.73, chi2=83.29, theta=0.92389), verdict='FAIL')

def model(xi):
    return MultiFieldBG(1, F=lambda f, N: F0 + xi * f[0], F_i=lambda f, N: np.array([xi]))

def verdict(o):
    if o is None: return 'FAIL (no solution)'
    if o['DR1']['gof_sigma'] < 2.0 and o['chi2_theta'] < 4.0: return 'PASS'
    if o['DR1']['delta_chi2_tot_vs_LCDM'] > 0: return 'PARTIAL'
    return 'FAIL'

print("=" * 72)
print("SIM135 — Phase 3: bi-scalar (Ψ frozen, φ curvature-sourced)   [audited]")
print("=" * 72)
rows = []
for xi in [0.01, 0.03, 0.05, 0.08, 0.10, 0.15]:
    m = model(xi)
    rec = dict(xi=xi)
    for mode, bg in (('faithful', m.calibrated([0.0], [0.0], h_target)),
                     ('theta_matched', m.theta_matched([0.0], [0.0], h_target))):
        o = drv.evaluate_mf(m, bg, 1, rsd=False)
        if o is not None:
            o['local_gravity'] = drv.local_gravity(o['F_eff_today'], [xi],
                                                   Fdot_over_HF=xi * o['y0'][0] / o['F_eff_today'], H0=o['H0'])
        rec[mode] = o; rec[f'verdict_{mode}'] = verdict(o)
        if o is None:
            print(f"  ξ={xi:<5} {mode:13s} no solution"); continue
        lg = o['local_gravity']
        print(f"  ξ={xi:<5} {mode:13s} H0={o['H0']:6.2f} 100θ*={o['theta_star_100']:.4f} φ0={o['f0'][0]:.4f} "
              f"F0={o['F_eff_today']:.4f} χ²DR1={o['DR1']['chi2_DESI']:7.2f} Δχ²DR1={o['DR1']['delta_chi2_tot_vs_LCDM']:+8.2f} "
              f"Δχ²DR2={o['DR2']['delta_chi2_tot_vs_LCDM']:+8.2f} Cassini×{lg['cassini_ratio']:.0f} → {rec[f'verdict_{mode}']}")
    rows.append(rec)

def overall(mode):
    vs = [r[f'verdict_{mode}'] for r in rows]
    return 'PASS' if 'PASS' in vs else 'PARTIAL' if 'PARTIAL' in vs else 'FAIL'
ov = dict(faithful=overall('faithful'), theta_matched=overall('theta_matched'))
print("=" * 72)
print(f"SIM135 RESULT (audited): faithful = {ov['faithful']}, θ*-matched = {ov['theta_matched']}   (published: FAIL, 'tension' 3.73σ)")
print("  Note: φ's coupling ξ is unscreened today; see Cassini ratios above.")
print("=" * 72)
drv.save_json(os.path.join(OUT, 'sim135_results.json'), dict(
    sim='SIM135', phase='Phase 3', verdict=ov['theta_matched'], verdict_modes=ov, F0=F0,
    action='S=∫d⁴x√g[(F₀+ξφ)R − ½(∂φ)²]+S_SM+Λ_bare, F₀=0.5206 (Ψ frozen), φ_ini=0',
    lcdm_reference={rel: obs.lcdm_reference(rel) for rel in ('DR1', 'DR2')},
    scan=rows, published_pre_audit=PUBLISHED, audit_note='Audited 2026-10-22; see AUDIT_NOTES.md.'))
print("SIM135 complete.")
