#!/usr/bin/env python3
"""
SIM144 — Mixed-source scalar φ: curvature (ξ_R) + matter (β_m) coupling
=========================================================================
Mechanism-completeness probe for Paper I §3.4.

Action (as published):
  S = ∫d⁴x √-g [(1+2Λ₀Ψ²)/2 · R − ½(∂Ψ)² − ½(∂φ)² + ξ_R φ R + 2β_m φ ρ_m] + S_SM + Λ_bare
  F_eff = ½ + Λ₀Ψ² + ξ_R φ
The 2β_mφρ_m term is not a covariant Lagrangian (ρ_m is not a field); it is treated
here, as in the published φ equation, as an effective interaction with ρ_m the
physical matter density (GR units, 3H² = ρ):
  energy density −2β_mφρ_m  (= −L_int),  force on φ  +2β_mρ_m.
The published Friedmann added +2β_mφρ_m (wrong sign w.r.t. its own φ equation)
and used ω_m a⁻³ in the φ equation (not the physical density). Both fixed.

Corrected equations (audit; Planck-calibrated ρ̂ = 2F₀ω; ρ_m = 3·2F₀·ω_m a⁻³):
  H²[2F + 2(F_ΨΨ' + ξ_Rφ') − Ψ'²/6 − φ'²/6] = 2F₀(ω_m a⁻³ + ω_r a⁻⁴) + Λ_b − 2β_mφρ_m/3
  Ψ'' + (3−ε)Ψ' = F_Ψ·6(2−ε)
  φ'' + (3−ε)φ' = ξ_R·6(2−ε) + 2β_mρ_m/H²
Growth: G_eff of the scalar-tensor part × (1 + 8β_m²) for the matter-coupled φ exchange.
Local gravity: α² from the F-couplings plus 8β_m² (fifth-force fraction of the direct coupling).
Grid: ξ_R ∈ {0, 0.01, 0.1, 1}, β_m ∈ {0, 0.01, 0.1, 1}, φ_ini = 0, Ψ_ini = 2.62.
AUDIT 2026-10-22: see AUDIT_NOTES.md; faithful (H₀ = 67.59) + θ*-matched modes.
"""

import json, os, sys, warnings
from datetime import datetime
import numpy as np
warnings.filterwarnings('ignore')

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', 'common'))
import cmstg_obs as obs
from mf_solver import MultiFieldBG
import audit_driver as drv

h_target, Lambda0, PSI_BAR = 0.6759, 0.003, 2.62
F0 = 0.5 + Lambda0 * PSI_BAR**2
RS = 2.0 * F0

def model(xi_R, beta_m):
    rho_m = lambda N: 3.0 * RS * obs.OMEGA_M * np.exp(-3.0 * N)
    return MultiFieldBG(2,
        F=lambda f, N: 0.5 + Lambda0 * f[0]**2 + xi_R * f[1],
        F_i=lambda f, N: np.array([2 * Lambda0 * f[0], xi_R]),
        V_int=lambda f, N: -2.0 * beta_m * f[1] * rho_m(N),
        V_int_i=lambda f, N: np.array([0.0, -2.0 * beta_m * rho_m(N)]),
        rho_scale=RS, fifth_force=lambda f, N: 1.0 + 8.0 * beta_m**2)

def verdict(o):
    if o is None: return 'FAIL_NO_SOLUTION'
    if o['DR1']['gof_sigma'] >= 2.0: return 'FAIL_DESI'
    if o['chi2_theta'] >= 4.0: return 'FAIL_CMB'
    if o['rsd_chi2N'] >= 1.5: return 'FAIL_RSD'
    return 'PASS'

print("=" * 72)
print("SIM144: Mixed-source φ (ξ_R, β_m) completeness probe   [audited]")
print("=" * 72)
results = []
for xi_R in [0.0, 0.01, 0.1, 1.0]:
    for beta_m in [0.0, 0.01, 0.1, 1.0]:
        m = model(xi_R, beta_m)
        rec = dict(xi_R=xi_R, beta_m=beta_m)
        for mode in ('faithful', 'theta_matched'):
            bg = m.calibrated([PSI_BAR, 0.0], [0, 0], h_target) if mode == 'faithful' else \
                 m.theta_matched([PSI_BAR, 0.0], [0, 0], h_target)
            o = drv.evaluate_mf(m, bg, (xi_R > 0) + (beta_m > 0))
            if o is not None:
                Fi = np.array([2 * Lambda0 * o['f0'][0], xi_R])
                Fdot = (Fi @ np.array(o['y0'])) / o['F_eff_today']
                o['local_gravity'] = drv.local_gravity(o['F_eff_today'], Fi, Fdot_over_HF=Fdot, H0=o['H0'],
                                                       alpha2_extra=8.0 * beta_m**2)
            rec[mode] = o; rec[f'verdict_{mode}'] = verdict(o)
            if o is None:
                print(f"  ξ_R={xi_R:<5} β_m={beta_m:<5} {mode:13s} no solution"); continue
            lg = o['local_gravity']
            print(f"  ξ_R={xi_R:<5} β_m={beta_m:<5} {mode:13s} H0={o['H0']:6.2f} 100θ*={o['theta_star_100']:.4f} "
                  f"χ²DR1={o['DR1']['chi2_DESI']:7.2f} Δχ²DR1={o['DR1']['delta_chi2_tot_vs_LCDM']:+8.2f} "
                  f"Δχ²DR2={o['DR2']['delta_chi2_tot_vs_LCDM']:+8.2f} RSD={o['rsd_chi2N']:5.2f} "
                  f"Cassini×{lg['cassini_ratio']:.0f} → {rec[f'verdict_{mode}']}")
        results.append(rec)

def overall(mode):
    return 'PASS' if any(r[f'verdict_{mode}'] == 'PASS' for r in results) else 'FAIL'
ov = dict(faithful=overall('faithful'), theta_matched=overall('theta_matched'))
print("=" * 72)
print(f"SIM144 RESULT (audited): faithful = {ov['faithful']}, θ*-matched = {ov['theta_matched']}   (published: FAIL/STRUCTURAL)")
print("=" * 72)
with open(os.path.join(HERE, 'SIM144_results.json'), 'w') as f:
    json.dump(dict(sim_id='SIM144', timestamp=datetime.utcnow().isoformat() + 'Z', verdict=ov['theta_matched'],
                   verdict_modes=ov, lcdm_reference={rel: obs.lcdm_reference(rel) for rel in ('DR1', 'DR2')},
                   scan_results=results, published_pre_audit=dict(verdict='FAIL/STRUCTURAL'),
                   audit_note='Audited 2026-10-22; see AUDIT_NOTES.md.'), f, indent=2, default=float)
print("SIM144 complete.")
