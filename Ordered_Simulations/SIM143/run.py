#!/usr/bin/env python3
"""
SIM143 — Bi-scalar: Phase 1 Ψ + decoupled quintessence φ

Action:
  L = (½+Λ₀Ψ²)R − ½(∂Ψ)² − ½(∂φ)² − U(φ) + L_matter + Λ_bare
  φ minimally coupled (not sourced by R).
Potentials (as published; U₀ ∈ {0.05, 0.20, 0.50}):
  (A) exponential U₀e^{−λφ}, λ ∈ {0.5, 1, 2}, φ_ini = 0.01
  (B) power law  U₀/φⁿ, n ∈ {1, 2},         φ_ini = 1
  (C) hilltop    U₀(1 − φ²/μ²)², μ = 1,      φ_ini = 0.01
Units: the published Friedmann effectively added U/(3F₀) to H², i.e. U in units where
ρ_φ = 2U on the GR scale 3H² = ρ. That energy weight is kept (V = 2U in ρ-units); the
φ equation then has force −V_φ/H² = −2U_φ/H² (published −U_φ/H², inconsistent by 2).

Corrected equations (audit; Planck-calibrated ρ̂ = 2F₀·ω):
  H²(2F + 2F_ΨΨ' − Ψ'²/6 − φ'²/6) = 2F₀(ω_m a⁻³ + ω_r a⁻⁴) + Λ_b + V(φ)/3
  Ψ'' + (3−ε)Ψ' = F_Ψ·6(2−ε),   φ'' + (3−ε)φ' = −V_φ/H²
r_s: integrated to z = ∞ (the published "R_S_CANON = 144.7 + Δr_s" anchoring is removed).

PASS (as published): |Δθ*| < 2σ, DESI GOF < 2σ, RSD χ²/N < 1.5, Ω_φ(z_drag) < 0.05,
Λ_b ≥ −0.05, |Ψ − Ψ̄| < 0.01 is replaced by "Ψ evolves as in SIM130" (Ψ is never frozen).
Cases with Δχ² vs Phase 1 below 1 are "TRIVIAL" (φ does nothing).
AUDIT 2026-10-22: see AUDIT_NOTES.md.
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

def U_pair(form, p, U0):
    if form == 'exp':
        lam = p['lam']
        return (lambda x: U0 * np.exp(-lam * x)), (lambda x: -lam * U0 * np.exp(-lam * x))
    if form == 'power':
        n = p['n']
        return (lambda x: U0 / max(abs(x), 1e-8)**n), (lambda x: -n * U0 / max(abs(x), 1e-8)**(n + 1) * np.sign(x))
    mu = p['mu']
    return (lambda x: U0 * max(1 - (x / mu)**2, 0.0)**2), (lambda x: -4 * U0 * x / mu**2 * max(1 - (x / mu)**2, 0.0))

def model(form, p, U0):
    U, Up = U_pair(form, p, U0)
    return MultiFieldBG(2, F=lambda f, N: 0.5 + Lambda0 * f[0]**2,
                        F_i=lambda f, N: np.array([2 * Lambda0 * f[0], 0.0]),
                        V=lambda f: 2.0 * U(f[1]), V_i=lambda f: np.array([0.0, 2.0 * Up(f[1])]),
                        rho_scale=2.0 * F0)

def omega_phi_drag(m, bg):
    zd = obs.Z_DRAG
    k = np.argmin(abs(bg['z'] - zd))
    f, y, H2 = bg['f'][:, k], bg['y'][:, k], (bg['H'][k] / 100)**2
    rho_phi = (0.5 * H2 * y[1]**2 + m.V(f)) / 3.0               # ω units
    return rho_phi / H2 / (2 * bg['F'][k])

cases = []
for lam in [0.5, 1.0, 2.0]:
    for U0 in [0.05, 0.20, 0.50]:
        cases.append(('exp', {'lam': lam}, f'exp λ={lam} U₀={U0:.2f}', 0.01, U0))
for n in [1, 2]:
    for U0 in [0.05, 0.20, 0.50]:
        cases.append(('power', {'n': n}, f'power n={n} U₀={U0:.2f}', 1.0, U0))
for U0 in [0.05, 0.20, 0.50]:
    cases.append(('hilltop', {'mu': 1.0}, f'hilltop U₀={U0:.2f}', 0.01, U0))

print("=" * 72)
print("SIM143: Bi-scalar Ψ + quintessence φ   [audited]")
print("=" * 72)
ref = {}
m0 = model('exp', {'lam': 0.5}, 0.0)
for mode, bg in (('faithful', m0.calibrated([PSI_BAR, 0.01], [0, 0], h_target)),
                 ('theta_matched', m0.theta_matched([PSI_BAR, 0.01], [0, 0], h_target))):
    ref[mode] = drv.evaluate_mf(m0, bg, 0)
    print(f"  Phase 1 reference (U₀=0) {mode:13s}: H0={ref[mode]['H0']:.2f} χ²DR1={ref[mode]['DR1']['chi2_DESI']:.2f} "
          f"Δχ²DR1 vs ΛCDM={ref[mode]['DR1']['delta_chi2_tot_vs_LCDM']:+.2f}")

def verdict(o, ede, mode):
    if o is None: return 'FAIL_NO_SOLUTION'
    if o['Lb'] < -0.05: return 'FAIL_NEGATIVE_LB'
    if ede > 0.05: return 'FAIL_EDE'
    if o['DR1']['gof_sigma'] >= 2.0: return 'FAIL_DESI'
    if o['chi2_theta'] >= 4.0: return 'FAIL_CMB'
    if o['rsd_chi2N'] >= 1.5: return 'FAIL_RSD'
    gain = o['DR1']['delta_chi2_tot_vs_LCDM'] - ref[mode]['DR1']['delta_chi2_tot_vs_LCDM']
    return 'TRIVIAL_PASS' if gain < 1.0 else 'PASS_PHENOMENOLOGICAL'

results = []
for form, p, label, phi_ini, U0 in cases:
    m = model(form, p, U0)
    rec = dict(form=form, params=p, label=label, U0=U0, phi_ini=phi_ini)
    for mode in ('faithful', 'theta_matched'):
        bg = m.calibrated([PSI_BAR, phi_ini], [0, 0], h_target) if mode == 'faithful' else \
             m.theta_matched([PSI_BAR, phi_ini], [0, 0], h_target)
        o = drv.evaluate_mf(m, bg, 2)
        ede = omega_phi_drag(m, bg) if bg is not None else np.nan
        if o is not None:
            o['Omega_phi_drag'] = float(ede)
        rec[mode] = o; rec[f'verdict_{mode}'] = verdict(o, ede, mode)
        if o is None:
            print(f"  {label:24s} {mode:13s} no solution"); continue
        print(f"  {label:24s} {mode:13s} H0={o['H0']:6.2f} 100θ*={o['theta_star_100']:.4f} χ²DR1={o['DR1']['chi2_DESI']:6.2f} "
              f"Δχ²DR1={o['DR1']['delta_chi2_tot_vs_LCDM']:+7.2f} Δχ²DR2={o['DR2']['delta_chi2_tot_vs_LCDM']:+7.2f} "
              f"RSD={o['rsd_chi2N']:.2f} Λb={o['Lb']:+.3f} → {rec[f'verdict_{mode}']}")
    results.append(rec)

def overall(mode):
    vs = [r[f'verdict_{mode}'] for r in results]
    if 'PASS_PHENOMENOLOGICAL' in vs: return 'PASS_PHENOMENOLOGICAL'
    if 'TRIVIAL_PASS' in vs: return 'TRIVIAL_PASS'
    return 'FAIL'
ov = dict(faithful=overall('faithful'), theta_matched=overall('theta_matched'))
print("=" * 72)
print(f"SIM143 RESULT (audited): faithful = {ov['faithful']}, θ*-matched = {ov['theta_matched']}   (published: FAIL/STRUCTURAL)")
print("  Any pass still inherits Ψ's unscreened coupling (Cassini ~50×, see SIM107/SIM130).")
print("=" * 72)
with open(os.path.join(HERE, 'output.json'), 'w') as f:
    json.dump(dict(sim_id='SIM143', timestamp=datetime.utcnow().isoformat() + 'Z', verdict=ov['theta_matched'],
                   verdict_modes=ov, phase1_reference=ref,
                   lcdm_reference={rel: obs.lcdm_reference(rel) for rel in ('DR1', 'DR2')},
                   scan_results=results, published_pre_audit=dict(verdict='FAIL/STRUCTURAL'),
                   audit_note='Audited 2026-10-22; see AUDIT_NOTES.md.'), f, indent=2, default=float)
print("SIM143 complete.")
