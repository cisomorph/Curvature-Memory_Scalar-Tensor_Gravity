#!/usr/bin/env python3
"""
SIM142 — P4-A: Galileon G₃(Ψ)□Ψ sector

Adds G₃(Ψ)□Ψ to the Phase 1 canonical action:
  L = (½+Λ₀Ψ²)R − ½(∂Ψ)² + G₃(Ψ)□Ψ + L_matter
  (a) linear G₃ = c₃Ψ,  (b) quadratic G₃ = c₃Ψ²;  c₃ ∈ {1e-4, 1e-3, 1e-2, 1e-1}
G₃(Ψ) with no X-dependence is a total derivative plus a kinetic rescaling
(as the published note says): G₃□Ψ ≅ −G₃_Ψ(∂Ψ)² ⇒ −½K(Ψ)(∂Ψ)²,  K = 1 − 2G₃_Ψ
(Horndeski check: E₂+E₃ = X(1 − 2G₃_φ)).

Corrected equations (audit; Planck-calibrated normalisation ρ̂ = 2F₀·ω):
  H²(2F + 2F_ΨΨ' − KΨ'²/6) = 2F₀(ω_m a⁻³ + ω_r a⁻⁴) + Λ_b
  Ψ'' = [F_Ψ·6(2−ε) − ½K_ΨΨ'²]/K − (3−ε)Ψ',  ε constraint-consistent
Ghost-free requires K > 0. For the quadratic form K = 1 − 4c₃Ψ, so c₃ = 0.1 with
Ψ ≈ 2.6–2.9 gives K < 0: a ghost. The published "α_K > 0, no ghost" was wrong there.
(published Friedmann: coef = 3F + 6F_ΨΨ' − ½Ψ'² − 3G₃_ΨΨ' — the G₃ term must be
quadratic in Ψ', and F_ΨΨ' and kinetic terms were doubled; the EOM lacked the 1/K.)

PASS (as published): stability (K > 0), |Δθ*| < 2σ, DESI GOF < 2σ, RSD χ²/N < 1.5, c_T = c.
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

def model(c3, form):
    if form == 'linear':
        K, Kp = (lambda f: np.array([1.0 - 2.0 * c3])), (lambda f: np.array([0.0]))
    else:
        K, Kp = (lambda f: np.array([1.0 - 4.0 * c3 * f[0]])), (lambda f: np.array([-4.0 * c3]))
    return MultiFieldBG(1, F=lambda f, N: 0.5 + Lambda0 * f[0]**2,
                        F_i=lambda f, N: np.array([2.0 * Lambda0 * f[0]]), K=K, K_i=Kp, rho_scale=2.0 * F0)

def verdict(o, ghost):
    if ghost: return 'FAIL_GHOST'
    if o is None: return 'FAIL_NO_SOLUTION'
    if o['DR1']['gof_sigma'] >= 2.0: return 'FAIL_DESI'
    if o['chi2_theta'] >= 4.0: return 'FAIL_CMB'
    if o['rsd_chi2N'] >= 1.5: return 'FAIL_RSD'
    return 'PASS'

print("=" * 72)
print("SIM142: Galileon G₃(Ψ)□Ψ (kinetic rescaling)   [audited]")
print("=" * 72)
results = []
for form in ('linear', 'quadratic'):
    for c3 in [1e-4, 1e-3, 1e-2, 1e-1]:
        m = model(c3, form)
        K_ini = m.K(np.array([PSI_BAR]))[0]
        rec = dict(form=form, c3=c3, K_at_Psi_bar=float(K_ini))
        ghost = K_ini <= 0.0
        for mode in ('faithful', 'theta_matched'):
            bg = None
            if not ghost:
                bg = m.calibrated([PSI_BAR], [0.0], h_target) if mode == 'faithful' else m.theta_matched([PSI_BAR], [0.0], h_target)
                if bg is not None:
                    Kmin = min(m.K(bg['f'][:, k])[0] for k in range(bg['f'].shape[1]))
                    ghost = Kmin <= 0
            o = None if (ghost or bg is None) else drv.evaluate_mf(m, bg, 1)
            if o is not None:
                Fi = 2 * Lambda0 * o['f0'][0]
                o['local_gravity'] = drv.local_gravity(o['F_eff_today'], [Fi], K0=[m.K(np.array(o['f0']))[0]],
                                                       Fdot_over_HF=Fi * o['y0'][0] / o['F_eff_today'], H0=o['H0'])
            rec[mode] = o; rec[f'verdict_{mode}'] = verdict(o, ghost)
            if o is None:
                print(f"  {form:9s} c₃={c3:.0e} {mode:13s} K(Ψ̄)={K_ini:+.3f} → {rec[f'verdict_{mode}']}"); continue
            print(f"  {form:9s} c₃={c3:.0e} {mode:13s} H0={o['H0']:6.2f} 100θ*={o['theta_star_100']:.4f} "
                  f"χ²DR1={o['DR1']['chi2_DESI']:6.2f} Δχ²DR1={o['DR1']['delta_chi2_tot_vs_LCDM']:+7.2f} "
                  f"Δχ²DR2={o['DR2']['delta_chi2_tot_vs_LCDM']:+7.2f} RSD={o['rsd_chi2N']:.2f} "
                  f"Cassini×{o['local_gravity']['cassini_ratio']:.0f} → {rec[f'verdict_{mode}']}")
        results.append(rec)

def overall(mode):
    vs = [r[f'verdict_{mode}'] for r in results]
    return 'PASS' if 'PASS' in vs else 'FAIL'
ov = dict(faithful=overall('faithful'), theta_matched=overall('theta_matched'))
print("=" * 72)
print(f"SIM142 RESULT (audited): faithful = {ov['faithful']}, θ*-matched = {ov['theta_matched']}   (published: FAIL/STRUCTURAL)")
print("  Cosmology passes here are inherited from the Phase 1 unlocked background (SIM130);")
print("  Cassini (×~50) still excludes them — see local_gravity in output.json.")
print("=" * 72)
with open(os.path.join(HERE, 'output.json'), 'w') as f:
    json.dump(dict(sim_id='SIM142', timestamp=datetime.utcnow().isoformat() + 'Z', verdict=ov['theta_matched'],
                   verdict_modes=ov, lcdm_reference={rel: obs.lcdm_reference(rel) for rel in ('DR1', 'DR2')},
                   scan_results=results, published_pre_audit=dict(verdict='FAIL/STRUCTURAL'),
                   audit_note='Audited 2026-10-22; see AUDIT_NOTES.md.'), f, indent=2, default=float)
print("SIM142 complete.")
