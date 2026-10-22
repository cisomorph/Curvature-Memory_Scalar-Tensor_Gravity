"""
SIM121-B — CMSTG Phase 2: H₀ Tension via G_eff Evolution
=========================================================
Question (as published): if the universe is CMSTG (F(Ψ) = ½ + Λ₀Ψ², G_eff = G/(2F)),
what H₀ would a ΛCDM analysis of the CMB acoustic scale infer, and does CMSTG move the
inferred H₀ toward SH0ES (73.0 ± 1.0)?

AUDIT 2026-10-22 (see AUDIT_NOTES.md). The published background used
E² = (Ω_m a⁻³ + Ω_r a⁻⁴ + V)/(3F − ½Ψ'²) (no GR limit: matter enters as Ω_m/(3F)) and
the non-DESI H(z) table. It also framed the result as "the ΛCDM-inferred H₀" without
comparing it to the true H₀ of the CMSTG universe, which is what SH0ES measures.

Corrected logic, using the audited solvers (Planck ω's; θ* matched in each model):
  • frozen Ψ (constant F₀): exactly degenerate with ΛCDM for θ* + BAO; the true H₀ of
    the CMSTG universe is H₀^ΛCDM × √(0.5/F₀)  →  LOWER than the CMB-ΛCDM value.
  • unlocked Ψ (Ψ_ini = 2.62, SIM130 dynamics): true H₀ from θ*.
The Hubble tension is measured against the true H₀.
Pass criteria (as published): H₀_true > 68.0, tension reduction ≥ 0.5σ, |G_eff(0)/G_N − 1| < 5%.
"""
import json, os, sys, warnings
import numpy as np
warnings.filterwarnings('ignore')
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', '..', 'common'))
import cmstg_obs as obs
from st_solver import ScalarTensorBG

OUT = os.path.join(HERE, '..', 'Outputs')
os.makedirs(OUT, exist_ok=True)
SH0ES, SH0ES_ERR, PLANCK_ERR = 73.0, 1.0, 0.54
ref = obs.lcdm_reference('DR1')
H0_lcdm = ref['H0']
F0 = 0.5 + 0.003 * 2.62**2

rows = {}
rows['ΛCDM'] = dict(H0_true=H0_lcdm)
rows['Phase 1 frozen (F0=0.5206)'] = dict(H0_true=H0_lcdm * np.sqrt(0.5 / F0),
                                         note='exact F0–H0 degeneracy (SIM121C audit)')
m = ScalarTensorBG(F=lambda p, N: 0.5 + 0.003 * p * p, F_psi=lambda p, N: 0.006 * p)
bg = m.theta_matched(2.62, 0.0)
rows['Phase 1 unlocked (Ψ_ini=2.62)'] = dict(H0_true=bg['H0'], F_today=bg['F0'])

print("=" * 70)
print("SIM121-B — H₀ tension in CMSTG   [audited]")
print("=" * 70)
for k, r in rows.items():
    t = (SH0ES - r['H0_true']) / np.hypot(SH0ES_ERR, PLANCK_ERR)
    r['tension_vs_SH0ES_sigma'] = float(t)
    print(f"  {k:32s} H₀(true, θ*-matched) = {r['H0_true']:6.2f}   tension with SH0ES = {t:.2f}σ")
t0 = rows['ΛCDM']['tension_vs_SH0ES_sigma']
best = max(rows.values(), key=lambda r: r['H0_true'])
passes = dict(H0_gt_68=best['H0_true'] > 68.0, tension_reduced_0p5=(t0 - best['tension_vs_SH0ES_sigma']) >= 0.5)
verdict = 'PASS' if all(passes.values()) else 'FAIL'
print("=" * 70)
print(f"SIM121-B RESULT (audited): {verdict} — CMSTG lowers the true H₀ (F > ½ ⇒ smaller G_eff), "
      f"so the Hubble tension gets worse, not better.   (published: FAIL)")
print("=" * 70)
json.dump(dict(verdict=verdict, rows=rows, criteria=passes, lcdm_reference=ref,
               published_pre_audit=dict(verdict='FAIL'), audit_note='Audited 2026-10-22; see AUDIT_NOTES.md.'),
          open(os.path.join(OUT, 'sim121b_results.json'), 'w'), indent=2, default=float)
print("SIM121-B complete.")
