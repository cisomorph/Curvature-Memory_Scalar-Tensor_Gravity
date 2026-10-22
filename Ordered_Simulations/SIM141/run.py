#!/usr/bin/env python3
"""
SIM141 — P4-B: Running Λ₀(a) (Brans-Dicke analog)

Phenomenological extension: Λ₀ carries cosmological-scale a-dependence on top of
its RG value (0.003 at the CMB epoch), decreasing after recombination.
  F(Ψ,a) = ½ + Λ₀(a)·Ψ²  →  F_eff weakens at late times
Forms scanned:
  (A) Linear:      Λ₀(a) = Λ₀_CMB · (1 − β(a − a_CMB))
  (B) Exponential: Λ₀(a) = Λ₀_CMB · exp(−β(a − a_CMB))
  (C) Tanh:        Λ₀(a) = Λ₀_CMB · [1 − γ·½(1 + tanh((a−a_trans)/σ))]
Λ₀(a) is NOT derived from the action — phenomenological (flagged).

Corrected equations (audit; ρ in ω units with the published "Planck-calibrated"
normalisation ρ̂ = 2F₀·ω so that F = F₀ reproduces ΛCDM-Planck at early times):
  H²[2F + 2(F_Ψ Ψ' + F_N) − Ψ'²/6] = 2F₀(ω_m a⁻³ + ω_r a⁻⁴) + Λ_b,   F_N = Ψ²·dΛ₀/dN
  Ψ'' + (3−ε)Ψ' = F_Ψ·6(2−ε),   ε constraint-consistent
(published: coef = 3F + 6F_ΨΨ' + 6F_N − ½Ψ'² with rhs 3F₀ω — F_ΨΨ', F_N and kinetic
terms each 2× too large relative to 3F.)

PASS criteria (as published):
  (a) CMB: |100θ* − Planck| < 2σ
  (b) DESI GOF < 2σ           (published: "tension" = √(χ²/N) < 2)
  (c) SIM105 RG: Λ₀(a) ≥ 0, perturbative; tanh form UV-consistent
  (d) RSD fσ₈ χ²/N < 1.5      (growth now includes G_eff(z) of the scalar-tensor theory)
  (e) F_eff > 0, no ghost
AUDIT 2026-10-22 (AUDIT_NOTES.md): corrected Friedmann/ε, official DESI likelihood, r_s
computed (published: hardcoded r_s = 144.7 Mpc, so θ* could not respond to the model),
θ* = r_s(z*)/D_M(z*), proper statistics, faithful (H₀ = 67.59) + θ*-matched modes.
The published "θ* 63σ" structural anti-correlation came from holding BOTH H₀ and r_s fixed.
"""

import json, os, sys, warnings
from datetime import datetime
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
warnings.filterwarnings('ignore')

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', 'common'))
import cmstg_obs as obs
from mf_solver import MultiFieldBG

FIG_DIR = os.path.join(HERE, 'figures')
os.makedirs(FIG_DIR, exist_ok=True)

h_target    = 0.6759
Lambda0_CMB = 0.003
PSI_BAR     = 2.62
F0          = 0.5 + Lambda0_CMB * PSI_BAR**2
a_CMB       = 1.0 / (1.0 + 1089.8)

# RSD fσ₈ compilation (as published; not re-verified in the audit)
RSD_Z   = np.array([0.02, 0.067, 0.10, 0.17, 0.22, 0.25, 0.30, 0.37, 0.41, 0.57, 0.60, 0.77, 0.80, 1.40])
RSD_FS8 = np.array([0.428, 0.423, 0.370, 0.510, 0.420, 0.351, 0.407, 0.460, 0.450, 0.427, 0.480, 0.490, 0.470, 0.482])
RSD_ERR = np.array([0.047, 0.055, 0.130, 0.060, 0.070, 0.058, 0.055, 0.038, 0.040, 0.023, 0.100, 0.080, 0.080, 0.116])

def L0_func(a, form, p):
    if form == 'constant':
        return Lambda0_CMB
    if form == 'linear':
        return max(Lambda0_CMB * (1.0 - p['beta'] * (a - a_CMB)), 0.0)
    if form == 'exponential':
        return Lambda0_CMB * np.exp(-p['beta'] * (a - a_CMB))
    return Lambda0_CMB * (1.0 - p['gamma'] * 0.5 * (1.0 + np.tanh((a - p['a_trans']) / p['sigma'])))

def dL0_dN(a, form, p):
    if form == 'constant':
        return 0.0
    if form == 'linear':
        return 0.0 if L0_func(a, form, p) <= 0 else a * (-Lambda0_CMB * p['beta'])
    if form == 'exponential':
        return a * (-Lambda0_CMB * p['beta'] * np.exp(-p['beta'] * (a - a_CMB)))
    arg = (a - p['a_trans']) / p['sigma']
    return a * (-Lambda0_CMB * p['gamma'] * 0.5 / np.cosh(arg)**2 / p['sigma'])

def model(form, p):
    return MultiFieldBG(1,
        F=lambda f, N: 0.5 + L0_func(np.exp(N), form, p) * f[0]**2,
        F_i=lambda f, N: np.array([2.0 * L0_func(np.exp(N), form, p) * f[0]]),
        F_N=lambda f, N: dL0_dN(np.exp(N), form, p) * f[0]**2,
        rho_scale=2.0 * F0)

def sim105_check(form, p, L0_today):
    L0c = L0_func(a_CMB, form, p)
    frac = (L0c - L0_today) / L0c if L0c > 0 else 0.0
    return dict(uv_consistent=(form in ('tanh', 'constant')), always_positive=bool(L0_today >= 0),
                perturbative=bool(frac < 1.0), change_frac=float(frac), L0_today=float(L0_today))

def evaluate(m, bg, form, p, n_params):
    if bg is None:
        return None
    Hf = m.H_func(bg); rd = obs.r_drag(Hf); th = obs.theta_star_100(Hf)
    o = dict(H0=bg['H0'], Psi0=float(bg['f0'][0]), F_eff_today=bg['F0'], theta_star_100=th,
             chi2_theta=obs.chi2_theta(th), r_drag=rd, L0_today=L0_func(1.0, form, p))
    for rel in ('DR1', 'DR2'):
        lik = obs.DesiBAO(rel); c2 = lik.chi2(Hf, rd); ref = obs.lcdm_reference(rel)
        s = obs.summarize_fit(c2, lik.n, n_params, o['chi2_theta'])
        s['delta_chi2_tot_vs_LCDM'] = (ref['chi2_DESI'] + ref['chi2_theta']) - (c2 + o['chi2_theta'])
        o[rel] = s
    fs8, s8 = m.fsigma8(bg, RSD_Z)
    o['sigma8'] = float(s8); o['rsd_chi2N'] = float(np.sum(((fs8 - RSD_FS8) / RSD_ERR)**2) / len(RSD_Z))
    return o

def verdict(o, rg):
    if o is None:
        return 'FAIL_NO_SOLUTION'
    if not rg['always_positive']:
        return 'FAIL_UNPHYSICAL'
    if o['DR1']['gof_sigma'] >= 2.0:
        return 'FAIL_DESI'
    if o['chi2_theta'] >= 4.0:
        return 'FAIL_CMB'
    if o['rsd_chi2N'] >= 1.5:
        return 'FAIL_RSD'
    if not rg['uv_consistent']:
        return 'PARTIAL_PHENOM'
    return 'PASS'

print("=" * 72)
print("SIM141: Running Λ₀(a) — P4-B Brans-Dicke analog   [audited]")
print("=" * 72)
for rel in ('DR1', 'DR2'):
    r = obs.lcdm_reference(rel)
    print(f"  ΛCDM reference, DESI {rel}: χ² = {r['chi2_DESI']:.2f}/{r['N_DESI']} (GOF {r['gof_sigma']:.2f}σ), H₀ = {r['H0']:.2f}")
print("  ⚠ Λ₀(a) is PHENOMENOLOGICAL — not derived from the action.\n")

cases = [('constant', {}, 'Phase 1 (const Λ₀)')]
for beta in [0.10, 0.30, 0.60, 1.00]:
    cases.append(('linear', {'beta': beta}, f'linear β={beta:.2f}'))
for beta in [0.10, 0.50, 1.00, 5.00]:
    cases.append(('exponential', {'beta': beta}, f'exp β={beta:.2f}'))
for a_tr in [0.33, 0.50, 0.67]:
    for gam in [0.50, 1.00]:
        cases.append(('tanh', {'a_trans': a_tr, 'sigma': 0.20, 'gamma': gam}, f'tanh a_tr={a_tr:.2f} γ={gam:.2f}'))

results = []
hdr = f"{'case':>22} | {'mode':>7} {'H0':>6} {'100θ*':>7} {'χ²DR1':>7} {'GOF':>5} {'Δχ²DR1':>8} {'Δχ²DR2':>8} {'RSD':>5} {'verdict':>16}"
print(hdr); print('-' * len(hdr))
for form, p, label in cases:
    m = model(form, p)
    k = 0 if form == 'constant' else (1 if form != 'tanh' else 2)
    bf = m.calibrated([PSI_BAR], [0.0], h_target)
    bm = m.theta_matched([PSI_BAR], [0.0], h_guess=h_target)
    rec = dict(form=form, params=p, label=label, n_params=k)
    for mode, bg in (('faithful', bf), ('theta_matched', bm)):
        o = evaluate(m, bg, form, p, k)
        rg = sim105_check(form, p, L0_func(1.0, form, p))
        v = verdict(o, rg)
        rec[mode] = o; rec[f'verdict_{mode}'] = v; rec['rg_consistency'] = rg
        if o is None:
            print(f"{label:>22} | {mode[:7]:>7} no solution"); continue
        print(f"{label:>22} | {mode[:7]:>7} {o['H0']:6.2f} {o['theta_star_100']:7.4f} {o['DR1']['chi2_DESI']:7.2f} "
              f"{o['DR1']['gof_sigma']:5.2f} {o['DR1']['delta_chi2_tot_vs_LCDM']:+8.2f} "
              f"{o['DR2']['delta_chi2_tot_vs_LCDM']:+8.2f} {o['rsd_chi2N']:5.2f} {v:>16}")
    results.append(rec)

def overall(mode):
    vs = [r[f'verdict_{mode}'] for r in results if r['form'] != 'constant']
    if 'PASS' in vs: return 'PASS'
    if 'PARTIAL_PHENOM' in vs: return 'PARTIAL'
    return 'FAIL'
ov = dict(faithful=overall('faithful'), theta_matched=overall('theta_matched'))
best = min((r for r in results if r['theta_matched']), key=lambda r: r['theta_matched']['DR1']['chi2_DESI'])
print("\n" + "=" * 72)
print(f"SIM141 RESULT (audited): faithful = {ov['faithful']}, θ*-matched = {ov['theta_matched']}   (published: PARTIAL, 'θ* 63σ')")
o = best['theta_matched']
print(f"  θ*-matched best DESI: {best['label']}: H0 = {o['H0']:.2f}, χ²(DR1) = {o['DR1']['chi2_DESI']:.2f}, "
      f"Δχ² vs ΛCDM = {o['DR1']['delta_chi2_tot_vs_LCDM']:+.2f} (DR1) / {o['DR2']['delta_chi2_tot_vs_LCDM']:+.2f} (DR2), "
      f"RSD χ²/N = {o['rsd_chi2N']:.2f}, verdict {best['verdict_theta_matched']}")
print("=" * 72)

# Figure: Δχ² vs ΛCDM per case, both modes
fig, ax = plt.subplots(figsize=(11, 5))
xs = np.arange(len(results))
for mode, mk, col in (('faithful', 'o', '#d73027'), ('theta_matched', 's', '#1a9850')):
    ys = [r[mode]['DR1']['delta_chi2_tot_vs_LCDM'] if r[mode] else np.nan for r in results]
    ax.plot(xs, ys, mk, color=col, label=mode.replace('_', ' '))
ax.axhline(0, color='k', ls=':')
ax.set_xticks(xs); ax.set_xticklabels([r['label'] for r in results], rotation=60, ha='right', fontsize=8)
ax.set_ylabel('Δχ²_tot vs ΛCDM (DESI DR1 + θ*; > 0 better)'); ax.set_yscale('symlog'); ax.legend(); ax.grid(alpha=.3)
ax.set_title('SIM141 (audited) — running Λ₀(a)')
fig.tight_layout()
for ext in ('pdf', 'png'):
    fig.savefig(os.path.join(FIG_DIR, f'sim141_running_L0.{ext}'), bbox_inches='tight', dpi=150)
plt.close(fig)

with open(os.path.join(HERE, 'output.json'), 'w') as f:
    json.dump(dict(sim_id='SIM141', timestamp=datetime.utcnow().isoformat() + 'Z',
                   action_spec='PHENOMENOLOGICAL: L = (½+Λ₀(a)Ψ²)R − ½(∂Ψ)², Λ₀(a) decreasing after CMB',
                   phenomenological_flag=True, verdict=ov['theta_matched'], verdict_modes=ov,
                   lcdm_reference={rel: obs.lcdm_reference(rel) for rel in ('DR1', 'DR2')},
                   scan_results=results,
                   published_pre_audit=dict(verdict='PARTIAL', best_desi_case='tanh a_tr=0.50 γ=1.00',
                                            tension_old_metric=0.547, theta_star_100=1.05939,
                                            note='θ* "63σ" with H0 and r_s both held fixed'),
                   audit_note='Audited 2026-10-22; see AUDIT_NOTES.md.'), f, indent=2, default=float)
print("SIM141 complete.")
