"""
SIM131 — CMSTG Phase 3: Curvature-Memory Scalar (Conformal Coupling)
====================================================================
Phase 3 opens a first-principles redesign of the CMSTG action. Direction D:
replace V(Ψ) with a linear curvature-memory coupling ξΨR.

Motivation (as published):
  - Linear ξΨR term sources Ψ from R even at Ψ=0 (breaks the frozen attractor)
  - Ψ = □⁻¹R: causal integral of curvature over past light cone — spacetime memory
  - ξ = 1/6 (conformal value) — no new free parameter

Phase 3 action (Jordan frame, 8πG = 1, GR is F = ½):
  S = ∫d⁴x√(-g) [ F_eff(Ψ)·R − ½(∂Ψ)² ] + S_SM + bare Λ
  F_eff(Ψ) = (1 + 2Λ₀Ψ²)/2 + ξΨ

Equations (FRW, N = ln a; densities ω = Ωh², GR: H² = Σω) — corrected in the audit:
  Friedmann: H²(2F + 2F'Ψ' − Ψ'²/6) = ω_m a⁻³ + ω_r a⁻⁴ + Λ_bare
  Field:     Ψ'' + (3−ε)Ψ' = F'·R/H²,   R/H² = 6(2−ε)
  ε from the constraint along the trajectory (equals the Raychaudhuri closed form).

Structure:
  A. Reference run: ξ=1/6, Λ₀=0.003, Ψ_ini=0
  B. Λ₀ scan: ξ=1/6, Λ₀ ∈ {0.001, 0.002, 0.003, 0.005, 0.008, 0.010}
  C. ξ scan: Λ₀=0.003, ξ ∈ {1/20, 1/12, 1/6, 1/4, 1/3}
  D. Best joint fit and verdict

AUDIT 2026-10-22 (see AUDIT_NOTES.md):
  1. Friedmann equation: published H²(3F + 6F'Ψ' − ½Ψ'²) = ω… has no GR limit
     (F = ½ gives H² = ω/1.5, i.e. ΛCDM with Ω_m = 0.21) and doubles the F'Ψ' and
     kinetic terms. Replaced by the form derived from the action (above).
  2. ε_H was ∂_N ln H² at fixed (Ψ, Ψ'); now includes the field's evolution.
  3. Observables: official DESI likelihood (not the hardcoded H(z) table), r_s to
     z = ∞, θ* = r_s(z*)/D_M(z*), Planck 100θ* = 1.04110 ± 0.00031.
  4. "tension" = √(χ²/N) replaced by χ², p-value, GOF σ, and Δχ² vs ΛCDM run through
     the same pipeline (the hardcoded "SIM121C χ² = 18.26 / 2.77σ" baseline is gone;
     the corrected Phase 1 frozen model is identical to ΛCDM for BAO+θ*).
  5. Each case is run faithfully (H₀ = 67.4 fixed, as published) and θ*-matched
     (H₀ solved from θ*), because a fixed H₀ makes θ* fail for any change of H(z).
Verdict thresholds unchanged: PASS = GOF < 2σ and |Δθ*| < 2σ; PARTIAL = beats ΛCDM.
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
plt.rcParams.update({'font.family': 'serif', 'font.size': 11, 'axes.labelsize': 12, 'legend.fontsize': 9})

h_target     = 0.674
XI_CONFORMAL = 1.0 / 6.0
PSI0_phase1  = 2.62
F0_phase1    = 0.5 + 0.003 * PSI0_phase1**2

PUBLISHED = dict(reference=dict(psi0=2.9090, Feff0=1.0102, theta=0.79661, chi2=213.92, tension=5.97),
                 map=dict(xi=0.05, Lambda0=0.003, theta=0.90242, chi2=91.33, tension=3.90),
                 verdict='FAIL', readme='FAIL — θ* violation 5.97σ')


def model(xi, L0):
    return ScalarTensorBG(F=lambda p, N: 0.5 + L0 * p * p + xi * p,
                          F_psi=lambda p, N: 2.0 * L0 * p + xi)


print("=" * 72)
print("SIM131 — CMSTG Phase 3: Curvature-Memory Scalar (ξΨR coupling)   [audited]")
print("=" * 72)
for rel in ('DR1', 'DR2'):
    r = obs.lcdm_reference(rel)
    print(f"  ΛCDM reference (same pipeline), DESI {rel}: χ² = {r['chi2_DESI']:.2f}/{r['N_DESI']} "
          f"(GOF {r['gof_sigma']:.2f}σ), H₀ = {r['H0']:.2f}")
print()

# ── Part A ───────────────────────────────────────────────────────────────────
print("Part A: Reference run (ξ=1/6, Λ₀=0.003, Ψ_ini=0)")
A = drv.run_case(model(XI_CONFORMAL, 0.003), 0.0, 0.0, h_target, 0, 'ξ=1/6, Λ₀=0.003')
print(drv.line(A))
bgA = A['bg_faithful']
if bgA is not None:
    print(f"    Ψ(z*) = {np.interp(1089.92, bgA['z'], bgA['psi']):.4f}, Ψ₀ = {bgA['psi0']:.4f}, "
          f"F(z*) = {np.interp(1089.92, bgA['z'], bgA['F']):.5f}, F₀ = {bgA['F0']:.5f}  (Phase 1 F₀ = {F0_phase1:.5f})")

# ── Part B ───────────────────────────────────────────────────────────────────
print("\nPart B: Λ₀ scan (ξ=1/6 fixed, Ψ_ini=0)")
B = []
for L0 in [0.001, 0.002, 0.003, 0.005, 0.008, 0.010]:
    c = drv.run_case(model(XI_CONFORMAL, L0), 0.0, 0.0, h_target, 1, f'Λ₀={L0}')
    c['Lambda0'], c['xi'] = L0, XI_CONFORMAL
    B.append(c); print(drv.line(c))

# ── Part C ───────────────────────────────────────────────────────────────────
print("\nPart C: ξ scan (Λ₀=0.003 fixed, Ψ_ini=0)")
C = []
for xi, lbl in zip([0.05, 1 / 12, 1 / 6, 1 / 4, 1 / 3], ['1/20', '1/12', '1/6', '1/4', '1/3']):
    c = drv.run_case(model(xi, 0.003), 0.0, 0.0, h_target, 1, f'ξ={lbl}')
    c['Lambda0'], c['xi'] = 0.003, xi
    C.append(c); print(drv.line(c))

# ── Part D ───────────────────────────────────────────────────────────────────
allc = B + C
bf = drv.best(allc, 'faithful')
bm = drv.best(allc, 'theta_matched')
v = verdict_from(bf['faithful'] if bf else None, bm['theta_matched'] if bm else None, 2)
print("\nPart D: best points")
if bf: print("  faithful best    :", drv.line(bf))
if bm: print("  θ*-matched best  :", drv.line(bm))
for rel in ('DR1', 'DR2'):
    if bm:
        s = bm['theta_matched'][rel]
        print(f"  DESI {rel} (θ*-matched best): χ² = {s['chi2_DESI']:.2f}, GOF {s['gof_sigma']:.2f}σ, "
              f"Δχ² vs ΛCDM = {s['delta_chi2_tot_vs_LCDM']:+.2f}")
verdict = v['theta_matched']
print("\n" + "=" * 72)
print(f"SIM131 RESULT (audited): faithful = {v['faithful']},  θ*-matched = {v['theta_matched']}")
print(f"  Published: {PUBLISHED['verdict']} (reference 'tension' {PUBLISHED['reference']['tension']}σ, "
      f"100θ* {PUBLISHED['reference']['theta']})")
print("=" * 72)

# ── Figures ──────────────────────────────────────────────────────────────────
if bgA is not None:
    fig, ax = plt.subplots(1, 2, figsize=(13, 5))
    m = bgA['z'] < 3000
    ax[0].semilogx(bgA['z'][m] + 1, bgA['psi'][m], color='#d73027', lw=2, label='Ψ(z), ξ=1/6, Λ₀=0.003')
    ax[0].axhline(PSI0_phase1, color='#2166ac', ls='--', label='Phase 1 Ψ̄ = 2.62')
    ax[0].set_xlabel('1+z'); ax[0].set_ylabel('Ψ [M_Pl]'); ax[0].legend(); ax[0].grid(alpha=.3)
    ax[1].semilogx(bgA['z'][m] + 1, bgA['F'][m], color='#d73027', lw=2, label='F_eff(z)')
    ax[1].axhline(0.5, color='gray', ls=':', label='GR (F=½)')
    ax[1].set_xlabel('1+z'); ax[1].set_ylabel('F_eff'); ax[1].legend(); ax[1].grid(alpha=.3)
    fig.suptitle('SIM131 (audited) — curvature-memory Ψ from Ψ_ini = 0')
    fig.tight_layout()
    for ext in ('pdf', 'png'):
        fig.savefig(os.path.join(OUT, f'sim131_psi_evolution.{ext}'), dpi=150, bbox_inches='tight')
    plt.close(fig)

fig, ax = plt.subplots(1, 2, figsize=(13, 5))
for axi, cases, key, xl in ((ax[0], B, 'Lambda0', 'Λ₀ (ξ = 1/6)'), (ax[1], C, 'xi', 'ξ (Λ₀ = 0.003)')):
    for mode, mk, col in (('faithful', 'o-', '#d73027'), ('theta_matched', 's--', '#1a9850')):
        xs = [c[key] for c in cases if c[mode]]
        ys = [c[mode]['DR1']['delta_chi2_tot_vs_LCDM'] for c in cases if c[mode]]
        axi.plot(xs, ys, mk, color=col, label=mode.replace('_', ' '))
    axi.axhline(0, color='k', lw=1, ls=':', label='ΛCDM')
    axi.set_xlabel(xl); axi.set_ylabel('Δχ²_tot vs ΛCDM  (> 0 is better)'); axi.set_yscale('symlog')
    axi.legend(); axi.grid(alpha=.3)
fig.suptitle('SIM131 (audited) — DESI DR1 + θ* relative to ΛCDM')
fig.tight_layout()
for ext in ('pdf', 'png'):
    fig.savefig(os.path.join(OUT, f'sim131_tension_landscape.{ext}'), dpi=150, bbox_inches='tight')
plt.close(fig)

if bm is not None:
    mdl = model(bm['xi'], bm['Lambda0'])
    fig, ax = plt.subplots(figsize=(7, 5))
    drv.plot_desi_ratio(ax, [(f"best ξ={bm['xi']:.3f}, Λ₀={bm['Lambda0']}", mdl.H_func(bm['bg_matched']))])
    ax.set_title('SIM131 (audited) — θ*-matched best point vs DESI DR1')
    fig.tight_layout()
    for ext in ('pdf', 'png'):
        fig.savefig(os.path.join(OUT, f'sim131_Hz_weff.{ext}'), dpi=150, bbox_inches='tight')
    plt.close(fig)

drv.save_json(os.path.join(OUT, 'sim131_results.json'), dict(
    sim='SIM131', phase='Phase 3', verdict=verdict, verdict_modes=v,
    action='S = ∫d⁴x√g [F_eff(Ψ)·R − ½(∂Ψ)²] + S_SM + Λ_bare, F_eff = (1+2Λ₀Ψ²)/2 + ξΨ',
    lcdm_reference={rel: obs.lcdm_reference(rel) for rel in ('DR1', 'DR2')},
    reference_run=drv.strip(A), scan_Lambda0=[drv.strip(c) for c in B], scan_xi=[drv.strip(c) for c in C],
    best_faithful=drv.strip(bf) if bf else None, best_theta_matched=drv.strip(bm) if bm else None,
    published_pre_audit=PUBLISHED,
    audit_note='Audited 2026-10-22: corrected Friedmann/ε, official DESI likelihood, r_s to ∞, θ* at z*, '
               'proper statistics, faithful + θ*-matched modes. See AUDIT_NOTES.md.'))
print(f"\nAll outputs in: {OUT}\nSIM131 complete.")
