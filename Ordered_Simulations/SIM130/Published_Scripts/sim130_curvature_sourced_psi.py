"""
SIM130 — CMSTG Phase 1: Curvature-Sourced Ψ Evolution (Unlocked Field)
======================================================================
Question: does the Phase 1 action itself, when Ψ is UNLOCKED, generate
curvature-sourced Ψ evolution that produces measurable H(z) deviation at DESI
redshifts, without any new parameters?

Phase 1 action (Jordan frame, 8πG = 1, GR is F = ½; README normalisation):
  S = ∫d⁴x√(-g) [F(Ψ) R − ½(∂Ψ)² − Λ_bare] + S_m,   F(Ψ) = ½ + Λ₀Ψ²
(The published docstring wrote "F(Ψ)/2 R", inconsistent with F = ½ + Λ₀Ψ²; the
README form F·R is used, which is also what the published code integrated.)

Equations (FRW, N = ln a, densities ω = Ωh², GR: H² = Σω) — corrected in the audit:
  Friedmann: H²(2F + 2F'Ψ' − Ψ'²/6) = ω_m a⁻³ + ω_r a⁻⁴ + Λ_bare
  Field:     Ψ'' + (3−ε)Ψ' = 2Λ₀Ψ·R/H²,  R/H² = 6(2−ε),  ε constraint-consistent

Scan: Ψ_ini ∈ [0.1, 3.5] M_Pl at z_ini = 10⁵, Ψ'_ini = 0; Λ₀ = 0.003 (locked).

AUDIT 2026-10-22 (see AUDIT_NOTES.md): corrected Friedmann (published form had no GR
limit: H² = ω/(3F)), constraint-consistent ε, official DESI likelihood, r_s to ∞,
θ* = r_s(z*)/D_M(z*), proper statistics, ΛCDM baseline in the same pipeline (replaces
the hardcoded "SIM121C χ² = 18.26 / 2.77σ"), faithful (H₀ = 67.4) + θ*-matched modes.
Verdict thresholds unchanged.
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

Lambda0  = 0.003
h_target = 0.674
PSI0_ref = 2.62
PUBLISHED = dict(reference=dict(psi_ini=2.62, psi0=2.886, H0=67.37, chi2=84.26, tension=3.75, theta=0.92280),
                 map=dict(psi_ini=3.5, theta=0.92452, chi2=91.68, tension=3.91), verdict='FAIL')

model = ScalarTensorBG(F=lambda p, N: 0.5 + Lambda0 * p * p, F_psi=lambda p, N: 2.0 * Lambda0 * p)

print("=" * 72)
print("SIM130 — CMSTG Phase 1: Curvature-Sourced Ψ Evolution (Unlocked)   [audited]")
print("=" * 72)
for rel in ('DR1', 'DR2'):
    r = obs.lcdm_reference(rel)
    print(f"  ΛCDM reference, DESI {rel}: χ² = {r['chi2_DESI']:.2f}/{r['N_DESI']} (GOF {r['gof_sigma']:.2f}σ), H₀ = {r['H0']:.2f}")

print("\nPart A: Reference trajectory (Ψ_ini = 2.62)")
A = drv.run_case(model, PSI0_ref, 0.0, h_target, 0, 'Ψ_ini=2.62')
print(drv.line(A))
bgA = A['bg_faithful']
if bgA is not None:
    psi_cmb = np.interp(1089.92, bgA['z'], bgA['psi'])
    print(f"    Ψ(z*) = {psi_cmb:.5f}, Ψ₀ = {bgA['psi0']:.5f} ({100*(bgA['psi0']/PSI0_ref-1):+.2f}% from Ψ_ini), "
          f"F(z*) = {np.interp(1089.92, bgA['z'], bgA['F']):.5f}, F₀ = {bgA['F0']:.5f}")

print("\nPart C: Ψ_ini scan (Λ₀ = 0.003)")
grid = [0.1, 0.3, 0.5, 0.8, 1.0, 1.5, 2.0, 2.3, 2.62, 2.8, 3.0, 3.2, 3.5]
C = []
for p in grid:
    c = drv.run_case(model, p, 0.0, h_target, 1, f'Ψ_ini={p}')
    c['psi_ini'] = p
    C.append(c); print(drv.line(c))

bf, bm = drv.best(C, 'faithful'), drv.best(C, 'theta_matched')
v = verdict_from(bf['faithful'] if bf else None, bm['theta_matched'] if bm else None, 1)
print("\nPart D: best points")
if bf: print("  faithful best   :", drv.line(bf))
if bm: print("  θ*-matched best :", drv.line(bm))
if bm:
    for rel in ('DR1', 'DR2'):
        s = bm['theta_matched'][rel]
        print(f"  DESI {rel} (θ*-matched best): χ² = {s['chi2_DESI']:.2f}, GOF {s['gof_sigma']:.2f}σ, "
              f"Δχ² vs ΛCDM = {s['delta_chi2_tot_vs_LCDM']:+.2f}")
print("\n" + "=" * 72)
print(f"SIM130 RESULT (audited): faithful = {v['faithful']},  θ*-matched = {v['theta_matched']}")
print(f"  Published: {PUBLISHED['verdict']} (reference 100θ* = {PUBLISHED['reference']['theta']}, "
      f"'tension' {PUBLISHED['reference']['tension']}σ)")
print("=" * 72)

# Figures
if bgA is not None:
    fig, ax = plt.subplots(1, 2, figsize=(13, 5))
    m = bgA['z'] < 3000
    ax[0].semilogx(bgA['z'][m] + 1, bgA['psi'][m], color='#d73027', lw=2, label='Ψ(z), Ψ_ini = 2.62')
    ax[0].set_xlabel('1+z'); ax[0].set_ylabel('Ψ [M_Pl]'); ax[0].legend(); ax[0].grid(alpha=.3)
    ax[1].semilogx(bgA['z'][m] + 1, bgA['F'][m], color='#d73027', lw=2, label='F(z)')
    ax[1].axhline(0.5, color='gray', ls=':', label='GR')
    ax[1].set_xlabel('1+z'); ax[1].set_ylabel('F'); ax[1].legend(); ax[1].grid(alpha=.3)
    fig.suptitle('SIM130 (audited) — unlocked Phase 1 Ψ')
    fig.tight_layout()
    for ext in ('pdf', 'png'):
        fig.savefig(os.path.join(OUT, f'sim130_psi_evolution.{ext}'), dpi=150, bbox_inches='tight')
    plt.close(fig)

fig, ax = plt.subplots(figsize=(7, 5))
for mode, mk, col in (('faithful', 'o-', '#d73027'), ('theta_matched', 's--', '#1a9850')):
    xs = [c['psi_ini'] for c in C if c[mode]]
    ys = [c[mode]['DR1']['delta_chi2_tot_vs_LCDM'] for c in C if c[mode]]
    ax.plot(xs, ys, mk, color=col, label=mode.replace('_', ' '))
ax.axhline(0, color='k', ls=':', label='ΛCDM')
ax.set_xlabel('Ψ_ini'); ax.set_ylabel('Δχ²_tot vs ΛCDM (> 0 better)'); ax.set_yscale('symlog'); ax.legend(); ax.grid(alpha=.3)
ax.set_title('SIM130 (audited) — DESI DR1 + θ*')
fig.tight_layout()
for ext in ('pdf', 'png'):
    fig.savefig(os.path.join(OUT, f'sim130_tension_scan.{ext}'), dpi=150, bbox_inches='tight')
plt.close(fig)

drv.save_json(os.path.join(OUT, 'sim130_results.json'), dict(
    sim='SIM130', verdict=v['theta_matched'], verdict_modes=v,
    description='Phase 1 CMSTG unlocked — curvature-sourced Psi evolution, no SSB, no new knobs',
    fixed_params=dict(Lambda0=Lambda0, h_target=h_target),
    lcdm_reference={rel: obs.lcdm_reference(rel) for rel in ('DR1', 'DR2')},
    reference_trajectory=drv.strip(A), scan=[drv.strip(c) for c in C],
    best_faithful=drv.strip(bf) if bf else None, best_theta_matched=drv.strip(bm) if bm else None,
    published_pre_audit=PUBLISHED,
    audit_note='Audited 2026-10-22; see AUDIT_NOTES.md.'))
print(f"\nAll outputs in: {OUT}\nSIM130 complete.")
