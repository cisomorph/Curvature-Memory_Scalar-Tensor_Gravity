"""
SIM123 — CMSTG Phase 3: Running F(z) — Full Background ODE Evolution
====================================================================
Hilltop (SSB) potential as the dark energy, no bare Λ:
  S = ∫√-g [ F(Ψ)R − ½(∂Ψ)² − V(Ψ) ] + S_m,  F = ½ + Λ₀Ψ²,  V = λ(Ψ² − v²)²,  v = 13.16 M_Pl
The field starts near the hilltop at z_ini = 1e5 and rolls; F(z) = ½ + Λ₀Ψ²(z) evolves.
Scan (as published): Λ₀ ∈ linspace(0.001, 0.008, 10), Ψ_ini ∈ {0.001 … 2.0} (10 values).

AUDIT 2026-10-22 (see AUDIT_NOTES.md). Published equations:
  3F H² = ρ_m + ρ_r + V   (ρ = ω: no GR limit; H² = ω/1.5 at F = ½)
  3H² dΨ/dN ≈ −V' + 24Λ₀ΨH²   (first-order slow roll with R = 12H², the de Sitter
                               value; R/H² = 6(2−ε) is 3 in matter and ~0 in radiation)
Corrected (full second-order, ρ in GR units, densities ω):
  H²(2F + 2F'Ψ' − Ψ'²/6) = ω_m a⁻³ + ω_r a⁻⁴ + V/3
  Ψ'' + (3−ε)Ψ' = F'·6(2−ε) − V'/H²,   ε constraint-consistent
λ is calibrated so that H(z=0) = 100h (faithful: h = 0.674, ≈ the published Ω_DE = 0.685
normalisation) or with h solved from θ* (θ*-matched). Observables: official DESI
likelihood, r_s to ∞, θ* = r_s(z*)/D_M(z*), proper statistics. (The published Part C fine
grid is not repeated; the coarse grid already spans the published ranges.)
"""
import numpy as np
import os, sys, warnings
warnings.filterwarnings('ignore')

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', '..', 'common'))
import cmstg_obs as obs
import cmstg_bg as bgm
from st_solver import ScalarTensorBG, verdict_from
import audit_driver as drv

OUT = os.path.join(HERE, '..', 'Outputs')
os.makedirs(OUT, exist_ok=True)
h_target, v = 0.674, 13.16
PUBLISHED = dict(map=dict(L0=0.001, psi_ini=0.001, H0=56.59, theta=0.88865), verdict='FAIL')


class HilltopBG(ScalarTensorBG):
    def __init__(self, L0):
        self.lam = 1e-6
        super().__init__(F=lambda p, N: 0.5 + L0 * p * p, F_psi=lambda p, N: 2 * L0 * p,
                         V=lambda p: self.lam * (p * p - v * v)**2,
                         V_psi=lambda p: 4 * self.lam * (p * p - v * v) * p)

    def calibrated(self, psi_ini, y_ini, h, n_iter=12, tol=1e-9):
        """No bare Λ: tune λ so that H(z=0) = 100h."""
        def resid(lmb):
            self.lam = lmb
            bg = self.integrate(psi_ini, y_ini, 0.0)
            return np.nan if bg is None else (bg['H0'] / 100.0)**2 - h * h
        lam0 = 3 * (2 * self.F(psi_ini, 0) * h * h - self.om - self.orr) / v**4
        lam = bgm.secant_match(resid, lam0, lam0 * 1.02, tol=1e-12 * lam0)
        if not np.isfinite(lam):
            return None
        self.lam = lam
        bg = self.integrate(psi_ini, y_ini, 0.0)
        if bg is not None:
            bg['h_target'] = h; bg['lambda'] = lam
        return bg


print("=" * 72)
print("SIM123 — Running F(z) with hilltop DE (no bare Λ)   [audited]")
print("=" * 72)
A = drv.run_case(HilltopBG(0.003), 0.01, 0.0, h_target, 0, 'ref Λ₀=0.003, Ψ_ini=0.01')
print(drv.line(A))
cases = []
for L0 in np.linspace(0.001, 0.008, 10):
    for psi_ini in [0.001, 0.01, 0.05, 0.1, 0.3, 0.5, 0.8, 1.0, 1.5, 2.0]:
        c = drv.run_case(HilltopBG(L0), psi_ini, 0.0, h_target, 2, f'Λ₀={L0:.4f} Ψi={psi_ini}')
        c['Lambda0'], c['psi_ini'] = float(L0), psi_ini
        cases.append(c); print(drv.line(c), flush=True)

bf, bm = drv.best(cases, 'faithful'), drv.best(cases, 'theta_matched')
v_ = verdict_from(bf['faithful'] if bf else None, bm['theta_matched'] if bm else None, 2)
print("=" * 72)
print(f"SIM123 RESULT (audited): faithful = {v_['faithful']}, θ*-matched = {v_['theta_matched']}   (published: FAIL)")
if bm:
    o = bm['theta_matched']
    print(f"  θ*-matched best: {bm['label']}: H0 = {o['H0']:.2f}, Ψ0 = {o['psi0']:.3f}, F0 = {o['F0']:.4f}, "
          f"Δχ² = {o['DR1']['delta_chi2_tot_vs_LCDM']:+.2f} (DR1) / {o['DR2']['delta_chi2_tot_vs_LCDM']:+.2f} (DR2)")
    lg = drv.local_gravity(o['F0'], [2 * bm['Lambda0'] * o['psi0']], H0=o['H0'])
    print(f"  local gravity at best: |γ−1| = {lg['gamma_minus_1']:.2e} ({lg['cassini_ratio']:.1f}× Cassini)")
print("=" * 72)
drv.save_json(os.path.join(OUT, 'sim123_results.json'), dict(
    verdict=v_['theta_matched'], verdict_modes=v_, reference=drv.strip(A), scan=[drv.strip(c) for c in cases],
    best_faithful=drv.strip(bf) if bf else None, best_theta_matched=drv.strip(bm) if bm else None,
    lcdm_reference={rel: obs.lcdm_reference(rel) for rel in ('DR1', 'DR2')},
    published_pre_audit=PUBLISHED, audit_note='Audited 2026-10-22; see AUDIT_NOTES.md.'))
print("SIM123 complete.")
