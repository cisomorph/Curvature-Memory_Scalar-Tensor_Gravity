#!/usr/bin/env python3
"""
SIM124 — CMSTG Phase 4: Non-canonical kinetic term K(Ψ) = exp(-αΨ²)
Question (as published): with K(Ψ) = exp(−αΨ²) < 1 (reduced Hubble friction), can the
hilltop field roll from near Ψ ≈ 0 at z_ini = 1500 to Ψ₀ = 2.62 today without fine-tuning?
Action: S = ∫√-g [F(Ψ)R − ½K(Ψ)(∂Ψ)² − V(Ψ)] + S_m,  F = ½ + Λ₀Ψ²,  V = λ(Ψ² − v²)², v = 13.16
No bare Λ: λ is calibrated so that H(z=0) = H₀_ref (67.0, as published) or from θ*.

AUDIT 2026-10-22 (see AUDIT_NOTES.md). Published bugs:
  • Friedmann H² = (ρ + V)/(3F − ½Ku²) with ρ = ω (no GR limit)
  • the F'R curvature source was missing from the field equation although F(Ψ)R is in
    the action
  • θ*: r_s integrated from z = 0 to 1090 (the sound horizon is ∫_{z*}^∞)
  • yet another variant of the non-DESI H(z) table (different error bars)
Corrected: H²(2F + 2F'Ψ' − KΨ'²/6) = ω_m a⁻³ + ω_r a⁻⁴ + V/3;
           Ψ'' = [F'·6(2−ε) − V'/H² − ½K'Ψ'²]/K − (3−ε)Ψ', ε constraint-consistent;
           Ψ'_ini = 0 at z_ini = 1500 (published: slow-roll velocity, quickly forgotten).
Parts: A — Ψ(z=0) over (α, Ψ_ini); B — Ψ_ini(α) with Ψ(z=0) = 2.62; C — observables.
"""
import json, os, sys, warnings
import numpy as np
from scipy.optimize import brentq
warnings.filterwarnings('ignore')

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', '..', 'common'))
import cmstg_obs as obs
import cmstg_bg as bgm
from mf_solver import MultiFieldBG
import audit_driver as drv

OUT = os.path.join(HERE, '..', 'Outputs')
os.makedirs(OUT, exist_ok=True)
L0, v, PSI0_REF, H0_REF = 0.003, 13.16, 2.62, 67.0


class KHilltop(MultiFieldBG):
    def __init__(self, alpha):
        self.lam = 1e-6
        K = lambda f: np.array([np.exp(-alpha * f[0]**2)])
        super().__init__(1, F=lambda f, N: 0.5 + L0 * f[0]**2, F_i=lambda f, N: np.array([2 * L0 * f[0]]),
                         K=K, K_i=lambda f: np.array([-2 * alpha * f[0] * K(f)[0]]),
                         V=lambda f: self.lam * (f[0]**2 - v * v)**2,
                         V_i=lambda f: np.array([4 * self.lam * (f[0]**2 - v * v) * f[0]]),
                         z_ini=1500.0, n_points=1500)

    def calibrated(self, f_ini, y_ini, h, **kw):
        def resid(lmb):
            self.lam = lmb
            bg = self.integrate(f_ini, y_ini, 0.0)
            return np.nan if bg is None else (bg['H0'] / 100.0)**2 - h * h
        lam0 = 3 * (2 * self.F(np.asarray(f_ini, float), 0) * h * h - self.om - self.orr) / v**4
        lam = bgm.secant_match(resid, lam0, lam0 * 1.02, tol=1e-12 * lam0)
        if not np.isfinite(lam):
            return None
        self.lam = lam
        bg = self.integrate(f_ini, y_ini, 0.0)
        if bg is not None:
            bg['h_target'] = h; bg['lambda'] = lam
        return bg


import signal

class _Timeout(Exception):
    pass

def _alarm(signum, frame):
    raise _Timeout()

signal.signal(signal.SIGALRM, _alarm)
UNRESOLVED = []

def psi_today(alpha, psi_ini, h=H0_REF / 100, budget=20):
    """Ψ(z=0) for one case; cases that exceed the time budget are recorded as unresolved."""
    signal.alarm(budget)
    try:
        bg = KHilltop(alpha).calibrated([psi_ini], [0.0], h)
        return np.nan if bg is None else float(bg['f0'][0])
    except _Timeout:
        UNRESOLVED.append(dict(alpha=alpha, psi_ini=psi_ini))
        return np.nan
    finally:
        signal.alarm(0)

print("=" * 70)
print("SIM124 — K(Ψ) = exp(−αΨ²) kinetic modification   [audited]")
print("=" * 70)
print("\nPart A: Ψ(z=0) over (α, Ψ_ini)")
scan_A = []
alphas = [0.0, 0.1, 0.3, 0.5, 1.0, 1.5, 2.0, 3.0, 5.0]
psis = [0.001, 0.01, 0.05, 0.1, 0.3, 0.5, 1.0, 2.0, 2.5]
print("  α \\ Ψ_ini " + " ".join(f"{p:>7}" for p in psis))
for a in alphas:
    row = []
    for p in psis:
        row.append(psi_today(a, p))
    scan_A += [dict(alpha=a, psi_ini=p, psi_z0=r) for p, r in zip(psis, row)]
    print(f"  {a:9.2f} " + " ".join(f"{r:7.3f}" for r in row), flush=True)

print("\nPart B: Ψ_ini(α) such that Ψ(z=0) = 2.62")
viable = []
for a in [0.0, 0.1, 0.2, 0.3, 0.5, 0.7, 1.0, 1.5, 2.0, 3.0, 5.0]:
    tries = np.logspace(-3, np.log10(PSI0_REF * 0.99), 18)
    vals = np.array([psi_today(a, p) for p in tries])
    found = False
    for i in range(len(tries) - 1):
        if np.isfinite(vals[i]) and np.isfinite(vals[i + 1]) and (vals[i] - PSI0_REF) * (vals[i + 1] - PSI0_REF) < 0:
            try:
                ps = brentq(lambda p: psi_today(a, p) - PSI0_REF, tries[i], tries[i + 1], xtol=1e-4, maxiter=30)
            except ValueError:
                continue
            viable.append(dict(alpha=a, psi_ini_sol=float(ps))); found = True
            print(f"  α={a:.2f}: Ψ_ini = {ps:.4f} → Ψ(z=0) = 2.62  (ΔΨ = {PSI0_REF - ps:.3f})", flush=True)
            break
    if not found:
        print(f"  α={a:.2f}: no crossing (max Ψ(z=0) = {np.nanmax(vals):.3f})", flush=True)

print("\nPart C: observables for the viable solutions (θ*-matched; λ recalibrated)")
results_C = []
for sol in viable:
    m = KHilltop(sol['alpha'])
    signal.alarm(600)
    try:
        bg = m.theta_matched([sol['psi_ini_sol']], [0.0], h_guess=H0_REF / 100)
    except _Timeout:
        UNRESOLVED.append(dict(alpha=sol['alpha'], psi_ini=sol['psi_ini_sol'], stage='C')); bg = None
    finally:
        signal.alarm(0)
    o = drv.evaluate_mf(m, bg, 2, rsd=False)
    if o is None:
        print(f"  α={sol['alpha']:.2f}: no θ* solution"); continue
    lg = drv.local_gravity(o['F_eff_today'], [2 * L0 * o['f0'][0]], K0=[np.exp(-sol['alpha'] * o['f0'][0]**2)],
                           Fdot_over_HF=2 * L0 * o['f0'][0] * o['y0'][0] / o['F_eff_today'], H0=o['H0'])
    o['local_gravity'] = lg
    v_ = 'PASS' if (o['DR1']['gof_sigma'] < 2 and o['chi2_theta'] < 4) else ('PARTIAL' if o['DR1']['delta_chi2_tot_vs_LCDM'] > 0 else 'FAIL')
    results_C.append(dict(**sol, obs=o, verdict=v_, fine_tuning_dPsi=PSI0_REF - sol['psi_ini_sol']))
    print(f"  α={sol['alpha']:.2f} Ψ_ini={sol['psi_ini_sol']:.3f}: H0={o['H0']:.2f} Ψ0={o['f0'][0]:.3f} "
          f"χ²DR1={o['DR1']['chi2_DESI']:.2f} Δχ²={o['DR1']['delta_chi2_tot_vs_LCDM']:+.2f} (DR2 {o['DR2']['delta_chi2_tot_vs_LCDM']:+.2f}) "
          f"Cassini×{lg['cassini_ratio']:.0f} → {v_}", flush=True)

near_hilltop = [r for r in results_C if r['psi_ini_sol'] < 0.5]
verdict = 'FAIL' if not near_hilltop else 'PARTIAL'
print("=" * 70)
print(f"  unresolved (time budget exceeded): {UNRESOLVED}")
print(f"SIM124 RESULT (audited): {verdict} — viable roll-ups need Ψ_ini ∈ "
      f"[{min([r['psi_ini_sol'] for r in results_C], default=np.nan):.2f}, {max([r['psi_ini_sol'] for r in results_C], default=np.nan):.2f}] "
      f"(fine-tuning criterion: Ψ_ini < 0.5 counts as genuine roll from the hilltop)   (published: FAIL)")
print("=" * 70)
drv.save_json(os.path.join(OUT, 'sim124_results.json'), dict(
    verdict=verdict, scan_A=scan_A, viable=viable, results_C=results_C, unresolved_cases=UNRESOLVED,
    published_pre_audit=dict(verdict='FAIL', note="'viable' solutions needed Ψ_ini ∈ [1.84, 2.57]"),
    audit_note='Audited 2026-10-22; see AUDIT_NOTES.md.'))
print("SIM124 complete.")
