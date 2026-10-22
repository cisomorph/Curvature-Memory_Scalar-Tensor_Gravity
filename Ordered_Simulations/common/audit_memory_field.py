#!/usr/bin/env python3
"""
audit_memory_field.py — corrected observables for the memory-field toy models SIM126 and SIM129
(audit 2026-10-22; see AUDIT_NOTES.md).

Both published scripts compared to the non-DESI H(z) table (a variant with errors
[2.9, 3.8, 4.5, 5.3, 7.3, 15.7]). SIM129 also integrated r_s from z = 0 to z* (must be
z* → ∞) and used Ω_b where ω_b is needed in the baryon loading (R off by 1/h² ≈ 2.2).
The model equations are reproduced unchanged (phenomenological; no action):
  SIM126: E² = Ω_m a⁻³ + Ω_r a⁻⁴ + Ω_Λ + ε M²,  dM/da = αC/(aE) − βM/a
  SIM129: E² = Ω_m a⁻³ + Ω_r a⁻⁴ + Ω_Λ(1 − γM̂(a)),  M̂ from the ΛCDM-H solution, M̂(1) = 1
with H₀ = 67.4 fixed as published (note: E(1) ≠ 1 when ε, γ ≠ 0, so H(z=0) ≠ 67.4).
Above a = 0.01 (z ≈ 99) both models are ΛCDM; r_s uses that extension.
Reported: official DESI DR1/DR2 χ² with the model's own r_d, 100θ*, and Δχ² relative to the
same script's ε = 0 / γ = 0 case (ΛCDM at H₀ = 67.4, Ω_m = 0.315).
"""
import json, os, sys
import numpy as np
from scipy.integrate import solve_ivp
from scipy.interpolate import interp1d

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import cmstg_obs as obs

H0, Om, Or, OL = 67.4, 0.315, 9.4e-5, 0.685
rho_m0, rho_L = 3 * Om, 3 * OL
alpha, beta, xi, sigma0, a_eq, a_nl, n_sh = 1.0, 1.5, 10.0, 0.1, 3.0e-4, 0.3, 2
a_arr = np.logspace(-2, 0, 800)

def C_func(a):
    rm = rho_m0 * a**-3
    return rm * sigma0 * (a / a_eq)**n_sh * np.exp(-a / a_nl) / (rm + rho_L / xi)

def H_from_E(Ea):
    z = 1 / a_arr - 1
    return obs.H_from_table(z, H0 * Ea, omega_m=Om * (H0 / 100)**2, omega_r=Or * (H0 / 100)**2)

def evaluate(Hf):
    rd = obs.r_drag(Hf)
    out = dict(theta_star_100=obs.theta_star_100(Hf), r_drag=rd, H_z0=float(Hf(np.array([0.0]))[0]))
    for rel in ('DR1', 'DR2'):
        lik = obs.DesiBAO(rel); c2 = lik.chi2(Hf, rd)
        out[rel] = obs.summarize_fit(c2, lik.n, 1)
    return out

def sim126():
    rows = []
    for eps in [0.0, 2e-11, 2e-10, 2e-9, 2e-8]:
        def ode(a, y):
            E = np.sqrt(max(Om * a**-3 + Or * a**-4 + OL + eps * y[0]**2, 1e-30))
            return [alpha * C_func(a) / (a * E) - beta * y[0] / a]
        M = solve_ivp(ode, [a_arr[0], 1.0], [0.0], t_eval=a_arr, rtol=1e-9, atol=1e-12).y[0]
        E = np.sqrt(Om * a_arr**-3 + Or * a_arr**-4 + OL + eps * M**2)
        rows.append(dict(eps=eps, **evaluate(H_from_E(E))))
    return rows

def sim129():
    E_l = lambda a: np.sqrt(Om * a**-3 + Or * a**-4 + OL)
    M = solve_ivp(lambda a, y: [alpha * C_func(a) / (a * E_l(a)) - beta * y[0] / a], [a_arr[0], 1.0], [0.0],
                  t_eval=a_arr, rtol=1e-10, atol=1e-13).y[0]
    Mhat = M / M[-1]
    rows = []
    for g in [0.00, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.40, 0.50]:
        E = np.sqrt(Om * a_arr**-3 + Or * a_arr**-4 + OL * (1 - g * Mhat))
        rows.append(dict(gamma=g, **evaluate(H_from_E(E))))
    return rows

if __name__ == '__main__':
    root = os.path.join(HERE, '..')
    for name, fn, key in (('SIM126', sim126, 'eps'), ('SIM129', sim129, 'gamma')):
        rows = fn()
        base = rows[0]
        print(f"{name} (audited observables; H0 = 67.4 fixed as published)")
        for r in rows:
            for rel in ('DR1', 'DR2'):
                r[rel]['delta_chi2_vs_param0'] = base[rel]['chi2_DESI'] - r[rel]['chi2_DESI']
            print(f"  {key}={r[key]:<8g} H(z=0)={r['H_z0']:6.2f} 100θ*={r['theta_star_100']:.4f} "
                  f"χ²DR1={r['DR1']['chi2_DESI']:7.2f} (Δ vs {key}=0: {r['DR1']['delta_chi2_vs_param0']:+6.2f})  "
                  f"χ²DR2={r['DR2']['chi2_DESI']:7.2f} (Δ {r['DR2']['delta_chi2_vs_param0']:+6.2f})")
        best = max(rows, key=lambda r: r['DR1']['delta_chi2_vs_param0'])
        print(f"  best {key} = {best[key]}: Δχ² (DR1) = {best['DR1']['delta_chi2_vs_param0']:+.2f}, "
              f"100θ* = {best['theta_star_100']:.4f} (Planck 1.04110 ± 0.00031)\n")
        outdir = os.path.join(root, name, 'Outputs'); os.makedirs(outdir, exist_ok=True)
        json.dump(dict(sim=name, rows=rows, audit_note='Audited observables 2026-10-22; model unchanged; see AUDIT_NOTES.md.'),
                  open(os.path.join(outdir, f'{name.lower()}_audit.json'), 'w'), indent=2, default=float)
