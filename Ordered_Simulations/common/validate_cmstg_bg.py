#!/usr/bin/env python3
"""
Validation of the corrected background equations (common/cmstg_bg.py).
  1. GR limit: F = 1/2 reproduces flat ΛCDM H(z) exactly.
  2. Constraint-consistent ε equals the closed-form Raychaudhuri ε along a strongly
     evolving trajectory (SIM131 reference: F = ½ + Λ₀Ψ² + ξΨ, ξ = 1/6, Ψ_ini = 0).
  3. Size of the published-equation error on the same trajectory.
Writes common/validation_bg_output.json.
"""
import json, os, sys
import numpy as np
from scipy.integrate import solve_ivp

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import cmstg_obs as obs
import cmstg_bg as bgm

om, orr, h = obs.OMEGA_M, obs.OMEGA_R, 0.674

def model(L0, xi):
    F   = lambda p: 0.5 + L0 * p**2 + xi * p
    Fp  = lambda p: 2 * L0 * p + xi
    Fpp = lambda p: 2 * L0
    return F, Fp, Fpp

def H2_new(F, Fp, Lb):
    def H2(p, y, N):
        a = np.exp(N)
        return (om / a**3 + orr / a**4 + Lb) / (2 * F(p) + 2 * Fp(p) * y - y**2 / 6)
    return H2

def H2_old(F, Fp, Lb):
    def H2(p, y, N):
        a = np.exp(N)
        coef = 3 * F(p) + 6 * Fp(p) * y - 0.5 * y**2
        if coef < 0.3:
            coef = 3 * F(p)
        return (om / a**3 + orr / a**4 + Lb) / coef
    return H2

def solve(L0, xi, Lb, equations='new', z_ini=1e5, n=4000):
    F, Fp, Fpp = model(L0, xi)
    H2 = (H2_new if equations == 'new' else H2_old)(F, Fp, Lb)
    def rhs(N, s):
        p, y = s
        if equations == 'new':
            A, B = 12 * Fp(p) - 3 * y, y - 6 * Fp(p)
            eps = bgm.eps_from_constraint(H2, p, y, N, A, B)
        else:   # published: ∂_N ln H² at fixed (Ψ, y)
            d = 5e-4
            eps = -0.5 * (H2(p, y, N + d) - H2(p, y, N - d)) / (2 * d * H2(p, y, N))
        return [y, Fp(p) * 6 * (2 - eps) - (3 - eps) * y]
    Ni = np.log(1 / (1 + z_ini))
    sol = solve_ivp(rhs, (Ni, 0), [0.0, 0.0], t_eval=np.linspace(Ni, 0, n),
                    method='RK45', rtol=1e-9, atol=1e-12, max_step=0.05)
    p, y, N = sol.y[0], sol.y[1], sol.t
    Hh = np.array([np.sqrt(H2(p[i], y[i], N[i])) for i in range(len(N))])
    return dict(N=N, z=np.exp(-N) - 1, psi=p, y=y, H=100 * Hh, F=F(p), Fp=Fp(p), Fpp=np.full_like(p, Fpp(p)), H2=Hh**2)

def calibrated(L0, xi, equations='new'):
    """Iterate Λ_bare so that H(z=0) = 100 h (as the sims do)."""
    F, Fp, _ = model(L0, xi)
    Lb = (2 * 0.5 if equations == 'new' else 3 * 0.5) * h * h - om - orr
    for _ in range(8):
        r = solve(L0, xi, Lb, equations)
        p0, y0 = r['psi'][-1], r['y'][-1]
        coef = (2 * F(p0) + 2 * Fp(p0) * y0 - y0**2 / 6) if equations == 'new' \
               else (3 * F(p0) + 6 * Fp(p0) * y0 - 0.5 * y0**2)
        Lb_new = h * h * coef - om - orr
        if abs(Lb_new - Lb) < 1e-12:
            break
        Lb = Lb_new
    r = solve(L0, xi, Lb, equations); r['Lb'] = Lb
    return r

out = {}
print("=" * 72); print("cmstg_bg validation"); print("=" * 72)

# 1. GR limit
r = calibrated(0.0, 0.0, 'new')
H_ref = obs.H_lcdm(h)
dev = np.max(np.abs(r['H'] / H_ref(r['z']) - 1))
r_old = calibrated(0.0, 0.0, 'old')
dev_old = np.max(np.abs(r_old['H'] / H_ref(r_old['z']) - 1))
print(f"\n[1] GR limit (Λ0 = ξ = 0): max |H/H_ΛCDM − 1| = {dev:.2e} (corrected),  {dev_old:.3f} (published equations)")
Hf_old = obs.H_from_table(r_old['z'], r_old['H'])
print(f"    published GR limit is ΛCDM with effective Ω_m = {om/1.5/h**2:.3f} (not {om/h**2:.3f}); "
      f"100θ* = {obs.theta_star_100(Hf_old):.4f}, DESI DR1 χ² = "
      f"{obs.DesiBAO('DR1').chi2(Hf_old, obs.r_drag(Hf_old)):.1f}")
out['gr_limit'] = dict(max_dev_corrected=dev, max_dev_published=dev_old)

# 2. ε consistency on the SIM131 reference trajectory
r = calibrated(0.003, 1 / 6, 'new')
eps_c, eps_k = [], []
H2f = H2_new(*model(0.003, 1 / 6)[:2], r['Lb'])
for i in range(0, len(r['N']), 40):
    p, y, N = r['psi'][i], r['y'][i], r['N'][i]
    a = np.exp(N)
    Fp = r['Fp'][i]
    eps_c.append(bgm.eps_from_constraint(H2f, p, y, N, 12 * Fp - 3 * y, y - 6 * Fp))
    rho_p = 3 * (om / a**3 + (4 / 3) * orr / a**4)
    eps_k.append(bgm.eps_closed_form(r['F'][i], Fp, r['Fpp'][i], y, r['H2'][i], rho_p))
eps_c, eps_k = np.array(eps_c), np.array(eps_k)
print(f"\n[2] SIM131 reference (ξ=1/6, Λ0=0.003): Ψ0 = {r['psi'][-1]:.4f}, F0 = {r['F'][-1]:.4f}")
print(f"    constraint-ε vs Raychaudhuri closed-form ε: max |Δε| = {np.max(np.abs(eps_c-eps_k)):.2e} "
      f"over {len(eps_c)} points (ε range {eps_k.min():.3f}..{eps_k.max():.3f})")
out['eps_consistency'] = dict(max_abs_diff=float(np.max(np.abs(eps_c - eps_k))))

# 3. published vs corrected on the same model
ro = calibrated(0.003, 1 / 6, 'old')
for lab, rr in (('published', ro), ('corrected', r)):
    Hf = obs.H_from_table(rr['z'], rr['H'])
    th = obs.theta_star_100(Hf)
    c2 = obs.DesiBAO('DR1').chi2(Hf, obs.r_drag(Hf))
    print(f"    {lab:9s} equations: Ψ0 = {rr['psi'][-1]:.4f}, F0 = {rr['F'][-1]:.4f}, "
          f"100θ* = {th:.4f}, DESI DR1 χ² = {c2:.1f}  (both with the corrected observables)")
    out[f'sim131_ref_{lab}'] = dict(psi0=float(rr['psi'][-1]), F0=float(rr['F'][-1]), theta=th, chi2_DR1=c2)

with open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'validation_bg_output.json'), 'w') as f:
    json.dump(out, f, indent=2, default=float)
print("\nSaved common/validation_bg_output.json")
