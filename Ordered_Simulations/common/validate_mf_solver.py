#!/usr/bin/env python3
"""Validation of mf_solver.py: single-field agreement with st_solver, 2-field GR limit, ΛCDM growth."""
import numpy as np, json, os
import cmstg_obs as obs
from mf_solver import MultiFieldBG
from st_solver import ScalarTensorBG
out = {}
# 1. single field vs st_solver (SIM130, Ψ_ini = 2.62)
L0 = 0.003
mf = MultiFieldBG(1, F=lambda f,N: 0.5+L0*f[0]**2, F_i=lambda f,N: np.array([2*L0*f[0]]))
st = ScalarTensorBG(F=lambda p,N: 0.5+L0*p*p, F_psi=lambda p,N: 2*L0*p)
bm = mf.theta_matched([2.62],[0.0]); bs = st.theta_matched(2.62,0.0)
lik = obs.DesiBAO('DR1')
c_mf = lik.chi2(mf.H_func(bm), obs.r_drag(mf.H_func(bm))); c_st = lik.chi2(st.H_func(bs), obs.r_drag(st.H_func(bs)))
print(f"[1] SIM130 Ψ_ini=2.62 θ*-matched: multi-field H0={bm['H0']:.4f} χ²={c_mf:.4f} | single-field H0={bs['H0']:.4f} χ²={c_st:.4f}")
out['single_vs_multi'] = dict(H0=[bm['H0'], bs['H0']], chi2=[c_mf, c_st])
# 2. two-field GR limit
gr = MultiFieldBG(2, F=lambda f,N: 0.5, F_i=lambda f,N: np.zeros(2))
b = gr.calibrated([0.3, 0.0],[0.0,0.0], 0.674)
dev = np.max(np.abs(b['H']/obs.H_lcdm(0.674)(b['z'])-1))
print(f"[2] two-field GR limit: max|H/H_ΛCDM−1| = {dev:.1e}")
out['gr_limit_2field'] = dev
# 3. growth in the GR limit (θ*-matched → ΛCDM reference)
b = gr.theta_matched([0.0,0.0],[0.0,0.0])
fs8, s8 = gr.fsigma8(b, np.array([0.0, 0.5, 1.0]))
Om = obs.OMEGA_M/(b['H0']/100)**2
print(f"[3] GR growth: σ8 = {s8:.4f} (input 0.811), fσ8(z=0,0.5,1) = {np.round(fs8,4)}, Ω_m = {Om:.4f}, "
      f"f(0) ≈ Ω_m^0.55 = {Om**0.55:.4f} vs {fs8[0]/s8:.4f}")
out['growth_gr'] = dict(s8=s8, fs8=list(fs8))
json.dump(out, open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'validation_mf_output.json'), 'w'), indent=2, default=float)
