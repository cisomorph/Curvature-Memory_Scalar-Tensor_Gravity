#!/usr/bin/env python3
"""
Validation of common/cmstg_obs.py against published Planck 2018 and DESI numbers.
Run:  python3 Ordered_Simulations/common/validate_cmstg_obs.py
Writes common/validation_output.json.
"""
import json, os, sys
import numpy as np
from scipy.optimize import minimize
from scipy.integrate import quad

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import cmstg_obs as obs

out = {}
print("=" * 72)
print("cmstg_obs validation")
print("=" * 72)

# 1. Sound horizon: calibration size and an independent explanation of it
H_ref = obs.H_lcdm(obs.H0_PLANCK / 100)
rs_star_raw = obs.sound_horizon(H_ref, obs.Z_STAR, calibrated=False)
rs_drag_raw = obs.sound_horizon(H_ref, obs.Z_DRAG, calibrated=False)
print(f"\n[1] Planck-2018 reference LCDM (h=0.6736, omega_m=0.1430, simplified radiation)")
print(f"    raw r_s(z*)    = {rs_star_raw:.3f} Mpc  (Planck 144.43)  -> K_STAR = {obs.K_STAR:.5f}")
print(f"    raw r_s(z_d)   = {rs_drag_raw:.3f} Mpc  (Planck 147.09)  -> K_DRAG = {obs.K_DRAG:.5f}")

# Independent check: proper massive-neutrino background (1 x 0.06 eV + 2 massless,
# N_eff = 3.046) instead of counting the massive species as both radiation and matter.
Neff, m_nu = 3.046, 0.06
w_nu1 = (7 / 8) * (4 / 11)**(4 / 3) * (Neff / 3) * obs.OMEGA_GAMMA       # one relativistic species
T_nu0_eV = 8.617333e-5 * 2.7255 * (4 / 11)**(1 / 3) * (Neff / 3)**0.25
I0 = 7 * np.pi**4 / 120
def I(y):
    return quad(lambda q: q * q * np.sqrt(q * q + y * y) / (np.exp(q) + 1), 0, 60)[0]
lna_tab = np.linspace(np.log(1e-9), 0, 400)
ratio_tab = np.array([I(m_nu * np.exp(x) / T_nu0_eV) / I0 for x in lna_tab])
w_nu_massive_today = w_nu1 * ratio_tab[-1]
w_cb = obs.OMEGA_M - w_nu_massive_today
h = obs.H0_PLANCK / 100
w_de = h * h - w_cb - obs.OMEGA_GAMMA - 2 * w_nu1 - w_nu_massive_today
def H_nu(z):
    z = np.asarray(z, float); a = 1 / (1 + z)
    r = np.interp(np.log(a), lna_tab, ratio_tab)
    return 100 * np.sqrt(w_cb / a**3 + (obs.OMEGA_GAMMA + 2 * w_nu1) / a**4 + w_nu1 * r / a**4 + w_de)
rs_star_nu = obs.sound_horizon(H_nu, obs.Z_STAR, calibrated=False)
rs_drag_nu = obs.sound_horizon(H_nu, obs.Z_DRAG, calibrated=False)
print(f"    with proper massive-nu background: r_s(z*) = {rs_star_nu:.3f}, r_s(z_d) = {rs_drag_nu:.3f}")
print(f"    -> the calibration (~{100*(obs.K_STAR-1):.2f}%) is the neutrino double-count plus the")
print(f"       finite-visibility effects the simple integral omits (residual "
      f"{100*(144.43/rs_star_nu-1):+.2f}% / {100*(147.09/rs_drag_nu-1):+.2f}%)")
out['sound_horizon'] = dict(raw_r_star=rs_star_raw, raw_r_drag=rs_drag_raw,
                            K_STAR=obs.K_STAR, K_DRAG=obs.K_DRAG,
                            massive_nu_r_star=rs_star_nu, massive_nu_r_drag=rs_drag_nu)

# truncated-integral bug, for the record
def rs_truncated(H, z1, z2):
    h_ = obs.H0_PLANCK / 100
    Og = 2.469e-5 / h_**2
    f = lambda z: (obs.C_KMS / np.sqrt(3 * (1 + (3 * obs.OMEGA_B / h_**2) / (4 * Og * (1 + z))))) / H(np.array([z]))[0]
    return quad(f, z1, z2, limit=200, epsrel=1e-5)[0]
for zmax in (1e4, 1e5):
    print(f"    old code r_s(z_drag) truncated at z={zmax:.0e}: {rs_truncated(H_ref, 1059.6, zmax):.2f} Mpc")
out['old_truncated_rs'] = {f'{zmax:.0e}': rs_truncated(H_ref, 1059.6, zmax) for zmax in (1e4, 1e5)}

# 2. theta*
th = obs.theta_star_100(H_ref)
h_th = obs.h_from_theta(lambda hh: obs.H_lcdm(hh))
print(f"\n[2] 100theta*(h=0.6736) = {th:.5f}  (Planck {obs.THETA_STAR_100} +- {obs.THETA_STAR_100_ERR};"
      f" diff {(th-obs.THETA_STAR_100)/obs.THETA_STAR_100_ERR:+.2f} sigma)")
print(f"    h solved from theta*  = {h_th:.5f}  (H0 = {100*h_th:.2f}; Planck 67.36 +- 0.54)")
out['theta'] = dict(theta_at_planck_h=th, h_from_theta=h_th)

# 3. DESI likelihood: reproduce DESI's own LCDM BAO-only fits
print("\n[3] DESI BAO-only flat LCDM fit (free Omega_m, h*r_d)")
published = {'DR1': (0.295, 0.015, 101.8, 1.3), 'DR2': (0.2975, 0.0086, 101.54, 0.73)}
out['desi_bao_only'] = {}
for rel in ('DR1', 'DR2'):
    lik = obs.DesiBAO(rel)
    def chi2(p):
        Om, hrd = p
        if not (0.05 < Om < 0.7):
            return 1e10
        H = lambda z: 100 * np.sqrt(Om * (1 + np.asarray(z))**3 + 1 - Om)   # h = 1 shape
        return lik.chi2(H, hrd)
    best = minimize(chi2, [0.3, 100.0], method='Nelder-Mead', options=dict(xatol=1e-6, fatol=1e-8))
    Om, hrd = best.x
    # 1-sigma on Omega_m from the profile (delta chi2 = 1)
    prof = lambda om: minimize(lambda q: chi2([om, q[0]]), [hrd], method='Nelder-Mead').fun - best.fun - 1
    from scipy.optimize import brentq
    lo, hi = brentq(prof, Om - 0.1, Om), brentq(prof, Om, Om + 0.1)
    pOm, pOme, phrd, phrde = published[rel]
    print(f"    {rel}: Omega_m = {Om:.4f} (+{hi-Om:.4f}/-{Om-lo:.4f}), h*r_d = {hrd:.2f} Mpc, "
          f"chi2 = {best.fun:.2f} / {lik.n - 2} dof   | DESI paper: {pOm} +- {pOme}, {phrd} +- {phrde}")
    out['desi_bao_only'][rel] = dict(Omega_m=Om, Omega_m_err=(hi - lo) / 2, h_rd=hrd, chi2=best.fun,
                                     dof=lik.n - 2, published=published[rel])

# 4. LCDM reference in the sims' setup, and w0wa preference
print("\n[4] Sims' setup (omega_m, omega_b fixed to Planck; h from theta*): LCDM vs w0waCDM")
out['pipeline_reference'] = {}
for rel in ('DR1', 'DR2'):
    ref = obs.lcdm_reference(rel)
    lik = obs.DesiBAO(rel)
    def c2w(p):
        w0, wa = p
        if w0 + wa >= 0:
            return 1e10
        hh = obs.h_from_theta(lambda x: obs.H_lcdm(x, w0=w0, wa=wa))
        if not np.isfinite(hh):
            return 1e10
        H = obs.H_lcdm(hh, w0=w0, wa=wa)
        return lik.chi2(H, obs.r_drag(H))
    bw = minimize(c2w, [-0.8, -0.6], method='Nelder-Mead', options=dict(xatol=1e-4, fatol=1e-6))
    d = ref['chi2_DESI'] - bw.fun
    sig = obs.delta_chi2_to_sigma(d, 2)
    print(f"    {rel}: LCDM H0={ref['H0']:.2f}, r_d={ref['r_drag']:.2f}, chi2={ref['chi2_DESI']:.2f}/{ref['N_DESI']} "
          f"(p={ref['p_value']:.3f}, {ref['gof_sigma']:.2f} sigma GOF)")
    print(f"          w0waCDM best w0={bw.x[0]:.3f}, wa={bw.x[1]:.3f}, chi2={bw.fun:.2f}; "
          f"dchi2={d:.2f} -> {sig:.2f} sigma preference over LCDM")
    out['pipeline_reference'][rel] = dict(lcdm=ref, w0wa=dict(w0=bw.x[0], wa=bw.x[1], chi2=bw.fun,
                                          dchi2_vs_lcdm=d, sigma=sig))

# 5. Statistics helpers
print("\n[5] significance helper")
for c2, dof in ((41.49, 6), (41.49, 3), (213.92, 6), (228390.0, 1)):
    p, s = obs.chi2_to_sigma(c2, dof)
    print(f"    chi2={c2:>9.2f}, dof={dof}: p={p:.3e} -> {s:.2f} sigma   (old metric sqrt(chi2/N) = {np.sqrt(c2/6):.2f})")

with open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'validation_output.json'), 'w') as f:
    json.dump(out, f, indent=2, default=float)
print("\nSaved common/validation_output.json")
