"""
audit_driver.py — shared run/report helpers for the audited Phase 3/4 sims (2026-10-22).

Every case is evaluated twice:
  faithful      : H0 fixed to the sim's published h_target (Λ_bare calibrated), as published
  theta-matched : H0 solved so that 100θ* equals Planck (profile over H0; the fair test)
Both use the official DESI likelihood, r_s to z = ∞, θ* = r_s(z*)/D_M(z*), and the
ΛCDM reference computed in the same pipeline.
"""
import json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

import cmstg_obs as obs


def run_case(model, psi_ini, y_ini, h_target, n_params=0, label='', match=True):
    bf = model.calibrated(psi_ini, y_ini, h_target)
    of = model.observe(bf, n_params) if bf is not None else None
    bm = model.theta_matched(psi_ini, y_ini, h_guess=h_target) if match else None
    om = model.observe(bm, n_params) if bm is not None else None
    return dict(label=label, psi_ini=psi_ini, y_ini=y_ini, faithful=of, theta_matched=om,
                bg_faithful=bf, bg_matched=bm)


def line(case, rel='DR1'):
    def fmt(o):
        if o is None:
            return 'no solution'.ljust(52)
        s = o[rel]
        return (f"H0={o['H0']:6.2f} 100θ*={o['theta_star_100']:.4f} χ²={s['chi2_DESI']:7.2f} "
                f"GOF={s['gof_sigma']:5.2f}σ Δχ²={s['delta_chi2_tot_vs_LCDM']:+8.2f}")
    return f"{case['label']:>24} | faithful: {fmt(case['faithful'])} | θ*-matched: {fmt(case['theta_matched'])}"


def strip(case):
    """JSON-safe copy without the background arrays."""
    return {k: v for k, v in case.items() if not k.startswith('bg_')}


def best(cases, mode='theta_matched', rel='DR1'):
    ok = [c for c in cases if c[mode] is not None]
    if not ok:
        return None
    return min(ok, key=lambda c: c[mode][rel]['chi2_DESI'] + c[mode]['chi2_theta'])


def plot_desi_ratio(ax, model_curves, rel='DR1'):
    """DESI points and model curves relative to the ΛCDM reference (θ*-matched ΛCDM)."""
    lik = obs.DesiBAO(rel)
    ref = obs.lcdm_reference(rel)
    H_l = obs.H_lcdm(ref['h'])
    rd_l = obs.r_drag(H_l)
    m_l = lik.model(H_l, rd_l)
    zz = np.linspace(0.05, 2.5, 100)

    def curves(Hf, rd):
        DM = np.array([obs.comoving_distance(Hf, z, n=300) for z in zz])
        DH = obs.C_KMS / Hf(zz)
        return {'DM_over_rs': DM / rd, 'DH_over_rs': DH / rd, 'DV_over_rs': (zz * DM**2 * DH)**(1 / 3) / rd}
    cl = curves(H_l, rd_l)
    cols = {'DM_over_rs': '#d6604d', 'DH_over_rs': '#2166ac', 'DV_over_rs': '#4dac26'}
    for kind, mk in (('DM_over_rs', 'o'), ('DH_over_rs', 's'), ('DV_over_rs', '^')):
        sel = [i for i, k in enumerate(lik.kind) if k == kind]
        ax.errorbar(lik.z[sel], lik.value[sel] / m_l[sel], yerr=lik.sigma[sel] / m_l[sel], fmt=mk,
                    color=cols[kind], ms=5, capsize=3, label=f'DESI {rel} {kind.split("_")[0]}')
    styles = ['-', '--', ':', '-.']
    for j, (lab, Hf) in enumerate(model_curves):
        cm = curves(Hf, obs.r_drag(Hf))
        for kind in ('DM_over_rs', 'DH_over_rs'):
            ax.plot(zz, cm[kind] / cl[kind], color=cols[kind], ls=styles[j % 4], lw=1.6,
                    label=f'{lab} {kind.split("_")[0]}')
    ax.axhline(1, color='k', lw=0.8, ls=':')
    ax.set_xlabel('Redshift $z$')
    ax.set_ylabel(r'$(D/r_d)/(D/r_d)_{\Lambda{\rm CDM}}$')
    ax.legend(fontsize=7, ncol=2)


# RSD fσ₈ compilation used by the published Phase 4 sims (kept as published; not re-verified)
RSD_Z   = np.array([0.02, 0.067, 0.10, 0.17, 0.22, 0.25, 0.30, 0.37, 0.41, 0.57, 0.60, 0.77, 0.80, 1.40])
RSD_FS8 = np.array([0.428, 0.423, 0.370, 0.510, 0.420, 0.351, 0.407, 0.460, 0.450, 0.427, 0.480, 0.490, 0.470, 0.482])
RSD_ERR = np.array([0.047, 0.055, 0.130, 0.060, 0.070, 0.058, 0.055, 0.038, 0.040, 0.023, 0.100, 0.080, 0.080, 0.116])


def evaluate_mf(m, bg, n_params, rsd=True):
    """Observables for a MultiFieldBG background: θ*, DESI DR1/DR2, Δχ² vs ΛCDM, RSD."""
    if bg is None:
        return None
    Hf = m.H_func(bg); rd = obs.r_drag(Hf); th = obs.theta_star_100(Hf)
    o = dict(H0=bg['H0'], f0=[float(x) for x in bg['f0']], y0=[float(x) for x in bg['y0']],
             F_eff_today=bg['F0'], theta_star_100=th, chi2_theta=obs.chi2_theta(th), r_drag=rd, Lb=bg['Lb'])
    for rel in ('DR1', 'DR2'):
        lik = obs.DesiBAO(rel); c2 = lik.chi2(Hf, rd); ref = obs.lcdm_reference(rel)
        s = obs.summarize_fit(c2, lik.n, n_params, o['chi2_theta'])
        s['delta_chi2_tot_vs_LCDM'] = (ref['chi2_DESI'] + ref['chi2_theta']) - (c2 + o['chi2_theta'])
        o[rel] = s
    if rsd:
        fs8, s8 = m.fsigma8(bg, RSD_Z)
        o['sigma8'] = float(s8)
        o['rsd_chi2N'] = float(np.sum(((fs8 - RSD_FS8) / RSD_ERR)**2) / len(RSD_Z))
    return o


def local_gravity(F0, F_i0, K0=None, Fdot_over_HF=0.0, H0=67.4, alpha2_extra=0.0):
    """
    Solar-System checks at z = 0 for massless scalars (no screening):
      PPN |γ−1| = 2α²/(1+α²), α² = S/(3S + 2F̃·…) with F̃ = 2F (DEF), S = Σ F̃_i²/K_i
      (+ alpha2_extra for direct matter couplings), and Ġ/G ≈ −Ḟ/F.
    """
    F_i0 = np.atleast_1d(F_i0)
    K0 = np.ones_like(F_i0) if K0 is None else np.atleast_1d(K0)
    Ft = 2.0 * F0
    S = float(np.sum((2.0 * F_i0)**2 / K0))
    a2 = S / (3.0 * S + 2.0 * Ft) + alpha2_extra
    g1 = 2.0 * a2 / (1.0 + a2)
    H0yr = H0 / 3.0857e19 * 3.156e7
    gdot = -Fdot_over_HF * H0yr
    return dict(alpha2=a2, gamma_minus_1=g1, cassini_ratio=g1 / 2.3e-5,
                Gdot_over_G_per_yr=gdot, llr_sigma=(gdot - 7.1e-14) / 7.6e-14,
                cassini_ok=bool(g1 < 2.3e-5), llr_ok=bool(abs((gdot - 7.1e-14) / 7.6e-14) < 3.0))


def save_json(path, payload):
    with open(path, 'w') as f:
        json.dump(payload, f, indent=2, default=lambda x: float(x) if isinstance(x, (np.floating, np.integer)) else str(x))
