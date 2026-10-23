"""Background solvers must recover GR exactly and be internally consistent (ROADMAP 1.1)."""
import numpy as np

import cmstg_obs as obs
import cmstg_bg as bgm
from st_solver import ScalarTensorBG
from mf_solver import MultiFieldBG


def test_gr_limit_single_field():
    m = ScalarTensorBG(F=lambda p, N: 0.5, F_psi=lambda p, N: 0.0)
    bg = m.calibrated(0.0, 0.0, 0.674)
    assert np.max(np.abs(bg['H'] / obs.H_lcdm(0.674)(bg['z']) - 1)) < 1e-10


def test_gr_limit_two_fields():
    m = MultiFieldBG(2, F=lambda f, N: 0.5, F_i=lambda f, N: np.zeros(2))
    bg = m.calibrated([0.3, 0.0], [0.0, 0.0], 0.674)
    assert np.max(np.abs(bg['H'] / obs.H_lcdm(0.674)(bg['z']) - 1)) < 1e-10


def test_theta_matched_gr_equals_lcdm_reference():
    m = ScalarTensorBG(F=lambda p, N: 0.5, F_psi=lambda p, N: 0.0)
    o = m.observe(m.theta_matched(0.0, 0.0))
    assert abs(o['DR1']['delta_chi2_tot_vs_LCDM']) < 1e-3
    assert abs(o['DR2']['delta_chi2_tot_vs_LCDM']) < 1e-3


def test_constraint_eps_equals_raychaudhuri():
    """SIM131 trajectory (ξ = 1/6): ε from the constraint vs the closed form."""
    L0, xi = 0.003, 1 / 6
    F = lambda p: 0.5 + L0 * p * p + xi * p
    Fp = lambda p: 2 * L0 * p + xi
    m = ScalarTensorBG(F=lambda p, N: F(p), F_psi=lambda p, N: Fp(p))
    bg = m.calibrated(0.0, 0.0, 0.674)
    H2f = lambda p, y, N: m.H2(p, y, N, bg['Lb'])
    worst = 0.0
    for i in range(0, len(bg['N']), 200):
        p, y, N = bg['psi'][i], bg['y'][i], bg['N'][i]
        a = np.exp(N)
        e_c = bgm.eps_from_constraint(H2f, p, y, N, 12 * Fp(p) - 3 * y, y - 6 * Fp(p))
        rho_p = 3 * (obs.OMEGA_M / a**3 + 4 / 3 * obs.OMEGA_R / a**4)
        e_k = bgm.eps_closed_form(F(p), Fp(p), 2 * L0, y, H2f(p, y, N), rho_p)
        worst = max(worst, abs(e_c - e_k))
    assert worst < 1e-8


def test_single_and_multi_field_agree():
    st = ScalarTensorBG(F=lambda p, N: 0.5 + 0.003 * p * p, F_psi=lambda p, N: 0.006 * p)
    mf = MultiFieldBG(1, F=lambda f, N: 0.5 + 0.003 * f[0]**2, F_i=lambda f, N: np.array([0.006 * f[0]]))
    a, b = st.theta_matched(2.62, 0.0), mf.theta_matched([2.62], [0.0])
    assert abs(a['H0'] - b['H0']) < 1e-6


def test_frozen_F0_is_degenerate_with_lcdm():
    """Bug-finding B11: with Ψ frozen, F0 only rescales H0 for BAO + θ*."""
    lik = obs.DesiBAO('DR1')
    out = []
    for F0 in (0.5, 0.5206):
        g = 0.5 / F0
        def make(h, g=g, F0=F0):
            Om, Or = obs.OMEGA_M / h**2, obs.OMEGA_R / h**2
            return lambda z: 100 * h * np.sqrt(g * (Om * (1 + np.asarray(z))**3 + Or * (1 + np.asarray(z))**4
                                                     + (2 * F0 - Om - Or)))
        h = obs.h_from_theta(make)
        H = make(h)
        out.append(lik.chi2(H, obs.r_drag(H)))
    assert abs(out[0] - out[1]) < 1e-6
