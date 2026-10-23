"""Observational pipeline must reproduce published DESI and Planck results (ROADMAP 1.1)."""
import numpy as np
import pytest
from scipy.optimize import brentq, minimize

import cmstg_obs as obs

# DESI BAO-only flat ΛCDM: (Ω_m, σ, h·r_d [Mpc], σ)
PUBLISHED = {'DR1': (0.295, 0.015, 101.8, 1.3),     # DESI 2024 VI
             'DR2': (0.2975, 0.0086, 101.54, 0.73)}  # DESI 2025 DR2 II


def _bao_only_fit(release):
    lik = obs.DesiBAO(release)

    def chi2(p):
        om, hrd = p
        if not 0.05 < om < 0.7:
            return 1e10
        return lik.chi2(lambda z: 100 * np.sqrt(om * (1 + np.asarray(z))**3 + 1 - om), hrd)

    best = minimize(chi2, [0.3, 100.0], method='Nelder-Mead', options=dict(xatol=1e-7, fatol=1e-9))
    return best.x, best.fun, chi2


@pytest.mark.parametrize('release', ['DR1', 'DR2'])
def test_desi_bao_only_reproduces_published(release):
    (om, hrd), _, _ = _bao_only_fit(release)
    p_om, s_om, p_hrd, s_hrd = PUBLISHED[release]
    assert abs(om - p_om) < 0.25 * s_om
    assert abs(hrd - p_hrd) < 0.25 * s_hrd


def test_desi_data_shapes():
    for rel, n in (('DR1', 12), ('DR2', 13)):
        lik = obs.DesiBAO(rel)
        assert lik.n == n and lik.cov.shape == (n, n)
        assert np.all(np.linalg.eigvalsh(lik.cov) > 0)


def test_h0_from_theta_matches_planck():
    h = obs.h_from_theta(lambda hh: obs.H_lcdm(hh))
    assert abs(h - 0.6736) < 0.0054          # Planck 2018 TT,TE,EE+lowE+lensing


def test_sound_horizon_calibration():
    H = obs.H_lcdm(obs.H0_PLANCK / 100)
    raw_star = obs.sound_horizon(H, obs.Z_STAR, calibrated=False)
    raw_drag = obs.sound_horizon(H, obs.Z_DRAG, calibrated=False)
    assert abs(raw_star / obs.R_STAR_PLANCK - 1) < 2e-3      # simple integral within 0.2%
    assert abs(raw_drag / obs.R_DRAG_PLANCK - 1) < 2e-3
    assert abs(obs.r_drag(H) - obs.R_DRAG_PLANCK) < 1e-6


def test_sound_horizon_not_truncated():
    """Regression for bug B3: the tail above z = 1e4 is ~27 Mpc and must be included."""
    H = obs.H_lcdm(obs.H0_PLANCK / 100)
    assert obs.r_drag(H) > 146.0


def test_significance_helpers():
    p, s = obs.chi2_to_sigma(41.49, 6)
    assert abs(s - 5.17) < 0.01
    p, s = obs.chi2_to_sigma(228390.0, 1)
    assert np.isfinite(s) and 470 < s < 490
    assert obs.delta_chi2_to_sigma(-3.0, 2) == 0.0
    assert abs(obs.delta_chi2_to_sigma(1.0, 1) - 1.0) < 1e-6
