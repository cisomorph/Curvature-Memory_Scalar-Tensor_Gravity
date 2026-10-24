"""
SIM156 background likelihood: compressed Planck 2018 + DESI DR2 BAO + Pantheon+ (SPEC.md).

    lik = BgLikelihood();  chi2 = lik.chi2(H_of_z, omega_b, omega_cdm, h)   -> dict of parts

- **Planck 2018 distance priors** (Chen, Huang & Wang 2019, Table I, ΛCDM column, TT,TE,EE+lowE). R, l_A and ω_b
  come with the 3×3 covariance built from σ = (0.0046, 0.090, 0.00015) and correlations R–l_A +0.46,
  R–ω_b −0.66, l_A–ω_b −0.33.
  - z* follows Hu & Sugiyama, as in the paper. R = √ω_m · 100 · D_M(z*)/c and l_A = π D_M(z*)/r_s(z*).
  - r_s is the integral in common/cmstg_obs.
  - One constant factor each for R and l_A is fixed at the Planck 2018 TT,TE,EE+lowE ΛCDM means
    (ω_b 0.02236, ω_cdm 0.1202, h 0.6727). It removes the simple-integral offset (validation V1).
- **DESI DR2**: common/cmstg_obs.DesiBAO('DR2'). r_d comes from the Aubourg et al. (2015) formula
  (0.02% accuracy; standard early universe, true for every SIM156 model).
- **Pantheon+**: exactly as Cobaya's sn.pantheonplus. That means m_b_corr, z_HD as z_cmb, z_HEL, the cut z_HD > 0.01,
  the STAT+SYS covariance, and the absolute magnitude marginalised by projecting the constant mode out of C⁻¹.
  The distance modulus is 5 log10[(1+z_hel) D_M(z_HD)/Mpc].
"""
import os, sys
import numpy as np
from scipy.integrate import cumulative_trapezoid
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', 'common'))
import cmstg_obs as obs

C = obs.C_KMS
OMEGA_NU = 0.06 / 93.14              # one 0.06 eV neutrino, counted as matter at late times
SN_DIR = '/home/aion/cobaya_packages/data/sn_data/PantheonPlus'

# Chen, Huang & Wang (2019), Table I, ΛCDM, Planck 2018 TT,TE,EE+lowE
DP_MEAN = np.array([1.7502, 301.471, 0.02236])
DP_SIG = np.array([0.0046, 0.090, 0.00015])
DP_CORR = np.array([[1.0, 0.46, -0.66], [0.46, 1.0, -0.33], [-0.66, -0.33, 1.0]])
DP_ICOV = np.linalg.inv(DP_CORR * np.outer(DP_SIG, DP_SIG))
PLANCK_MEANS = dict(omega_b=0.02236, omega_cdm=0.1202, h=0.6727)


def z_star_hs(omega_b, omega_m):
    """Hu & Sugiyama (1996) decoupling redshift, eqs. (8)-(10) of Chen et al."""
    g1 = 0.0738 * omega_b**-0.238 / (1 + 39.5 * omega_b**0.763)
    g2 = 0.560 / (1 + 21.1 * omega_b**1.81)
    return 1048.0 * (1 + 0.00124 * omega_b**-0.738) * (1 + g1 * omega_m**g2)


def r_drag_aubourg(omega_b, omega_cdm, omega_nu=OMEGA_NU):
    """Aubourg et al. (2015), eq. (16): r_d in Mpc for N_eff = 3.046."""
    return 55.154 * np.exp(-72.3 * (omega_nu + 0.0006)**2) / ((omega_b + omega_cdm)**0.25351 * omega_b**0.12807)


def _raw_compressed(H, omega_b, omega_m):
    zs = z_star_hs(omega_b, omega_m)
    DM = obs.comoving_distance(H, zs)
    rs = obs.sound_horizon(H, zs, omega_b, calibrated=False)
    return np.sqrt(omega_m) * 100.0 * DM / C, np.pi * DM / rs


def _calibrate():
    p = PLANCK_MEANS
    om = p['omega_b'] + p['omega_cdm'] + OMEGA_NU
    R, lA = _raw_compressed(obs.H_lcdm(p['h'], om, obs.OMEGA_R), p['omega_b'], om)
    return DP_MEAN[0] / R, DP_MEAN[1] / lA, R, lA


K_R, K_LA, R_RAW_PLANCK, LA_RAW_PLANCK = _calibrate()


class Pantheon:
    def __init__(self):
        with open(os.path.join(SN_DIR, 'Pantheon+SH0ES.dat')) as f:
            cols = [c.strip().lower() for c in f.readline().lstrip('#').split()]
        d = np.genfromtxt(os.path.join(SN_DIR, 'Pantheon+SH0ES.dat'), skip_header=1, dtype=float,
                          usecols=[cols.index('m_b_corr'), cols.index('zhd'), cols.index('zhel')])
        mag, zcmb, zhel = d[:, 0], d[:, 1], d[:, 2]
        cov = np.loadtxt(os.path.join(SN_DIR, 'Pantheon+SH0ES_STAT+SYS.cov'))
        n = int(cov[0]); cov = cov[1:].reshape(n, n)
        m = zcmb > 0.01
        self.mag, self.zcmb, self.zhel = mag[m], zcmb[m], zhel[m]
        icov = np.linalg.inv(cov[np.ix_(m, m)])
        one = np.ones_like(self.mag)[:, None]
        u = icov @ one
        self.icov = icov - u @ u.T / float(one.T @ u)          # analytic marginalisation over M (Cobaya)
        self.n = int(m.sum())
        self.zgrid = np.linspace(0.0, self.zcmb.max() * 1.001, 3000)

    def chi2(self, H):
        Dc = cumulative_trapezoid(C / H(self.zgrid), self.zgrid, initial=0.0)    # D_M (flat), Mpc
        DM = np.interp(self.zcmb, self.zgrid, Dc)
        diff = self.mag - 5.0 * np.log10((1.0 + self.zhel) * DM)
        return float(diff @ self.icov @ diff)


class BgLikelihood:
    def __init__(self):
        self.desi = obs.DesiBAO('DR2')
        self.sn = Pantheon()

    def compressed(self, H, omega_b, omega_cdm):
        om = omega_b + omega_cdm + OMEGA_NU
        R, lA = _raw_compressed(H, omega_b, om)
        x = np.array([R * K_R, lA * K_LA, omega_b])
        d = x - DP_MEAN
        return float(d @ DP_ICOV @ d), x

    def chi2(self, H, omega_b, omega_cdm, h):
        c_cmb, x = self.compressed(H, omega_b, omega_cdm)
        c_bao = float(self.desi.chi2(H, r_drag_aubourg(omega_b, omega_cdm)))
        c_sn = self.sn.chi2(H)
        return dict(cmb=c_cmb, bao=c_bao, sn=c_sn, total=c_cmb + c_bao + c_sn, R=float(x[0]), l_A=float(x[1]))
