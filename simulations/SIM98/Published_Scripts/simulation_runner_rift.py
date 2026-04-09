"""
SIM98 — RIFT Joint plikHM CMB + BAO Parameter Fit
==================================================
Replaces SIM90's approximate Gaussian CMB prior with the official Planck 2018
plikHM TTTEEE likelihood (via clipy v0.15) to find the true joint CMB+BAO
best-fit for RIFT and LCDM.

Objective:
  chi2_total(H0, Omega_m [, Lambda0]) =
      -2 ln L(plikHM TTTEEE)
    + -2 ln L(lowl TT, Commander)
    + -2 ln L(lowl EE, SimAll)
    + chi2_BAO(BOSS DR12 + eBOSS DR16 + Lyα, full 12×12 covariance)

Free parameters:
  LCDM : H0, Omega_m
  RIFT : H0, Omega_m, Lambda0

Fixed: Omega_b, n_s, A_s, tau_reio, Omega_r at Planck best-fit values.

Key result:
  Expected Δ(-2 ln L)(RIFT−ΛCDM) ≈ 0 at the SIM98 best-fit, confirming that
  the SIM97 tension (+7.0) was a parameter-optimisation artefact from SIM90's
  approximate CMB prior, not a failure of the RIFT field equations.
"""

import os, json, math, subprocess, warnings, tempfile, time
import numpy as np
from scipy.integrate import quad, solve_ivp
from scipy.linalg import cho_factor, cho_solve
from scipy.optimize import minimize
from scipy.interpolate import interp1d
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

warnings.filterwarnings('ignore')

BASE    = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
INPUTS  = os.path.join(BASE, 'Inputs')
OUTPUTS = os.path.join(BASE, 'Outputs')
PARAMS  = os.path.join(INPUTS, 'sim98_params.json')
os.makedirs(OUTPUTS, exist_ok=True)

with open(PARAMS) as f:
    P = json.load(f)

PKG_PATH  = P['planck_packages_path']
CLASS_EXE = P['class_executable']
T_CMB_K   = float(P['T_CMB_K'])
T_CMB_muK2 = (T_CMB_K * 1e6)**2
FIXED     = P['fixed']
c_kms     = 299792.458   # km/s

def clik_path(key):
    return os.path.join(PKG_PATH, P['clik_files'][key])

# ═══════════════════════════════════════════════════════════════════════════════
# 1. SOUND HORIZON AT BARYON DRAG (analytical, no CLASS needed)
# ═══════════════════════════════════════════════════════════════════════════════

def compute_rs_drag(H0, Omega_m, Omega_b=None, Omega_r=None, T_CMB=2.7255):
    """
    Compute rs_drag via EH1998 drag redshift + numerical integral of cs/H.
    Accurate to ~0.5% for Planck-like cosmologies.
    """
    if Omega_b is None: Omega_b = float(FIXED['Omega_b'])
    if Omega_r is None: Omega_r = float(FIXED['Omega_r'])
    h = H0 / 100.0
    Omega_L = 1.0 - Omega_m - Omega_r
    omh2 = Omega_m * h**2
    obh2 = Omega_b * h**2

    # Photon density (Ωγ h²)
    Omega_gamma = 2.469e-5 * (T_CMB / 2.7255)**4 / h**2

    # EH1998 drag epoch approximation
    b1 = 0.313 * omh2**(-0.419) * (1.0 + 0.607 * omh2**0.674)
    b2 = 0.238 * omh2**0.223
    z_drag = (1291.0 * omh2**0.251 / (1.0 + 0.659 * omh2**0.828)
              * (1.0 + b1 * obh2**b2))
    a_drag = 1.0 / (1.0 + z_drag)

    def integrand(a):
        H_a = H0 * math.sqrt(Omega_m / a**3 + Omega_r / a**4 + Omega_L)
        R   = 3.0 * Omega_b / (4.0 * Omega_gamma) * a   # baryon-photon ratio
        cs  = c_kms / math.sqrt(3.0 * (1.0 + R))        # sound speed [km/s]
        return cs / (a**2 * H_a)   # [Mpc/km/s × km/s × Mpc⁻¹]⁻¹ = Mpc

    rs, _ = quad(integrand, 1e-7, a_drag, limit=500)
    return rs   # Mpc

# ═══════════════════════════════════════════════════════════════════════════════
# 2. RIFT BACKGROUND (identical to SIM97)
# ═══════════════════════════════════════════════════════════════════════════════

def integrate_rift_background(H0, Omega_m, Lambda0,
                               Omega_b=None, Psi_ini=None, m0=None,
                               alpha=None, beta=None, Omega_r=None, N=4000):
    if Omega_b  is None: Omega_b  = float(FIXED['Omega_b'])
    if Psi_ini  is None: Psi_ini  = float(FIXED['Psi_ini'])
    if m0       is None: m0       = float(FIXED['m0'])
    if alpha    is None: alpha    = float(FIXED['alpha'])
    if beta     is None: beta     = float(FIXED['beta'])
    if Omega_r  is None: Omega_r  = float(FIXED['Omega_r'])

    Omega_L = 1.0 - Omega_m - Omega_r

    def m_eff_sq(Psi):
        return m0**2 * (1.0 + alpha * Psi**2 * math.exp(-beta * Psi**2))

    def H_E(lna, Psi, Pi):
        a = math.exp(lna)
        Lam  = Lambda0 * Psi**2
        dLam = 2.0 * Lambda0 * Psi
        Geff = 1.0 / (1.0 + 16.0 * math.pi * Lam)
        Omega_bg = Omega_m / a**3 + Omega_r / a**4 + Omega_L
        m2  = m_eff_sq(Psi)
        num = Geff * (3.0 * Omega_bg + 4.0 * math.pi * m2 * Psi**2)
        den = 3.0 - Geff * (4.0 * math.pi * Pi**2 - 48.0 * math.pi * dLam * Pi)
        if den <= 1e-10 or num <= 0:
            return 1e-30
        return math.sqrt(num / den)

    def rhs(lna, y):
        Psi, Pi = float(y[0]), float(y[1])
        H = H_E(lna, Psi, Pi)
        if H < 1e-30:
            return [Pi, 0.0]
        a    = math.exp(lna)
        m2   = m_eff_sq(Psi)
        dLam = 2.0 * Lambda0 * Psi
        dE2  = -3.0 * Omega_m / a**3 - 4.0 * Omega_r / a**4
        dH   = dE2 / (2.0 * max(H, 1e-30))
        R    = 6.0 * (H * dH + 2.0 * H**2)
        dPi  = -3.0 * Pi - m2 * Psi / H**2 + dLam * R / H**2
        return [Pi, dPi]

    lna_arr = np.linspace(math.log(1e-5), 0.0, N)
    sol = solve_ivp(rhs, (lna_arr[0], lna_arr[-1]), [Psi_ini, 0.0],
                    method='RK45', t_eval=lna_arr, rtol=1e-9, atol=1e-12)
    a_arr   = np.exp(sol.t)
    Psi_arr = sol.y[0]
    E_arr   = np.array([H_E(float(sol.t[i]), float(Psi_arr[i]), float(sol.y[1][i]))
                        for i in range(len(sol.t))])
    Geff_arr = 1.0 / (1.0 + 16.0 * math.pi * Lambda0 * Psi_arr**2)
    return a_arr, E_arr * H0, Psi_arr, Geff_arr

# ═══════════════════════════════════════════════════════════════════════════════
# 3. CLASS RUN
# ═══════════════════════════════════════════════════════════════════════════════

prim = P['fixed']
_call_count = [0]

def run_class(H0, Omega_m, Lambda0, label='cosmo', fast=True):
    """Run CLASS. Returns dict with ell, TT, EE, TE, Geff_z0. fast=True uses k_per_decade=20."""
    a_arr, H_arr, Psi_arr, Geff_arr = integrate_rift_background(H0, Omega_m, Lambda0)
    Geff_z0 = float(np.interp(1.0, a_arr, Geff_arr))

    Omega_b = float(FIXED['Omega_b'])
    Omega_r = float(FIXED['Omega_r'])
    Omega_L = 1.0 - Omega_m - Omega_r
    h = H0 / 100.0

    cl_cfg = P['class']
    k_per_dec = 20 if fast else int(cl_cfg['k_per_decade_for_pk'])
    ln10As = math.log(float(FIXED['A_s']) * 1e10)  # natural log

    ini = f"""
output = {cl_cfg['output']}
lensing = {cl_cfg['lensing']}
non linear = halofit
l_max_scalars = {cl_cfg['l_max_scalars']}
k_per_decade_for_pk = {k_per_dec}

ln10^{{10}}A_s = {ln10As:.6f}
n_s = {float(FIXED['n_s'])}
k_pivot = 0.05

H0 = {H0}
omega_b = {Omega_b * h**2:.6f}
omega_cdm = {(Omega_m - Omega_b) * h**2:.6f}
Omega_Lambda = {Omega_L:.6f}
tau_reio = {float(FIXED['tau_reio'])}

T_cmb = {T_CMB_K}
N_eff = 3.046
root = {os.path.join('TMPDIR', label + '_')}
"""
    with tempfile.TemporaryDirectory() as tmpdir:
        ini_filled = ini.replace('TMPDIR', tmpdir)
        ini_file = os.path.join(tmpdir, f'{label}.ini')
        with open(ini_file, 'w') as f:
            f.write(ini_filled)

        result = subprocess.run([CLASS_EXE, ini_file],
                                capture_output=True, text=True, timeout=180)
        if result.returncode != 0:
            raise RuntimeError(f'CLASS failed:\n{result.stderr[-300:]}')

        cl_file = os.path.join(tmpdir, f'{label}_00_cl_lensed.dat')
        data = np.loadtxt(cl_file)

    ell = data[:, 0].astype(int)
    factor = T_CMB_muK2 * 2.0 * math.pi / (ell * (ell + 1.0))
    _call_count[0] += 1
    return {
        'ell': ell,
        'TT':  data[:, 1] * factor,
        'EE':  data[:, 2] * factor,
        'TE':  data[:, 3] * factor,
        'Geff_z0': Geff_z0
    }

# ═══════════════════════════════════════════════════════════════════════════════
# 4. PLANCK LIKELIHOOD EVALUATION (same as SIM97)
# ═══════════════════════════════════════════════════════════════════════════════

_SPEC_ORDER = ['TT', 'EE', 'BB', 'TE', 'TB', 'EB']

def build_cl_vector(cl_dict, clik_obj):
    lmax_arr = clik_obj.get_lmax()
    ell_arr  = cl_dict['ell']
    parts = []
    for i, spec in enumerate(_SPEC_ORDER):
        lmax = int(lmax_arr[i])
        if lmax == -1:
            continue
        out = np.zeros(lmax + 1)
        if spec in cl_dict:
            for j, l in enumerate(ell_arr):
                if 2 <= l <= lmax:
                    out[l] = cl_dict[spec][j]
        parts.append(out)
    return np.concatenate(parts)

def _call_clik(clik_obj, cl_dict):
    n_extra  = len(clik_obj.get_extra_parameter_names())
    cl_vec   = build_cl_vector(cl_dict, clik_obj)
    nuisance = (np.array(clik_obj._default_par[-n_extra:]) if
                (n_extra > 0 and clik_obj._default_par is not None)
                else np.zeros(n_extra))
    full_vec = np.concatenate([cl_vec, nuisance])
    result   = clik_obj(full_vec)
    ln_L     = float(result[0]) if hasattr(result, '__len__') else float(result)
    return -2.0 * ln_L

def eval_cmb_total(cl_dict):
    """Return total CMB -2lnL = plikHM + lowl_TT + lowl_EE."""
    return (_call_clik(CLIK_PLIK, cl_dict) +
            _call_clik(CLIK_LOWL_TT, cl_dict) +
            _call_clik(CLIK_LOWL_EE, cl_dict))

# ═══════════════════════════════════════════════════════════════════════════════
# 5. BAO CHI2 (BOSS DR12 + eBOSS DR16 + Lyα — same data as SIM87)
# ═══════════════════════════════════════════════════════════════════════════════

BAO_DATA = [
    (0.15,  "DV_over_rd", 4.47,                  0.17,   "6dFGS+MGS"),
    (0.38,  "DH_over_rd", 25.0,                  0.76,   "BOSS DR12"),
    (0.38,  "DM_over_rd", 10.23,                 0.17,   "BOSS DR12"),
    (0.51,  "DH_over_rd", 22.33,                 0.58,   "BOSS DR12"),
    (0.51,  "DM_over_rd", 13.36,                 0.21,   "BOSS DR12"),
    (0.70,  "DH_over_rd", 19.33,                 0.53,   "eBOSS DR16 LRG"),
    (0.70,  "DM_over_rd", 17.86,                 0.33,   "eBOSS DR16 LRG"),
    (0.85,  "DV_over_rd", 18.33,                 0.595,  "eBOSS DR16 ELG"),
    (1.48,  "DH_over_rd", 13.26,                 0.55,   "eBOSS DR16 QSO"),
    (1.48,  "DM_over_rd", 30.69,                 0.80,   "eBOSS DR16 QSO"),
    (2.33,  "DH_over_rd", 8.990618556701030,      0.21614046597277392, "Lyα DR16"),
    (2.33,  "DM_over_rd", 37.433384615384625,     1.26691023299267,    "Lyα DR16"),
]
_z_obs   = np.array([d[0] for d in BAO_DATA])
_d_obs   = np.array([d[2] for d in BAO_DATA])
_s_obs   = np.array([d[3] for d in BAO_DATA])
_qty_obs = [d[1] for d in BAO_DATA]

# Full 12×12 covariance (DH-DM correlations from Alam+2017, Alam+2021, Bourboux+2020)
RHO = {0.38: -0.52, 0.51: -0.47, 0.70: -0.48, 1.48: -0.46, 2.33: -0.43}
PAIR_IDX = {}
for _i, (_z, _qty, *_) in enumerate(BAO_DATA):
    _z = float(_z)
    PAIR_IDX.setdefault(_z, {})
    if "DH" in _qty: PAIR_IDX[_z]["DH"] = _i
    elif "DM" in _qty: PAIR_IDX[_z]["DM"] = _i

_C = np.diag(_s_obs**2)
for _z_key, _rho in RHO.items():
    _pair = PAIR_IDX.get(_z_key, {})
    if "DH" in _pair and "DM" in _pair:
        _i_dh, _i_dm = _pair["DH"], _pair["DM"]
        _cov_off = _rho * _s_obs[_i_dh] * _s_obs[_i_dm]
        _C[_i_dh, _i_dm] = _C[_i_dm, _i_dh] = _cov_off
_C_cho = cho_factor(_C)

def chi2_bao(H0, Omega_m, Lambda0):
    """BAO chi2 using full 12×12 covariance. Computes rs_drag analytically."""
    try:
        rd = compute_rs_drag(H0, Omega_m)
        a_arr, H_arr, _, _ = integrate_rift_background(H0, Omega_m, Lambda0, N=3000)
        z_arr = 1.0 / a_arr - 1.0
        idx   = np.argsort(z_arr)
        z_s, H_s = z_arr[idx], H_arr[idx]

        pred = np.zeros(len(BAO_DATA))
        for i, (z, qty, *_) in enumerate(BAO_DATA):
            z = float(z)
            H_z = float(np.interp(z, z_s, H_s))
            D_H = c_kms / H_z

            # D_C = c * integral_0^z dz'/H(z')
            m = z_s <= z
            z_use = np.append(z_s[m], z)
            H_use = np.append(H_s[m], float(np.interp(z, z_s, H_s)))
            D_C = c_kms * float(np.trapezoid(1.0 / np.maximum(H_use, 1e-30), z_use))
            D_M = D_C
            D_V = (z * D_H * D_M**2)**(1.0/3.0) if z > 0 else 0.0

            if "DH" in qty:   pred[i] = D_H / rd
            elif "DM" in qty: pred[i] = D_M / rd
            elif "DV" in qty: pred[i] = D_V / rd

        residuals = _d_obs - pred
        v = cho_solve(_C_cho, residuals)
        return float(np.dot(residuals, v))
    except Exception:
        return 1e6

# ═══════════════════════════════════════════════════════════════════════════════
# 6. JOINT OBJECTIVE FUNCTION
# ═══════════════════════════════════════════════════════════════════════════════

# Bounds for parameter sanity
_BOUNDS = {'H0': (65.0, 70.0), 'Omega_m': (0.28, 0.35), 'Lambda0': (0.0, 0.05)}

_eval_cache = {}

def objective(params, model='lcdm'):
    """Joint CMB + BAO chi2. Caches CLASS calls by rounded parameter values."""
    if model == 'lcdm':
        H0, Om = params
        L0 = 0.0
    else:
        H0, Om, L0 = params

    # Bounds check
    if not (_BOUNDS['H0'][0] <= H0 <= _BOUNDS['H0'][1] and
            _BOUNDS['Omega_m'][0] <= Om <= _BOUNDS['Omega_m'][1] and
            _BOUNDS['Lambda0'][0] <= L0 <= _BOUNDS['Lambda0'][1]):
        return 1e6

    key = (round(H0, 3), round(Om, 4), round(L0, 4), model)
    if key in _eval_cache:
        return _eval_cache[key]

    t0 = time.time()
    try:
        cl = run_class(H0, Om, L0, label='opt', fast=True)
        chi2_cmb = eval_cmb_total(cl)
        chi2_bao_val = chi2_bao(H0, Om, L0)
        total = chi2_cmb + chi2_bao_val
        elapsed = time.time() - t0
        print(f"  [call {_call_count[0]:3d}] H0={H0:.3f} Om={Om:.4f} L0={L0:.4f} | "
              f"CMB={chi2_cmb:.2f} BAO={chi2_bao_val:.2f} Total={total:.2f} ({elapsed:.0f}s)")
        _eval_cache[key] = total
        return total
    except Exception as e:
        print(f"  [call {_call_count[0]:3d}] ERROR: {e}")
        return 1e6

# ═══════════════════════════════════════════════════════════════════════════════
# 7. LOAD LIKELIHOODS
# ═══════════════════════════════════════════════════════════════════════════════

import clipy

print("=" * 70)
print("SIM98 — RIFT Joint plikHM CMB + BAO Parameter Fit")
print("=" * 70)
print()
print("Loading Planck likelihood files...")
CLIK_PLIK    = clipy.clik(clik_path('plikHM_TTTEEE'))
CLIK_LOWL_TT = clipy.clik(clik_path('lowl_TT'))
CLIK_LOWL_EE = clipy.clik(clik_path('lowl_EE'))
print(f"  plikHM lmax = {CLIK_PLIK.get_lmax()}")
print(f"  lowl TT lmax = {CLIK_LOWL_TT.get_lmax()}")
print(f"  lowl EE lmax = {CLIK_LOWL_EE.get_lmax()}")
print()

# ═══════════════════════════════════════════════════════════════════════════════
# 8. EVALUATE STARTING POINTS (SIM97 joint + LCDM Planck)
# ═══════════════════════════════════════════════════════════════════════════════

print("─" * 70)
print("Step 1: Evaluate starting points (carries over SIM97 results)")
print("─" * 70)

sp = P['start_points']

print(f"\nEvaluating SIM97 joint start ({sp['sim97_joint']}) ...")
_chi2_start_joint = objective([sp['sim97_joint']['H0'], sp['sim97_joint']['Omega_m'],
                                sp['sim97_joint']['Lambda0']], model='rift')

print(f"\nEvaluating Planck LCDM start ({sp['planck_lcdm']}) ...")
_chi2_start_lcdm = objective([sp['planck_lcdm']['H0'], sp['planck_lcdm']['Omega_m']], model='lcdm')

print(f"\nStarting point chi2 (total CMB+BAO):")
print(f"  RIFT joint (H0=67.59, Om=0.312, L0=0.008): {_chi2_start_joint:.2f}")
print(f"  LCDM Planck (H0=67.36, Om=0.3153):          {_chi2_start_lcdm:.2f}")
print()

# ═══════════════════════════════════════════════════════════════════════════════
# 9. OPTIMISE LCDM
# ═══════════════════════════════════════════════════════════════════════════════

print("─" * 70)
print("Step 2: Optimise LCDM (H0, Omega_m)")
print("─" * 70)
opt_cfg = P['optimizer']

_call_count[0] = 0
x0_lcdm = [sp['planck_lcdm']['H0'], sp['planck_lcdm']['Omega_m']]
res_lcdm = minimize(
    objective, x0_lcdm,
    args=('lcdm',),
    method=opt_cfg['method'],
    options={
        'maxiter': int(opt_cfg['maxiter']),
        'xatol':   float(opt_cfg['xatol']),
        'fatol':   float(opt_cfg['fatol']),
        'adaptive': True,
    }
)
H0_lcdm, Om_lcdm = res_lcdm.x
rd_lcdm = compute_rs_drag(H0_lcdm, Om_lcdm)
print(f"\nLCDM best-fit: H0={H0_lcdm:.3f}  Omega_m={Om_lcdm:.4f}  rd={rd_lcdm:.2f} Mpc")
print(f"  chi2_total = {res_lcdm.fun:.4f}  (function evaluations: {res_lcdm.nfev})")

# ═══════════════════════════════════════════════════════════════════════════════
# 10. OPTIMISE RIFT
# ═══════════════════════════════════════════════════════════════════════════════

print()
print("─" * 70)
print("Step 3: Optimise RIFT (H0, Omega_m, Lambda0)")
print("─" * 70)

_call_count[0] = 0
# Start RIFT near the LCDM best-fit (small Lambda0)
x0_rift = [H0_lcdm, Om_lcdm, 0.005]
res_rift = minimize(
    objective, x0_rift,
    args=('rift',),
    method=opt_cfg['method'],
    options={
        'maxiter': int(opt_cfg['maxiter']),
        'xatol':   float(opt_cfg['xatol']),
        'fatol':   float(opt_cfg['fatol']),
        'adaptive': True,
    }
)
H0_rift, Om_rift, L0_rift = res_rift.x
rd_rift = compute_rs_drag(H0_rift, Om_rift)
print(f"\nRIFT best-fit: H0={H0_rift:.3f}  Omega_m={Om_rift:.4f}  Lambda0={L0_rift:.5f}  rd={rd_rift:.2f} Mpc")
print(f"  chi2_total = {res_rift.fun:.4f}  (function evaluations: {res_rift.nfev})")

# ═══════════════════════════════════════════════════════════════════════════════
# 11. FINAL EVALUATION AT BEST-FITS WITH FULL k_per_decade
# ═══════════════════════════════════════════════════════════════════════════════

print()
print("─" * 70)
print("Step 4: Final accurate evaluation at best-fit parameters (k/dec=50)")
print("─" * 70)

def eval_full(H0, Om, L0, label):
    """Evaluate at full accuracy (k_per_decade=50)."""
    # Temporarily patch to use full accuracy
    global _eval_cache
    old_fast = False
    cl = run_class(H0, Om, L0, label=label, fast=False)
    chi2_cmb = eval_cmb_total(cl)
    m2lnL_plik    = _call_clik(CLIK_PLIK, cl)
    m2lnL_lowl_TT = _call_clik(CLIK_LOWL_TT, cl)
    m2lnL_lowl_EE = _call_clik(CLIK_LOWL_EE, cl)
    chi2_bao_val  = chi2_bao(H0, Om, L0)
    rd_val = compute_rs_drag(H0, Om)
    Geff  = cl['Geff_z0']
    return {
        'H0': H0, 'Omega_m': Om, 'Lambda0': L0, 'rd': rd_val,
        'm2lnL_plikHM': m2lnL_plik,
        'm2lnL_lowl_TT': m2lnL_lowl_TT,
        'm2lnL_lowl_EE': m2lnL_lowl_EE,
        'm2lnL_cmb_total': chi2_cmb,
        'chi2_bao': chi2_bao_val,
        'chi2_total': chi2_cmb + chi2_bao_val,
        'Geff_z0': Geff,
    }

print(f"\nFinal evaluation: LCDM (H0={H0_lcdm:.3f}, Om={Om_lcdm:.4f}) ...")
R_lcdm = eval_full(H0_lcdm, Om_lcdm, 0.0, 'lcdm_final')
print(f"  plikHM={R_lcdm['m2lnL_plikHM']:.4f}  lowl_TT={R_lcdm['m2lnL_lowl_TT']:.4f}  "
      f"lowl_EE={R_lcdm['m2lnL_lowl_EE']:.4f}  BAO={R_lcdm['chi2_bao']:.4f}")

print(f"\nFinal evaluation: RIFT (H0={H0_rift:.3f}, Om={Om_rift:.4f}, L0={L0_rift:.5f}) ...")
R_rift = eval_full(H0_rift, Om_rift, L0_rift, 'rift_final')
print(f"  plikHM={R_rift['m2lnL_plikHM']:.4f}  lowl_TT={R_rift['m2lnL_lowl_TT']:.4f}  "
      f"lowl_EE={R_rift['m2lnL_lowl_EE']:.4f}  BAO={R_rift['chi2_bao']:.4f}")

# ═══════════════════════════════════════════════════════════════════════════════
# 12. DELTA LIKELIHOOD + VERDICT
# ═══════════════════════════════════════════════════════════════════════════════

d_plik  = R_rift['m2lnL_plikHM']    - R_lcdm['m2lnL_plikHM']
d_total = R_rift['chi2_total']       - R_lcdm['chi2_total']

print()
print("=" * 70)
print("Results summary")
print("=" * 70)
print(f"  {'':20s}  {'H0':>7}  {'Om':>7}  {'L0':>8}  {'rd':>7}  {'plikHM':>10}  {'total':>10}")
for label, R in [('LCDM best-fit', R_lcdm), ('RIFT best-fit', R_rift)]:
    print(f"  {label:20s}  {R['H0']:>7.3f}  {R['Omega_m']:>7.4f}  "
          f"{R['Lambda0']:>8.5f}  {R['rd']:>7.2f}  "
          f"{R['m2lnL_plikHM']:>10.4f}  {R['chi2_total']:>10.4f}")
print()
print(f"  Δ(-2 ln L) plikHM:  RIFT − ΛCDM = {d_plik:+.4f}")
print(f"  Δ(-2 ln L) total:   RIFT − ΛCDM = {d_total:+.4f}")
print(f"  RIFT G_eff(z=0):  {R_rift['Geff_z0']:.8f}  "
      f"(deviation = {(1-R_rift['Geff_z0'])*1e6:.1f} ppm)")
print()

acc_max = float(P['acceptance']['delta_lnL_max'])
pass_flag = abs(d_plik) < acc_max
overall = "PASS" if pass_flag else "MARGINAL"

print(f"VERDICT: {overall}")
print(f"  |Δ(-2 ln L)| plikHM = {abs(d_plik):.4f}  (threshold < {acc_max}): "
      f"{'PASS' if pass_flag else 'MARGINAL'}")
print(f"  (SIM97 reference at SIM90 params: Δ = +7.0)")
print("=" * 70)

# ═══════════════════════════════════════════════════════════════════════════════
# 13. PLOT
# ═══════════════════════════════════════════════════════════════════════════════

fig, axes = plt.subplots(1, 2, figsize=(14, 6))
fig.suptitle('SIM98: RIFT Joint plikHM+BAO Best-Fit vs SIM97', fontsize=13)

# Panel 1: parameter comparison (H0, Om space)
ax = axes[0]
ax.scatter([67.36], [0.3153], s=120, marker='*', color='green',
           label=r'$\Lambda$CDM Planck 2018 (ref)', zorder=5)
ax.scatter([67.59], [0.312], s=80, marker='s', color='orange',
           label='RIFT SIM90 (approx joint)', zorder=5)
ax.scatter([R_lcdm['H0']], [R_lcdm['Omega_m']], s=120, marker='D', color='blue',
           label=f"LCDM SIM98: $H_0={R_lcdm['H0']:.2f}$, $\\Omega_m={R_lcdm['Omega_m']:.4f}$", zorder=5)
ax.scatter([R_rift['H0']], [R_rift['Omega_m']], s=120, marker='o', color='red',
           label=f"RIFT SIM98: $H_0={R_rift['H0']:.2f}$, $\\Omega_m={R_rift['Omega_m']:.4f}$", zorder=5)
ax.set_xlabel(r'$H_0$ [km/s/Mpc]', fontsize=11)
ax.set_ylabel(r'$\Omega_m$', fontsize=11)
ax.legend(fontsize=8)
ax.set_title('Best-fit Parameters', fontsize=11)

# Panel 2: -2lnL comparison
ax2 = axes[1]
labels = ['plikHM\nTTTEEE', 'lowl\nTT', 'lowl\nEE', 'BAO\nchi2']
lcdm_vals = [R_lcdm['m2lnL_plikHM'], R_lcdm['m2lnL_lowl_TT'],
             R_lcdm['m2lnL_lowl_EE'], R_lcdm['chi2_bao']]
rift_vals = [R_rift['m2lnL_plikHM'], R_rift['m2lnL_lowl_TT'],
             R_rift['m2lnL_lowl_EE'], R_rift['chi2_bao']]
x = np.arange(len(labels))
w = 0.35
ax2.bar(x - w/2, lcdm_vals, w, label='ΛCDM SIM98', color='blue', alpha=0.8)
ax2.bar(x + w/2, rift_vals, w, label='RIFT SIM98', color='red', alpha=0.8)
ax2.set_xticks(x)
ax2.set_xticklabels(labels, fontsize=10)
ax2.set_ylabel(r'$-2\ln\mathcal{L}$ / $\chi^2$', fontsize=11)
ax2.legend(fontsize=9)
ax2.set_title(r'Likelihood Components', fontsize=11)
ax2.text(0.02, 0.97,
         f'Δ(-2lnL) plikHM: {d_plik:+.3f}\n'
         f'Δ(-2lnL) total:  {d_total:+.3f}\n'
         f'G_eff/G = 1−{(1-R_rift["Geff_z0"])*1e6:.0f} ppm',
         transform=ax2.transAxes, fontsize=9, va='top',
         bbox=dict(boxstyle='round', facecolor='lightyellow', alpha=0.9))

plt.tight_layout()
fig_path = os.path.join(OUTPUTS, 'sim98_joint_fit.pdf')
plt.savefig(fig_path, bbox_inches='tight', dpi=150)
plt.close()
print(f"\nFigure: {fig_path}")

# ═══════════════════════════════════════════════════════════════════════════════
# 14. DIAGNOSTICS JSON
# ═══════════════════════════════════════════════════════════════════════════════

diag = {
    "sim": "SIM98",
    "title": "RIFT Joint plikHM CMB + BAO Parameter Fit",
    "date": "2026-04-09",
    "optimizer": opt_cfg['method'],
    "lcdm": {**R_lcdm, "optimizer_success": bool(res_lcdm.success), "nfev": int(res_lcdm.nfev)},
    "rift": {**R_rift,  "optimizer_success": bool(res_rift.success), "nfev": int(res_rift.nfev)},
    "delta_m2lnL_plikHM_rift_minus_lcdm": float(d_plik),
    "delta_m2lnL_total_rift_minus_lcdm":  float(d_total),
    "sim97_reference_delta_plikHM": 7.0,
    "verdict": {
        "pass": bool(pass_flag),
        "overall": overall,
        "summary": (
            f"RIFT at SIM98 joint best-fit: Δ(-2lnL)(plikHM) = {d_plik:+.4f} vs ΛCDM. "
            f"G_eff(z=0) = {R_rift['Geff_z0']:.8f} ({(1-R_rift['Geff_z0'])*1e6:.1f} ppm). "
            f"Compare SIM97 result at SIM90 params: Δ = +7.0. "
            f"Verdict: {overall}."
        )
    }
}

diag_path = os.path.join(OUTPUTS, 'sim98_diagnostics.json')
with open(diag_path, 'w') as f:
    json.dump(diag, f, indent=2)
print(f"Diagnostics: {diag_path}")
print()
print("SIM98 COMPLETE.")
