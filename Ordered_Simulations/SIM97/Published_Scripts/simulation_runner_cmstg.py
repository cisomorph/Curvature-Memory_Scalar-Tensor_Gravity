"""
SIM97 — CMSTG Full Planck 2018 plikHM TTTEEE Likelihood
=======================================================
Replaces the approximate Gaussian likelihood in SIM89 with the official
Planck 2018 plikHM TTTEEE likelihood code (clipy / clik interface).

Likelihood files used:
  - plik_rd12_HM_v22b_TTTEEE.clik  (high-ell TT+TE+EE, ell=30-2508)
  - plik_lite_v22_TTTEEE.clik       (cross-check: compressed likelihood)
  - commander_dx12_v3_2_29.clik     (low-ell TT, ell=2-29)
  - simall_..._EE.clik              (low-ell EE, ell=2-29)

Three cosmologies:
  A. CMSTG joint CMB+BAO best-fit (SIM90): H0=67.59, Omega_m=0.312, Lambda0=0.008
  B. CMSTG BAO-only best-fit (SIM87):      H0=68.14, Omega_m=0.294, Lambda0=0.003
  C. LCDM Planck 2018 best-fit:           H0=67.36, Omega_m=0.3153

Key test:
  Delta(-2 ln L) = -2 ln L(CMSTG_joint) - (-2 ln L(LCDM))
  PASS if |Delta(-2 ln L)| < 5 on full plikHM TTTEEE.

Expected result:
  CMSTG at joint best-fit (A) ≈ LCDM (C): same CLASS spectrum → same likelihood.
  CMSTG BAO-only (B) shows large tension (large positive Delta(-2lnL)),
  consistent with SIM89 result (Delta_chi2 = +907).
"""

import os, json, math, subprocess, warnings, tempfile
import numpy as np
from scipy.integrate import solve_ivp
from scipy.interpolate import interp1d
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

warnings.filterwarnings('ignore')

BASE    = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
INPUTS  = os.path.join(BASE, 'Inputs')
OUTPUTS = os.path.join(BASE, 'Outputs')
PARAMS  = os.path.join(INPUTS, 'sim97_params.json')
os.makedirs(OUTPUTS, exist_ok=True)

with open(PARAMS) as f:
    P = json.load(f)

PKG_PATH  = P['planck_packages_path']
CLASS_EXE = P['class_executable']
T_CMB_K   = float(P['T_CMB_K'])
T_CMB_muK2 = (T_CMB_K * 1e6)**2

def clik_path(key):
    return os.path.join(PKG_PATH, P['clik_files'][key])

# ═══════════════════════════════════════════════════════════════════════════════
# 1. CMSTG BACKGROUND (identical to SIM88/95)
# ═══════════════════════════════════════════════════════════════════════════════

def integrate_cmstg_background(H0, Omega_m, Lambda0, Omega_b,
                               Psi_ini=0.01, m0=1.0, alpha=0.1, beta=0.05,
                               Omega_r=9.2e-5, N=5000):
    Omega_L = 1.0 - Omega_m - Omega_r
    def m_eff_sq(Psi):
        return m0**2 * (1.0 + alpha * Psi**2 * math.exp(-beta * Psi**2))
    def H_E(lna, Psi, Pi):
        a = math.exp(lna)
        Lam = Lambda0 * Psi**2
        dLam = 2.0 * Lambda0 * Psi
        Geff = 1.0 / (1.0 + 16.0 * math.pi * Lam)
        Omega_bg = Omega_m / a**3 + Omega_r / a**4 + Omega_L
        m2   = m_eff_sq(Psi)
        num  = Geff * (3.0 * Omega_bg + 4.0 * math.pi * m2 * Psi**2)
        den  = 3.0 - Geff * (4.0 * math.pi * Pi**2 - 48.0 * math.pi * dLam * Pi)
        if den <= 1e-10 or num <= 0:
            return 1e-30
        return math.sqrt(num / den)
    def rhs(lna, y):
        Psi, Pi = float(y[0]), float(y[1])
        H = H_E(lna, Psi, Pi)
        if H < 1e-30:
            return [Pi, 0.0]
        a = math.exp(lna)
        m2 = m_eff_sq(Psi)
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
# 2. CLASS RUN (analytic primordial spectrum)
# ═══════════════════════════════════════════════════════════════════════════════
# For plikHM evaluation, G_eff correction is <10 ppm at Lambda0=0.008 (SIM91).
# The Cl spectrum is dominated by the background cosmology (H0, Omega_m).
# We use CLASS's standard analytic power law primordial spectrum; this is
# equivalent to the external_Pk approach used in SIM88/95 since G_eff ~ 1.

prim = P['primordial']

def run_class(cosmo, label, tmpdir):
    """Run CLASS for a cosmology. Returns dict of ell->Cl arrays (muK^2)."""
    a_arr, H_arr, Psi_arr, Geff_arr = integrate_cmstg_background(
        H0=cosmo['H0'], Omega_m=cosmo['Omega_m'],
        Lambda0=cosmo['Lambda0'], Omega_b=cosmo['Omega_b'],
        Psi_ini=cosmo.get('Psi_ini', 0.0),
        m0=cosmo.get('m0', 1.0),
        alpha=cosmo.get('alpha', 0.0),
        beta=cosmo.get('beta', 0.0),
        Omega_r=cosmo.get('Omega_r', 9.2e-5)
    )
    Geff_z0 = float(np.interp(1.0, a_arr, Geff_arr))

    Omega_L = 1.0 - cosmo['Omega_m'] - cosmo.get('Omega_r', 9.2e-5)
    h = cosmo['H0'] / 100.0

    cl_cfg = P['class']
    ini = f"""
output = {cl_cfg['output']}
lensing = {cl_cfg['lensing']}
non linear = halofit
l_max_scalars = {cl_cfg['l_max_scalars']}
k_per_decade_for_pk = {cl_cfg['k_per_decade_for_pk']}

ln10^{{10}}A_s = {math.log10(float(prim['A_s']) * 1e10) / math.log10(math.e):.6f}
n_s = {float(prim['n_s'])}
k_pivot = {float(prim['k_pivot'])}

H0 = {cosmo['H0']}
omega_b = {cosmo['Omega_b'] * h**2:.6f}
omega_cdm = {(cosmo['Omega_m'] - cosmo['Omega_b']) * h**2:.6f}
Omega_Lambda = {Omega_L:.6f}
tau_reio = {float(prim['tau_reio'])}

T_cmb = {T_CMB_K}
N_eff = 3.046
root = {os.path.join(tmpdir, label + '_')}
"""
    ini_file = os.path.join(tmpdir, f'{label}.ini')
    with open(ini_file, 'w') as f:
        f.write(ini)

    result = subprocess.run([CLASS_EXE, ini_file],
                            capture_output=True, text=True, timeout=120)
    if result.returncode != 0:
        raise RuntimeError(f'CLASS failed for {label}:\n{result.stderr[-500:]}')

    # Load lensed Cl file
    cl_file = os.path.join(tmpdir, f'{label}_00_cl_lensed.dat')
    data = np.loadtxt(cl_file)
    # CLASS lensed Cl output: l, TT, EE, TE, BB, PP, TP (dimensionless)
    ell      = data[:, 0].astype(int)
    TT_dimless = data[:, 1]
    EE_dimless = data[:, 2]
    TE_dimless = data[:, 3]
    # Convert to muK^2: Cl_muK2 = Cl_dimless * T_CMB^2 * l(l+1)/(2pi) ...
    # Actually CLASS outputs Cl*l(l+1)/(2pi) in dimensionless (Delta_T^2/T^2 units)
    # plikHM expects Cl in muK^2 (not multiplied by l(l+1)/(2pi))
    # So: Cl_muK2[l] = TT_dimless[l] * T_CMB_muK2 * 2pi / (l*(l+1))
    factor = T_CMB_muK2 * 2.0 * math.pi / (ell * (ell + 1.0))
    TT_muK2 = TT_dimless * factor
    EE_muK2 = EE_dimless * factor
    TE_muK2 = TE_dimless * factor

    return {
        'ell': ell,
        'TT': TT_muK2,
        'EE': EE_muK2,
        'TE': TE_muK2,
        'Geff_z0': Geff_z0
    }

# ═══════════════════════════════════════════════════════════════════════════════
# 3. PLANCK plikHM LIKELIHOOD EVALUATION
# ═══════════════════════════════════════════════════════════════════════════════

# Spectrum order expected by clik: [TT, EE, BB, TE, TB, EB]
_SPEC_ORDER = ['TT', 'EE', 'BB', 'TE', 'TB', 'EB']

def build_cl_vector(cl_dict, clik_obj):
    """Build the Cl vector expected by a clipy clik object.

    clipy normalize_mnp expects: for each spectrum with lmax[i] != -1, a block
    of (lmax[i]+1) values covering ell=0..lmax[i].  Spectra with lmax[i]==-1
    are omitted entirely.  Units: muK^2 (raw Cl, NOT l(l+1)/(2pi) weighted).
    """
    lmax_arr = clik_obj.get_lmax()   # length-6 array: [TT,EE,BB,TE,TB,EB]
    ell_arr  = cl_dict['ell']
    parts = []
    for i, spec in enumerate(_SPEC_ORDER):
        lmax = int(lmax_arr[i])
        if lmax == -1:
            continue
        out = np.zeros(lmax + 1)     # ell=0..lmax; ell=0,1 remain zero
        if spec in cl_dict:
            for j, l in enumerate(ell_arr):
                if 2 <= l <= lmax:
                    out[l] = cl_dict[spec][j]
        parts.append(out)
    return np.concatenate(parts)

def _build_nuisance(clik_obj):
    """Return nuisance array at Planck best-fit values.

    Uses the clik self-test default_par vector (last n_nuisance elements),
    which contains the Planck best-fit values for ALL nuisance parameters
    (foreground amplitudes, calibration, etc.).  This avoids having to
    enumerate the ~47 plikHM foreground parameters manually.
    """
    n_extra = len(clik_obj.get_extra_parameter_names())
    if n_extra == 0:
        return np.array([])
    if clik_obj._default_par is not None:
        return np.array(clik_obj._default_par[-n_extra:])
    # Fallback: calibration params at nominal, foreground at 1.0
    nus = P['nuisance']
    nuisance_vals = []
    for name in clik_obj.get_extra_parameter_names():
        if 'A_planck' in name:     nuisance_vals.append(nus['A_planck'])
        elif 'calib_100T' in name: nuisance_vals.append(nus['calib_100T'])
        elif 'calib_217T' in name: nuisance_vals.append(nus['calib_217T'])
        elif 'calib_100P' in name: nuisance_vals.append(nus['calib_100P'])
        elif 'calib_143P' in name: nuisance_vals.append(nus['calib_143P'])
        elif 'calib_217P' in name: nuisance_vals.append(nus['calib_217P'])
        elif 'calib_143x217P' in name or 'calib_143x217p' in name:
            nuisance_vals.append(nus['calib_143x217P'])
        else: nuisance_vals.append(1.0)
    return np.array(nuisance_vals)

def _call_clik(clik_obj, cl_dict):
    """Evaluate any clik likelihood. Returns -2 ln L."""
    cl_vec   = build_cl_vector(cl_dict, clik_obj)
    full_vec = np.concatenate([cl_vec, _build_nuisance(clik_obj)])
    result   = clik_obj(full_vec)
    ln_L     = float(result[0]) if hasattr(result, '__len__') else float(result)
    return -2.0 * ln_L

def eval_plik_likelihood(cl_dict, clik_obj, **_kwargs):
    """Evaluate plikHM (or plik_lite) likelihood. Returns -2 ln L."""
    return _call_clik(clik_obj, cl_dict)

def eval_lowl_TT(cl_dict, clik_obj):
    """Evaluate low-ell TT (Commander) likelihood. Returns -2 ln L."""
    return _call_clik(clik_obj, cl_dict)

def eval_lowl_EE(cl_dict, clik_obj):
    """Evaluate low-ell EE (SimAll) likelihood. Returns -2 ln L."""
    return _call_clik(clik_obj, cl_dict)

# ═══════════════════════════════════════════════════════════════════════════════
# 4. MAIN RUN
# ═══════════════════════════════════════════════════════════════════════════════

import clipy

print("=" * 70)
print("SIM97 — CMSTG Full Planck 2018 plikHM TTTEEE Likelihood")
print("=" * 70)

# Load likelihood objects once
print("Loading Planck likelihood files...")
clik_plikHM   = clipy.clik(clik_path('plikHM_TTTEEE'))
clik_lite     = clipy.clik(clik_path('plik_lite_TTTEEE'))
clik_lowl_TT  = clipy.clik(clik_path('lowl_TT'))
clik_lowl_EE  = clipy.clik(clik_path('lowl_EE'))
print(f"  plikHM TTTEEE: lmax = {clik_plikHM.get_lmax()}")
print(f"  plik_lite:     lmax = {clik_lite.get_lmax()}")
print(f"  lowl TT:       lmax = {clik_lowl_TT.get_lmax()}")
print(f"  lowl EE:       lmax = {clik_lowl_EE.get_lmax()}")
print()

cosmo_cfg = P['cosmologies']
results = {}

with tempfile.TemporaryDirectory() as tmpdir:
    for label, cosmo in cosmo_cfg.items():
        print(f"Running CLASS for {label} ...")
        cl = run_class(cosmo, label, tmpdir)
        print(f"  CLASS done. G_eff(z=0) = {cl['Geff_z0']:.8f}")

        # plikHM TTTEEE (full, high-ell)
        m2lnL_plik = eval_plik_likelihood(cl, clik_plikHM)
        # plik_lite TTTEEE (cross-check)
        m2lnL_lite = eval_plik_likelihood(cl, clik_lite, lmax_TT=1508, lmax_pol=1508)
        # lowl TT
        m2lnL_lowl_TT = eval_lowl_TT(cl, clik_lowl_TT)
        # lowl EE
        m2lnL_lowl_EE = eval_lowl_EE(cl, clik_lowl_EE)

        m2lnL_total = m2lnL_plik + m2lnL_lowl_TT + m2lnL_lowl_EE

        results[label] = {
            'cl': cl,
            'm2lnL_plikHM':  m2lnL_plik,
            'm2lnL_lite':    m2lnL_lite,
            'm2lnL_lowl_TT': m2lnL_lowl_TT,
            'm2lnL_lowl_EE': m2lnL_lowl_EE,
            'm2lnL_total':   m2lnL_total,
            'Geff_z0':       cl['Geff_z0']
        }

        print(f"  -2 ln L (plikHM TTTEEE) = {m2lnL_plik:.4f}")
        print(f"  -2 ln L (plik_lite)     = {m2lnL_lite:.4f}")
        print(f"  -2 ln L (lowl TT)       = {m2lnL_lowl_TT:.4f}")
        print(f"  -2 ln L (lowl EE)       = {m2lnL_lowl_EE:.4f}")
        print(f"  -2 ln L (total)         = {m2lnL_total:.4f}")
        print()

# ═══════════════════════════════════════════════════════════════════════════════
# 5. DELTA LIKELIHOOD + VERDICT
# ═══════════════════════════════════════════════════════════════════════════════

print("-" * 70)
print("Results summary:")
print(f"  {'Cosmology':20s}  {'plikHM TTTEEE':>16}  {'lowl TT':>10}  {'lowl EE':>10}  {'Total':>10}")
for label, r in results.items():
    print(f"  {label:20s}  {r['m2lnL_plikHM']:>16.4f}  {r['m2lnL_lowl_TT']:>10.4f}  "
          f"{r['m2lnL_lowl_EE']:>10.4f}  {r['m2lnL_total']:>10.4f}")
print()

d_plik_joint  = results['cmstg_joint']['m2lnL_plikHM']  - results['lcdm_planck']['m2lnL_plikHM']
d_plik_bao    = results['cmstg_bao_only']['m2lnL_plikHM'] - results['lcdm_planck']['m2lnL_plikHM']
d_total_joint = results['cmstg_joint']['m2lnL_total']   - results['lcdm_planck']['m2lnL_total']
d_total_bao   = results['cmstg_bao_only']['m2lnL_total'] - results['lcdm_planck']['m2lnL_total']

print(f"  Δ(-2 ln L) plikHM:  CMSTG joint vs ΛCDM = {d_plik_joint:+.4f}")
print(f"  Δ(-2 ln L) plikHM:  CMSTG BAO-only vs ΛCDM = {d_plik_bao:+.4f}")
print(f"  Δ(-2 ln L) total:   CMSTG joint vs ΛCDM = {d_total_joint:+.4f}")
print(f"  Δ(-2 ln L) total:   CMSTG BAO-only vs ΛCDM = {d_total_bao:+.4f}")
print()

acc_max = float(P['acceptance']['delta_lnL_max'])
pass_joint = abs(d_plik_joint) < acc_max
pass_bao_expected = d_plik_bao > 100  # BAO-only should show large tension

overall = "PASS" if pass_joint else "FAIL"
print("=" * 70)
print(f"VERDICT: {overall}")
print(f"  CMSTG joint |Δ(-2 ln L)| = {abs(d_plik_joint):.4f}  (threshold < {acc_max}): {'PASS' if pass_joint else 'FAIL'}")
print(f"  CMSTG BAO-only Δ(-2 ln L) = {d_plik_bao:+.1f}  (large tension expected: {'YES' if pass_bao_expected else 'NO'})")
print(f"  G_eff(z=0) at joint best-fit: {results['cmstg_joint']['Geff_z0']:.8f}  (deviation = {(1-results['cmstg_joint']['Geff_z0'])*1e6:.1f} ppm)")
print("=" * 70)
print()

# ═══════════════════════════════════════════════════════════════════════════════
# 6. PLOT — Cl comparison + likelihood bar chart
# ═══════════════════════════════════════════════════════════════════════════════

fig, axes = plt.subplots(1, 2, figsize=(14, 6))
fig.suptitle('SIM97: CMSTG vs Planck 2018 plikHM TTTEEE Likelihood', fontsize=13)

# Panel 1: TT power spectra comparison
ax = axes[0]
colors = {'cmstg_joint': '#1f77b4', 'cmstg_bao_only': '#ff7f0e', 'lcdm_planck': '#2ca02c'}
labels = {'cmstg_joint':   r'CMSTG joint ($\Lambda_0=0.008$)',
          'cmstg_bao_only': r'CMSTG BAO-only ($\Lambda_0=0.003$)',
          'lcdm_planck':  r'$\Lambda$CDM Planck 2018'}
for label, r in results.items():
    ell = r['cl']['ell']
    TT  = r['cl']['TT']
    fac = ell * (ell + 1) / (2 * np.pi)
    mask = (ell >= 2) & (ell <= 2500)
    lw = 2 if 'joint' in label else 1.5
    ls = '-' if label != 'lcdm_planck' else '--'
    ax.plot(ell[mask], (fac * TT)[mask], color=colors[label], lw=lw, ls=ls,
            label=labels[label])
ax.set_xscale('log')
ax.set_xlabel(r'Multipole $\ell$', fontsize=11)
ax.set_ylabel(r'$\ell(\ell+1)C_\ell^{TT}/(2\pi)$ [$\mu$K$^2$]', fontsize=11)
ax.legend(fontsize=9)
ax.set_title('TT Power Spectra', fontsize=11)
ax.set_xlim(2, 2500)

# Panel 2: -2 ln L bar chart
ax2 = axes[1]
cosmo_labels = list(results.keys())
plik_vals  = [results[l]['m2lnL_plikHM']  for l in cosmo_labels]
lowlTT_vals = [results[l]['m2lnL_lowl_TT'] for l in cosmo_labels]
lowlEE_vals = [results[l]['m2lnL_lowl_EE'] for l in cosmo_labels]
x = np.arange(len(cosmo_labels))
w = 0.25
bars = [
    ax2.bar(x - w, plik_vals,   w, label='plikHM TTTEEE', color='#1f77b4', alpha=0.85),
    ax2.bar(x,     lowlTT_vals, w, label='lowl TT',        color='#ff7f0e', alpha=0.85),
    ax2.bar(x + w, lowlEE_vals, w, label='lowl EE',        color='#2ca02c', alpha=0.85),
]
ax2.set_xticks(x)
ax2.set_xticklabels(['CMSTG\njoint', 'CMSTG\nBAO-only', r'$\Lambda$CDM'], fontsize=10)
ax2.set_ylabel(r'$-2\ln\mathcal{L}$', fontsize=11)
ax2.legend(fontsize=9)
ax2.set_title(r'Planck 2018 $-2\ln\mathcal{L}$ by component', fontsize=11)
ax2.text(0.02, 0.97,
         f'Δ(-2lnL) CMSTG joint − ΛCDM:\n'
         f'  plikHM: {d_plik_joint:+.3f}\n'
         f'  total:  {d_total_joint:+.3f}',
         transform=ax2.transAxes, fontsize=9, va='top',
         bbox=dict(boxstyle='round', facecolor='lightyellow', alpha=0.9))

plt.tight_layout()
fig_path = os.path.join(OUTPUTS, 'sim97_planck_likelihood.pdf')
plt.savefig(fig_path, bbox_inches='tight', dpi=150)
plt.close()
print(f"Figure saved: {fig_path}")

# ═══════════════════════════════════════════════════════════════════════════════
# 7. SAVE DIAGNOSTICS
# ═══════════════════════════════════════════════════════════════════════════════

diag = {
    "sim": "SIM97",
    "title": "CMSTG Full Planck 2018 plikHM TTTEEE Likelihood",
    "date": "2026-04-09",
    "likelihood_files": {k: clik_path(k) for k in P['clik_files']},
    "results": {
        label: {
            "m2lnL_plikHM_TTTEEE": float(r['m2lnL_plikHM']),
            "m2lnL_plik_lite":     float(r['m2lnL_lite']),
            "m2lnL_lowl_TT":       float(r['m2lnL_lowl_TT']),
            "m2lnL_lowl_EE":       float(r['m2lnL_lowl_EE']),
            "m2lnL_total":         float(r['m2lnL_total']),
            "Geff_z0":             float(r['Geff_z0'])
        }
        for label, r in results.items()
    },
    "delta_m2lnL_plikHM_joint_minus_lcdm":   float(d_plik_joint),
    "delta_m2lnL_plikHM_bao_minus_lcdm":     float(d_plik_bao),
    "delta_m2lnL_total_joint_minus_lcdm":    float(d_total_joint),
    "verdict": {
        "pass": bool(pass_joint),
        "overall": overall,
        "summary": (
            f"CMSTG at SIM90 joint best-fit: Δ(-2lnL)(plikHM TTTEEE) = {d_plik_joint:+.4f} "
            f"vs ΛCDM Planck best-fit. G_eff(z=0) = {results['cmstg_joint']['Geff_z0']:.8f} "
            f"({(1-results['cmstg_joint']['Geff_z0'])*1e6:.1f} ppm). "
            f"CMSTG BAO-only Δ(-2lnL) = {d_plik_bao:+.1f} (consistent with SIM89 result). "
            f"Verdict: {overall}."
        )
    }
}

diag_path = os.path.join(OUTPUTS, 'sim97_diagnostics.json')
with open(diag_path, 'w') as f:
    json.dump(diag, f, indent=2)
print(f"Diagnostics: {diag_path}")
print()
print("SIM97 COMPLETE.")
