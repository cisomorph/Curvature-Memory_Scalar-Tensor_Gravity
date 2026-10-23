"""
cobaya_hiclass.py — run hi_class (Python wrapper `hiclassy`, built with the cmstg_attractor patch in
hi_class_model/) through Cobaya's CLASS interface. ROADMAP Phase 3.

Usage in a Cobaya info dict:
    'theory': {'hiclass': {'external': HiClassCMSTG, 'cmstg_model': True, 'cmstg_family': 0, 'extra_args': {...}}}
Sampled CMSTG parameters: cmstg_L0, cmstg_xi. They are intercepted here and turned into
hi_class's `parameters_smg = p0, Lambda0, xi, family, phi_ini (0), phi_prime_ini (0 = attractor IC)`.

Helium (pre-registered, ROADMAP 3.4): YHe = 0.2454 + 0.438(ω_b − 0.02238) + 0.16(S − 1), S = √(2F),
with F the value at the attractor, F* = ½ + ξ²/(4Λ0) (family A). The pre-registered S uses F0 today;
F* = F0 once the field has reached the attractor, which is checked at every best fit. ΛCDM and w0wa
use the same formula with S = 1.
"""
import importlib
import types
import numpy as np
import hiclassy

# the package re-exports the class under the module's name, so fetch the module itself
_cc = importlib.import_module('cobaya.theories.classy.classy')

_SHIM = types.SimpleNamespace(Class=hiclassy.HiClass,
                              CosmoComputationError=hiclassy.CosmoComputationError,
                              CosmoSevereError=hiclassy.CosmoSevereError,
                              __version__='v3.4.0')
_orig_loader = _cc.load_external_module


def _loader(module_name, *args, **kwargs):
    if module_name == 'classy':
        return _SHIM
    return _orig_loader(module_name, *args, **kwargs)


_cc.load_external_module = _loader

# Settings validated in hi_class_model/README.md
SMG_SETTINGS = dict(Omega_Lambda=0., Omega_fld=0., Omega_smg=-1., gravity_model='cmstg_attractor',
                    kineticity_safe_smg=1e-6, a_ini_over_a_today_default=1e-10, a_ini_test_qs_smg=1e-10,
                    pert_ic_ini_z_ref_smg=1e9, cs2_safe_smg=1e-5)

# Planck 2018 baseline neutrinos (one massive 0.06 eV state, N_eff = 3.046)
NEUTRINOS = dict(N_ur=2.0328, N_ncdm=1, m_ncdm=0.06, T_ncdm=0.71611)


def yhe(omega_b, S=1.0):
    return 0.2454 + 0.438 * (omega_b - 0.02238) + 0.16 * (S - 1.0)


def F_attractor(L0, xi, family=0):
    if family == 0:
        return 0.5 + xi * xi / (4.0 * L0)
    from scipy.optimize import brentq
    Fp = lambda p: -2 * L0 * p * np.exp(-2 * L0 * p * p) + xi
    if xi >= np.sqrt(L0) * np.exp(-0.5):
        return np.nan                      # no attractor in family B
    ps = brentq(Fp, 1e-12, 0.5 / np.sqrt(L0))
    return 0.5 * np.exp(-2 * L0 * ps * ps) + xi * ps


class HiClassCMSTG(_cc.classy):
    cmstg_model: bool = False   # True: CMSTG attractor; False: plain CLASS (ΛCDM, or w0wa via the fluid)
    cmstg_family: int = 0       # 0 subtractive (SIM132), 1 exponential (SIM134)

    def initialize(self):
        super().initialize()
        for k, v in NEUTRINOS.items():
            self.extra_args.setdefault(k, v)
        if self.cmstg_model:
            for k, v in SMG_SETTINGS.items():
                self.extra_args.setdefault(k, v)

    def set(self, params_values_dict):
        p = dict(params_values_dict)
        L0 = p.pop('cmstg_L0', None)
        xi = p.pop('cmstg_xi', None)
        S = 1.0
        if self.cmstg_model:
            self.extra_args['parameters_smg'] = f"0.7, {L0:.12g}, {xi:.12g}, {int(self.cmstg_family)}, 0., 0."
            Fs = F_attractor(L0, xi, int(self.cmstg_family))
            S = float(np.sqrt(2.0 * Fs)) if np.isfinite(Fs) and Fs > 0 else 1.0
        self.extra_args['YHe'] = yhe(p['omega_b'], S)
        super().set(p)
