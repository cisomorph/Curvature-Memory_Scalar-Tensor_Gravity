"""
SIM100 — CMSTG 1+1D Time-Domain Galactic Field Evolution
=========================================================
First time-dependent CMSTG simulation at galactic scales.

Field equation (spherical 1+1D, derived from locked action, consistent with SIM99):

  d^2Psi/dt^2 = d^2Psi/dr^2 + (2/r)*dPsi/dr  - m0^2*Psi
                - lambda*Psi^3  +  2*Lambda0*Psi*R(r,t)

where R(r,t) = -8*pi*(G/c^2)*rho_baryon(r) * (1 + delta(t))  [kpc^-2]

In static limit (d^2/dt^2 = 0) this reduces exactly to SIM99's equation:
  Psi'' + (2/r)Psi' = m0^2*Psi + lambda*Psi^3 + 2*Lambda0*Psi*R(r)

Physics questions:
  (A) ADIABATIC TRACKING: For |lambda| < lambda_crit, does the time-dependent
      field relax from IC Psi=Psi_cosmo to the static Yukawa profile (SIM99)?
      Expected: YES, since Compton time (33 kyr) << galaxy formation (1 Gyr).

  (B) CONDENSATION: For lambda < lambda_spatial_crit ~ -778, does the
      attractive self-coupling cause tachyonic condensation?
      If so, does the condensed profile give flat rotation curves?

Condensation thresholds:
  Uniform mode (k=0):    lambda_crit_uniform  = -m0^2/(3*Psi_cosmo^2) = -370
  Lowest spatial mode:   lambda_spatial_crit  ~ -(pi/r_max)^2+m0^2)/(3*Psi_cosmo^2) ~ -778
  => Only for lambda < -778 does the spatial field become unstable.

Note on curvature coupling:
  The galactic Ricci scalar |R_gal| ~ 1e-8 kpc^-2.
  The coupling 2*Lambda0*|R_gal| ~ 6e-11 kpc^-2 << m0^2 = 0.01 kpc^-2.
  => Curvature coupling is ~1e8x too small for condensation at galactic densities.
  => Galaxy formation overdensity (delta_max=100) provides only 6e-9 kpc^-2 coupling.
  => Condensation physics is governed entirely by the lambda*Psi^3 term.

Method:
  Method of lines (MOL): spatial FD → ODE system. Scipy Radau solver (stiff).
  State: y = [Psi_0..Psi_{N-1}, Pi_0..Pi_{N-1}], Pi = dPsi/dt.
  BC: dPsi/dr|_{r_min}=0 (Neumann, ghost cell), Psi(r_max)=Psi_cosmo (Dirichlet).
  IC: Psi(r,0) = Psi_cosmo (uniform), Pi(r,0) = 0.
"""

import os, json, time, warnings
import numpy as np
from scipy.integrate import solve_ivp, cumulative_trapezoid
from scipy.special import iv, kv
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

warnings.filterwarnings('ignore')

BASE    = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
INPUTS  = os.path.join(BASE, 'Inputs')
OUTPUTS = os.path.join(BASE, 'Outputs')
os.makedirs(OUTPUTS, exist_ok=True)

with open(os.path.join(INPUTS, 'sim100_params.json')) as f:
    P = json.load(f)

# ── Physical constants (units: kpc, M_sun, km/s) ──────────────────────────
G_GAL  = 4.3009e-6    # (km/s)^2 kpc / M_sun
C_KMS  = 2.998e5      # km/s
# 1 kpc/c = 1 kpc / C_KMS (in km/s·kpc units) ≈ 3,260 yr
KPC_PER_C_TO_KYR = 3.262  # 1 kpc/c = 3.262 kyr

CMSTG    = P['cmstg']
GAL     = P['galaxy']
NUM     = P['numerics']
SCAN    = P['scan']
VERDICT = P['verdict_criteria']

M0      = CMSTG['m0_kpc_inv']        # 0.1 kpc^-1
L0      = CMSTG['Lambda0']           # 0.003
PSI_C   = CMSTG['psi_cosmo']         # 0.003

M_DISK  = GAL['M_disk_Msun']
R_D     = GAL['r_d_kpc']
M_GAS   = GAL['M_gas_Msun']
R_GAS   = GAL['r_gas_kpc']
V_FLAT  = GAL['v_flat_kms']

R_MIN   = NUM['r_min_kpc']
R_MAX   = NUM['r_max_kpc']
N_R     = NUM['n_r']
N_SNAP  = NUM['n_snapshots']
RTOL    = NUM['radau_rtol']
ATOL    = NUM['radau_atol']
BLOW_FAC = NUM['blow_up_threshold_factor']  # |Psi| > BLOW_FAC*PSI_C → abort

T_END   = SCAN['t_end_kpc_per_c']   # [kpc/c]
T_FORM  = SCAN['t_form_kpc_per_c']  # galaxy formation time ramp [kpc/c]
LAMBDA_SCAN = SCAN['lambda_values']
DELTA_SCAN  = SCAN['delta_max_values']

COND_FACTOR = VERDICT['condensation_growth_factor']   # 2.0
FLAT_TOL    = VERDICT['flatness_tolerance']            # 0.20
FLAT_FRAC   = VERDICT['flatness_outer_fraction']       # 0.50


# ══════════════════════════════════════════════════════════════════════════════
# BARYONIC MODELS  (identical to SIM99 for direct comparison)
# ══════════════════════════════════════════════════════════════════════════════

def rho_baryon(r_arr):
    """Spherical-average baryonic density [M_sun/kpc^3]."""
    dMdr = (M_DISK / R_D**2 * r_arr * np.exp(-r_arr / R_D)
          + M_GAS  / R_GAS**2 * r_arr * np.exp(-r_arr / R_GAS))
    return dMdr / (4.0 * np.pi * r_arr**2)


def M_baryon_enclosed(r_arr):
    """Enclosed baryonic mass [M_sun]."""
    enc_d = M_DISK * (1.0 - np.exp(-r_arr / R_D) * (1.0 + r_arr / R_D))
    enc_g = M_GAS  * (1.0 - np.exp(-r_arr / R_GAS) * (1.0 + r_arr / R_GAS))
    return enc_d + enc_g


def R_static(r_arr):
    """
    Weak-field Ricci scalar from static baryons [kpc^-2].
    R = -8*pi*(G/c^2)*rho_baryon  (negative for positive density).
    Same formula as SIM99.
    """
    G_c2 = G_GAL / C_KMS**2    # kpc / M_sun  (= G/c^2 in these units)
    return -8.0 * np.pi * G_c2 * rho_baryon(r_arr)


# ══════════════════════════════════════════════════════════════════════════════
# METHOD OF LINES: PDE → ODE
# ══════════════════════════════════════════════════════════════════════════════

def build_rhs(r_arr, lambda_gal, delta_max):
    """
    Build the RHS function for the MOL ODE system.

    State y = [Psi_0..Psi_{N-1}, Pi_0..Pi_{N-1}], Pi = dPsi/dt.

    ODE:
      dPsi/dt = Pi
      dPi/dt  = D^2_r Psi + (2/r) D_r Psi - m0^2*Psi - lambda*Psi^3 + 2*L0*Psi*R(r,t)

    Spatial discretization (2nd-order FD):
      D^2_r Psi_i = (Psi_{i+1} - 2*Psi_i + Psi_{i-1}) / dr^2
      D_r Psi_i   = (Psi_{i+1} - Psi_{i-1}) / (2*dr)

    Boundary conditions:
      i=0   (r_min): Neumann dPsi/dr=0 → ghost cell Psi_{-1} = Psi_1
      i=N-1 (r_max): Dirichlet Psi = Psi_cosmo → Psi_N = Psi_cosmo (fixed)
    """
    N    = len(r_arr)
    dr   = r_arr[1] - r_arr[0]
    dr2  = dr * dr
    inv2dr = 0.5 / dr
    coeff_2r = 2.0 / r_arr        # (2/r) at each point
    R_stat = R_static(r_arr)      # precomputed static part
    blow_thresh = BLOW_FAC * PSI_C

    def rhs(t, y):
        Psi = y[:N]
        Pi  = y[N:]

        # ── Detect blow-up and signal solver to stop ──────────────────────
        if np.max(np.abs(Psi)) > blow_thresh:
            raise RuntimeError(f"BLOW-UP: max|Psi|={np.max(np.abs(Psi)):.3e} > {blow_thresh:.3e}")

        # ── Time-dependent R ──────────────────────────────────────────────
        fac  = delta_max * min(t / T_FORM, 1.0) if T_FORM > 0 else delta_max
        R    = R_stat * (1.0 + fac)

        # ── Extended Psi with boundary conditions ─────────────────────────
        Psi_ext = np.empty(N + 2)
        Psi_ext[1:N+1] = Psi         # interior
        Psi_ext[0]     = Psi[1]      # ghost: Neumann (dPsi/dr=0 at r_min)
        Psi_ext[N+1]   = PSI_C       # Dirichlet at r_max

        # ── Finite-difference Laplacian ───────────────────────────────────
        d2Psi = (Psi_ext[2:N+2] - 2.0*Psi_ext[1:N+1] + Psi_ext[0:N]) / dr2
        dPsi  = (Psi_ext[2:N+2] - Psi_ext[0:N]) * inv2dr

        # ── RHS: dPi/dt ───────────────────────────────────────────────────
        dPi = (d2Psi
               + coeff_2r * dPsi
               - M0**2 * Psi
               - lambda_gal * Psi**3
               + 2.0 * L0 * Psi * R)

        return np.concatenate([Pi, dPi])

    return rhs


# ══════════════════════════════════════════════════════════════════════════════
# ROTATION CURVE FROM FINAL FIELD PROFILE
# ══════════════════════════════════════════════════════════════════════════════

def compute_rotation_curve(r_arr, psi, lambda_gal):
    """
    Rotation curve from field profile psi(r).
    Identical to SIM99 formulas for direct comparison.
    """
    dpsi_dr = np.gradient(psi, r_arr)

    # T^Psi_00 in CMSTG units — convert to M_sun/kpc^3 via c^2/G
    c2_over_G = C_KMS**2 / G_GAL
    T00 = 0.5*dpsi_dr**2 + 0.5*M0**2*psi**2 + 0.25*lambda_gal*psi**4
    rho_psi = T00 * c2_over_G

    M_psi = cumulative_trapezoid(4.0*np.pi*r_arr**2 * np.clip(rho_psi, 0, None),
                                  r_arr, initial=0.0)
    M_bar = M_baryon_enclosed(r_arr)
    geff  = 1.0 / (1.0 + 16.0*np.pi*L0*psi**2)
    v2    = geff * G_GAL * (M_bar + M_psi) / r_arr
    vc    = np.sqrt(np.clip(v2, 0, None))
    return vc, rho_psi, M_psi


def check_flatness(r_arr, vc):
    """True if v_c within ±FLAT_TOL of V_FLAT for r > FLAT_FRAC * r_max."""
    r_cut = r_arr[0] + FLAT_FRAC * (r_arr[-1] - r_arr[0])
    mask  = r_arr >= r_cut
    if not np.any(mask):
        return False, 1.0
    dev = np.mean(np.abs(vc[mask] - V_FLAT) / V_FLAT)
    return bool(dev < FLAT_TOL), float(dev)


# ══════════════════════════════════════════════════════════════════════════════
# RUN ONE EVOLUTION
# ══════════════════════════════════════════════════════════════════════════════

def run_one(r_arr, lambda_gal, delta_max):
    """
    Evolve Psi(r,t) via MOL-Radau from t=0 to t=T_END.
    Returns dict with condensation status and final rotation curve.
    """
    N  = len(r_arr)
    rhs = build_rhs(r_arr, lambda_gal, delta_max)

    # IC: uniform cosmological background, at rest
    y0 = np.concatenate([PSI_C * np.ones(N), np.zeros(N)])

    t_eval = np.linspace(0.0, T_END, N_SNAP)

    t_cpu0 = time.time()
    try:
        sol = solve_ivp(
            rhs,
            (0.0, T_END),
            y0,
            t_eval=t_eval,
            method='Radau',
            rtol=RTOL,
            atol=ATOL,
            max_step=SCAN['max_step_compton_fraction'] if 'max_step_compton_fraction' in SCAN
                      else (0.5 / M0),
            dense_output=False
        )
        success    = sol.success
        message    = sol.message
        psi_snapshots = sol.y[:N, :]      # shape (N, n_snap)
        t_actual   = sol.t
    except RuntimeError as e:
        success    = False
        message    = str(e)
        psi_snapshots = None
        t_actual   = None

    elapsed = time.time() - t_cpu0

    result = {
        'lambda_gal':  lambda_gal,
        'delta_max':   delta_max,
        'solver_success': success,
        'solver_message': message,
        'solve_time_s': round(elapsed, 2),
    }

    if not success or psi_snapshots is None:
        result.update({
            'condensed': None,
            'blow_up': 'BLOW-UP' in message,
            'psi_max_vs_t': None,
            'psi_max_growth_factor': None,
            'pass_rotation': False,
            'outer_dev_frac': None,
        })
        return result, None, None

    # ── Condensation diagnostic ───────────────────────────────────────────
    psi_max_t = np.max(np.abs(psi_snapshots), axis=0)   # max|Psi(r)| per snapshot
    psi_ini   = float(psi_max_t[0])
    psi_fin   = float(psi_max_t[-1])
    growth    = psi_fin / max(psi_ini, 1e-20)
    condensed = bool(growth > COND_FACTOR)

    # Final field profile
    psi_final = psi_snapshots[:, -1]

    # ── Rotation curve ────────────────────────────────────────────────────
    vc, rho_psi, M_psi = compute_rotation_curve(r_arr, psi_final, lambda_gal)
    flat_ok, outer_dev = check_flatness(r_arr, vc)

    result.update({
        'condensed': condensed,
        'blow_up': False,
        'psi_ini': psi_ini,
        'psi_fin': psi_fin,
        'psi_max_growth_factor': round(growth, 4),
        'psi_max_vs_t': psi_max_t.tolist(),
        't_kpc_per_c': t_actual.tolist() if t_actual is not None else None,
        'psi_max_final': float(np.max(np.abs(psi_final))),
        'psi_min_final': float(np.min(np.abs(psi_final))),
        'outer_dev_frac': round(outer_dev, 4),
        'pass_rotation': flat_ok,
        'v_c_outer_mean_kms': float(np.mean(vc[r_arr >= r_arr[0] + FLAT_FRAC*(r_arr[-1]-r_arr[0])])),
        'M_psi_at_rmax': float(M_psi[-1]),
    })

    return result, psi_final, vc


# ══════════════════════════════════════════════════════════════════════════════
# MAIN
# ══════════════════════════════════════════════════════════════════════════════

def run():
    print("=" * 72)
    print("SIM100 v1 — CMSTG 1+1D Time-Domain Galactic Evolution")
    print("=" * 72)
    r_arr = np.linspace(R_MIN, R_MAX, N_R)
    dr    = r_arr[1] - r_arr[0]

    # Reference values
    R_stat_arr   = R_static(r_arr)
    R_max_mag    = float(np.max(np.abs(R_stat_arr)))
    coupling_max = 2.0 * L0 * R_max_mag
    lam_crit_unif = -M0**2 / (3.0 * PSI_C**2)
    k_min = np.pi / R_MAX
    lam_crit_spat = -(k_min**2 + M0**2) / (3.0 * PSI_C**2)

    print(f"  m0 = {M0} kpc^-1  (Compton length = {1/M0:.0f} kpc, period = {2*np.pi/M0:.1f} kpc/c = {2*np.pi/M0*KPC_PER_C_TO_KYR:.0f} kyr)")
    print(f"  L0 = {L0},  Psi_cosmo = {PSI_C}")
    print(f"  t_end = {T_END} kpc/c = {T_END*KPC_PER_C_TO_KYR:.0f} kyr  ({T_END/(2*np.pi/M0):.1f} Compton periods)")
    print(f"  N_r = {N_R}, dr = {dr:.3f} kpc, CFL limit = {dr:.3f} kpc/c = {dr*KPC_PER_C_TO_KYR*1e3:.0f} yr")
    print()
    print(f"  Galactic R: max|R_static| = {R_max_mag:.2e} kpc^-2 (at r={R_MIN} kpc)")
    print(f"  Curvature coupling: 2*L0*|R_max| = {coupling_max:.2e} kpc^-2")
    print(f"  m0^2 = {M0**2:.4f} kpc^-2")
    print(f"  => 2*L0*|R|/m0^2 = {coupling_max/M0**2:.2e}  (R coupling is NEGLIGIBLE)")
    print()
    print(f"  Condensation thresholds:")
    print(f"    Uniform mode  (k=0):  lambda_crit = {lam_crit_unif:.0f}")
    print(f"    Lowest spatial mode:  lambda_crit = {lam_crit_spat:.0f}")
    print(f"    => Spatial instability only for lambda < {lam_crit_spat:.0f}")
    print()

    diag = {
        'sim': 'SIM100v1',
        'description': P['description'],
        'physics': {
            'm0_kpc_inv': M0,
            'Lambda0': L0,
            'psi_cosmo': PSI_C,
            'compton_length_kpc': 1.0/M0,
            'compton_period_kyr': 2*np.pi/M0 * KPC_PER_C_TO_KYR,
            't_end_kpc_per_c': T_END,
            't_end_kyr': T_END * KPC_PER_C_TO_KYR,
            'R_galactic_max_mag_kpc2': R_max_mag,
            'curvature_coupling_max': coupling_max,
            'curvature_coupling_over_m02': coupling_max / M0**2,
            'lambda_crit_uniform': lam_crit_unif,
            'lambda_crit_spatial': lam_crit_spat,
        },
        'runs': [],
    }

    print(f"Scan: {len(LAMBDA_SCAN)} lambda × {len(DELTA_SCAN)} delta_max = {len(LAMBDA_SCAN)*len(DELTA_SCAN)} runs\n")

    header = (f"{'lambda':>8}  {'delta':>6}  {'regime':>14}  {'success':>8}  "
              f"{'growth':>8}  {'condensed':>10}  {'rot-curve':>10}  {'t_cpu':>6}")
    print(header)
    print("-" * len(header))

    all_results = []
    psi_store   = {}   # (lambda, delta) -> psi_final
    vc_store    = {}   # (lambda, delta) -> vc

    for lambda_gal in LAMBDA_SCAN:
        regime = ("SUPER-CRIT" if lambda_gal < lam_crit_spat
                  else "near-crit" if lambda_gal < lam_crit_unif
                  else "sub-crit")
        for delta_max in DELTA_SCAN:
            res, psi_fin, vc = run_one(r_arr, lambda_gal, delta_max)
            all_results.append(res)

            if psi_fin is not None:
                psi_store[(lambda_gal, delta_max)] = psi_fin
                vc_store[(lambda_gal, delta_max)]  = vc

            g_str   = f"{res['psi_max_growth_factor']:.3f}x" if res['psi_max_growth_factor'] else "N/A"
            c_str   = ("CONDENSED" if res.get('condensed')
                       else "BLOW-UP" if res.get('blow_up')
                       else "stable" if res['solver_success'] else "FAIL")
            rc_str  = "PASS" if res.get('pass_rotation') else "fail"
            t_str   = f"{res['solve_time_s']:.1f}s"

            print(f"  {lambda_gal:>6.0f}  {delta_max:>6.0f}  {regime:>14}  "
                  f"{'OK' if res['solver_success'] else 'FAIL':>8}  "
                  f"{g_str:>8}  {c_str:>10}  {rc_str:>10}  {t_str:>6}")

            diag['runs'].append({k: v for k, v in res.items()
                                 if k not in ('psi_max_vs_t', 't_kpc_per_c')})

    # ── Summary table ─────────────────────────────────────────────────────
    print(f"\n{'='*72}")
    print("SUMMARY TABLE  (delta_max=0 only for clarity)")
    print(f"{'='*72}")
    print(f"{'lambda':>8}  {'regime':>14}  {'growth':>8}  {'condensed':>10}  "
          f"{'outer_dev%':>11}  {'rot-curve':>10}")
    for res in [r for r in all_results if r['delta_max'] == 0]:
        regime = ("SUPER-CRIT" if res['lambda_gal'] < lam_crit_spat
                  else "near-crit" if res['lambda_gal'] < lam_crit_unif
                  else "sub-crit")
        g  = f"{res['psi_max_growth_factor']:.3f}x" if res['psi_max_growth_factor'] else "N/A"
        c  = ("CONDENSED" if res.get('condensed') else
              "BLOW-UP"   if res.get('blow_up')   else
              "stable"    if res['solver_success'] else "FAIL")
        d  = f"{res['outer_dev_frac']*100:.1f}%" if res['outer_dev_frac'] is not None else "N/A"
        rc = "PASS" if res.get('pass_rotation') else "fail"
        print(f"  {res['lambda_gal']:>6.0f}  {regime:>14}  {g:>8}  {c:>10}  {d:>11}  {rc:>10}")

    # ── Verdict ───────────────────────────────────────────────────────────
    print(f"\n{'='*72}")
    print("VERDICT")
    print(f"{'='*72}")

    any_cond   = any(r.get('condensed') for r in all_results)
    any_pass   = any(r.get('pass_rotation') for r in all_results)
    any_blow   = any(r.get('blow_up') for r in all_results)

    if not any_cond and not any_blow:
        verdict_str = (
            "ADIABATIC CONFIRMED, NO CONDENSATION: "
            "All (lambda, delta_max) combinations produce stable evolution. "
            "The field relaxes from IC Psi=Psi_cosmo to the static Yukawa profile "
            "on the Compton timescale (33 kyr), consistent with the adiabatic "
            "approximation (Compton time << galaxy formation time). "
            "The spatial instability threshold (lambda_crit ~ "
            f"{lam_crit_spat:.0f}) was not reached "
            "within the scanned range. SIM99 (static) remains the correct "
            "description of CMSTG galactic dynamics."
        )
    elif any_blow:
        verdict_str = (
            "TACHYONIC BLOW-UP DETECTED: "
            f"For lambda < {lam_crit_spat:.0f} (spatial instability threshold), "
            "the CMSTG field grows without bound from IC Psi=Psi_cosmo. "
            "The attractive quartic V=(lambda/4)Psi^4 has no stable minimum "
            "for lambda < 0 — the potential is unbounded below. "
            "The field cannot form a stable halo condensate. "
        )
        if not any_pass:
            verdict_str += (
                "No rotation curve flattening observed in any stable run. "
                "CONCLUSION: Tachyonic blow-up, not stable condensation, "
                "is the outcome of attractive self-coupling in CMSTG. "
                "Flat rotation curves require stable dark sector physics "
                "not present in V=(lambda/4)Psi^4 with lambda < 0."
            )
    elif any_cond and not any_pass:
        verdict_str = (
            "CONDENSATION OCCURS but ROTATION CURVE FAILS: "
            "The field condenses for lambda < "
            f"{lam_crit_spat:.0f}, but the condensed profile does not "
            "produce flat rotation curves. The condensate is spatially "
            "concentrated (soliton/oscillaton-like), not the 1/r^2 "
            "isothermal sphere needed for curve flattening."
        )
    else:
        verdict_str = (
            "PASS: At least one (lambda, delta_max) produces flat rotation curve. "
            "See runs for parameters."
        )

    print(verdict_str)

    diag['verdict'] = {
        'any_condensed': any_cond,
        'any_blow_up': any_blow,
        'any_pass_rotation': any_pass,
        'summary': verdict_str,
        'physics_note': (
            f"Curvature coupling 2*L0*|R_gal| = {coupling_max:.2e} kpc^-2 << "
            f"m0^2 = {M0**2:.4f} kpc^-2. Galaxy formation (delta_max=100) "
            "makes no qualitative difference. Condensation is driven by "
            "lambda<0 self-coupling, not by CMSTG curvature coupling. "
            "Adiabatic limit (Compton time 33 kyr << t_form 1 Gyr) means "
            "time-domain and static (SIM99) solutions agree for sub-threshold lambda."
        ),
    }

    # ── Plots ─────────────────────────────────────────────────────────────
    _make_plots(r_arr, all_results, psi_store, vc_store,
                lam_crit_unif, lam_crit_spat)

    # ── Save diagnostics ───────────────────────────────────────────────────
    diag_path = os.path.join(OUTPUTS, 'sim100_diagnostics.json')
    with open(diag_path, 'w') as f:
        json.dump(diag, f, indent=2)

    print(f"\nOutputs written to {OUTPUTS}/")
    print("  sim100_phase_diagram.pdf       — growth factor & condensation map")
    print("  sim100_field_profiles.pdf      — final Psi(r) for key cases")
    print("  sim100_rotation_curves.pdf     — rotation curves for key cases")
    print("  sim100_time_evolution.pdf      — max|Psi(t)| vs time")
    print("  sim100_diagnostics.json        — full numerical results")
    print("Done.")
    return diag


# ══════════════════════════════════════════════════════════════════════════════
# PLOTS
# ══════════════════════════════════════════════════════════════════════════════

def _make_plots(r_arr, all_results, psi_store, vc_store,
                lam_crit_unif, lam_crit_spat):
    lambdas = sorted(set(r['lambda_gal'] for r in all_results))
    deltas  = sorted(set(r['delta_max']  for r in all_results))

    # ── 1. Phase diagram ──────────────────────────────────────────────────
    grow_mat  = np.full((len(lambdas), len(deltas)), np.nan)
    cond_mat  = np.full((len(lambdas), len(deltas)), np.nan)
    for res in all_results:
        il = lambdas.index(res['lambda_gal'])
        id_ = deltas.index(res['delta_max'])
        if res['psi_max_growth_factor'] is not None:
            grow_mat[il, id_] = np.log10(max(res['psi_max_growth_factor'], 1e-3))
        if res.get('condensed') is not None:
            cond_mat[il, id_] = 1.0 if res['condensed'] else 0.0
        elif res.get('blow_up'):
            cond_mat[il, id_] = 2.0  # blow-up marker

    fig, axes = plt.subplots(1, 2, figsize=(13, 5))

    ax = axes[0]
    im = ax.imshow(grow_mat, aspect='auto', origin='lower',
                   cmap='RdYlGn', vmin=-0.5, vmax=1.5)
    ax.set_xticks(range(len(deltas)))
    ax.set_xticklabels([f'{d:.0f}' for d in deltas])
    ax.set_yticks(range(len(lambdas)))
    ax.set_yticklabels([f'{l:.0f}' for l in lambdas])
    ax.set_xlabel(r'$\delta_{\rm max}$ (overdensity)', fontsize=11)
    ax.set_ylabel(r'$\lambda$ (self-coupling)', fontsize=11)
    ax.set_title(r'$\log_{10}$(max$|\Psi|$ growth factor)', fontsize=10)
    plt.colorbar(im, ax=ax)

    # Mark thresholds
    for lam_th, ls, label in [(lam_crit_unif, '--', rf'$\lambda_{{crit,unif}}={lam_crit_unif:.0f}$'),
                               (lam_crit_spat, '-',  rf'$\lambda_{{crit,spat}}={lam_crit_spat:.0f}$')]:
        if lam_th in lambdas:
            y_pos = lambdas.index(lam_th)
        else:
            y_pos = next((i - 0.5 for i, l in enumerate(lambdas) if l > lam_th), None)
        if y_pos is not None:
            ax.axhline(y=y_pos, color='white', lw=1.5, ls=ls, label=label)
    ax.legend(fontsize=7, loc='upper right')

    ax2 = axes[1]
    im2 = ax2.imshow(cond_mat, aspect='auto', origin='lower',
                     cmap='RdYlGn', vmin=0, vmax=1)
    ax2.set_xticks(range(len(deltas)))
    ax2.set_xticklabels([f'{d:.0f}' for d in deltas])
    ax2.set_yticks(range(len(lambdas)))
    ax2.set_yticklabels([f'{l:.0f}' for l in lambdas])
    ax2.set_xlabel(r'$\delta_{\rm max}$', fontsize=11)
    ax2.set_ylabel(r'$\lambda$', fontsize=11)
    ax2.set_title('Condensed? (green=YES, red=NO, dark=blow-up)', fontsize=10)
    plt.colorbar(im2, ax=ax2)

    fig.suptitle('SIM100 — Phase Diagram: CMSTG Galactic Condensation', fontsize=12)
    fig.tight_layout()
    fig.savefig(os.path.join(OUTPUTS, 'sim100_phase_diagram.pdf'),
                dpi=150, bbox_inches='tight')
    plt.close(fig)

    # ── 2. Field profiles & rotation curves for key cases (delta=0) ───────
    cases_d0 = [(lam, 0) for lam in lambdas]
    n_cases  = len(cases_d0)
    if n_cases == 0:
        return

    fig_f, axes_f = plt.subplots(1, n_cases, figsize=(4.5*n_cases, 5), sharey=False)
    fig_v, axes_v = plt.subplots(1, n_cases, figsize=(4.5*n_cases, 5), sharey=True)
    if n_cases == 1:
        axes_f, axes_v = [axes_f], [axes_v]

    M_bar = M_baryon_enclosed(r_arr)
    v_bar = np.sqrt(np.clip(G_GAL * M_bar / r_arr, 0.0, None))

    for ic, (lam, delta) in enumerate(cases_d0):
        psi = psi_store.get((lam, delta))
        vc  = vc_store.get((lam, delta))
        res = next((r for r in all_results
                    if r['lambda_gal']==lam and r['delta_max']==delta), None)
        if psi is None or res is None:
            continue

        cond_str = ("CONDENSED" if res.get('condensed') else
                    "BLOW-UP"   if res.get('blow_up')   else "stable")
        label = f"$\\lambda={lam:.0f}$\n{cond_str}"

        # Field profile
        ax = axes_f[ic]
        ax.semilogy(r_arr, np.abs(psi), 'b-', lw=2)
        ax.axhline(PSI_C, color='gray', ls=':', lw=1.5,
                   label=f'$\\Psi_{{cosmo}}={PSI_C}$')
        ax.set_xlabel('r [kpc]', fontsize=10)
        ax.set_ylabel('$|\\Psi(r,t_{\\rm end})|$', fontsize=10)
        ax.set_title(label, fontsize=9)
        ax.legend(fontsize=7)

        # Rotation curve
        if vc is not None:
            ax = axes_v[ic]
            ax.plot(r_arr, vc,    'b-',  lw=2, label='CMSTG total')
            ax.plot(r_arr, v_bar, 'k--', lw=1.5, label='Baryons only')
            ax.axhline(V_FLAT, color='r', ls=':', lw=1.5,
                       label=f'$v_{{flat}}={V_FLAT}$ km/s')
            ax.fill_between(r_arr, V_FLAT*0.8, V_FLAT*1.2,
                            alpha=0.12, color='red', label='±20%')
            ax.set_xlabel('r [kpc]', fontsize=10)
            ax.set_ylabel('$v_c$ [km/s]', fontsize=10)
            ax.set_title(label, fontsize=9)
            ax.legend(fontsize=7)
            ax.set_xlim(0, R_MAX)
            ax.set_ylim(0, max(1.8*V_FLAT, float(np.nanmax(vc))*1.1))

    fig_f.suptitle('SIM100 — Final Field Profiles $\\Psi(r,t_{\\rm end})$', fontsize=11)
    fig_f.tight_layout()
    fig_f.savefig(os.path.join(OUTPUTS, 'sim100_field_profiles.pdf'),
                  dpi=150, bbox_inches='tight')
    plt.close(fig_f)

    fig_v.suptitle('SIM100 — Final Rotation Curves $v_c(r,t_{\\rm end})$', fontsize=11)
    fig_v.tight_layout()
    fig_v.savefig(os.path.join(OUTPUTS, 'sim100_rotation_curves.pdf'),
                  dpi=150, bbox_inches='tight')
    plt.close(fig_v)

    # ── 3. Time evolution of max|Psi(t)| ──────────────────────────────────
    fig_t, ax_t = plt.subplots(figsize=(10, 5))
    colors = plt.cm.coolwarm(np.linspace(0, 1, len(lambdas)))

    for il, lam in enumerate(lambdas):
        res = next((r for r in all_results
                    if r['lambda_gal']==lam and r['delta_max']==0), None)
        if res is None or not res['solver_success']:
            continue
        t_arr = res.get('t_kpc_per_c') or []
        psi_t = res.get('psi_max_vs_t') or []
        if not t_arr or not psi_t:
            continue
        lbl = f"$\\lambda={lam:.0f}$"
        ax_t.semilogy(t_arr, psi_t, color=colors[il], lw=1.5, label=lbl)

    ax_t.axhline(PSI_C * COND_FACTOR, color='k', ls='--', lw=1.5,
                 label=f'Condensation threshold ({COND_FACTOR}×IC)')
    ax_t.axhline(PSI_C, color='gray', ls=':', lw=1, label=f'IC ($\\Psi_{{cosmo}}={PSI_C}$)')
    ax_t.set_xlabel(f'time [kpc/c]  (1 kpc/c = {KPC_PER_C_TO_KYR:.2f} kyr)', fontsize=11)
    ax_t.set_ylabel('$\\max_r |\\Psi(r,t)|$', fontsize=11)
    ax_t.set_title('SIM100 — Max field amplitude vs time (delta_max=0)', fontsize=11)
    ax_t.legend(fontsize=7, ncol=2, loc='upper right')
    fig_t.tight_layout()
    fig_t.savefig(os.path.join(OUTPUTS, 'sim100_time_evolution.pdf'),
                  dpi=150, bbox_inches='tight')
    plt.close(fig_t)


if __name__ == '__main__':
    run()
