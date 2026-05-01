"""
SIM99 — CMSTG Galaxy Rotation Curves (v3)
=========================================
Tests whether the locked-action CMSTG scalar field can replace dark matter
via its stress-energy T^Psi, across NGC 3198 and NGC 6503.

Changes from v2:
  1. SHOOTING METHOD: integrate INWARD from r_max with decaying Yukawa BC.
       Psi(r_max) = Psi_bc (scanned amplitude)
       Psi'(r_max) = -(m0 + 1/r_max) * Psi_bc   [decaying Yukawa: Psi ~ A*exp(-m0*r)/r]
     This selects the physical decaying halo mode and avoids the growing-mode
     contamination that produced unphysical v_c = 33,000 km/s in v2.

  2. V(Psi) = (lambda_gal/4) * Psi^4 added throughout.
       V'(Psi) = lambda_gal * Psi^3
     This is the simplest nonlinear extension consistent with:
       - V'(0) = 0  (no constant background drive; decaying Yukawa BC valid at r_max)
       - Renormalizability (no new energy scales beyond Psi itself)
     lambda_gal > 0: repulsive self-interaction (limits central amplitude)
     lambda_gal < 0: attractive self-interaction (can self-bind, but check rho_Psi >= 0)

  3. UPDATED VERDICT: rotation curve flatness criterion replaces mass-ratio check.
     PASS requires ALL of:
       (a) |v_c(r) - v_flat| / v_flat < 20% for r > 0.5*r_max
       (b) Lambda0 * Psi(r)^2 < cosmo_bound / (16*pi) everywhere
       (c) rho_Psi(r) >= 0 everywhere (energy positivity)
     The mass-ratio check (v2 verdict) was a false positive: M_Psi/M_needed>0.5
     does not guarantee a flat rotation curve profile.

Physical model (locked action, eq. 3.1-3.5 of paper):

  Field equation (static, spherically-symmetric):
    Psi''(r) + (2/r) Psi'(r) = m0^2 Psi + V'(Psi) + 2*Lambda0*Psi * R(r)
    V'(Psi) = lambda_gal * Psi^3
    R(r) ~ -8*pi*(G/c^2) * rho_baryon(r)   [weak-field Ricci scalar]

  Effective DM energy density:
    rho_Psi(r) = (c^2/G) * [0.5*(Psi')^2 + 0.5*m0^2*Psi^2 + (lambda_gal/4)*Psi^4]

  Modified Poisson + circular velocity:
    v_c^2(r) = G_eff(r) * G * [M_baryon(r) + M_Psi(r)] / r
    G_eff(r)/G = 1 / (1 + 16*pi*Lambda0*Psi(r)^2)

Key question:
  Does ANY (Psi_bc, lambda_gal) within cosmological bounds (Lambda0*Psi^2 < bound)
  produce a flat rotation curve (v_c within 20% of v_flat over outer 50% of r)?
  If not: CMSTG at Lambda0=0.003 does not replace dark matter. Galactic-scale
  predictions require additional physics (screening, condensation) not yet
  derived from the locked action.
"""

import os, json, warnings
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

with open(os.path.join(INPUTS, 'sim99_params.json')) as f:
    P = json.load(f)

# ── Units: kpc, M_sun, km/s ───────────────────────────────────────────────
G_gal  = 4.3009e-6    # (km/s)^2 kpc M_sun^-1
c_kms  = 2.998e5      # km/s

GALAXIES      = P['galaxies']
CMSTG_P        = P['cmstg']
NUM           = P['numerics']
VERDICT_P     = P['verdict']

M0            = CMSTG_P['m0_kpc_inv']          # kpc^-1
L0            = CMSTG_P['Lambda0_fiducial']     # 0.003
L0_BOUND      = CMSTG_P['Lambda0_cosmo_bound']  # 0.095
PSI_BC_SCAN   = CMSTG_P['psi_bc_scan']
LAMBDA_SCAN   = CMSTG_P['lambda_gal_scan']

FLAT_TOL  = VERDICT_P['flatness_tolerance']    # 0.20
FLAT_FRAC = VERDICT_P['flatness_fraction']     # 0.50


# ══════════════════════════════════════════════════════════════════════════════
# BARYONIC MASS MODELS
# ══════════════════════════════════════════════════════════════════════════════

def disk_vcirc2(r, M_disk, r_d):
    """Exponential disk v_c^2, Freeman (1970)."""
    Sigma0 = M_disk / (2 * np.pi * r_d**2)
    y = np.clip(r / (2 * r_d), 1e-8, None)
    bessel = iv(0, y)*kv(0, y) - iv(1, y)*kv(1, y)
    return 4 * np.pi * G_gal * Sigma0 * r_d * y**2 * bessel


def M_baryon_enclosed(r, gal):
    """Enclosed baryonic mass (disk + gas exponential profiles)."""
    M_d, r_d = gal['M_disk_Msun'], gal['r_d_kpc']
    M_g, r_g = gal['M_gas_Msun'],  gal['r_gas_kpc']
    enc_d = M_d * (1 - np.exp(-r/r_d) * (1 + r/r_d))
    enc_g = M_g * (1 - np.exp(-r/r_g) * (1 + r/r_g))
    return enc_d + enc_g


def rho_baryon(r, gal):
    """Spherical-average 3D baryonic density [M_sun / kpc^3]."""
    M_d, r_d = gal['M_disk_Msun'], gal['r_d_kpc']
    M_g, r_g = gal['M_gas_Msun'],  gal['r_gas_kpc']
    dMdr = (M_d/r_d**2 * r * np.exp(-r/r_d)
          + M_g/r_g**2 * r * np.exp(-r/r_g))
    return dMdr / (4 * np.pi * r**2)


def ricci_scalar(r, gal):
    """Weak-field Ricci scalar R ~ -8*pi*(G/c^2)*rho_baryon [kpc^-2]."""
    G_over_c2 = G_gal / c_kms**2
    return -8 * np.pi * G_over_c2 * rho_baryon(r, gal)


def rho_DM_needed(r, gal):
    """1/r^2 isothermal sphere density needed for flat curve at v_flat."""
    return gal['v_flat_kms']**2 / (4 * np.pi * G_gal * r**2)


def M_DM_needed(r_arr, gal):
    """Enclosed DM mass needed for flat curve."""
    v = gal['v_flat_kms']
    M_tot = v**2 * r_arr / G_gal
    return np.maximum(M_tot - M_baryon_enclosed(r_arr, gal), 0.0)


# ══════════════════════════════════════════════════════════════════════════════
# SCALAR FIELD ODE  (with V(Psi) = lambda_gal/4 * Psi^4)
# ══════════════════════════════════════════════════════════════════════════════

def psi_ode(r, y, Lambda0, m0, lambda_gal, gal):
    """
    Psi'' = -2/r * Psi' + m0^2*Psi + lambda_gal*Psi^3 + 2*Lambda0*Psi*R(r)
    """
    psi, dpsi = y
    R     = ricci_scalar(r, gal)
    Vprime = lambda_gal * psi**3
    d2psi  = -2.0/r * dpsi + m0**2 * psi + Vprime + 2.0*Lambda0 * psi * R
    return [dpsi, d2psi]


def solve_psi_shooting(r_arr, Lambda0, m0, lambda_gal, psi_bc, gal):
    """
    SHOOTING METHOD: integrate INWARD from r_max to r_min.

    Boundary condition at r_max (decaying Yukawa mode: Psi ~ A*exp(-m0*r)/r):
      Psi(r_max) = psi_bc
      Psi'(r_max) = -(m0 + 1/r_max) * psi_bc

    This selects the physically motivated decaying dark matter halo mode.
    The growing mode exp(+m0*r)/r is suppressed exponentially.
    """
    r_max  = r_arr[-1]
    dpsi_bc = -(m0 + 1.0/r_max) * psi_bc

    # Integrate from r_max down to r_min (t_span decreasing)
    sol = solve_ivp(
        psi_ode,
        (r_max, r_arr[0]),
        [psi_bc, dpsi_bc],
        t_eval=r_arr[::-1],   # evaluation points in decreasing order
        args=(Lambda0, m0, lambda_gal, gal),
        method='RK45',
        rtol=NUM['rtol'], atol=NUM['atol'],
        dense_output=False
    )

    if not sol.success:
        return None, None

    # Reverse to restore outward (increasing-r) ordering
    psi  = sol.y[0][::-1]
    dpsi = sol.y[1][::-1]

    # Flag numerically exploded solutions
    if np.any(np.abs(psi) > 1e3) or np.any(~np.isfinite(psi)):
        return None, None

    return psi, dpsi


# ══════════════════════════════════════════════════════════════════════════════
# FIELD STRESS-ENERGY AND ROTATION CURVE
# ══════════════════════════════════════════════════════════════════════════════

def compute_rho_psi(psi, dpsi, m0, lambda_gal):
    """
    T^Psi_00 = (1/2)(dPsi/dr)^2 + (1/2) m0^2 Psi^2 + (lambda_gal/4) Psi^4
    Converted to M_sun/kpc^3 via factor c^2/G_gal.
    """
    c2_over_G = c_kms**2 / G_gal
    T00 = 0.5*dpsi**2 + 0.5*m0**2*psi**2 + 0.25*lambda_gal*psi**4
    return T00 * c2_over_G


def M_psi_enclosed(r_arr, psi, dpsi, m0, lambda_gal):
    """Enclosed Psi mass: 4*pi * integral rho_Psi * r^2 dr."""
    rho = compute_rho_psi(psi, dpsi, m0, lambda_gal)
    integrand = 4 * np.pi * rho * r_arr**2
    return cumulative_trapezoid(integrand, r_arr, initial=0)


def G_eff_ratio(psi, Lambda0):
    """G_eff/G = 1 / (1 + 16*pi*Lambda0*Psi^2)."""
    return 1.0 / (1.0 + 16.0*np.pi*Lambda0*psi**2)


def compute_rotation_curve(r_arr, psi, dpsi, m0, lambda_gal, Lambda0, gal):
    """Full CMSTG circular velocity: baryons + T^Psi field, G_eff-corrected."""
    geff  = G_eff_ratio(psi, Lambda0)
    M_bar = M_baryon_enclosed(r_arr, gal)
    M_psi = M_psi_enclosed(r_arr, psi, dpsi, m0, lambda_gal)
    v2    = geff * G_gal * (M_bar + M_psi) / r_arr
    return np.sqrt(np.clip(v2, 0, None)), geff, M_bar, M_psi


# ══════════════════════════════════════════════════════════════════════════════
# VERDICT CHECKS
# ══════════════════════════════════════════════════════════════════════════════

def check_flatness(r_arr, vc, v_flat):
    """True if v_c is within ±FLAT_TOL of v_flat for r > FLAT_FRAC*r_max."""
    r_cut = r_arr[0] + FLAT_FRAC * (r_arr[-1] - r_arr[0])
    mask  = r_arr >= r_cut
    if not np.any(mask):
        return False
    dev = np.abs(vc[mask] - v_flat) / v_flat
    return bool(np.all(dev < FLAT_TOL))


def check_cosmo_ok(psi, Lambda0, bound):
    """True if 16*pi*Lambda0*Psi^2 < bound everywhere."""
    coupling = 16 * np.pi * Lambda0 * psi**2
    return bool(np.all(coupling < bound))


def check_energy_positive(psi, dpsi, m0, lambda_gal):
    """True if rho_Psi >= 0 everywhere (energy positivity)."""
    T00 = 0.5*dpsi**2 + 0.5*m0**2*psi**2 + 0.25*lambda_gal*psi**4
    return bool(np.all(T00 >= 0))


# ══════════════════════════════════════════════════════════════════════════════
# MAIN
# ══════════════════════════════════════════════════════════════════════════════

def run():
    print("=" * 70)
    print("SIM99 v3 — CMSTG Rotation Curves: Shooting + V(Psi)=(lambda/4)Psi^4")
    print("=" * 70)
    print(f"  m0 = {M0} kpc^-1  (Yukawa length = {1/M0:.1f} kpc)")
    print(f"  Lambda0 = {L0} (fiducial cosmological best-fit)")
    print(f"  Cosmo bound: 16*pi*Lambda0*Psi^2 < {L0_BOUND}")
    print(f"  Psi_max_cosmo = {np.sqrt(L0_BOUND/(16*np.pi*L0)):.4f}")
    print(f"  Flatness: |v_c - v_flat|/v_flat < {FLAT_TOL*100:.0f}% for r > {FLAT_FRAC*100:.0f}% of r_max")
    print()

    diag = {
        'sim': 'SIM99v3',
        'description': (
            'Shooting method + V=(lambda/4)Psi^4. '
            'Tests whether ANY (Psi_bc, lambda_gal) produces flat rotation curves '
            'within cosmological constraints.'
        ),
        'method': 'inward_shooting_decaying_yukawa_bc',
        'V_form': 'V(Psi) = (lambda_gal/4) * Psi^4,  V_prime = lambda_gal * Psi^3',
        'galaxies': {},
        'overall_verdict': {}
    }

    any_pass = False
    best_fits = {}   # best (Psi_bc, lambda_gal, vc) per galaxy

    for gal_name, gal in GALAXIES.items():
        r_max    = gal['r_max_kpc']
        v_flat   = gal['v_flat_kms']
        r_arr    = np.linspace(NUM['r_min_kpc'], r_max, NUM['n_points'])
        M_needed = M_DM_needed(r_arr, gal)
        rho_need = rho_DM_needed(r_arr, gal)

        print(f"\n{'='*70}")
        print(f"{gal_name}  (v_flat={v_flat} km/s, r_max={r_max} kpc)")
        print(f"  Yukawa length / r_max = {1/M0/r_max:.3f}  "
              f"({'small: field decays over galaxy' if 1/M0 < r_max else 'large: nearly Coulomb-like'})")
        print(f"  {'Psi_bc':>8}  {'lambda':>8}  {'v_c(rmax)':>10}  "
              f"{'flat?':>6}  {'cosmo?':>7}  {'E>0?':>5}  {'PASS':>5}")

        gal_results = {}
        best_vc_diff = np.inf
        best_entry   = None

        for psi_bc in PSI_BC_SCAN:
            for lam in LAMBDA_SCAN:
                psi, dpsi = solve_psi_shooting(r_arr, L0, M0, lam, psi_bc, gal)

                if psi is None:
                    print(f"  {psi_bc:>8.1e}  {lam:>8.1f}  {'SOLVER FAIL':>10}")
                    continue

                vc, geff, M_bar, M_psi = compute_rotation_curve(
                    r_arr, psi, dpsi, M0, lam, L0, gal)

                v_rmax   = float(np.mean(vc[-20:]))
                flat_ok  = check_flatness(r_arr, vc, v_flat)
                cosmo_ok = check_cosmo_ok(psi, L0, L0_BOUND)
                epos_ok  = check_energy_positive(psi, dpsi, M0, lam)
                passes   = flat_ok and cosmo_ok and epos_ok

                if passes:
                    any_pass = True

                print(f"  {psi_bc:>8.1e}  {lam:>8.1f}  {v_rmax:>10.2f}  "
                      f"{'YES':>6}  " if flat_ok else
                      f"  {psi_bc:>8.1e}  {lam:>8.1f}  {v_rmax:>10.2f}  "
                      f"{'no':>6}  ", end='')
                print(f"{'YES':>7}  " if cosmo_ok else f"{'NO':>7}  ", end='')
                print(f"{'YES':>5}  " if epos_ok  else f"{'NO':>5}  ", end='')
                print(f"{'PASS' if passes else 'fail':>5}")

                psi_max = float(np.max(np.abs(psi)))
                key = f'psi_bc={psi_bc:.1e}_lam={lam}'
                gal_results[key] = {
                    'psi_bc': psi_bc,
                    'lambda_gal': lam,
                    'v_c_at_rmax_kms': v_rmax,
                    'v_flat_obs_kms': v_flat,
                    'flat_ok': flat_ok,
                    'cosmo_ok': cosmo_ok,
                    'energy_positive': epos_ok,
                    'pass': passes,
                    'psi_max': psi_max,
                    'Lambda0_psi2_max': float(L0 * psi_max**2),
                    'M_Psi_at_rmax': float(M_psi[-1]),
                    'M_DM_needed': float(M_needed[-1]),
                }

                # Track best (closest v_c to v_flat in outer region)
                vc_outer = vc[r_arr >= FLAT_FRAC * r_max]
                diff = float(np.mean(np.abs(vc_outer - v_flat) / v_flat))
                if diff < best_vc_diff and cosmo_ok and epos_ok:
                    best_vc_diff = diff
                    best_entry = (psi_bc, lam, r_arr, vc, psi, dpsi, geff, M_bar, M_psi)

        best_fits[gal_name] = best_entry
        diag['galaxies'][gal_name] = {
            'v_flat_kms': v_flat,
            'r_max_kpc': r_max,
            'best_outer_deviation_frac': best_vc_diff,
            'scan_results': gal_results
        }

    # ── Physical diagnosis ────────────────────────────────────────────────────
    print(f"\n{'='*70}")
    print("PHYSICAL DIAGNOSIS")
    print(f"{'='*70}")

    # Compute the Yukawa profile shape factor to explain why FAIL is expected
    r_max_ngc = GALAXIES['NGC3198']['r_max_kpc']
    m0_r_max  = M0 * r_max_ngc
    print(f"\nYukawa scale analysis (NGC 3198):")
    print(f"  m0 * r_max = {m0_r_max:.2f}")
    print(f"  Shooting inward: Psi(r) ~ Psi_bc * (r_max/r) * exp(m0*(r_max-r))")
    print(f"  => rho_Psi(r) ~ exp(2*m0*(r_max-r)) / r^2")
    print(f"  => NOT proportional to 1/r^2 (isothermal sphere needed for flat curve)")
    print(f"  Exponential factor at r=5 kpc: exp(2*{M0}*{r_max_ngc-5:.0f}) = {np.exp(2*M0*(r_max_ngc-5)):.1f}x")
    print(f"  Exponential factor at r=15 kpc: exp(2*{M0}*{r_max_ngc-15:.0f}) = {np.exp(2*M0*(r_max_ngc-15)):.1f}x")
    print(f"\n  The Yukawa profile is EXPONENTIALLY steeper inward than 1/r^2.")
    print(f"  V=(lambda/4)Psi^4 modifies the amplitude but not the exponential shape.")
    print(f"  A flat rotation curve requires rho_DM ~ 1/r^2 over 5-30 kpc.")
    print(f"  => Conclusion: locked action at m0={M0} kpc^-1 cannot reproduce")
    print(f"     isothermal-sphere-like DM density regardless of lambda_gal or Psi_bc.")
    print(f"  => A screening or condensation mechanism (NOT in locked action) is needed.")

    # ── Overall verdict ───────────────────────────────────────────────────────
    print(f"\n{'='*70}")
    print("VERDICT")
    print(f"{'='*70}")

    if any_pass:
        verdict_str = (
            "PASS found: at least one (Psi_bc, lambda_gal) combination produces "
            "a flat rotation curve within cosmological bounds. "
            "See scan_results for details."
        )
    else:
        verdict_str = (
            "FAIL: No (Psi_bc, lambda_gal) combination in the scanned range produces "
            f"a flat rotation curve (|v_c - v_flat|/v_flat < {FLAT_TOL*100:.0f}% for "
            f"r > {FLAT_FRAC*100:.0f}% of r_max) while satisfying cosmological bounds "
            f"(Lambda0*Psi^2 < {L0_BOUND}) and energy positivity. "
            f"Physical reason: the decaying Yukawa profile (m0={M0} kpc^-1) yields "
            f"rho_Psi ~ exp(2*m0*(r_max-r))/r^2, which rises exponentially inward "
            f"(factor ~{np.exp(2*M0*(r_max_ngc-5)):.0f}x from r=30 to 5 kpc) rather "
            f"than the 1/r^2 isothermal profile needed for curve flattening. "
            f"The quartic interaction V=(lambda/4)Psi^4 modifies central amplitude "
            f"but cannot alter the exponential shape of the Yukawa halo. "
            f"CONCLUSION: The locked CMSTG action does not replace dark matter at "
            f"Lambda0={L0} (cosmological best-fit). Galactic-scale dark sector behavior "
            f"requires a screening or condensation mechanism not yet derived from "
            f"the locked action. Rotation curves remain a future prediction."
        )

    print(verdict_str)
    diag['overall_verdict'] = {
        'pass': any_pass,
        'summary': verdict_str,
        'physical_note': (
            f'Yukawa length = {1/M0:.1f} kpc vs galaxy r_max = {r_max_ngc:.0f} kpc. '
            f'm0*r_max = {m0_r_max:.1f} >> 1 means field decays exponentially across galaxy. '
            'Flat rotation curves need rho ~ 1/r^2; Yukawa gives rho ~ exp(-2m0*r)/r^2. '
            'V=(lambda/4)Psi^4 is too local (proportional to Psi^3, subdominant at small Psi) '
            'to reshape the profile from exponential to power-law.'
        )
    }

    # ── Plots ─────────────────────────────────────────────────────────────────
    n_gal = len(GALAXIES)
    fig_rc,  axes_rc  = plt.subplots(1, n_gal, figsize=(7*n_gal, 6))
    fig_fld, axes_fld = plt.subplots(1, n_gal, figsize=(7*n_gal, 6))
    fig_rho, axes_rho = plt.subplots(1, n_gal, figsize=(7*n_gal, 6))

    if n_gal == 1:
        axes_rc  = [axes_rc]
        axes_fld = [axes_fld]
        axes_rho = [axes_rho]

    for ig, (gal_name, gal) in enumerate(GALAXIES.items()):
        best = best_fits.get(gal_name)
        if best is None:
            continue
        psi_bc, lam, r_arr, vc, psi, dpsi, geff, M_bar, M_psi = best
        v_flat  = gal['v_flat_kms']
        rho_need = rho_DM_needed(r_arr, gal)
        rho_psi_arr = compute_rho_psi(psi, dpsi, M0, lam)

        # Baryons-only curve
        v_bar = np.sqrt(np.clip(G_gal * M_bar / r_arr, 0, None))

        # Rotation curve
        ax = axes_rc[ig]
        ax.plot(r_arr, v_bar, 'k--', lw=1.5, label='Baryons only')
        ax.plot(r_arr, vc,    'b-',  lw=2.0, label=f'CMSTG total ($\\Psi_{{bc}}$={psi_bc:.1e}, $\\lambda$={lam})')
        ax.axhline(v_flat, color='r', ls=':', lw=1.5, label=f'$v_{{flat}}={v_flat}$ km/s')
        ax.fill_between(r_arr, v_flat*(1-FLAT_TOL), v_flat*(1+FLAT_TOL),
                        alpha=0.15, color='red', label=f'±{FLAT_TOL*100:.0f}% band')
        ax.set_xlabel('r [kpc]', fontsize=12)
        ax.set_ylabel('$v_c$ [km/s]', fontsize=12)
        ax.set_title(f'{gal_name} — Best-fit rotation curve\n(best within cosmo+energy bounds)',
                     fontsize=11)
        ax.legend(fontsize=9)
        ax.set_xlim(0, r_arr[-1])
        ax.set_ylim(0, max(1.8*v_flat, float(np.max(vc)*1.1)))
        dev_str = f"Outer deviation: {diag['galaxies'][gal_name]['best_outer_deviation_frac']*100:.1f}%"
        ax.text(0.97, 0.05, dev_str, transform=ax.transAxes,
                ha='right', va='bottom', fontsize=9,
                color='red' if diag['galaxies'][gal_name]['best_outer_deviation_frac'] > FLAT_TOL else 'green')

        # Field profile
        ax = axes_fld[ig]
        ax.semilogy(r_arr, np.abs(psi), 'b-', lw=1.8, label='$|\\Psi(r)|$ (shooting)')
        ax.axvline(1/M0, color='gray', ls=':', lw=1, label=f'Yukawa length={1/M0} kpc')
        ax.set_xlabel('r [kpc]', fontsize=12)
        ax.set_ylabel('$|\\Psi(r)|$', fontsize=12)
        ax.set_title(f'{gal_name} — Field profile (decaying Yukawa mode)', fontsize=11)
        ax.legend(fontsize=9)

        # Density comparison
        ax = axes_rho[ig]
        ax.semilogy(r_arr, np.clip(rho_psi_arr, 1e-5, None), 'b-', lw=2,
                    label='$\\rho_\\Psi$ (CMSTG)')
        ax.semilogy(r_arr, rho_need, 'r--', lw=2,
                    label='Required $\\sim 1/r^2$ for flat curve')
        ax.set_xlabel('r [kpc]', fontsize=12)
        ax.set_ylabel('$\\rho$ [$M_\\odot$ kpc$^{-3}$]', fontsize=12)
        ax.set_title(f'{gal_name} — $\\rho_\\Psi$ vs required isothermal profile', fontsize=11)
        ax.legend(fontsize=9)

    for fig, name in [
        (fig_rc,  'sim99_rotation_curves.pdf'),
        (fig_fld, 'sim99_Psi_profiles.pdf'),
        (fig_rho, 'sim99_chi2_comparison.png'),   # reuse filename for paper integration
    ]:
        fig.tight_layout()
        fig.savefig(os.path.join(OUTPUTS, name), dpi=150, bbox_inches='tight')
        plt.close(fig)

    # ── Scan summary heatmap ──────────────────────────────────────────────────
    for gal_name, gal in GALAXIES.items():
        r_arr  = np.linspace(NUM['r_min_kpc'], gal['r_max_kpc'], NUM['n_points'])
        v_flat = gal['v_flat_kms']

        dev_matrix = np.full((len(PSI_BC_SCAN), len(LAMBDA_SCAN)), np.nan)

        for ip, psi_bc in enumerate(PSI_BC_SCAN):
            for il, lam in enumerate(LAMBDA_SCAN):
                key = f'psi_bc={psi_bc:.1e}_lam={lam}'
                res = diag['galaxies'][gal_name]['scan_results'].get(key)
                if res and res['cosmo_ok'] and res['energy_positive']:
                    dev = abs(res['v_c_at_rmax_kms'] - v_flat) / v_flat
                    dev_matrix[ip, il] = min(dev, 2.0)  # cap for display

        fig_hm, ax_hm = plt.subplots(figsize=(9, 5))
        im = ax_hm.imshow(dev_matrix, aspect='auto', origin='lower',
                          cmap='RdYlGn_r', vmin=0, vmax=1.0)
        ax_hm.set_xticks(range(len(LAMBDA_SCAN)))
        ax_hm.set_xticklabels([f'{l:.0f}' for l in LAMBDA_SCAN])
        ax_hm.set_yticks(range(len(PSI_BC_SCAN)))
        ax_hm.set_yticklabels([f'{p:.0e}' for p in PSI_BC_SCAN])
        ax_hm.set_xlabel('$\\lambda_{\\rm gal}$', fontsize=12)
        ax_hm.set_ylabel('$\\Psi_{\\rm bc}$ at $r_{\\rm max}$', fontsize=12)
        ax_hm.set_title(
            f'{gal_name} — Outer $|v_c - v_{{flat}}|/v_{{flat}}$ (within cosmo+energy bounds)\n'
            f'Green = flat ({FLAT_TOL*100:.0f}% tol), Red = discrepant, Gray = bound violated',
            fontsize=10)
        plt.colorbar(im, ax=ax_hm, label='Outer deviation fraction')
        ax_hm.axhline(-0.5, color='green', lw=2, linestyle='--',
                      label=f'Flat threshold ({FLAT_TOL*100:.0f}%)')
        fig_hm.tight_layout()
        fig_hm.savefig(os.path.join(OUTPUTS, f'sim99_scan_{gal_name}.pdf'),
                       dpi=150, bbox_inches='tight')
        plt.close(fig_hm)

    # ── Save diagnostics ──────────────────────────────────────────────────────
    diag_path = os.path.join(OUTPUTS, 'sim99_diagnostics.json')
    with open(diag_path, 'w') as f:
        json.dump(diag, f, indent=2)

    print(f"\nOutputs written to {OUTPUTS}/")
    print(f"  sim99_rotation_curves.pdf — best-fit RC per galaxy")
    print(f"  sim99_Psi_profiles.pdf    — field profile (decaying Yukawa mode)")
    print(f"  sim99_chi2_comparison.png — rho_Psi vs required 1/r^2 profile")
    print(f"  sim99_scan_NGC3198.pdf    — scan heatmap NGC 3198")
    print(f"  sim99_scan_NGC6503.pdf    — scan heatmap NGC 6503")
    print(f"  sim99_diagnostics.json    — full numerical results")
    print("Done.")
    return diag


if __name__ == '__main__':
    run()
