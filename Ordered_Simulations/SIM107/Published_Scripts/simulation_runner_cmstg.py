#!/usr/bin/env python3
"""
SIM107: Effective Matter-Psi Coupling and Solar System Viability
CMSTG: Curvature-Memory Scalar-Tensor Gravity

Derives the Psi-matter coupling from the linearised CMSTG field equation,
computes the resulting fifth-force strength, and tests against Solar System
constraints (Cassini PPN bound, Lunar Laser Ranging).

Key physics:
  - Linearise Psi equation around background Psi0 in presence of matter source
  - Effective Psi mass: m_eff^2(rho) = m0^2 + 2*Lambda0*|R| ~ m0^2 + 2*Lambda0*rho
    (R ~ rho in units 8piG=1, chameleon direction: heavier in dense regions)
  - Yukawa fifth-force potential V5(r) = -alpha^2 * G*M/r * exp(-m_eff*r)
  - PPN parameter: gamma - 1 = -2*alpha^2/(1+alpha^2)
  - Thin-shell chameleon screening inside dense bodies
  - Cassini bound: |gamma-1| < 2.3e-5

Units: H0=1, 8piG=1 (=> G_N=1/(8pi), M_Pl=1/sqrt(8pi)).
"""

import numpy as np
from scipy.integrate import solve_ivp, quad
from scipy.optimize import brentq
import json, os

# ─── CMSTG Parameters (observationally constrained) ────────────────────────────
Lambda0 = 0.003       # scalar-gravity coupling (dimensionless)
m0      = 0.001       # bare Psi mass [H0 units]
H0      = 1.0         # Hubble constant (sets units)
G_N     = 1.0 / (8.0 * np.pi)   # from 8piG=1

# Physical conversion factors
hubble_Mpc  = 4282.0        # 1/H0 in Mpc
rho_crit_SI = 9.47e-27      # kg/m^3 (Planck 2018)
rho_crit_today = 3.0 * H0**2   # = 3.0 in simulation units

def rho_SI_to_sim(rho_kg_m3):
    """Convert SI density to simulation units (rho_crit_today=3)."""
    return (rho_kg_m3 / rho_crit_SI) * rho_crit_today

# Cassini PPN bound (Bertotti et al. 2003)
cassini_gamma_bound = 2.3e-5     # |gamma - 1|
alpha2_cassini      = cassini_gamma_bound / 2.0   # alpha^2 for small alpha
# LLR bound (Williams et al.)
alpha2_llr = 1.0e-4

# ─── Environment Densities ────────────────────────────────────────────────────
environments = {
    'Deep void'          : rho_SI_to_sim(1.0e-30),
    'Cosmic mean (z=0)'  : 0.315 * rho_crit_today,
    'Galaxy cluster'     : rho_SI_to_sim(1.0e-24),
    'Milky Way disk'     : rho_SI_to_sim(1.0e-21),
    'Interstellar (SS)'  : rho_SI_to_sim(1.0e-22),   # Solar neighbourhood
    'Solar corona'       : rho_SI_to_sim(1.0e-12),
    'Solar interior'     : rho_SI_to_sim(1.4e+03),   # ~1400 kg/m^3 mean
    'Earth surface'      : rho_SI_to_sim(5.5e+03),
}

# ─── 1. Effective Mass ─────────────────────────────────────────────────────────
def m_eff_squared(rho):
    """
    m_eff^2(rho) = m0^2 + 2*Lambda0*|R|
    In units 8piG=1: trace of perturbed Einstein gives |R| = rho for pressureless
    matter (taking the magnitude; the chameleon direction adds mass in dense regions).
    """
    return m0**2 + 2.0 * Lambda0 * rho

def m_eff(rho):
    ms2 = m_eff_squared(rho)
    return np.sqrt(ms2) if ms2 >= 0 else None

def compton_wavelength_Mpc(rho):
    """Compton wavelength lambda = 1/m_eff in Mpc."""
    meff = m_eff(rho)
    return hubble_Mpc / meff if meff and meff > 0 else np.inf

# ─── 2. Scalar Coupling Strength (from linearised CMSTG action) ────────────────
def alpha_coupling(Psi0):
    """
    Effective dimensionless scalar-matter coupling alpha for PPN.

    From varying CMSTG action: the Psi-matter vertex comes from
    delta(Lambda0*Psi^2*R)/delta(Psi) = 2*Lambda0*Psi*R at background Psi=Psi0.

    In Brans-Dicke / scalar-tensor language, the conformal coupling is
    A^2(Psi) = 1 + 2*Lambda0*Psi^2  (using 8piG=1, so 16piG=2)
    alpha = d(ln A)/d(Psi)|_{Psi0} = Lambda0*Psi0 / (1 + 2*Lambda0*Psi0^2)
    """
    return Lambda0 * Psi0 / (1.0 + 2.0 * Lambda0 * Psi0**2)

def alpha_squared(Psi0):
    a = alpha_coupling(Psi0)
    return a**2

def gamma_ppn(Psi0):
    """PPN gamma = (1 - alpha^2)/(1 + alpha^2)."""
    a2 = alpha_squared(Psi0)
    return (1.0 - a2) / (1.0 + a2)

def gamma_minus_1(Psi0):
    return gamma_ppn(Psi0) - 1.0

# ─── 3. Fifth-Force Ratio (with Yukawa suppression) ──────────────────────────
def fifth_force_ratio(r_sim, Psi0, rho_env=0.0):
    """
    |F_fifth(r)| / |F_grav(r)|

    For Yukawa force F5/Fg = 2*alpha^2 * exp(-m_eff*r) * (1 + m_eff*r)
    The factor (1+m_eff*r) accounts for the gradient of the Yukawa potential.
    For m_eff*r << 1 this -> 2*alpha^2 (long-range limit).
    """
    a2   = alpha_squared(Psi0)
    meff = m_eff(rho_env) or 0.0
    x    = meff * r_sim
    return 2.0 * a2 * np.exp(-x) * (1.0 + x)

# ─── 4. Thin-Shell Chameleon Suppression ─────────────────────────────────────
def thin_shell_fraction(M_sun_sim, R_sun_sim, rho_sun_interior, rho_exterior):
    """
    Khoury-Weltman thin-shell parameter chi = 3*Phi_N * (m_eff_out^2 / alpha^2)^(-1)
    where Phi_N = G_N*M/R is the Newtonian potential.

    If chi << 1: strong chameleon suppression, effective alpha_eff = alpha * chi
    If chi >= 1: no suppression (full fifth force)
    """
    # Gravitational potential at surface (dimensionless)
    Phi_N = G_N * M_sun_sim / R_sun_sim

    # Mass scales inside vs outside
    meff_in_sq  = m_eff_squared(rho_sun_interior)
    meff_out_sq = m_eff_squared(rho_exterior)

    delta_m2 = meff_in_sq - meff_out_sq
    if delta_m2 <= 0:
        return 1.0   # no thin shell possible

    # chi = 6 * alpha * Phi_N / (delta_phi_env)
    # Simplified: chi ~ 3 * Phi_N / (Lambda0 * Psi0 * delta_Psi_min / M_Pl)
    # For the chameleon: delta_Psi_min ~ (m_eff_out - m_eff_in) / (Lambda0)
    # Standard result: chi = (m_eff_out^2 - m_eff_in^2)^{-1} * 6*Phi_N * m_eff_in^2
    # See Khoury & Weltman 2004 Eq.(39)
    chi = 6.0 * Phi_N * meff_in_sq / delta_m2
    return min(chi, 1.0)

# ─── 5. Critical Psi0 from Cassini ────────────────────────────────────────────
def find_Psi0_critical():
    """Find largest Psi0 satisfying |gamma-1| < cassini_gamma_bound."""
    # gamma-1 is negative and |gamma-1| is increasing with Psi0
    # Find root of |gamma_minus_1(Psi0)| - cassini_gamma_bound = 0
    f = lambda p: abs(gamma_minus_1(p)) - cassini_gamma_bound
    try:
        # Check brackets
        if f(0.001) > 0:
            return 0.0    # already violated at tiny Psi0
        if f(100.0) < 0:
            return 100.0  # always satisfied
        Psi0_crit = brentq(f, 0.001, 100.0, xtol=1e-6)
        return Psi0_crit
    except ValueError:
        return None

# ─── MAIN ─────────────────────────────────────────────────────────────────────
print("="*70)
print("SIM107: CMSTG Matter-Psi Coupling and Solar System Viability")
print("="*70)
print(f"Lambda0 = {Lambda0},  m0 = {m0} H0,  G_N = {G_N:.5f},  8piG = 1")
print(f"Cassini bound: |gamma-1| < {cassini_gamma_bound:.1e}\n")

# ── Table 1: Effective mass environment scan ───────────────────────────────
print("─── Table 1: Effective Psi Mass Across Environments ───")
header = f"{'Environment':<22} {'rho_sim':>11} {'m_eff/H0':>12} {'lambda_Psi [Mpc]':>18} {'Chameleon active':>17}"
print(header)
print("-" * 83)
for name, rho in environments.items():
    meff    = m_eff(rho)
    lam_Mpc = compton_wavelength_Mpc(rho)
    dm2     = 2.0 * Lambda0 * rho
    active  = "YES" if dm2 > m0**2 else "no"
    meff_s  = f"{meff:.4e}" if meff else "tachyonic"
    lam_s   = f"{lam_Mpc:.3e}" if lam_Mpc < 1e14 else ">Hubble"
    print(f"{name:<22} {rho:>11.3e} {meff_s:>12} {lam_s:>18} {active:>17}")

print()

# ── Table 2: Coupling strength vs Psi0 ────────────────────────────────────
print("─── Table 2: Fifth-Force Coupling vs Background Scalar Value ───")
header2 = f"{'Psi_0':>8} {'alpha':>12} {'alpha^2':>12} {'|gamma-1|':>12} {'Cassini':>9} {'LLR':>9}"
print(header2)
print("-" * 67)
Psi0_scan = [0.001, 0.01, 0.1, 0.5, 1.0, 1.16, 2.0, 5.0, 10.0]
for p in Psi0_scan:
    a  = alpha_coupling(p)
    a2 = a**2
    g1 = abs(gamma_minus_1(p))
    c  = "PASS" if g1 < cassini_gamma_bound else "FAIL"
    l  = "PASS" if 2.0*a2 < alpha2_llr else "FAIL"
    print(f"{p:>8.3f} {a:>12.4e} {a2:>12.4e} {g1:>12.4e} {c:>9} {l:>9}")

print()

# ── Table 3: Fifth-force ratio at 1 AU ────────────────────────────────────
r_1AU_Mpc = 4.848e-9   # 1 AU in Mpc
r_1AU_sim = r_1AU_Mpc / hubble_Mpc
meff_vac   = m_eff(0.0)   # vacuum (interplanetary space)
x_1AU      = meff_vac * r_1AU_sim   # m_eff * r

print("─── Table 3: Fifth-Force / Gravity Ratio at 1 AU (Cassini geometry) ───")
print(f"  r_1AU = {r_1AU_Mpc:.3e} Mpc | m0*r_1AU = {x_1AU:.3e} (Yukawa suppression factor)")
header3 = f"{'Psi_0':>8} {'F5/Fg (no shell)':>18} {'F5/Fg (thin shell Sun)':>24} {'Cassini':>9}"
print(header3)
print("-" * 65)

# Sun parameters in simulation units
M_sun_SI  = 1.989e30      # kg
M_sun_sim = G_N * M_sun_SI / rho_crit_SI / hubble_Mpc**3  # rough conversion
# Actually: G_N*M_sun in length units...
# Use: G_N*M_sun ~ 1.5 km ~ 1.5e-17 Mpc / hubble_Mpc -> sim units
G_M_sun_Mpc = 1.476e-3 / 1e3 / 1e6   # G*M_sun ~ 1476 m, convert to Mpc: 1.476e3 / 3.086e22
G_M_sun_sim = G_M_sun_Mpc / hubble_Mpc
R_sun_Mpc   = 6.957e8 / 3.086e22   # R_sun in Mpc
R_sun_sim   = R_sun_Mpc / hubble_Mpc
rho_sun_sim = rho_SI_to_sim(1.4e3)  # mean solar density ~1400 kg/m^3
rho_ISM_sim = rho_SI_to_sim(1e-22)  # interplanetary medium

for p in [0.1, 0.5, 1.0, 1.16, 2.0, 5.0]:
    ratio_naive = fifth_force_ratio(r_1AU_sim, p, rho_env=0.0)
    chi   = thin_shell_fraction(G_M_sun_sim / G_N, R_sun_sim, rho_sun_sim, rho_ISM_sim)
    ratio_shell = ratio_naive * chi**2   # thin-shell reduces coupling by chi
    g1    = abs(gamma_minus_1(p)) * chi**2
    c     = "PASS" if g1 < cassini_gamma_bound else "FAIL"
    print(f"{p:>8.3f} {ratio_naive:>18.4e} {ratio_shell:>24.4e} {c:>9}")

print()

# ── Critical Psi0 ──────────────────────────────────────────────────────────
Psi0_crit = find_Psi0_critical()
print("─── Critical Constraint on Psi_0 ───")
if Psi0_crit is not None:
    a2_crit = alpha_squared(Psi0_crit)
    g1_crit = abs(gamma_minus_1(Psi0_crit))
    print(f"  Cassini requires: |Psi_0| < {Psi0_crit:.4f} (Planck units, 8piG=1)")
    print(f"  At Psi0_crit: alpha^2 = {a2_crit:.4e}, |gamma-1| = {g1_crit:.4e}")
else:
    print("  Constraint: always satisfied (alpha too small)")

print()

# ── Chameleon screening analysis ───────────────────────────────────────────
print("─── Chameleon Screening Detail ───")
rho_SS   = rho_SI_to_sim(1e-22)
delta_m2 = 2.0 * Lambda0 * rho_SS
frac_SS  = delta_m2 / m0**2
print(f"  Solar System (ISM) rho = {rho_SS:.3e}  ->  Delta(m_eff^2)/m0^2 = {frac_SS:.3e}")
print(f"  Compton wavelength (vacuum): {compton_wavelength_Mpc(0.0):.3e} Mpc  >> Solar System")
print(f"  Compton wavelength (ISM):    {compton_wavelength_Mpc(rho_SS):.3e} Mpc")
rho_sol  = rho_SI_to_sim(1.4e3)
print(f"  Compton wavelength (Sun):    {compton_wavelength_Mpc(rho_sol):.3e} Mpc  (screened inside!)")
chi_sun  = thin_shell_fraction(G_M_sun_sim / G_N, R_sun_sim, rho_sun_sim, rho_ISM_sim)
print(f"  Thin-shell fraction chi_Sun = {chi_sun:.4e}")
print(f"  Coupling suppression factor chi^2 = {chi_sun**2:.4e}")
print()

# ── RG-flow implication: Psi0 today ────────────────────────────────────────
print("─── Cosmological Constraint on Psi_0 ───")
print(f"  For m0 >> H0: Psi decays as a^(-3/2) (matter-like oscillations)")
print(f"  For m0 ~ H0 ({m0} H0): slow roll; Psi_0 ~ O(1) today")
print(f"  Cassini (no thin shell): |Psi_0| < {Psi0_crit:.3f} M_Pl")
print(f"  Cassini (with Sun thin-shell chi={chi_sun:.2e}): constraint relaxed by 1/chi")
print()

# ── SUMMARY ─────────────────────────────────────────────────────────────────
print("="*70)
print("SUMMARY")
print("="*70)
print(f"1. Chameleon mechanism: ACTIVE inside dense bodies (lambda_Sun << AU)")
print(f"   but INACTIVE in interplanetary space (lambda >> Solar System)")
print()
print(f"2. Naive fifth-force coupling (no thin shell):")
print(f"   alpha^2(Psi0=1) = {alpha_squared(1.0):.4e}")
print(f"   |gamma-1|(Psi0=1) = {abs(gamma_minus_1(1.0)):.4e}  [Cassini: {cassini_gamma_bound:.1e}]")
status_naive = "PASS" if abs(gamma_minus_1(1.0)) < cassini_gamma_bound else "FAIL"
print(f"   Status: {status_naive}")
print()
print(f"3. Thin-shell suppression from solar interior:")
print(f"   chi_Sun = {chi_sun:.4e}  ->  alpha_eff^2 = alpha^2 * chi^2")
a2_eff = alpha_squared(1.0) * chi_sun**2
g1_eff = abs(gamma_minus_1(1.0)) * chi_sun**2
print(f"   Effective |gamma-1|(Psi0=1) = {g1_eff:.4e}  [Cassini: {cassini_gamma_bound:.1e}]")
status_shell = "PASS" if g1_eff < cassini_gamma_bound else "FAIL"
print(f"   Status: {status_shell}")
print()
print(f"4. Critical Psi0 (no thin shell): {Psi0_crit:.4f} M_Pl")
print()
overall = "PASS" if status_naive == "PASS" or status_shell == "PASS" else "FAIL"
print(f"OVERALL VERDICT: {overall}")
if overall == "PASS":
    print("CMSTG satisfies Solar System fifth-force constraints.")
    if status_naive == "PASS":
        print(f"Reason: alpha^2 = O(Lambda0^2*Psi0^2) is small for Psi0 < {Psi0_crit:.2f}.")
    else:
        print(f"Reason: thin-shell chameleon screening inside the Sun.")
    print(f"The locked action requires no additional screening mechanism.")
else:
    print("CMSTG violates Solar System constraints.")
    print("Requires either Psi0 << 1 or an additional screening mechanism.")

# ── Save diagnostics ─────────────────────────────────────────────────────────
out_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'Outputs')
os.makedirs(out_dir, exist_ok=True)

diag = {
    'Lambda0'               : Lambda0,
    'm0'                    : m0,
    'Psi0_critical_Cassini' : float(Psi0_crit) if Psi0_crit else None,
    'alpha2_at_Psi0_1'      : float(alpha_squared(1.0)),
    'gamma_minus_1_Psi0_1'  : float(abs(gamma_minus_1(1.0))),
    'cassini_bound'         : cassini_gamma_bound,
    'status_naive'          : status_naive,
    'chi_sun'               : float(chi_sun),
    'status_thin_shell'     : status_shell,
    'overall_verdict'       : overall,
    'm_eff_environments'    : {
        k: {
            'rho_sim'       : float(v),
            'm_eff'         : float(m_eff(v)) if m_eff(v) else None,
            'lambda_Mpc'    : float(compton_wavelength_Mpc(v))
                               if compton_wavelength_Mpc(v) < 1e14 else None,
        }
        for k, v in environments.items()
    },
    'fifth_force_range_Mpc' : float(hubble_Mpc / m0),
}

with open(os.path.join(out_dir, 'sim107_diagnostics.json'), 'w') as f:
    json.dump(diag, f, indent=2)

print(f"\nDiagnostics saved to Outputs/sim107_diagnostics.json")
