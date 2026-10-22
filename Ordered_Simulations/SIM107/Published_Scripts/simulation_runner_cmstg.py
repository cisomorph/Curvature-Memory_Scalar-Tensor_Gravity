#!/usr/bin/env python3
"""
SIM107: Effective Matter-Psi Coupling and Solar System Viability
CMSTG: Curvature-Memory Scalar-Tensor Gravity

Derives the Psi-matter coupling from the CMSTG action, computes the resulting
fifth-force strength, and tests against Solar System constraints (Cassini PPN
bound, Lunar Laser Ranging) at the Phase 1 canonical field value.

Action (Jordan frame, signature -+++, units 8piG = 1 so GR is F = 1/2):
  S = ∫d⁴x √-g [F(Ψ) R − ½(∂Ψ)² − ½m0²Ψ²] + S_m,   F(Ψ) = ½ + Λ0 Ψ²

Key physics (derived):
  - DEF form: (1/16πG*)∫√-g[F̃R − Z(∂Ψ)²] with F̃ = 2F = 1 + 2Λ0Ψ², Z = 1.
    Einstein-frame coupling  α² = F̃'² / (3F̃'² + 2ZF̃),   F̃' = 4Λ0Ψ
    PPN                      γ − 1 = −2α²/(1+α²) = −F̃'²/(ZF̃ + 2F̃'²)
    (Brans–Dicke check: F̃ = Φ, Z = ω/Φ gives γ − 1 = −1/(ω+2).)
  - Matter coupling: substituting the trace equation into the scalar EOM gives
    □Ψ(1 + 3F'²/F) = V' − (F'/2F)ρ + ...,  i.e.  V_eff = V − (ρ/2) ln F.
    Effective mass  m_eff² = [m0² − (ρ/2)(F''/F − F'²/F²)] / (1 + 3F'²/F).
    The matter term is NEGATIVE for this F (anti-chameleon).
  - Screening needs a heavy field inside the body (m_in R ≫ 1) and a thin shell
    (Khoury–Weltman ΔR/R = |ΔΨ|/(6βΦ_N) ≪ 1, β = α/√2 for canonical Ψ).
  - Cassini bound: |γ − 1| < 2.3e-5.

AUDIT 2026-10-22 (see AUDIT_NOTES.md). Changes from the published version:
  1. The published run evaluated γ only at Ψ0 = 1.0; its own critical value was
     Ψ0 = 1.139, below the canonical Ψ̄ = 2.62. The verdict is now taken at Ψ̄ = 2.62.
  2. Coupling formula: published α = Λ0Ψ/(1+2Λ0Ψ²) (from "A² = 1+2Λ0Ψ²") was off by
     2 and lacked the canonical normalisation. Replaced by the DEF result above.
  3. m_eff² = m0² + 2Λ0|R| had its sign chosen by hand ("taking the magnitude").
     Replaced by the derived expression (negative matter contribution).
  4. Unit bugs: G·M_sun was 1.476e-12 Mpc (correct: 4.78e-20 Mpc), which made
     Φ_N(Sun) = 65 instead of 2.1e-6; 1 AU was 4.848e-9 Mpc (correct: 4.848e-12).
  5. The ad-hoc thin-shell factor is replaced by the two Khoury–Weltman
     conditions (m_in R ≫ 1 and ΔR/R ≪ 1).

Units: H0 = 1, 8piG = 1 (reduced Planck mass = 1); ρ_crit,0 = 3 H0².
"""

import numpy as np
from scipy.optimize import brentq
import json, os

# ─── CMSTG Parameters ──────────────────────────────────────────────────────────
Lambda0 = 0.003       # scalar-gravity coupling (locked)
PSI_BAR = 2.62        # Phase 1 canonical background value [reduced M_Pl]
m0      = 0.001       # bare Psi mass [H0 units]

# Physical conversion factors
MPC_M        = 3.0857e22     # metres per Mpc
hubble_Mpc   = 4282.0        # c/H0 in Mpc (h ≈ 0.70)
rho_crit_SI  = 9.47e-27      # kg/m^3
rho_crit_sim = 3.0           # ρ_crit,0 in sim units

def rho_SI_to_sim(rho_kg_m3):
    return (rho_kg_m3 / rho_crit_SI) * rho_crit_sim

# Bounds
cassini_gamma_bound = 2.3e-5         # |γ − 1| (Bertotti, Iess & Tortora 2003)
alpha2_llr          = 1.0e-4         # rough LLR bound on 2α² (as in the published sim)

# Sun
GM_sun_m   = 1476.6                  # G M_sun / c² in metres
R_sun_m    = 6.957e8
Phi_sun    = GM_sun_m / R_sun_m      # Newtonian potential at the surface (~2.1e-6)
rho_sun_SI = 1.41e3                  # mean density
AU_m       = 1.496e11

environments = {
    'Deep void'          : rho_SI_to_sim(1.0e-30),
    'Cosmic mean (z=0)'  : 0.315 * rho_crit_sim,
    'Galaxy cluster'     : rho_SI_to_sim(1.0e-24),
    'Milky Way disk'     : rho_SI_to_sim(1.0e-21),
    'Interplanetary 1 AU': rho_SI_to_sim(1.0e-20),
    'Solar corona'       : rho_SI_to_sim(1.0e-12),
    'Solar interior'     : rho_SI_to_sim(rho_sun_SI),
    'Earth surface'      : rho_SI_to_sim(5.5e+03),
}

# ─── Coupling functions ────────────────────────────────────────────────────────
def F(Psi):    return 0.5 + Lambda0 * Psi**2
def Fp(Psi):   return 2.0 * Lambda0 * Psi
def Fpp(Psi):  return 2.0 * Lambda0

def alpha_squared(Psi0):
    """Einstein-frame coupling α² (DEF normalisation), F̃ = 2F, Z = 1."""
    Ft, Ftp = 2.0 * F(Psi0), 2.0 * Fp(Psi0)
    return Ftp**2 / (3.0 * Ftp**2 + 2.0 * Ft)

def gamma_minus_1(Psi0):
    a2 = alpha_squared(Psi0)
    return -2.0 * a2 / (1.0 + a2)

def alpha_squared_published(Psi0):
    """The pre-audit formula, kept for the record."""
    return (Lambda0 * Psi0 / (1.0 + 2.0 * Lambda0 * Psi0**2))**2

def m_eff_squared(rho, Psi0=PSI_BAR):
    """Derived: [m0² − (ρ/2)(F''/F − F'²/F²)] / (1 + 3F'²/F)."""
    f, fp, fpp = F(Psi0), Fp(Psi0), Fpp(Psi0)
    return (m0**2 - 0.5 * rho * (fpp / f - fp**2 / f**2)) / (1.0 + 3.0 * fp**2 / f)

def m_eff_squared_published(rho):
    return m0**2 + 2.0 * Lambda0 * rho

def range_m(m2):
    """Compton wavelength 1/sqrt|m²| in metres (H0 units → Mpc → m)."""
    return hubble_Mpc * MPC_M / np.sqrt(abs(m2))

def find_Psi0_critical():
    f = lambda p: abs(gamma_minus_1(p)) - cassini_gamma_bound
    return brentq(f, 1e-4, 100.0, xtol=1e-10)

# ─── MAIN ─────────────────────────────────────────────────────────────────────
print("=" * 72)
print("SIM107: CMSTG Matter-Psi Coupling and Solar System Viability  [audited]")
print("=" * 72)
print(f"Lambda0 = {Lambda0}, Psi_bar = {PSI_BAR}, m0 = {m0} H0, F(Psi_bar) = {F(PSI_BAR):.5f}")
print(f"Cassini bound: |gamma-1| < {cassini_gamma_bound:.1e}\n")

# ── Table 1: coupling vs Ψ0 ────────────────────────────────────────────────
print("─── Table 1: Fifth-force coupling vs background scalar value ───")
print(f"{'Psi_0':>7} {'alpha^2':>11} {'|gamma-1|':>11} {'x Cassini':>10} {'Cassini':>8} "
      f"{'LLR':>6}   {'published |gamma-1|':>20}")
scan = [0.1, 0.4, 0.5, 1.0, 1.139, 2.0, PSI_BAR, 2.9, 5.0]
table1 = []
for p in scan:
    a2 = alpha_squared(p); g1 = abs(gamma_minus_1(p))
    a2p = alpha_squared_published(p); g1p = 2 * a2p / (1 + a2p)
    c = "PASS" if g1 < cassini_gamma_bound else "FAIL"
    l = "PASS" if 2 * a2 < alpha2_llr else "FAIL"
    tag = "  <- canonical" if p == PSI_BAR else ""
    print(f"{p:>7.3f} {a2:>11.3e} {g1:>11.3e} {g1/cassini_gamma_bound:>10.1f} {c:>8} {l:>6}   {g1p:>20.3e}{tag}")
    table1.append(dict(Psi0=p, alpha2=a2, gamma_minus_1=g1, cassini=c, llr=l,
                       published_formula_gamma_minus_1=g1p))
print()

# ── Table 2: effective mass across environments (at Ψ̄) ─────────────────────
print("─── Table 2: Effective Psi mass at Psi_bar (derived) vs published formula ───")
print(f"{'Environment':<21} {'rho_sim':>10} {'m_eff^2 [H0^2]':>15} {'range [m]':>11} "
      f"{'published m^2':>14}")
table2 = {}
for name, rho in environments.items():
    m2 = m_eff_squared(rho); m2p = m_eff_squared_published(rho)
    print(f"{name:<21} {rho:>10.3e} {m2:>15.4e} {range_m(m2):>11.3e} {m2p:>14.4e}")
    table2[name] = dict(rho_sim=rho, m_eff2_derived=m2, range_m=range_m(m2),
                        tachyonic=bool(m2 < 0), m_eff2_published=m2p)
print("  (m_eff² < 0: the matter coupling destabilises rather than screens; a tachyonic")
print("   mode only grows on scales longer than its range, far larger than any of these bodies.)\n")

# ── Screening analysis (Khoury–Weltman) ─────────────────────────────────────
a2_bar  = alpha_squared(PSI_BAR)
beta_kw = np.sqrt(a2_bar / 2.0)                     # canonical-field coupling
rho_sun = rho_SI_to_sim(rho_sun_SI)
m_in_R  = R_sun_m / range_m(m_eff_squared(rho_sun))
m_in_R_pub = R_sun_m / range_m(m_eff_squared_published(rho_sun))
dPsi    = PSI_BAR                                    # Δψ ~ Ψ̄ − Ψ_inside, Ψ_inside ∈ [0, Ψ̄]
thin_shell = dPsi / (6.0 * beta_kw * Phi_sun)
rho_1AU = rho_SI_to_sim(1.0e-20)
range_1AU = range_m(m_eff_squared(rho_1AU))
yukawa_1AU = np.exp(-AU_m / range_1AU)
print("─── Screening (Sun) ───")
print(f"  Phi_N(Sun) = {Phi_sun:.3e}   (published code used 65.5 from a unit error)")
print(f"  m_in * R_sun = {m_in_R:.3e} derived, {m_in_R_pub:.3e} with the published sign  (screening needs >> 1)")
print(f"  KW thin-shell parameter ΔR/R = ΔΨ/(6βΦ_N) = {thin_shell:.3e}  (screening needs << 1)")
print(f"  Range at 1 AU = {range_1AU:.3e} m  → Yukawa factor over 1 AU = {yukawa_1AU:.6f}")
screened = (m_in_R > 1.0) and (thin_shell < 1.0)
print(f"  Screened: {screened}\n")

# ── Critical Ψ0 and implied F0 ───────────────────────────────────────────────
Psi0_crit = find_Psi0_critical()
F0_max = F(Psi0_crit)
g1_bar = abs(gamma_minus_1(PSI_BAR))
print("─── Cassini constraint ───")
print(f"  |gamma-1|(Psi_bar = {PSI_BAR}) = {g1_bar:.3e}  = {g1_bar/cassini_gamma_bound:.1f} x bound")
print(f"  Cassini requires |Psi_0| < {Psi0_crit:.4f}  →  F0 = 1/2 + Λ0Ψ0² < {F0_max:.6f}"
      f"  (|F0 − 1/2| < {F0_max-0.5:.2e})")
print()

status = "PASS" if (g1_bar < cassini_gamma_bound or screened) else "FAIL"
print("=" * 72)
print(f"OVERALL VERDICT (canonical Psi_bar = {PSI_BAR}): {status}")
if status == "FAIL":
    print("  The unscreened Ψ-mediated force violates Cassini by "
          f"{g1_bar/cassini_gamma_bound:.0f}x. No screening operates: the field is")
    print("  effectively massless on Solar System scales and the Sun has no thin shell.")
    print(f"  Viability needs |Ψ0| < {Psi0_crit:.2f} today (F0 within {F0_max-0.5:.1e} of GR) or a")
    print("  screening mechanism that is not in the locked action.")
print("=" * 72)

# ── Save diagnostics ─────────────────────────────────────────────────────────
out_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'Outputs')
os.makedirs(out_dir, exist_ok=True)
diag = {
    'Lambda0': Lambda0, 'Psi_bar': PSI_BAR, 'm0': m0,
    'alpha2_at_Psi_bar': a2_bar,
    'gamma_minus_1_at_Psi_bar': g1_bar,
    'cassini_bound': cassini_gamma_bound,
    'excess_over_cassini': g1_bar / cassini_gamma_bound,
    'Psi0_critical_Cassini': Psi0_crit,
    'F0_max_Cassini': F0_max,
    'screening': dict(Phi_sun=Phi_sun, m_in_R_sun=m_in_R, m_in_R_sun_published_sign=m_in_R_pub,
                      thin_shell_dR_over_R=thin_shell, range_at_1AU_m=range_1AU,
                      yukawa_factor_1AU=yukawa_1AU, screened=bool(screened)),
    'coupling_scan': table1,
    'environments': table2,
    'overall_verdict': status,
    'published_pre_audit': {'overall_verdict': 'PASS', 'evaluated_at_Psi0': 1.0,
                            'Psi0_critical_Cassini': 1.1391967488087065,
                            'gamma_minus_1_Psi0_1': 1.778577039546203e-05, 'chi_sun': 1.0},
    'audit_note': 'Audited 2026-10-22: verdict at canonical Psi_bar; DEF coupling formula; '
                  'derived m_eff sign; unit fixes (GM_sun, AU); Khoury-Weltman screening test.',
}
with open(os.path.join(out_dir, 'sim107_diagnostics.json'), 'w') as f:
    json.dump(diag, f, indent=2, default=float)
print("\nDiagnostics saved to Outputs/sim107_diagnostics.json")
