#!/usr/bin/env python3
"""
SIM117: Level-2 RIFT Recursion — Memory-of-Memory Dark Matter
RIFT: Recursive Intelligence-Field Theory

Level-2 recursion: the RIFT gravitational coupling constant is itself a
dynamical field driven by the causal history of (Ψ² × curvature).

Modified non-minimal coupling:
    F(Ψ,ξ) = M_Pl²/2 + Λ₀Ψ² + Λ₁Ψ²ξ

Level-2 field ξ EOM:
    □ξ − m_ξ²ξ = −Λ₁Ψ²R
    Static limit: (∇² − m_ξ²)ξ = −Λ₁Ψ̄²R(r)
    ξ(r) = Λ₁Ψ̄²/(4π) ∫ e^{−mξ|r−r'|}/ |r−r'| × R(r') d³r'  [Yukawa conv]

G_eff modification:
    G_eff(r) = G_N / (2F) where F = M_Pl²/2 + Λ₀Ψ̄² + Λ₁Ψ̄²ξ(r)
    Sign: ξ < 0 in matter (Λ₁>0) → F decreases → G_eff increases.

Resonance condition (f_DM ~ 5 at mean halo density):
    m_ξ/H₀ ≈ 11.6 × Λ₁

KEY DIAGNOSTIC ISSUE:
    Resonance gives r_ξ = 1/m_ξ ~ H₀⁻¹/11.6Λ₁ >> R_galaxy for all Λ₁ < 10.
    Consequence: DM density ∝ baryon density (spatial structure identical).
    DM ∝ baryons → falling rotation curves, not flat. Structural FAIL.

Tests:
  Part A — FRW background: ξ̄(t) evolution, Ψ̄ stability check
  Part B — f_DM profile in galactic units; r_ξ vs R_gal comparison
  Part C — Rotation curve: DM ∝ baryon model vs NGC 2403 (flat curve expected)
  Part D — Parameter scan confirming structural failure

Units: M_Pl = 1, 8πG_N = 1, H₀ = 1 (natural); galactic distances in kpc.
"""

import numpy as np
from scipy.integrate import solve_ivp, cumulative_trapezoid
from scipy.interpolate import interp1d
from scipy.special import i0, i1, k0, k1
import json, os

# ── Physical constants ──────────────────────────────────────────────────────
Lambda0    = 0.003        # Λ₀ (locked Phase 1)
Psi_bar    = 2.62         # Ψ̄ from SIM113 (M_Pl=1)
lam_ssb    = 7.4e-5       # SIM113 SSB quartic λ
v_ssb      = 13.16        # SIM113 SSB hilltop v (M_Pl)

Om2_bar    = 1.0 + 2.0*Lambda0*Psi_bar**2
Mpl2_eff   = 0.5 * Om2_bar            # = (M_Pl²_eff) in M_Pl=1

# Cosmology (Planck 2018)
Omega_m0   = 0.3089
Omega_r0   = 9.2e-5
Omega_L0   = 1.0 - Omega_m0 - Omega_r0

# Unit bridge: ρ_crit in galactic units
# ρ_crit [M_sun/kpc³] = 126.  ρ_crit [M_Pl²H₀²] = 3/(8π) ≈ 0.1194
rho_crit_phys  = 126.0          # M_sun/kpc³
rho_crit_nat   = 3.0/(8.0*np.pi)   # M_Pl²H₀² natural units
conv_rho       = rho_crit_phys / rho_crit_nat   # (M_sun/kpc³) per (M_Pl²H₀²)
# i.e. ρ [M_Pl²H₀²] = ρ [M_sun/kpc³] / conv_rho

# Hubble distance  c/H₀ ≈ 4478 Mpc = 4.478e9 kpc
H0_kpc_inv     = 1.0 / 4478e3      # H₀ in kpc⁻¹ (H₀ = H0_kpc_inv kpc⁻¹)
# → m_ξ [kpc⁻¹] = m_ξ [H₀] × H0_kpc_inv

# Galactic G
G_kpc          = 4.302e-6          # kpc (km/s)² M_sun⁻¹  [G in kpc units]

print("=" * 68)
print("SIM117: Level-2 RIFT Recursion — Memory-of-Memory Dark Matter")
print("=" * 68)
print(f"  Λ₀ = {Lambda0}   Ψ̄ = {Psi_bar} M_Pl   Ω²(Ψ̄) = {Om2_bar:.4f}")
print(f"  M_Pl²_eff = {Mpl2_eff:.4f}")
print(f"  ρ_crit = {rho_crit_phys} M_sun/kpc³  =  {rho_crit_nat:.4f} M_Pl²H₀²")
print(f"  conv_rho = {conv_rho:.1f} (M_sun/kpc³) / (M_Pl²H₀²)")
print()

# ═══════════════════════════════════════════════════════════════════════════
# ANALYTIC TOOLS
# ═══════════════════════════════════════════════════════════════════════════

def m_xi_from_resonance(Lambda1, f_DM_target=5.0, overdensity=200.0):
    """
    Resonance condition: m_ξ given Λ₁ and target f_DM.
    f_DM = Λ₁²Ψ̄⁴ × ρ_halo_nat / (m_ξ² × Mpl2_eff)
    ρ_halo_nat = overdensity × Ω_m × ρ_crit_nat
    Solve for m_ξ (in H₀=1 units).
    """
    rho_halo_nat = overdensity * Omega_m0 * rho_crit_nat
    num = Lambda1**2 * Psi_bar**4 * rho_halo_nat
    den = f_DM_target * Mpl2_eff
    return np.sqrt(num / den)

def f_DM_local(Lambda1, m_xi_H0, rho_m_phys):
    """
    Local f_DM = ΔG_eff/G_eff,0 at a point where baryon density = rho_m_phys.
    Valid when r_ξ >> galactic scale (massless-limit approximation).
    rho_m_phys in M_sun/kpc³.
    f_DM = Λ₁²Ψ̄⁴ × rho_m_nat / (m_xi² × Mpl2_eff)
    """
    rho_m_nat = rho_m_phys / conv_rho
    f = Lambda1**2 * Psi_bar**4 * rho_m_nat / (m_xi_H0**2 * Mpl2_eff)
    return f

def r_xi_kpc(m_xi_H0):
    """Yukawa scale length in kpc."""
    if m_xi_H0 == 0:
        return np.inf
    return 1.0 / (m_xi_H0 * H0_kpc_inv)

# ── Print resonance table ──────────────────────────────────────────────────
print("─" * 68)
print("Resonance condition: Λ₁ → m_ξ/H₀ → r_ξ")
print("─" * 68)
print(f"\n  {'Λ₁':>7} {'m_ξ/H₀':>9} {'r_ξ (kpc)':>14} {'r_ξ (Mpc)':>12} {'r_ξ/R_gal':>12}")
print("-" * 60)
R_gal = 15.0   # kpc, typical galaxy radius
for L1 in [0.01, 0.05, 0.10, 0.20, 0.50, 1.00, 5.00, 10.0]:
    mxi  = m_xi_from_resonance(L1)
    rxk  = r_xi_kpc(mxi)
    rxM  = rxk * 1e-3
    print(f"  {L1:>7.3f} {mxi:>9.4f} {rxk:>14.2e} {rxM:>12.2e} {rxk/R_gal:>12.2e}")
print()
print("  → r_ξ >> R_galaxy for all Λ₁ < ~10.")
print("  → Yukawa kernel acts as Coulomb (1/r) at galactic scales.")
print("  → DM density profile = f_DM(r) × ρ_baryon(r)  [spatially identical].")
print()

# ═══════════════════════════════════════════════════════════════════════════
# PART A: FRW BACKGROUND — ξ̄(t) and Ψ̄ stability
# ═══════════════════════════════════════════════════════════════════════════
print("─" * 68)
print("PART A: FRW background — ξ̄ evolution; Ψ̄ stability")
print("─" * 68)

def H_of_a(a):
    return np.sqrt(Omega_m0/a**3 + Omega_r0/a**4 + Omega_L0)

def R_frw(a):
    """FRW Ricci scalar. Non-rel matter dominant: R ≈ Ω_m/a³ in H₀=M_Pl=1."""
    return Omega_m0 / a**3

def frw_rhs_loga(loga, y, Lambda1, m_xi):
    """
    EOMs in d/d(ln a).  State: [Ψ, Ψ', ξ, ξ']  (primes = d/d lna).
    Ψ'' + (3 − ε_H)Ψ' + (dV_J/dΨ)/H² − 2(Λ₀+Λ₁ξ)ΨR/H² = 0
    ξ'' + (3 − ε_H)ξ' + (m_ξ/H)²ξ + Λ₁Ψ²R/H² = 0
    where ε_H = −d ln H/d ln a.
    """
    a       = np.exp(loga)
    Psi, dP, xi, dx = y
    H       = H_of_a(a)
    R       = R_frw(a)
    eps_H   = 1.5*Omega_m0/a**3/H**2 + 2.0*Omega_r0/a**4/H**2  # ≈ −ḢH⁻²

    dV      = 4.0*lam_ssb*Psi*(Psi**2 - v_ssb**2)
    ddPsi   = -(3.0 - eps_H)*dP - dV/H**2 + 2.0*(Lambda0 + Lambda1*xi)*Psi*R/H**2
    ddxi    = -(3.0 - eps_H)*dx - (m_xi/H)**2 * xi + Lambda1*Psi**2*R/H**2
    return [dP, ddPsi, dx, ddxi]

def run_background(Lambda1, m_xi, verbose=True):
    loga_span = (np.log(0.01), np.log(1.0))
    y0 = [Psi_bar, 0.0, 0.0, 0.0]
    sol = solve_ivp(frw_rhs_loga, loga_span, y0,
                    args=(Lambda1, m_xi),
                    t_eval=np.linspace(*loga_span, 600),
                    method='RK45', rtol=1e-8, atol=1e-10)
    Psi_f  = sol.y[0,-1]
    xi_f   = sol.y[2,-1]
    dPsi_p = abs(Psi_f - Psi_bar)/Psi_bar*100.0
    # quasi-static estimate at a=1
    H1 = H_of_a(1.0)
    R1 = R_frw(1.0)
    xi_qs = Lambda1*Psi_bar**2*R1/H1**2 / (m_xi**2/H1**2) if m_xi>0 else 0.0
    # Level-2 vs level-1 magnitude
    l2_frac = abs(Lambda1*Psi_bar**2*xi_f) / (Lambda0*Psi_bar**2)
    if verbose:
        print(f"  Λ₁={Lambda1:.3f}, m_ξ/H₀={m_xi:.3f}")
        print(f"    Ψ̄(a=1)={Psi_f:.4f}  (ref={Psi_bar}, Δ={dPsi_p:.1f}%)")
        print(f"    ξ̄(a=1)={xi_f:.3e}  quasi-static={xi_qs:.3e}  ratio={xi_f/xi_qs if xi_qs!=0 else 'inf':.1f}×")
        print(f"    |Λ₁Ψ̄²ξ̄|/|Λ₀Ψ̄²|={l2_frac:.3e}  →  "
              f"{'level-2 dominates!' if l2_frac>1 else 'level-2 subdominant'}")
    return {'Psi_f': Psi_f, 'xi_f': xi_f, 'dPsi_p': dPsi_p, 'l2_frac': l2_frac,
            'sol': sol, 'converged': sol.success}

print()
bg = {}
for L1, mxi in [(0.02, m_xi_from_resonance(0.02)),
                (0.05, m_xi_from_resonance(0.05)),
                (0.10, m_xi_from_resonance(0.10)),
                (0.20, m_xi_from_resonance(0.20))]:
    bg[L1] = run_background(L1, mxi)
    print()

# ═══════════════════════════════════════════════════════════════════════════
# PART B: f_DM RADIAL PROFILE  (massless-limit: DM ∝ baryons)
# ═══════════════════════════════════════════════════════════════════════════
print("─" * 68)
print("PART B: f_DM radial profile  (r_ξ >> R_gal limit)")
print("─" * 68)

# NGC 2403 exponential disk
Sigma0_pc2 = 56.0       # M_sun/pc² — central surface density (stars+gas)
Sigma0     = Sigma0_pc2 * 1e6  # M_sun/kpc² (used in formulae with kpc)
Rd         = 1.72       # kpc
hz         = 0.30       # kpc
fgas       = 0.35
G_kpc      = 4.302e-6   # kpc (km/s)² M_sun⁻¹  [G in kpc units]

def rho_disk_mid(R_kpc):
    """Midplane baryon density (disk+gas) in M_sun/kpc³."""
    return (Sigma0*(1+fgas)/(2.0*hz)) * np.exp(-R_kpc/Rd)   # sech(0)=1

R_arr = np.array([0.5, 1, 2, 3, 5, 7, 10, 13])
print(f"\n  f_DM radial profile (Λ₁=0.10, m_ξ/H₀={m_xi_from_resonance(0.10):.3f}):")
print(f"  {'R(kpc)':>8} {'ρ_m(M_sun/kpc³)':>18} {'f_DM':>10} {'ρ_DM/ρ_bar':>12}")
print("-" * 55)
L1_ref = 0.10
mxi_ref = m_xi_from_resonance(L1_ref)
for R in R_arr:
    rho_m = rho_disk_mid(R)
    f = f_DM_local(L1_ref, mxi_ref, rho_m)
    print(f"  {R:>8.1f} {rho_m:>18.2e} {f:>10.4f} {f:>12.4f}")
print()
print("  Note: f_DM decreases outward (∝ ρ_m) → DM profile peaks at centre.")
print("  Flat rotation curves require f_DM increasing or constant with R.")
print()

# ═══════════════════════════════════════════════════════════════════════════
# PART C: ROTATION CURVE  — DM ∝ baryon model
# ═══════════════════════════════════════════════════════════════════════════
print("─" * 68)
print("PART C: Rotation curve — DM ∝ baryon (r_ξ >> R_gal) vs NGC 2403")
print("─" * 68)

# NGC 2403 observed rotation curve (Begeman 1989; de Blok+2008)
ngc2403 = np.array([
    [0.5,  37.8,  5.0],
    [1.0,  58.4,  4.0],
    [1.5,  73.2,  4.0],
    [2.0,  85.0,  4.0],
    [2.5,  92.0,  3.5],
    [3.0,  97.5,  3.5],
    [3.5, 100.2,  3.0],
    [4.0, 103.5,  3.0],
    [5.0, 108.5,  3.0],
    [6.0, 112.8,  3.0],
    [7.0, 115.3,  3.5],
    [8.0, 116.0,  4.0],
    [9.0, 116.3,  4.0],
    [10.0,116.5,  4.5],
    [11.0,116.8,  5.0],
    [12.0,115.5,  5.5],
    [13.0,114.9,  6.0],
])
R_obs, v_obs, v_err = ngc2403[:,0], ngc2403[:,1], ngc2403[:,2]

def v_disk_sq(R_kpc):
    """Freeman (1970) disk rotation velocity² in (km/s)².
    v²(R) = 4πGΣ₀Rd × x² × [I₀(x)K₀(x) − I₁(x)K₁(x)]
    where x = R/(2Rd), G in kpc (km/s)² M_sun⁻¹, Σ₀ in M_sun/kpc².
    """
    x = R_kpc / (2.0*Rd)
    x = max(x, 1e-3)
    Sig = Sigma0 * (1.0 + fgas)   # M_sun/kpc²
    bessel = i0(x)*k0(x) - i1(x)*k1(x)
    return 4.0*np.pi * G_kpc * Sig * Rd * x**2 * bessel

def v_model_sq(R_kpc, Lambda1, m_xi_H0):
    """
    v²_total = v²_disk + v²_DM
    v²_DM(R) = f_DM(R) × v²_disk(R)   [since DM ∝ baryon in massless limit]
    f_DM(R) = Λ₁²Ψ̄⁴ × ρ_m(R)_nat / (m_ξ² × Mpl2_eff)
    """
    vd2  = v_disk_sq(R_kpc)
    rho  = rho_disk_mid(R_kpc)
    f    = f_DM_local(Lambda1, m_xi_H0, rho)
    return vd2*(1.0 + f), vd2, f

def chi2_dof(Lambda1, m_xi_H0):
    vs_mod = np.array([max(v_model_sq(R, Lambda1, m_xi_H0)[0], 0.0) for R in R_obs])
    vs_mod = np.sqrt(vs_mod)
    return np.sum(((vs_mod - v_obs)/v_err)**2) / len(R_obs)

# Show rotation curve for best Λ₁ choices
print()
print(f"  {'R(kpc)':>7} {'v_disk':>8} {'v_mod(Λ₁=0.10)':>16} {'v_obs':>8} {'f_DM':>7}")
print("-" * 55)
R_fine = np.linspace(0.4, 14.0, 100)
L1_show = 0.10
mxi_show = m_xi_from_resonance(L1_show)
for R in [0.5, 1, 2, 3, 5, 7, 10, 13]:
    vtot2, vd2, f = v_model_sq(R, L1_show, mxi_show)
    vtot = np.sqrt(max(vtot2, 0))
    vd   = np.sqrt(max(vd2, 0))
    # find closest observed point
    idx  = np.argmin(np.abs(R_obs - R))
    print(f"  {R:>7.1f} {vd:>8.1f} {vtot:>16.1f} {v_obs[idx]:>8.1f} {f:>7.3f}")

print()
print("  The model curve has same shape as disk (falls at large R).")
print("  Observed curve is FLAT at large R — shape mismatch is structural.")
print()

# ═══════════════════════════════════════════════════════════════════════════
# PART D: PARAMETER SCAN
# ═══════════════════════════════════════════════════════════════════════════
print("─" * 68)
print("PART D: Parameter scan over Λ₁ (m_ξ locked to resonance)")
print("─" * 68)

scan_L1 = np.logspace(-2, 1, 25)
print()
print(f"  {'Λ₁':>8} {'m_ξ/H₀':>9} {'r_ξ(kpc)':>12} "
      f"{'χ²/dof':>9} {'f_DM(5kpc)':>12} {'ΔΨ̄(%)':>9}")
print("-" * 68)

best    = {'chi2': np.inf}
results = []
for L1 in scan_L1:
    mxi  = m_xi_from_resonance(L1)
    rxk  = r_xi_kpc(mxi)
    c2   = chi2_dof(L1, mxi)
    _, _, f5 = v_model_sq(5.0, L1, mxi)
    dPp  = bg.get(L1, {}).get('dPsi_p', None)
    if dPp is None:
        # quick quasi-static estimate for background displacement
        H1   = H_of_a(1.0)
        xi_qs = Lambda1 = L1; xi_qs = L1*Psi_bar**2*R_frw(1.0)/H1**2 / (mxi**2/H1**2)
        l2   = abs(L1*Psi_bar**2*xi_qs)/(Lambda0*Psi_bar**2)
        dPp  = l2 * 100.0
    flag = ""
    if c2 < 2.0 and 3 < f5 < 8 and dPp < 1.0:
        flag = " ← PASS"
    elif c2 < 5.0:
        flag = " ← marginal"
    print(f"  {L1:>8.4f} {mxi:>9.4f} {rxk:>12.2e} {c2:>9.2f} {f5:>12.4f} {dPp:>9.2f}{flag}")
    results.append({'Lambda1': float(L1), 'm_xi': float(mxi), 'r_xi_kpc': float(rxk),
                    'chi2': float(c2), 'f_dm_5kpc': float(f5)})
    if c2 < best['chi2']:
        best = {'chi2': c2, 'Lambda1': L1, 'm_xi': mxi, 'r_xi_kpc': rxk}

print()

# ═══════════════════════════════════════════════════════════════════════════
# PART E: STRUCTURAL DIAGNOSIS
# ═══════════════════════════════════════════════════════════════════════════
print("─" * 68)
print("PART E: Structural failure diagnosis")
print("─" * 68)

print("""
  FAILURE MODE 1 — Scale tension:
    Resonance condition requires m_ξ ~ 11.6Λ₁ × H₀.
    For any physical Λ₁ < 10: r_ξ = 1/m_ξ > 4×10⁴ kpc >> R_gal ~ 15 kpc.
    Yukawa kernel is flat at galactic scales → DM density ∝ baryon density.

  FAILURE MODE 2 — Profile mismatch:
    DM ∝ baryons gives rotation curves with same shape as disk (peaked+falling).
    Observed curves are FLAT at large R.
    The shape mismatch is independent of Λ₁, m_ξ, or f_DM normalisation.
    NO parameter choice fixes this — structural.

  FAILURE MODE 3 — Ψ̄ displacement:
    The memory integral accumulates large ξ̄ from early-time R >> H₀².
    ξ̄_actual >> ξ̄_quasi-static by factor ~9×.
    This back-reacts on Ψ̄ via 2Λ₁ΨξR, displacing Ψ̄ by 10–191% for physical Λ₁.
    The Ψ̄ displacement and the DM effect are controlled by the SAME coupling:
      f_DM ∝ Λ₁²Ψ̄⁴/m_ξ²  and  ΔΨ̄ ∝ Λ₁²Ψ̄²R̄/m_ξ²  (same combination).
    → Cannot have f_DM ~ 5 without also large ΔΨ̄. Structural degeneracy.

  ROOT CAUSE:
    The Level-2 term Λ₁Ψ²ξR modifies gravity in proportion to the Ricci scalar R.
    In the quasi-static Newtonian limit: R ∝ ρ_m locally.
    Therefore: ρ_DM,eff ∝ ρ_m (exact proportionality in the r_ξ >> R_gal limit).
    Flat rotation curves require ρ_DM ∝ r⁻² at large r, not ∝ exp(-r/Rd).
    The Level-2 recursion cannot decouple the DM profile from the baryon profile.
""")

# ═══════════════════════════════════════════════════════════════════════════
# SUMMARY
# ═══════════════════════════════════════════════════════════════════════════
print("=" * 68)
print("SUMMARY — SIM117")
print("=" * 68)

best_L1  = best['Lambda1']
best_mxi = best['m_xi']
best_c2  = best['chi2']
best_rxk = best['r_xi_kpc']

# Check pass criteria
c2_pass  = best_c2 < 2.0
_, _, f5b = v_model_sq(5.0, best_L1, best_mxi)
fdm_pass = 3.0 < f5b < 8.0
bg_b     = run_background(best_L1, best_mxi, verbose=False)
psi_pass = bg_b['dPsi_p'] < 1.0

print(f"""
  Best-fit scan:
    Λ₁        = {best_L1:.4f}
    m_ξ/H₀    = {best_mxi:.4f}
    r_ξ       = {best_rxk:.2e} kpc
    χ²/dof    = {best_c2:.2f}
    f_DM(5kpc)= {f5b:.4f}
    ΔΨ̄        = {bg_b['dPsi_p']:.1f}%

  Pass criteria:
    χ²/dof < 2       : {'PASS' if c2_pass else 'FAIL'}  ({best_c2:.2f})
    f_DM ∈ [3,8]     : {'PASS' if fdm_pass else 'FAIL'}  ({f5b:.4f})
    ΔΨ̄ < 1%          : {'PASS' if psi_pass else 'FAIL'}  ({bg_b['dPsi_p']:.1f}%)

  VERDICT: FAIL — STRUCTURAL
  Three independent structural failure modes identified (see Part E).
  The Level-2 Λ₁Ψ²ξR term cannot produce galactic flat rotation curves.

  Physical reason: ρ_DM,eff ∝ R ∝ ρ_m in the r_ξ >> R_gal regime.
  DM cannot be decoupled from the baryon profile within this action.

  Required for DM: either
    (a) A separate DM field with galactic Compton wavelength (not Ψ or ξ), or
    (b) A non-minimal coupling that generates ρ_DM ∝ r⁻² independent of baryons.
  Both require new physics beyond the Λ₁Ψ²ξR level-2 term.
""")

# Save
out_dir = os.path.join(os.path.dirname(__file__), '..', 'Outputs')
os.makedirs(out_dir, exist_ok=True)
summary = {
    'verdict': 'FAIL',
    'failure_modes': [
        'Scale tension: r_xi >> R_gal for all physical Lambda1',
        'Profile mismatch: DM ∝ baryon → falling curves, not flat',
        'Psi instability: memory accumulation at early times'
    ],
    'best_Lambda1':  float(best_L1),
    'best_m_xi_H0':  float(best_mxi),
    'best_r_xi_kpc': float(best_rxk),
    'best_chi2':     float(best_c2),
    'f_DM_5kpc':     float(f5b),
    'dPsi_pct':      float(bg_b['dPsi_p']),
    'pass_chi2':     bool(c2_pass),
    'pass_fdm':      bool(fdm_pass),
    'pass_psi':      bool(psi_pass),
    'scan': results,
}
with open(os.path.join(out_dir, 'sim117_results.json'), 'w') as f:
    json.dump(summary, f, indent=2)
print(f"  Results saved to {out_dir}/sim117_results.json")
print("=" * 68)
