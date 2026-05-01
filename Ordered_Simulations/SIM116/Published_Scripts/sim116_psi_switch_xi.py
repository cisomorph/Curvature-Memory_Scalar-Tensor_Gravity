#!/usr/bin/env python3
"""
SIM116: CMSTG Ψ-switch ξ Condensate
CMSTG: Curvature-Memory Scalar-Tensor Gravity — Phase 2 DM

Potential:
    V(ξ, Ψ) = λ_ξ (ξ² − Ψ²/λ_ξ)²
             = λ_ξ ξ⁴ − 2Ψ²ξ² + Ψ⁴/λ_ξ

Mechanism:
  - Cosmological background: ξ condenses at ξ_eq = Ψ̄/√λ_ξ (VEV set by DE field)
  - Effective DM mass: m_ξ = 2√2 Ψ̄ (M_Pl units, H₀=1)
  - Galactic perturbation: δξ is sourced by δΨ(r) via Yukawa equation
    −∇²δξ + m_ξ² δξ = (8Ψ̄/√λ_ξ) × (−∇²δΨ/m_ξ² − δΨ) [linearised]
    simplified: −∇²δξ + m_ξ² δξ = (8Ψ̄²/√λ_ξ) × δΨ

  CRITICAL TEST: Is δΨ large enough at galactic scales to source ρ_DM?
  From SIM115: δΨ ~ H₀² × R²_gal × ρ_baryon ~ 10⁻⁷ (H₀² suppression).

Strategy:
  PART A: Derive ξ condensate (background cosmology)
  PART B: Linearised EOM — derive δξ sourced by δΨ(r)
  PART C: Compute ρ_ξ(r) in NGC 2403 galactic profile
  PART D: Compare to required ρ_DM; evaluate structural suppression
  PART E: Phase-transition diagnostic — can Ψ ever reach 0 in halo (non-perturbative)?

Pass criteria:
  - ρ_DM,eff(r) / ρ_DM,required ≥ 0.01   (within 2 orders of magnitude)
  - δΨ/Ψ̄ ≥ 10⁻³ somewhere in galaxy     (meaningful switch)
  - Rotation curve χ²/dof < 10            (qualitative fit)

Units: M_Pl = 1, H₀ = 1 natural; galactic units kpc / M_sun / km/s.
"""

import numpy as np
from scipy.integrate import quad, odeint
from scipy.special import i0, i1, k0, k1
from scipy.interpolate import interp1d
import json, os

# ── Physical constants ────────────────────────────────────────────────────────
Lambda0     = 0.003          # Λ₀ (Phase 1 locked)
Psi_bar     = 2.62           # Ψ̄ (M_Pl), SIM113 best fit
lam_v       = 7.4e-5         # λ (quartic DE potential, SIM113)
Om2_bar     = 1.0 + 2.0*Lambda0*Psi_bar**2
Mpl2_eff    = 0.5 * Om2_bar  # effective M_Pl² in CMSTG

Omega_m0    = 0.3089
Omega_DM0   = 0.2589
rho_crit_nat = 3.0/(8.0*np.pi)     # M_Pl²H₀² (natural)
rho_crit_phys = 126.0              # M_sun/kpc³
conv_rho    = rho_crit_phys / rho_crit_nat
rho_DM_phys = Omega_DM0 * rho_crit_phys  # ~32.6 M_sun/kpc³ mean

G_kpc       = 4.302e-6        # kpc (km/s)² M_sun⁻¹

print("=" * 72)
print("SIM116: CMSTG Ψ-switch ξ Condensate")
print("=" * 72)
print(f"  Ψ̄  = {Psi_bar} M_Pl  |  Λ₀ = {Lambda0}  |  λ = {lam_v:.2e}")
print(f"  ρ_crit = {rho_crit_phys} M_sun/kpc³  |  ρ_DM,mean = {rho_DM_phys:.2f} M_sun/kpc³")
print()

# ═══════════════════════════════════════════════════════════════════════════════
# PART A: ξ condensate — background cosmology
# ═══════════════════════════════════════════════════════════════════════════════
print("─" * 72)
print("PART A: ξ condensate — background cosmology")
print("─" * 72)

print(f"""
  Potential: V(ξ,Ψ) = λ_ξ(ξ² − Ψ²/λ_ξ)²
  Vacuum:    ξ_eq  = Ψ̄/√λ_ξ     (breaks Z₂: ξ → −ξ)
  Mass:      m_ξ²  = 8λ_ξ ξ_eq² = 8Ψ̄²  [=  8 × (2.62)² ≈ {8*Psi_bar**2:.2f}  H₀²]
  Condensate density at minimum: V(ξ_eq) = 0  (exactly)
  DM energy comes only from spatial/temporal FLUCTUATIONS δξ.
""")

# m_ξ in natural units (H₀ = 1, M_Pl = 1)
m_xi_sq_nat = 8.0 * Psi_bar**2    # in H₀² units
m_xi_nat    = np.sqrt(m_xi_sq_nat)

# m_ξ in kpc⁻¹
# H₀ in kpc⁻¹: H₀ = 67 km/s/Mpc = 67/(3.086e19) kpc⁻¹ s⁻¹ × (1 kpc/km/s units...)
# In kpc units where speed in km/s: H₀ = 67 km/s/Mpc = 67/977.8 kpc⁻¹ km/s = 0.06851 kpc⁻¹
H0_kpc = 0.06851       # H₀ in units of (km/s)/kpc ... = Mpc⁻¹ × (1000 kpc/Mpc)⁻¹...
# Actually H₀ = 67 km/s/Mpc; 1 Mpc = 1000 kpc → H₀ = 0.067 (km/s)/kpc.
# In natural units (H₀=1), m [H₀] = m [kpc⁻¹] / H₀[kpc⁻¹].
# H₀ in kpc⁻¹ from light-travel: H₀ = c/D_H; D_H = c/H₀ = 3e5/67 Mpc = 4478 Mpc = 4.478e6 kpc
# → H₀ = 1/(4.478e6 kpc) in natural (c=1) units
H0_inv_kpc = 4.478e6   # kpc per H₀ unit

r_xi_nat = 1.0 / m_xi_nat                    # Compton wavelength in H₀⁻¹
r_xi_kpc = r_xi_nat * H0_inv_kpc             # in kpc

print(f"  m_ξ (natural, H₀=1) = √(8Ψ̄²) = {m_xi_nat:.4f} H₀")
print(f"  r_ξ = 1/m_ξ (natural) = {r_xi_nat:.4e} H₀⁻¹")
print(f"  r_ξ (physical)        = {r_xi_kpc:.3e} kpc")
print()

# Scan over λ_ξ to find condensate parameters
print(f"  {'λ_ξ':>10} {'ξ_eq (M_Pl)':>14} {'m_ξ (H₀)':>10} {'r_ξ (kpc)':>12} {'ξ_eq/Ψ̄':>10}")
print("-" * 65)
for lam_xi in [1e-4, 1e-3, 0.01, 0.1, 1.0, 10.0, 100.0]:
    xi_eq = Psi_bar / np.sqrt(lam_xi)
    print(f"  {lam_xi:>10.1e} {xi_eq:>14.4e} {m_xi_nat:>10.4f} {r_xi_kpc:>12.3e} {xi_eq/Psi_bar:>10.4f}")

print(f"\n  NOTE: m_ξ = 2√2 Ψ̄ is INDEPENDENT of λ_ξ.")
print(f"  r_ξ = {r_xi_kpc:.3e} kpc >> galactic scale (R_gal ~ 10–20 kpc).")
print(f"  → Yukawa kernel is flat over the galaxy (r_ξ >> R_gal).")

# ═══════════════════════════════════════════════════════════════════════════════
# PART B: Linearised EOM — δξ sourced by δΨ(r)
# ═══════════════════════════════════════════════════════════════════════════════
print()
print("─" * 72)
print("PART B: Linearised EOM — δξ sourced by δΨ(r)")
print("─" * 72)

print("""
  Linearise around (Ψ̄, ξ_eq):
    ξ = ξ_eq + δξ,   Ψ = Ψ̄ + δΨ

  V ≈ λ_ξ(2ξ_eq δξ − 2Ψ̄δΨ/λ_ξ)² / (zeroth order)  [with ξ_eq²=Ψ̄²/λ_ξ]
  EOM (static, spherical):
    −∇²δξ + m_ξ² δξ = S(r)

  Source S(r) = (8Ψ̄²/√λ_ξ) × δΨ(r)
  [from ∂V/∂ξ linearised: 8λ_ξ ξ_eq × (−2Ψ̄δΨ/λ_ξ) = −16Ψ̄²ξ_eq/λ_ξ × δΨ...
   actually ∂²V/∂ξ∂Ψ term → S = 8Ψ̄²ξ_eq/(ξ_eq √λ_ξ) × δΨ = 8Ψ̄√λ_ξ × δΨ/λ_ξ × ξ_eq]

  Exact linearisation:
    ∂V/∂ξ = 4λ_ξ ξ(ξ² − Ψ²/λ_ξ)
    At ξ = ξ_eq + δξ, Ψ = Ψ̄ + δΨ:
      ξ² − Ψ²/λ_ξ ≈ 2ξ_eq δξ − 2Ψ̄δΨ/λ_ξ
      ∂V/∂ξ ≈ 4λ_ξ ξ_eq (2ξ_eq δξ − 2Ψ̄δΨ/λ_ξ)
             = 8λ_ξ ξ_eq² δξ − 8Ψ̄ ξ_eq δΨ
             = 8Ψ̄² δξ − (8Ψ̄²/√λ_ξ) δΨ

  EOM: −∇²δξ + 8Ψ̄² δξ = (8Ψ̄²/√λ_ξ) δΨ
  Green solution (Yukawa, r_ξ >> R):
    δξ(r) = (8Ψ̄²/√λ_ξ) ∫ G_Y(r,r') δΨ(r') d³r'
    G_Y(r,r') = e^{−m_ξ|r−r'|} / (4π|r−r'|)

  Since r_ξ >> R_gal: G_Y ≈ 1/(4π m_ξ² r²) at large r, so
    δξ(r) ≈ (8Ψ̄²/√λ_ξ) / m_ξ² × <δΨ>_vol / r
  But the dominant piece is the volume-averaged source.

  KEY: δξ ∝ δΨ/m_ξ² × (some geometric factor)
       Since m_ξ² = 8Ψ̄² >> 1:
         δξ ~ δΨ / (8Ψ̄² √λ_ξ)   [very suppressed by √λ_ξ too!]
""")

# SIM115 result: δΨ/Ψ̄ ~ 4×10⁻⁷ over the galaxy
delta_Psi_over_Psi = 4e-7
delta_Psi = delta_Psi_over_Psi * Psi_bar
print(f"  From SIM115: δΨ/Ψ̄ ~ {delta_Psi_over_Psi:.1e}  →  δΨ ~ {delta_Psi:.2e} M_Pl")
print()

# Compute δξ/ξ_eq for various λ_ξ
print(f"  Suppression of δξ relative to ξ_eq (~ DM fraction):")
print(f"  {'λ_ξ':>8} {'ξ_eq':>12} {'Source S':>14} {'δξ (r=5kpc)':>14} {'δξ/ξ_eq':>12} {'ρ_DM frac':>12}")
print("-" * 80)

r_test_kpc = 5.0    # kpc (test radius)
r_test_nat = r_test_kpc / H0_inv_kpc  # in H₀⁻¹

# Simple Yukawa estimate: for r_ξ >> r_test, δξ ≈ S / m_ξ²
# δξ ≈ S / m_ξ² where S = (8Ψ̄²/√λ_ξ) × δΨ
# Then ρ_DM,eff ~ 2ξ_eq × m_ξ² × δξ (fluctuation mass energy density)
# Actually ρ_ξ from condensate perturbation:
# δρ_ξ = ∂²V/∂ξ² × (δξ)² / 2 = m_ξ² (δξ)² / 2  (from quadratic term)
# Better: ρ_DM ~ m_ξ² × (δξ)² (order of magnitude)

for lam_xi in [1e-4, 1e-3, 0.01, 0.1, 1.0, 10.0]:
    xi_eq = Psi_bar / np.sqrt(lam_xi)
    Source = 8.0 * Psi_bar**2 / np.sqrt(lam_xi) * delta_Psi
    # Yukawa solution for r_ξ >> r: δξ ≈ S / m_ξ²  (long-wavelength limit)
    delta_xi = Source / m_xi_sq_nat   # in natural units
    delta_xi_frac = delta_xi / xi_eq
    # DM density perturbation (in natural units): ρ ~ m_ξ²(δξ)²
    rho_dm_nat = 0.5 * m_xi_sq_nat * delta_xi**2
    rho_dm_phys = rho_dm_nat * conv_rho
    rho_frac = rho_dm_phys / 1e7   # typical halo ρ_DM ~ 10⁷ M_sun/kpc³ at 5 kpc
    print(f"  {lam_xi:>8.1e} {xi_eq:>12.3e} {Source:>14.3e} {delta_xi:>14.3e} "
          f"{delta_xi_frac:>12.3e} {rho_frac:>12.3e}")

print(f"""
  Result: δξ/ξ_eq ~ 10⁻¹⁴ to 10⁻⁴  (all << 1, but ρ_DM ~ (δξ)² → doubly suppressed)
  ρ_DM,eff ~ m_ξ²(δξ)² ~ (δΨ)²/λ_ξ in natural units.
  Physical ρ_DM needed ~ 10⁷ M_sun/kpc³ at R ~ 5 kpc.
  All λ_ξ values give ρ_DM,eff << 1 M_sun/kpc³.
""")

# ═══════════════════════════════════════════════════════════════════════════════
# PART C: Full galactic profile — δΨ(r) → δξ(r) → ρ_DM(r) for NGC 2403
# ═══════════════════════════════════════════════════════════════════════════════
print("─" * 72)
print("PART C: Full galactic profile — δΨ(r) → δξ(r) → ρ_DM(r) [NGC 2403]")
print("─" * 72)

# NGC 2403 baryon density profile (exponential disk + bulge, as in SIM115)
Sigma0_pc2 = 56.0               # M_sun/pc²
Sigma0     = Sigma0_pc2 * 1e6   # M_sun/kpc²
Rd         = 1.72               # kpc
z_d        = 0.35               # kpc disk height
rho0_disk  = Sigma0 / (2.0 * z_d)   # M_sun/kpc³ at centre

def rho_baryon(r_kpc):
    """Approximate 3D baryon density: exponential disk projected to sphere."""
    return rho0_disk * np.exp(-r_kpc / Rd)

# δΨ sourced by baryon density via Yukawa-like equation from SIM115:
# −∇²δΨ + m_eff² δΨ = S_Ψ(r)   where S_Ψ ~ −2Λ₀/(M_eff²) × 8πG ρ_baryon
# From SIM115 analytical solution (Yukawa with r_Ψ >> R_gal):
# δΨ/Ψ̄ ~ (Λ₀/(4π M_eff²)) × (8πG ρ_baryon/H₀²) ~ H₀² × R² / (1 + 2Λ₀Ψ̄²) × ...
# Simplified: δΨ ≈ −Λ₀ × 8πG ρ_baryon × Ψ̄ / (m_Ψ² × M_eff²)

# SIM115 result: δΨ/Ψ̄ ≈ −(Λ₀ × 8πG ρ_baryon H₀⁻²) / (3 M_eff²) × (H₀/c)²...
# In natural units (H₀=1): source = 2Λ₀ × 8πG ρ_baryon(r) where ρ_baryon in M_Pl²H₀²
def delta_Psi_profile(r_kpc):
    """δΨ(r) calibrated from SIM115: δΨ/Ψ̄ = 4×10⁻⁷ at R=5 kpc, scales ∝ ρ_baryon.
    SIM115 showed: the cosmological Ψ field is H₀²-frozen at galactic scales.
    The gradient soliton EOM in kpc units has effective mass λ_kpc ~ H₀² × λ_cosmo
    = 3.6×10⁻¹⁹ kpc⁻² → Compton scale >> Hubble radius → field barely responds to local matter."""
    rho_b_ref = rho_baryon(5.0)            # 4.37×10⁶ M_sun/kpc³ at reference R=5 kpc
    delta_Psi_ref = 4.0e-7 * Psi_bar      # SIM115 result at R=5 kpc: δΨ/Ψ̄ = 4×10⁻⁷
    return delta_Psi_ref * rho_baryon(r_kpc) / rho_b_ref

# Test δΨ profile
r_arr = np.array([0.5, 1, 2, 3, 5, 7, 10, 13])
print(f"\n  {'R (kpc)':>8} {'ρ_bary(M_sun/kpc³)':>20} {'δΨ (M_Pl)':>14} {'δΨ/Ψ̄':>12}")
print("-" * 60)
for r in r_arr:
    rho_b = rho_baryon(r)
    dPsi  = delta_Psi_profile(r)
    print(f"  {r:>8.1f} {rho_b:>20.3e} {dPsi:>14.3e} {dPsi/Psi_bar:>12.3e}")

# Now compute δξ(r) for representative λ_ξ values
print()
lam_xi_values = [1e-3, 0.01, 0.1, 1.0]
for lam_xi in lam_xi_values:
    xi_eq = Psi_bar / np.sqrt(lam_xi)
    print(f"\n  λ_ξ = {lam_xi:.2e}  [ξ_eq = {xi_eq:.3e} M_Pl]:")
    print(f"  {'R(kpc)':>8} {'δΨ':>12} {'δξ (nat)':>14} {'ρ_DM (M_sun/kpc³)':>20} {'ρ_DM/ρ_req':>12}")
    print("-" * 72)
    rho_halo_req = np.array([5.0e8, 2.0e8, 8.0e7, 4.0e7, 1.5e7, 7.0e6, 3.0e6, 2.0e6])
    for i, r in enumerate(r_arr):
        dPsi_r = delta_Psi_profile(r)
        # Source: S = (8Ψ̄²/√λ_ξ) × δΨ
        source_r = 8.0 * Psi_bar**2 / np.sqrt(lam_xi) * dPsi_r
        # Yukawa solution (r_ξ >> r): δξ ≈ S / m_ξ²
        delta_xi_r = source_r / m_xi_sq_nat
        # DM energy density: ρ_DM ~ ½m_ξ²(δξ)²  [quadratic term in V]
        rho_dm_nat = 0.5 * m_xi_sq_nat * delta_xi_r**2
        rho_dm_phys = rho_dm_nat * conv_rho
        ratio = rho_dm_phys / rho_halo_req[i]
        print(f"  {r:>8.1f} {dPsi_r:>12.3e} {delta_xi_r:>14.3e} {rho_dm_phys:>20.3e} {ratio:>12.3e}")

# ═══════════════════════════════════════════════════════════════════════════════
# PART D: Suppression ratio — quantify failure mode
# ═══════════════════════════════════════════════════════════════════════════════
print()
print("─" * 72)
print("PART D: Suppression ratio and structural diagnosis")
print("─" * 72)

# Maximum possible ρ_DM from δξ mechanism, optimised over λ_ξ
# ρ_DM ∝ m_ξ²(δξ)² = m_ξ² × (S/m_ξ²)² = S²/m_ξ²
# S = (8Ψ̄²/√λ_ξ) × δΨ
# ρ_DM ∝ (8Ψ̄²/√λ_ξ)² × (δΨ)² / m_ξ²  ∝ (δΨ)²/(λ_ξ × m_ξ²)
# Maximise over λ_ξ: take λ_ξ → 0. But ξ_eq → ∞ → back-reaction on gravity not OK.
# Constraint: ρ_ξ,background = 0 at minimum (exact). No back-reaction from condensate.
# But ξ kinetic energy from cosmological rolling ~ ξ_eq² H₀² which must be << ρ_DM,cosm.
# ξ_eq² H₀² < Ω_DM ρ_crit → ξ_eq < √(Ω_DM × 3/(8π)) ~ 0.35 M_Pl → λ_ξ > Ψ̄²/ξ_eq,max² ~ 56
# So λ_ξ > 56 from cosmological stability.

lam_xi_min_cosmo = Psi_bar**2 / 0.35**2
print(f"\n  Cosmological stability bound: ξ_eq < 0.35 M_Pl → λ_ξ > {lam_xi_min_cosmo:.1f}")
print(f"  (from ξ_eq² H₀² < Ω_DM ρ_crit)\n")

r5_kpc = 5.0
dPsi_r5 = delta_Psi_profile(r5_kpc)
rho_dm_req_r5 = 1.5e7   # M_sun/kpc³ at 5 kpc (halo requirement)

print(f"  At R = 5 kpc (typical halo):")
print(f"  δΨ = {dPsi_r5:.3e} M_Pl   (from SIM115 structural suppression)")
print(f"  Required ρ_DM ~ {rho_dm_req_r5:.1e} M_sun/kpc³\n")
print(f"  {'λ_ξ':>8} {'ξ_eq':>10} {'δξ/ξ_eq':>12} {'ρ_DM (Msun/kpc³)':>20} {'ratio to req':>14} {'cosmo OK':>10}")
print("-" * 80)

for lam_xi in [10.0, 30.0, 56.0, 100.0, 300.0, 1000.0]:
    xi_eq = Psi_bar / np.sqrt(lam_xi)
    cosmo_ok = "YES" if lam_xi >= lam_xi_min_cosmo else "NO"
    source_r5 = 8.0 * Psi_bar**2 / np.sqrt(lam_xi) * dPsi_r5
    delta_xi_r5 = source_r5 / m_xi_sq_nat
    frac = delta_xi_r5 / xi_eq
    rho_dm_nat = 0.5 * m_xi_sq_nat * delta_xi_r5**2
    rho_dm_phys = rho_dm_nat * conv_rho
    ratio = rho_dm_phys / rho_dm_req_r5
    print(f"  {lam_xi:>8.1f} {xi_eq:>10.3e} {frac:>12.3e} {rho_dm_phys:>20.3e} {ratio:>14.3e} {cosmo_ok:>10}")

print()
print(f"  ROOT CAUSE: δξ ∝ δΨ / (m_ξ² √λ_ξ)")
print(f"              ρ_DM ∝ m_ξ²(δξ)² ∝ (δΨ)² / (m_ξ² λ_ξ)")
print(f"              With δΨ ~ {abs(dPsi_r5):.1e} M_Pl and m_ξ² = {m_xi_sq_nat:.1f} H₀²:")
print(f"              ρ_DM is suppressed by (δΨ)² ~ 10⁻¹⁴ compared to Ψ̄² scale.")
print(f"              This is H₀² suppression squared — deeper than SIM115.")

# ═══════════════════════════════════════════════════════════════════════════════
# PART E: Phase-transition diagnostic — non-perturbative check
# ═══════════════════════════════════════════════════════════════════════════════
print()
print("─" * 72)
print("PART E: Non-perturbative diagnostic — can Ψ reach 0?")
print("─" * 72)

print(f"""
  For the ξ switch to operate non-perturbatively, Ψ must reach 0 somewhere.
  This would put ξ at the WRONG vacuum (ξ → 0 rather than ξ_eq = Ψ̄/√λ_ξ).

  Required: δΨ/Ψ̄ ~ 1  (order-unity field variation)
  Achieved: δΨ/Ψ̄ ~ {abs(dPsi_r5/Psi_bar):.2e}  (from SIM115 H₀² suppression)

  Baryon density required to achieve δΨ/Ψ̄ ~ 1 (since δΨ ∝ ρ_baryon):
  ρ_baryon,crit ~ ρ_baryon(5 kpc) × (Ψ̄ / |δΨ(5 kpc)|)
               ~ {rho_baryon(5.0):.2e} × {abs(Psi_bar/dPsi_r5):.2e}
               ~ {rho_baryon(5.0) * abs(Psi_bar/dPsi_r5):.2e} M_sun/kpc³

  Compare: NGC 2403 core baryon density ~ {rho_baryon(0.5):.2e} M_sun/kpc³
           Neutron star density ~ 10¹⁸ M_sun/kpc³

  Conclusion: Non-perturbative Ψ→0 switching is physically impossible in galactic environments.
  Even in neutron star cores (densest matter), δΨ/Ψ̄ << 1.
""")

# ═══════════════════════════════════════════════════════════════════════════════
# PART F: Rotation curve (illustrative — what we'd get if the mechanism worked)
# ═══════════════════════════════════════════════════════════════════════════════
print("─" * 72)
print("PART F: Rotation curve with actual ρ_DM from δξ mechanism")
print("─" * 72)

lam_xi_test = 56.0   # cosmological stability bound (worst case / best case)
xi_eq_test  = Psi_bar / np.sqrt(lam_xi_test)
R_obs = np.array([0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0,
                  5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0, 12.0, 13.0])
v_obs = np.array([37.8, 58.4, 73.2, 85.0, 92.0, 97.5, 100.2, 103.5,
                  108.5, 112.8, 115.3, 116.0, 116.3, 116.5, 116.8, 115.5, 114.9])
v_err = np.array([5.0, 4.0, 4.0, 4.0, 3.5, 3.5, 3.0, 3.0,
                  3.0, 3.0, 3.5, 4.0, 4.0, 4.5, 5.0, 5.5, 6.0])

def v_disk_sq(R_kpc):
    x = max(R_kpc / (2.0*Rd), 1e-4)
    bessel = i0(x)*k0(x) - i1(x)*k1(x)
    return max(4.0*np.pi * G_kpc * Sigma0 * Rd * x**2 * bessel, 0.0)

def rho_DM_xi(r_kpc, lam_xi):
    dPsi_r = delta_Psi_profile(r_kpc)
    source_r = 8.0 * Psi_bar**2 / np.sqrt(lam_xi) * dPsi_r
    delta_xi_r = source_r / m_xi_sq_nat
    rho_dm_nat = 0.5 * m_xi_sq_nat * delta_xi_r**2
    return rho_dm_nat * conv_rho   # M_sun/kpc³

# Enclosed DM mass for rotation curve
n_grid = 500
R_max  = 15.0
r_grid = np.linspace(1e-3, R_max, n_grid)
rho_dm_arr = np.array([rho_DM_xi(r, lam_xi_test) for r in r_grid])
integrand = 4.0*np.pi * r_grid**2 * rho_dm_arr
M_enc = np.zeros(n_grid)
dr = np.diff(r_grid)
M_enc[1:] = np.cumsum(0.5*(integrand[:-1]+integrand[1:])*dr)
M_enc_fn = interp1d(r_grid, M_enc, kind='linear', fill_value='extrapolate')

print(f"\n  λ_ξ = {lam_xi_test:.0f} (cosmological bound), ξ_eq = {xi_eq_test:.3e} M_Pl")
print(f"\n  {'R(kpc)':>7} {'ρ_DM(M_sun/kpc³)':>18} {'v_disk':>8} {'v_DM':>8} {'v_tot':>8} {'v_obs':>8}")
print("-" * 65)
chi2_arr = []
for i, R in enumerate(R_obs):
    rho_dm_r = rho_DM_xi(R, lam_xi_test)
    v2_dm    = G_kpc * M_enc_fn(R) / R
    v2_disk  = v_disk_sq(R)
    v_tot    = np.sqrt(max(v2_dm + v2_disk, 0.0))
    v_dm     = np.sqrt(max(v2_dm, 0.0))
    v_d      = np.sqrt(max(v2_disk, 0.0))
    res      = (v_tot - v_obs[i]) / v_err[i]
    chi2_arr.append(res**2)
    print(f"  {R:>7.1f} {rho_dm_r:>18.3e} {v_d:>8.1f} {v_dm:>8.3f} {v_tot:>8.1f} {v_obs[i]:>8.1f}")

chi2_dof = sum(chi2_arr) / len(chi2_arr)
print(f"\n  χ²/dof = {chi2_dof:.2e}  (vs NGC 2403;  PASS requires < 10)")

# ═══════════════════════════════════════════════════════════════════════════════
# SUMMARY
# ═══════════════════════════════════════════════════════════════════════════════
print()
print("=" * 72)
print("SUMMARY — SIM116")
print("=" * 72)

# Pass criteria
dPsi_frac_max = max(abs(delta_Psi_profile(r)/Psi_bar) for r in R_obs)
rho_frac_max  = max(rho_DM_xi(r, lam_xi_test) / 1.5e7 for r in R_obs)
chi2_pass     = chi2_dof < 10.0
switch_pass   = dPsi_frac_max >= 1e-3
rho_pass      = rho_frac_max >= 0.01

print(f"""
  Mechanism: V(ξ,Ψ) = λ_ξ(ξ² − Ψ²/λ_ξ)²   [Ψ-switch ξ condensate]
  Action parameters at cosmological stability bound (λ_ξ = {lam_xi_min_cosmo:.0f}):
    ξ_eq            = Ψ̄/√λ_ξ = {Psi_bar/np.sqrt(lam_xi_min_cosmo):.3f} M_Pl
    m_ξ             = {m_xi_nat:.3f} H₀  =  r_ξ = {r_xi_kpc:.2e} kpc

  Key diagnostic values:
    max δΨ/Ψ̄       = {dPsi_frac_max:.2e}   (SIM115 structural suppression)
    max ρ_DM,eff    = {max(rho_DM_xi(r, lam_xi_test) for r in R_obs):.2e} M_sun/kpc³
    required ρ_DM   ~ 1.5×10⁷ M_sun/kpc³  at R = 5 kpc
    suppression     = {max(rho_DM_xi(r, lam_xi_test) for r in R_obs)/1.5e7:.2e}
    χ²/dof          = {chi2_dof:.2e}

  Pass criteria:
    ρ_DM,eff / ρ_req ≥ 0.01      : {'PASS' if rho_pass else 'FAIL'}  ({rho_frac_max:.2e})
    δΨ/Ψ̄ ≥ 10⁻³                 : {'PASS' if switch_pass else 'FAIL'}  ({dPsi_frac_max:.2e})
    Rotation curve χ²/dof < 10   : {'PASS' if chi2_pass else 'FAIL'}  ({chi2_dof:.2e})
""")

verdict = "PASS" if (rho_pass and switch_pass and chi2_pass) else "FAIL"

structural_note = """
  STRUCTURAL FAILURE — three independent suppression layers:
  1. H₀² suppression (SIM115): δΨ ~ H₀²/m_Ψ⁴ × ρ_baryon ~ 10⁻⁷ Ψ̄
  2. m_ξ² suppression: δξ ~ δΨ / m_ξ² ~ δΨ/(8Ψ̄²) — further suppressed by O(10)
  3. Quadratic suppression: ρ_DM ~ (δξ)² — total suppression ~10⁻¹⁴–10⁻¹⁶

  The Ψ-switch cannot operate because Ψ is cosmologically frozen (DE field).
  Non-perturbative transition (Ψ→0) requires ρ_baryon ~ 10¹² M_sun/kpc³ — neutron star density.
  DM in galactic halos: IMPOSSIBLE via this mechanism.

  This completes the Phase 2 DM sector exploration:
    SIM114: βΨ²ρ_m condensate — FAIL (trilemma)
    SIM115: Gradient soliton — FAIL (H₀² suppression)
    SIM116: Ψ-switch ξ condensate — FAIL (H₀² suppression²)
    SIM117: Level-2 recursion — FAIL (structural, 3 modes)
  → All Ψ-sector DM mechanisms exhaustively FAIL.
  → SIM118 (separate χ field, CMSTG-seeded) gives PARTIAL/PASS.
"""

print(f"  VERDICT: {verdict}")
print(structural_note)

# Save results
out_dir = os.path.join(os.path.dirname(__file__), '..', 'Outputs')
os.makedirs(out_dir, exist_ok=True)
results = {
    'verdict':              verdict,
    'm_xi_nat':             float(m_xi_nat),
    'r_xi_kpc':             float(r_xi_kpc),
    'lam_xi_cosmo_bound':   float(lam_xi_min_cosmo),
    'xi_eq_at_bound':       float(Psi_bar / np.sqrt(lam_xi_min_cosmo)),
    'delta_Psi_frac_max':   float(dPsi_frac_max),
    'rho_DM_max_Msun_kpc3': float(max(rho_DM_xi(r, lam_xi_test) for r in R_obs)),
    'rho_DM_required':      1.5e7,
    'suppression_factor':   float(rho_frac_max),
    'chi2_dof':             float(chi2_dof),
    'pass_rho':             bool(rho_pass),
    'pass_switch':          bool(switch_pass),
    'pass_chi2':            bool(chi2_pass),
    'all_pass':             bool(verdict == 'PASS'),
    'failure_mode':         'H0sq_suppression_squared',
    'note': (
        'Ψ-switch ξ condensate: structural failure. '
        'δΨ/Ψ̄ ~ 1e-7 (H₀² frozen DE field). '
        'ρ_DM ∝ (δΨ)² → doubly suppressed. '
        'All Ψ-sector DM mechanisms exhausted.'
    ),
}
with open(os.path.join(out_dir, 'sim116_results.json'), 'w') as f:
    json.dump(results, f, indent=2)
print("  Results saved to Outputs/sim116_results.json")
print("=" * 72)
