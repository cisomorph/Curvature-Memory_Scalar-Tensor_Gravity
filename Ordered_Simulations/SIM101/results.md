# SIM101 — Results

**Phase:** Phase 1  
**Title:** CMSTG BAO Sound Horizon from First Principles  
**Verdict:** PASS (ΛCDM baseline) / CMSTG shifts sub-ppm at Λ₀=0.003

## What was tested

Sound horizon r_d by integrating the CMSTG FLRW+scalar system from the radiation era through recombination. Scanned Λ₀ ∈ {0, 0.001, 0.003, 0.008, 0.013, 0.05, 0.095}.

## Key results

- ΛCDM baseline: r_d = 146.67 Mpc (Planck ref: 147.09 Mpc) — PASS
- CMSTG shift at Λ₀=0.003: −0.679 ppm
- CMSTG shift at Λ₀=0.095: −21.5 ppm
- BAO observational precision: ~3000 ppm → CMSTG r_d is degenerate with ΛCDM across all physical Λ₀

## Bug fix (2026-04-26)

Original script used total radiation density Ω_r in R_b = 3ρ_b/(4ρ_γ), but R_b requires photon-only density. Fix: derive Ω_γ = Ω_r / (1 + N_eff × (7/8) × (4/11)^(4/3)) with N_eff=3.046 and use Ω_γ in R_b. H(z) continues to use Ω_r (total). ΛCDM baseline improved from 153.05 Mpc (FAIL) to 146.67 Mpc (PASS). CMSTG fractional shifts unchanged.

## Context

Phase 1 canonical: Λ₀=0.003, Ψ̄=2.62 M_Pl, F_eff(z_CMB)=0.521, H₀≈67.59 km/s/Mpc.
