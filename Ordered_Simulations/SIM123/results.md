# SIM123 — Results

**Phase:** Phase 2  
**Title:** Phase 3 Direction B: Inverse Ψ  
**Verdict:** FAIL

## What was tested

1/Ψ coupling extension; FAIL — GR recovery broken.

## Context

Phase 1 canonical: Λ₀=0.003, Ψ̄=2.62 M_Pl, F_eff(z_CMB)=0.521, H₀≈67.59 km/s/Mpc.


## Audit rerun (2026-10-22)

**Audited verdict:** faithful PASS / θ*-matched PASS (published: FAIL) — but Cassini 162×

- Corrected full equations; hilltop V as the dark energy; no bare Λ.
- Best: Λ0 = 0.0064, Ψ_ini = 2.0 → Δχ² = +9.30 (DR1) / +21.51 (DR2). All 100 grid points have Δχ² ≥ 0.
- |γ−1| = 3.7e-3 at the best point (162× Cassini).
- Data/statistics: official DESI DR1/DR2 likelihood, r_s to z = ∞, θ* = r_s(z*)/D_M(z*), χ²/p-value/GOF σ, and Δχ² vs ΛCDM in the same pipeline.
- Details: `AUDIT_NOTES.md` (repo root); log: `audit_logs/`.
