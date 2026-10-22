# SIM126 — Results

**Phase:** Phase 3  
**Title:** Exponential Memory Decay  
**Verdict:** N/A

## What was tested

e^{−βM} coupling; FAIL — insufficient H(z) modification.

## Context

Phase 1 canonical: Λ₀=0.003, Ψ̄=2.62 M_Pl, F_eff(z_CMB)=0.521, H₀≈67.59 km/s/Mpc.


## Audit rerun (2026-10-22)

**Audited verdict:** N/A (toy; observables corrected only)

- ε = 2e-10 improves on its own baseline by Δχ² = +13 (DR1), but θ* = 1.0398 and E(1) ≠ 1.
- See Outputs/sim126_audit.json.
- Data/statistics: official DESI DR1/DR2 likelihood, r_s to z = ∞, θ* = r_s(z*)/D_M(z*), χ²/p-value/GOF σ, and Δχ² vs ΛCDM in the same pipeline.
- Details: `AUDIT_NOTES.md` (repo root); log: `audit_logs/`.
