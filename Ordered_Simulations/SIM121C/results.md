# SIM121C — Results

**Phase:** Phase 2  
**Title:** Joint DESI+Planck MCMC  
**Verdict:** PARTIAL

## What was tested

Full joint DESI+Planck MCMC at Phase 2 canonical; PARTIAL — 2.77σ floor.

## Key numerical results

- DESI tension = 2.63

## Context

Phase 1 canonical: Λ₀=0.003, Ψ̄=2.62 M_Pl, F_eff(z_CMB)=0.521, H₀≈67.59 km/s/Mpc.


## Audit rerun (2026-10-22)

**Audited verdict:** PARTIAL (published: PARTIAL, 'tension 2.63σ')

- Official DESI DR1/DR2 likelihood; r_s to ∞; θ* = r_s(z*)/D_M(z*); H0 closure fixed; proper statistics.
- Frozen F0 is exactly degenerate with H0 for θ* + BAO: F0 adds Δχ² = 0.00. The Phase-1 frozen model has the same χ² as ΛCDM (22.04 DR1 / 32.30 DR2).
- All of the improvement is the phenomenological CPL (w0, wa): Δχ² = +10.53 (DR1), +23.78 (DR2). DR2 outputs are in Outputs/DR2/.
- Data/statistics: official DESI DR1/DR2 likelihood, r_s to z = ∞, θ* = r_s(z*)/D_M(z*), χ²/p-value/GOF σ, and Δχ² vs ΛCDM in the same pipeline.
- Details: `AUDIT_NOTES.md` (repo root); log: `audit_logs/`.
