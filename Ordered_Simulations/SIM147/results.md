# SIM147 Results

**Phase:** 5 (post-JCAP exploratory)  
**Topic:** Curvature-Memory Kernel Decay Characterization  
**Gate for:** SIM148 (Ψ_pre), SIM149 (Ψ_local)  
**Date:** 2026-05-14

## Gate: RED

Both Ψ_pre and Ψ_local conjectures are effectively dead at the locked k_m = 10 Mpc⁻¹.

## Key Numbers

| Quantity | Value |
|---|---|
| Kernel form | Gaussian: K(Δt) ∝ Δt · exp(−(k_m Δt)²/4) |
| τ_mem = 2/k_m | 0.2 Mpc = **2.05×10⁻⁴ Gyr** ≈ 205,000 yr |
| Δt_peak = √2/k_m | 0.141 Mpc = 1.45×10⁻⁴ Gyr |
| w(10 Gyr) | 0 (double underflow) |
| w(13.4 Gyr) | 0 (double underflow) |
| ℓ_coh = 1/k_m | 0.100 Mpc = R_vir (L*) |
| m₀/k_m | 2.25×10⁻⁸ (k_m dominates) |

## Handoff artifacts

- `sims/sim147_output/kernel_sim147.pkl` — kernel grid + formula string + unit conversions
- `sims/sim147_output/sim147_metadata.json` — τ_mem, gate, ℓ_coh, sweep table
- `sims/sim147_output/sim147_kernel_decay.pdf` — plots

See RESULT.md for full derivation and SIM148/149 loading instructions.
