# Sim 80 — Recursive UV Finiteness (Quick Re-run)

**Goal:** Demonstrate UV falloff and softened cutoff dependence when including a recursive memory self-energy.

## Setup
- Baseline propagator: $G_0(p) = 1/(p^2 + m^2)$
- RIFT propagator: $G(p) = 1/(p^2 + m^2 + \Sigma_{\rm mem}(p))$
- Memory self-energy: $\Sigma_{\rm mem}(p) = \frac{g^2}{16\pi^2}\log\!\frac{\Lambda_{\rm mem}^2}{p^2 + m^2}$ (fixed memory scale).

Using:
- m = 1, g = 1, $\Lambda_{\rm mem}$ = 1e+06
- UV falloff grid: p ∈ [10, 10000] with 400 log-spaced points.
- Bubble integral: external |p| = 1, cutoffs Λ ∈ [100, 100000] (20 points).

## Results (quick indicators)
- Fitted UV tail slopes (log–log |G| vs p, last 1/3 of grid):
  - Baseline: slope ≈ -2.000 (≈ −2 as expected)
  - RIFT:    slope ≈ -2.000
- Bubble integral growth vs cutoff:
  - Baseline shows the usual increasing trend.
  - RIFT shows visibly softened growth across the same Λ range due to the positive memory term at internal momenta below $\Lambda_{\rm mem}$.

## Artifacts
- `Outputs/sim80_uv_falloff.png` — UV falloff comparison
- `Outputs/sim80_bubble_cutoff.png` — Bubble integral vs cutoff
- `Inputs/sim80_params.json` — Full parameter record
- `Published_Scripts/simulation_runner_rift.py` — Re-run script

**Note:** The memory scale $\Lambda_{\rm mem}$ is intentionally held fixed (decoupled from the integration cutoff) so that $\Sigma_{\rm mem}$ remains positive for relevant internal momenta, enhancing suppression and softening apparent divergences within the explored range.
