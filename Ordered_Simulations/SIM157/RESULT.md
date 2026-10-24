# SIM157 — RESULT (2026-10-24)

**Verdict (pre-registered criteria, SPEC.md): FAIL** on criterion (a).

- **The model.** "Dark curvature memory" puts CMSTG's coupling on dark matter only. Local gravity is satisfied by construction.
- **Right signature, wrong shape.** It reproduces the *qualitative* DESI signature, an apparent phantom crossing in the right direction. But the shape is wrong, and it improves the background fit by only **Δχ² = +2.88**, against a pass line of 6. That is +0.11 over plain thawing quintessence; w0waCDM gets +7.68.
- **Growth gets worse.** It raises S8 by 6%: the dark fifth force wins over the mass-change drag.

## Joint best fits (compressed Planck 2018 + DESI DR2 + Pantheon+ + eBOSS DR16 fσ8)

| model | Δχ² joint | Δχ² background | Δχ² fσ8 | ΔAIC (joint) | S8 | dark-matter mass today |
|---|---|---|---|---|---|---|
| ΛCDM | 0 (χ² = 1427.79) | 0 | 0 | 0 | 0.809 | — |
| w0waCDM | **+8.01** | **+7.68** | +0.33 | −4.01 | 0.829 | — |
| thawing quintessence | +2.67 | +2.77 | −0.10 | −0.67 | 0.815 | — |
| **DPT (dark curvature memory)** | **+3.15** | **+2.88** | +0.26 | +0.85 | **0.857** | **+6.8%** |

DPT best fit: Λ0 = 3.07, v1 = 0.30. Full table: `results.md`; fits: `fits.json`; 63-point grid: `grid.json`;
figure: `figures/sim157_dark.png`.

## Pre-registered criteria
- **(a)** Δχ²_background ≥ 6: **FAIL** (+2.88).
- **(b)** Growth no worse than ΛCDM: pass (fσ8 χ² better by 0.26).
- **(c)** Coupling off before recombination: pass (|β(z ≥ 1000)| = 9e-10).
- **Competitive with w0waCDM** (within 3): no, it trails by 4.8.
- **Expectation stated in advance** ("TQ < DPT < w0wa"): formally true but marginal, +0.11 over TQ on the background data.

## What the model does: the right signature, the wrong shape
Apparent dark-energy equation of state for an observer who assumes ordinary dark matter (right panel of the figure):

| z | 0 | 0.3 | 0.5 | 1.0 | 1.5 |
|---|---|---|---|---|---|
| DPT best fit | −0.96 | −1.02 | −1.06 | −1.22 | −1.57 |
| w0waCDM best fit (w0 = −0.85, wa = −0.54) | −0.85 | −0.97 | −1.03 | −1.12 | −1.17 |

- **The crossing is genuine.** Dark matter gaining mass at late times does mimic a phantom crossing (w > −1 today, w < −1 in the past), as Khoury, Lin & Trodden (2025) found for their dark-axion model.
- **But the shape is wrong.** In CMSTG's form the crossing is too weak today (w0 + 1 = 0.04 vs 0.15) and too steep at z ≳ 1, where the apparent dark-energy density heads toward zero. The data prefer w0waCDM's gentler shape.

## Why: the attractor that protects the early universe also weakens the late-time effect
- **The coupling starts at zero.** Around CMSTG's attractor the dark-matter mass factor is quadratic: A ≈ 1 + Λ0δφ², so the coupling β = 2Λ0δφ vanishes where the field starts. That is what keeps the early universe standard: criterion (c) passes by eight orders of magnitude.
- **So the mass change is second order and late.** It goes as δφ², which compresses the apparent phantom phase into z ≲ 1.5.
- **Forcing it costs growth.** A large mass change needs strong coupling, and the fifth force between dark-matter particles is 2β² = 8Λ0·(A − 1). At the best fit it is ≈ 2× gravity today (β ≈ 1.0).
  - That raises S8 by 6%, and at stronger coupling far more (S8 up to 1.6 in the grid; middle panels).
  - The growth data then cut off the strong-coupling region. The joint best fit balances background gain against growth cost.
- **The S8 direction is opposite to Khoury–Lin–Trodden.** Their model lowers S8 through a delayed matter-radiation equality: dark matter lighter *since early times*. CMSTG's attractor forbids exactly that, since the coupling must vanish before recombination. So here the fifth force dominates and S8 goes up.

## Validation
- **V1–V3** (inherited from SIM156): passed.
- **V4.** DPT with v1 = 0 is identical to ΛCDM, and DPT with Λ0 → 0 is identical to TQ, to 1e-4 in χ². H(z) agrees with analytic ΛCDM to < 1e-9.
- **V5.** ΛCDM fσ8 matches hi_class's at z = 0, 0.5, 1.0 to 0.1%.

**Method fixes made before any result** (see RESEARCH_LOG):
- H(z) from the dense ODE solution instead of a table;
- a guard at the singular point of A (first written with a wrong threshold, caught in testing);
- grid points where every cosmology is unphysical are now recorded instead of crashing the run.

**After the fits:** the printed w_eff diagnostics had a z-ordering bug. It was fixed and the values recomputed; the criteria were unaffected (they do not use w_eff).

## Caveats
- **Background + linear growth**, not a Boltzmann treatment of the coupled perturbations. Since the verdict is FAIL, the full treatment (ROADMAP 9.3) is not needed.
- **Light-field fifth force.** At the best fit m/H ~ √(6Λ0Ω_c) ≈ 3, so the force acts on scales below ~1.5 Gpc, which includes the σ8 scale. The approximation is adequate there.
- **One coupling shape.** Only the family-A attractor shape A = (1 − 2Λ0δφ²)^{-1/2} with a linear tilt was tested. A coupling linear in δφ (non-attractor) behaves differently, but it would not be CMSTG's attractor structure, and it is active early unless something else pins it.

## Conclusion for Phase 9 (options 1 and 2)
- **Option 2** (hide Ġ with universal coupling) is closed in theory (SCREENING_NOTES §6).
- **Option 1** (dark-sector coupling) works qualitatively: local gravity is automatic, and the apparent phantom crossing goes in the right direction. Quantitatively it adds essentially nothing over thawing quintessence and makes S8 worse.
- **Across SIM155–157**, every tested form of CMSTG curvature memory tops out around Δχ² ≈ +3 against w0waCDM's +8 to +9. The form is set by where the coupling acts (universal or dark) and when (early or late). The limit comes from BBN (early, universal), LLR (late, universal) or growth (late, dark).
