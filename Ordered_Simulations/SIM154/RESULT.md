# SIM154 — RESULT (2026-10-22)

**Verdict:** SURVIVORS FOUND. Two of 70 grid points pass all pre-registered criteria C1–C4 (SPEC.md).

| family | Λ0 | ξ | H0 | F0 | Δχ² vs ΛCDM (DR1 / DR2) | DR2 GOF | \|γ−1\| | Ġ/G [1/yr] | Y_p (pull) | D/H pull | ΔN_eff equiv. |
|---|---|---|---|---|---|---|---|---|---|---|---|
| A (subtractive) | 0.15 | 1/12 | 67.75 | 0.5116 | +6.98 / +17.72 | 1.27σ | 1.1e-7 | −1.9e-17 | 0.2505 (+1.8σ) | +0.6σ | 0.14 |
| B (exponential) | 0.15 | 1/12 | 67.79 | 0.5117 | +7.14 / +18.08 | 1.22σ | 1.3e-7 | −1.9e-17 | 0.2506 (+1.8σ) | +0.6σ | 0.14 |

## What the scan shows
- **BBN is the binding constraint.** The total change of G since BBN is set roughly by F* − ½ = ξ²/(4Λ0) (family A). Helium allows only ≈ 2.5% (Y_p pull ≤ 2σ).
  - Even standard BBN sits at +1.3σ in Y_p with these fits and data, so the headroom is small.
- The DESI improvement grows with the size of the late-time change in G, so DESI pushes up against the BBN edge. The survivors sit on that edge (Y_p pull +1.8σ).
- Several "CB" points (pass Solar System and BBN) improve on ΛCDM but miss C1's GOF < 2σ:
  - (0.2, 1/12): +10.8, GOF 2.19σ
  - (0.3, 1/8): +8.3, GOF 2.49σ
  - (0.1, 1/24): +11.4, GOF 2.11σ
- The attractor keeps every point with Λ0 ≥ 0.1 Solar-System safe (|γ−1| ≤ 2e-5), except the largest ξ in family B.
- Small Λ0 (0.02–0.05) never reaches the attractor by today, so Cassini fails.
- Families A and B behave almost identically near the attractor, as expected; they differ only away from Ψ*.

## Caveats (why this goes to Phase 3, not to a paper)
- Background-level test only: θ* + BAO with Planck ω_b, ω_m fixed. In this setup the DESI preference is overstated: w0waCDM gets Δχ² = +23.8 on DR2 here, against about 3σ in DESI's full analysis.
- G changes by about 2% between BBN and today, mostly before recombination. That affects the CMB damping tail and early ISW, which θ* does not test. A full Boltzmann treatment is required.
- BBN uses linear fitting formulas (Kneller & Steigman 2004). A BBN code (PRyMordial / AlterBBN) should confirm the survivors.
- Grid is coarse; the survivor region along the band ξ²/Λ0 ≈ const has not been mapped finely.
- Dark energy is still a bare Λ.

## Next (ROADMAP Phase 3)
Implement families A and B in hi_class and fit Planck + DESI DR2 + supernovae, comparing against ΛCDM and w0waCDM.
