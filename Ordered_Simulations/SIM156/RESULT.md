# SIM156 — RESULT (2026-10-24)

**Verdict (pre-registered criteria, SPEC.md): FAIL**, as predicted in advance in SCREENING_NOTES §3.

- **The model.** "Pinned-thawing" CMSTG is curvature memory with a late trigger: the curvature R holds the field at the GR point early, and a tilt releases it at late times. Cassini is Vainshtein-screened.
- **Within LLR it is just quintessence.** Within the Lunar Laser Ranging bound on the running of the Planck mass (|α_M(0)| ≤ 0.002, which Vainshtein cannot screen), it improves on ΛCDM by only **Δχ² = +2.77**. That is exactly what minimally coupled thawing quintessence achieves: the fit switches the coupling off.
- **Without LLR it still falls short.** Even ignoring LLR it reaches only +4.95, needing a Planck mass that runs 56× faster than LLR allows. Both are below the pass line (6) and below w0waCDM (+7.68).

## Results (background level: compressed Planck 2018 + DESI DR2 + Pantheon+)

| model | Δχ² vs ΛCDM | ΔAIC | α_M(0) | notes |
|---|---|---|---|---|
| ΛCDM | 0 (χ² = 1421.39) | 0 | — | H0 = 68.37 (SIM155 full fit: 68.30) |
| w0waCDM | **+7.68** | −3.68 | — | w0 = −0.852, wa = −0.524 (SIM155: −0.850, −0.552) |
| thawing quintessence (TQ) | +2.77 | −0.77 | 0 | linear potential, minimal coupling |
| **PT within LLR** | **+2.77** | +1.23 | −0.0019 | Λ0 = 0.007: the coupling is off; identical to TQ |
| PT unconstrained | +4.95 | −0.95 | **−0.111** (56× LLR) | Λ0 = 0.88, v1 = 0.43; G today 5% stronger than early on |

Full table: `results.md`; best fits: `fits.json`; grid: `grid.json`; strong-coupling refits: `refine.json`;
figure: `figures/sim156_screened.png`.

## Pre-registered criteria
- **(a)** Δχ² ≥ 6 vs ΛCDM within |α_M(0)| ≤ 0.002: **FAIL** (+2.77).
- **(b)** BBN at that point: pass (Y_p pull +1.23, D/H −0.40; S − 1 = −2.5e-4).
- **Competitive with w0waCDM** (within 3): no, it trails by 4.9.
- **Prediction stated in advance** ("PT ≈ TQ < w0wa within LLR"): confirmed. PT − TQ = +0.00.

## Validation (all passed before the verdict)
- **V1.** The distance priors are reproduced at the Planck 2018 means after one constant calibration each: K_R = 0.999988, i.e. none needed, and K_lA = 0.998020, the known 0.2% offset of the simple r_s integral.
- **V2.** The Pantheon+ χ² at the SIM155 D1 ΛCDM best fit is 1405.58, identical to SIM155's hi_class + Cobaya value. DESI χ² there: 12.91 here vs 12.66.
- **V3.** Δχ²(w0wa vs ΛCDM) = +7.68, in the required [6.2, 12.2] around SIM155's full-likelihood +9.19.
- **Consistency.** TQ with zero tilt reproduces ΛCDM to 0.01 in χ².

## What the grid shows (left panel of the figure)
- **The coupling helps only by running G.** Improvements beyond thawing quintessence (+2.8) appear only where the Planck mass runs fast today: |α_M(0)| ≳ 0.01 (5× LLR) for +3…+4, and ≈ 0.1 (55×) for the maximum +4.9.
- **LLR-allowed points cannot exceed quintessence.** Every point that passes LLR has Δχ² ≤ 2.8. Strong pinning (Λ0 ≳ 3) keeps the field near GR (ΛCDM-like). Weak pinning (Λ0 ≲ 0.01) is quintessence.
- **The analytic picture is confirmed.** The CMSTG ingredient helps only by changing G over time today, which is exactly what LLR forbids (SCREENING_NOTES §3).
- **The free fits missed the strong-coupling region.** They started near Λ0 = 0.1 and converged to the weak-coupling optimum. The grid revealed the strong-coupling region, and three extra starts there (`refine.py`) give the unconstrained optimum +4.95 (+4.84 and +4.94 from the other two starts, another flat valley). Reporting the unconstrained optimum was required by the SPEC.

## Caveats
- **Background level.** The compressed CMB prior is valid here because the early universe is standard in all four models. Validation V3 shows it slightly *understates* w0waCDM's advantage (7.7 vs 9.2). A similar understatement for PT would not change the verdict: even +4.95 × 9.2/7.7 ≈ 5.9 is below 6, and that point violates LLR by 56×.
- **Assumed Vainshtein screening.** Screening by a small cubic Galileon term is assumed, not simulated. Its background effect (α_B ≲ 1e-3) is neglected; its nonlinear screening was estimated in SCREENING_NOTES §2.
- **One model family.** Only a quadratic coupling around a maximum with a linear tilt was tested. SCREENING_NOTES §3 argues the LLR limit is generic for coupling-driven late-time effects. Other shapes could change the details, not the |α_M| ≤ 0.002 ceiling.
- **SN luminosity and G.** A G-dependent SN absolute magnitude was neglected: ≲ 0.01 mag within LLR. In the unconstrained, LLR-violating region it would matter.

## Conclusion for the screening idea
1. **Thin-shell screening** (the natural screened form of the CMSTG attractor) cannot change cosmology (the WHK no-go, derived explicitly in SCREENING_NOTES §1).
2. **Vainshtein screening** frees Cassini, but LLR still caps the Planck-mass running. Within that cap, curvature memory adds nothing to quintessence (this test).
3. Together with SIM155 (early curvature memory is BBN-limited to about +3): **in universally coupled F(Ψ)R form, curvature memory cannot produce the DESI-preferred dynamics.** The remaining options are different theories: a dark-sector-only coupling, screening of Ġ, or accepting tuned bands (SCREENING_NOTES §4).
