# SIM157 — SPEC (pre-registered 2026-10-24, before any run)

**Question** (DARK_COUPLING_NOTES.md). If CMSTG's coupling acts only on dark matter ("dark curvature memory"),
can it reproduce the DESI-preferred dynamics? Local gravity is automatically satisfied: baryons see GR with a constant G.
And what does it do to structure growth?

## Models (background + linear growth; bare Λ by closure)
| label | definition | extra parameters vs ΛCDM |
|---|---|---|
| ΛCDM | ω_b, ω_c, h | 0 |
| w0wa | + w0, wa (CPL; growth as in GR with smooth DE) | 2 |
| TQ | thawing quintessence: V = 3 v1 φ, minimally coupled, φ = φ' = 0 at z = 3000 | 1 |
| **DPT** | dark pinned thawing: A(φ) = (1 − 2Λ0δφ²)^{-1/2} on dark matter only, V = 3 v1 δφ, δφ = δφ' = 0 at z = 3000 | 2 (Λ0, v1) |

ω_c is the early-time dark-matter density (the one the CMB measures); today's is ω_c A(φ0). Equations are in DARK_COUPLING_NOTES §2.

## Data
- **Background (as SIM156, same code):**
  - compressed Planck 2018 (R, l_A, ω_b; Chen, Huang & Wang 2019; R and z* evaluated with the early-time ω_m);
  - DESI DR2 BAO (r_d from Aubourg et al. with early ω_b, ω_c);
  - Pantheon+.
- **Growth:** eBOSS DR16 RSD-only fσ8 (Alam et al. 2021, Table III; 6 points; treated as independent, as in that table).
  - σ8(z) = σ8_ref · δ_m(z)/δ_m^ref(0), with σ8_ref = 0.8076 (the SIM155 hi_class ΛCDM D1 best fit).
  - δ_m is normalised to a common value at z = 100, where every model equals ΛCDM.

## Validation (before any verdict)
- **V1–V3.** Inherited from SIM156 (same likelihood code; all passed).
- **V4.** DPT with v1 = 0 reproduces ΛCDM. DPT with Λ0 = 0 reproduces TQ (to 0.01 in χ²).
- **V5.** The ΛCDM growth f(z) from this solver matches hi_class's ΛCDM best fit (SIM155 growth.json: fσ8 at z = 0, 0.5, 1.0) to 2%.

## Criteria (DPT passes if all hold at its best fit)
- **(a)** Δχ²_background vs ΛCDM ≥ 6.
- **(b)** Growth not worse than ΛCDM: χ²_RSD(DPT) − χ²_RSD(ΛCDM) ≤ 4.
- **(c)** Coupling inactive before recombination: |β(z ≥ 1000)| < 1e-3. This also keeps the compressed CMB prior valid.

**Always reported:**
- Δχ² vs w0wa ("competitive" = within 3) and vs TQ (what the dark coupling adds);
- S8/S8_ΛCDM;
- the dark-matter mass change A(φ0) − 1;
- w_eff(z) and whether it crosses −1.

## Expectation stated in advance (not a criterion)
- DPT beats TQ, unlike SIM156's universally coupled PT, because the dark coupling can mimic phantom crossing without an LLR limit.
- Its Δχ² is expected between TQ (+2.8) and w0wa (+7.7).
- S8 is expected slightly lower (mass-change drag) unless the fifth force dominates.

## Caveats known before running
- **Background + linear growth only**, not a Boltzmann treatment of the coupled perturbations. A DPT that passes would go to a full hi_class or CLASS implementation before any claim.
- **Light-field (all-scale) fifth force.** This overstates the fifth force for Λ0 ≳ 1.
- **σ8 normalisation** neglects ~1% transfer-function differences between best fits.
- **Correlations.** The eBOSS fσ8 points are treated as uncorrelated (BOSS z = 0.38 and 0.51 overlap).
