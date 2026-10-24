# SIM156 — SPEC (pre-registered 2026-10-24, before any run)

**Question** (SCREENING_NOTES.md §5). Can a screened version of curvature memory give a DESI-sized late-time effect?
The test case is the open branch: Vainshtein-screened Cassini, with Lunar Laser Ranging (LLR) applied to the
unscreened time variation of G.

## Models (all with bare Λ by closure; background level)
| label | definition | extra parameters vs ΛCDM |
|---|---|---|
| ΛCDM | ω_b, ω_cdm, h | 0 |
| w0wa | + w0, wa (CPL) | 2 |
| TQ | minimally coupled thawing quintessence: F = ½, V = 3 v1 Ψ, Ψ = Ψ' = 0 at z = 1e5 | 1 (v1) |
| PT | **pinned-thawing CMSTG**: F = ½ − Λ0 δΨ² (family-A attractor at δΨ = 0, F* = ½), V = 3 v1 δΨ, δΨ = δΨ' = 0 at z = 1e5 | 2 (Λ0, v1) |

- PT is curvature memory with a late trigger. The restoring force 2Λ0δΨR pins the field at the GR point while R is large; the tilt v1 releases it as R falls.
- F* = ½ makes the early universe exactly standard, so BBN and the CMB are untouched; S = √(F0/F*).
- Cassini is assumed screened by a cubic Galileon term with small c3 (α_B ≲ 1e-2; its background effect is neglected).

## Data (background level)
- **Compressed Planck 2018** (TT,TE,EE+lowE distance priors; Chen, Huang & Wang 2019, Table I, ΛCDM column):
  - R = 1.7502 ± 0.0046, l_A = 301.471 ± 0.090, ω_b = 0.02236 ± 0.00015;
  - correlations R–l_A +0.46, R–ω_b −0.66, l_A–ω_b −0.33;
  - z* from Hu & Sugiyama (as in the paper).
  - Valid for models whose early universe is standard, as all four are (the field is frozen until late times).
- **DESI DR2 BAO**: official Gaussian likelihood (`common/cmstg_obs.DesiBAO('DR2')`).
- **Pantheon+** (no SH0ES): as Cobaya `sn.pantheonplus` (z_HD > 0.01, STAT+SYS covariance, analytic marginalisation over M).

## Validation (must hold before any verdict; otherwise SIM156 is "inconclusive")
- **V1.** At the Planck 2018 TT,TE,EE+lowE ΛCDM means (ω_b 0.02236, ω_cdm 0.1202, h 0.6727), the pipeline's R and l_A are within 0.5σ of the priors. A constant multiplicative calibration of r_s and D_M is allowed and must be documented; it is then fixed for all models.
- **V2.** The Pantheon+ χ² of the SIM155 D1 ΛCDM best fit agrees with SIM155's value (1405.58) to within 1.0.
- **V3.** Δχ²(w0wa vs ΛCDM) is within ±3 of SIM155's full-likelihood value on D1 (+9.19), i.e. in [6.2, 12.2].

## Method
- Best fits by `common/fitcore.minimize` from 3 seeded starts per model (best reported).
- For PT also a grid in (Λ0, v1) to map where |α_M(0)| ≤ 0.002.

## Criteria for PT (the screened curvature-memory model)
PT **passes** if its best fit within |α_M(z=0)| ≤ 0.002 has:
- **(a)** Δχ² vs ΛCDM ≥ 6, and
- **(b)** passes BBN (C3: Y_p and D/H pulls ≤ 2 with S = √(F0/F*)).

**Always reported:**
- Δχ² and ΔAIC vs ΛCDM, TQ and w0wa ("competitive" = within 3 of w0wa);
- the unconstrained PT best fit and its α_M(0), i.e. how much LLR costs;
- α_M is d ln F/d ln a at z = 0. LLR (Burrage & Dombrowski 2020): |α_M| ≲ 0.002.

## Prediction stated in advance (SCREENING_NOTES §3)
Within |α_M| ≤ 0.002, PT does no better than TQ, and neither reaches w0wa.

## Caveats (known before running)
- **Background-level test.** It overstated gains in SIM154; the compressed-CMB and validation checks are there to limit that. A PT survivor would go to a full hi_class fit before any claim.
- **SN Ia luminosity and G.** If G varies (unscreened), the absolute magnitude of SNe may drift. With |α_M| ≤ 0.002 that drift is ≲ 0.01 mag over z < 1, and it is neglected.
