# CMSTG research log — post-audit roadmap

A detailed, chronological record of the work done under `ROADMAP.md` (the plan and its pre-registered criteria):
every step, problem, fix, decision and number. The audit itself is recorded in `AUDIT_NOTES.md`.

Related files:
- `THEORY_NOTES.md`: Phase 4 theory
- `PAPER_CORRECTIONS.md`: Phase 6
- `hi_class_model/README.md`: the hi_class model, its settings and its validation

---

## Environment (as of 2026-10-23)
- **Machine:** 12 cores, 78 GB RAM, about 25 GB free disk.
- **Python venv `/home/aion/cmstg-venv`** (outside the repo): Python 3.10, numpy 2.2.6, scipy 1.15.3, matplotlib, emcee 3.1.6, pytest 9.1.1, cobaya 3.6.2, cython 3.3.0, Py-BOBYQA, iminuit 2.33.0, hiclassy 3.4.0.1 (built from the patched source).
- **hi_class:** `/home/aion/cmstg-codes/hi_class`, upstream v3.4.0.1 (`ece81a57`) plus local commit `6861aea5` on branch `cmstg-attractor`. The patch is mirrored in `hi_class_model/cmstg_attractor.patch`.
- **Cobaya data:** `/home/aion/cobaya_packages` (250 MB):
  - Planck 2018 native low-ℓ TT and EE, plik-lite TTTEEE, lensing
  - DESI DR1 and DR2 BAO
  - Pantheon+, Union3, DES-Y5
- **Run directory for hi_class validation:** `/home/aion/cmstg-codes/runs` (the key `.ini` files are copied into `hi_class_model/`).
- **CI:** `.github/workflows/tests.yml` runs `pytest tests` on every push (green since commit `f2a1076`).

## Commits on `audit-fixes` for this phase

| commit | content |
|---|---|
| `f2a1076` | ROADMAP locked; 16 tests; CI |
| `2537370` | SIM154 SPEC (pre-registered) |
| `fc93dff` | SIM154 results (2 survivors) |
| `9dde0ef` | cmstg_attractor model in hi_class, validated |
| `566a9a2` | Phase 3 pre-registration (data, models, criteria) |
| `6c7c1ca` | Cobaya adapter validated; SIM155 v1 launched |
| `d709e27` | THEORY_NOTES, PAPER_CORRECTIONS, survivor predictions |
| `d1f9cc0` | SIM155 summary script |
| `00494b6` | turning-point c_s² tolerance fix; SIM155 v2 driver |

---

## Chronology

### 2026-10-22 — Phase 1: lock in the pipeline
- `tests/` has 16 pytest checks. All pass locally in about 8 s and on GitHub Actions. They cover:
  - DESI DR1/DR2 BAO-only fits reproduce the published Ω_m and h·r_d within 0.25σ
  - H0 from θ* matches Planck within 1σ
  - r_s calibration within 0.2%, plus a regression test against the truncated-r_s bug
  - exact GR limit of both solvers (1e-10)
  - constraint-ε equals Raychaudhuri ε (1e-8)
  - frozen-F0 degeneracy
  - Brans–Dicke and Cassini limits
  - statistics helpers

### 2026-10-22 — Phase 2: SIM154 (attractor scan + local gravity + BBN)
- SPEC written and committed **before** `run.py`.
- Grid: 2 families × 7 Λ0 × 5 ξ = 70 points, θ*-matched, DESI DR1 + DR2.
- **Result: 2 survivors**, family A and B at Λ0 = 0.15, ξ = 1/12:
  - Δχ² vs ΛCDM +7.0 / +17.7 (A) and +7.1 / +18.1 (B), DR1 / DR2
  - |γ−1| ≈ 1e-7; Ġ/G ≈ −2e-17 /yr
  - Y_p = 0.2505 (+1.8σ, against a 2σ cut)
- BBN is the binding constraint. The total G change since BBN is ≈ ξ²/(4Λ0) and must stay ≲ 2.5%.
- Log: `audit_logs/SIM154.log`. Figure: `Ordered_Simulations/SIM154/figures/sim154_scan.png`.

### 2026-10-22/23 — Phase 3.3: the model in hi_class

Implementation: covariant Horndeski model with G2 = X − V0 and G4 = F(φ), V0 tuned by hi_class's closure shooting. Three edits:
- `include/background.h` (enum)
- `gravity_smg/gravity_models_smg.c` (parser, G-functions, initial conditions, stdout)
- `source/background.c` (attractor-IC refinement)

Problems met and how each was solved:

1. **GR-limit run crashed in perturbations ("singular matrix").**
   - Cause: with zero scalar kinetic energy, hi_class's kinetic matrix D = 0.
   - Fix: `kineticity_safe_smg` regulator.
   - Check: κ = 1e-4 vs 1e-6 gave identical GR-limit results, and the survivor's C_ℓ converge (κ = 1e-5 vs 1e-6: ≤ 1e-5).
2. **GR-limit C_ℓ differ from ΛCDM by 2e-4 (TT), 4e-4 (EE), 9e-4 (TE) of peak.**
   - Not the regulator: identical at κ = 1e-4 and 1e-6.
   - Not hi_class's precision override: `perturbations_sampling_stepsize` 0.05 was tested.
   - Conclusion: hi_class scalar-sector numerics. Estimated Δχ² ≈ 0.1 on Planck. **Open: measure it directly in Cobaya.**
3. **Survivor background took more than 5 minutes even with no closure loop.**
   - Cause: started at rest, the field relaxes to a radiation-era attractor; the stiff solver resolved that at machine-noise level. hi_class documents the same problem for Brans–Dicke.
   - Fix: analytic attractor initial condition φ' = 3a²ρ_m F′ / [4F aH(1 + 3F′²/F)], refined by hi_class's own secant step (extended to our model), with `a_ini_over_a_today_default = 1e-10` (as in hi_class's BD example).
   - Background now about 1.4 s.
4. **hi_class's c_s² test failed with min c_s² = −1.2e-7 at a = 0.12.**
   - Diagnosis: the field's velocity changes sign at z ≈ 7.36 (it turns around near the attractor). There D → 0 and the c_s² formula is 0/0.
   - Physically c_s² = 1 for this theory class. Tolerance set to `cs2_safe_smg = 1e-5`, later 1e-2 (see SIM155 v1).
5. **A consistency check needed `a_ini_test_qs_smg ≥ a_ini`.** Set to 1e-10, with `pert_ic_ini_z_ref_smg = 1e9`.

Validation of the finished model:

| test | result |
|---|---|
| GR-limit background | identical 100θ_s = 1.043001; z_drag to 6e-9 |
| survivor background vs our `st_solver` (both started at z = 1e8) | H/H_ΛCDM to 2e-5; F(z) to 1e-5 for z = 0 … 1090 |
| survivor C_ℓ vs ΛCDM at equal H0 (physical effect) | 2.5e-3 (TT), 6.6e-3 (EE), 1.1e-2 (TE) |

### 2026-10-23 — Phase 3.1–3.2: Cobaya integration
- **Adapter `common/cobaya_hiclass.py`.** Cobaya's CLASS theory loads a module named `classy`; this hi_class version calls its wrapper `hiclassy`, so the adapter swaps in a shim with `Class = HiClass`. It also:
  - turns the sampled `cmstg_L0` and `cmstg_xi` into hi_class's `parameters_smg`;
  - adds the Planck neutrino setup (N_ur = 2.0328, one 0.06 eV state);
  - sets the pre-registered YHe formula.
- Two bugs on the way:
  - `import cobaya.theories.classy.classy` returned the *class*, because the package re-exports it under the module's name. Fixed with `importlib.import_module`.
  - The helper parameter `logA` was passed through to CLASS, which rejects unknown inputs. Fixed with `drop: True`.
- **Validation at the Planck 2018 ΛCDM point:** low-ℓ TT 23.45, low-ℓ EE 396.12, plik-lite 583.16, lensing 8.86. These match Planck's published values. DESI DR2 gave 29.96, Pantheon+ 1403.74.
- **Timing per evaluation:** ΛCDM 5.8 s (1 thread) / 1.7 s (4 threads); CMSTG 14.5 s / 4.0 s.

### 2026-10-23 — Phase 3.4 pre-registration (commit `566a9a2`, before any fit)
- Data D1/D2/D3, models, method and pass criteria (a)–(c) are fixed in ROADMAP.md and copied to SIM155/SPEC.md.

### 2026-10-23 — Phase 5 predictions at the SIM154 survivor (background level)
- **GW sirens:** |d_GW/d_EM − 1| < 1e-6 for z < 10. The field is frozen at the attractor since z ≈ 7, so sirens can't test it.
- **Growth:** fσ8 is −2.5% (z = 0), −1.3% (z = 0.5), −0.6% (z = 1) vs ΛCDM; σ8 is −0.55%. A DESI peculiar-velocity/RSD test; points toward the S8 tension. To redo at the SIM155 best fit.

### 2026-10-23 — SIM155 v1 (Cobaya `minimize`, default BOBYQA), stopped after D1

v1 best-fit χ² (Planck + DESI DR2 + Pantheon+, D1):

| model | start 1 | start 2 | start 3 |
|---|---|---|---|
| ΛCDM | 2433.18 | 2432.49 | 2434.84 |
| w0waCDM | 2429.95 | 2424.41 | 2424.61 |
| CMSTG A | crashed | crashed | crashed |

- The w0wa best fit (start 2) has w0 = −0.72, wa = −0.92, H0 = 66.7. Best-to-best: w0wa improves on ΛCDM by Δχ² ≈ 8.1, before polishing.
- **Problem 1: early termination.** Starting points disagree by up to 2.3 (ΛCDM) and 5.5 (w0wa) in χ².
  - The χ² noise floor, measured around the ΛCDM best point, is only about 0.005 (Δχ² for H0 ± 1e-4 … 1e-3).
  - So the disagreement is the minimizer quitting early (rhoend = 0.05 and slow-progress rules along the H0–ω_cdm degeneracy), not noise.
- **Problem 2: every CMSTG start crashed after about 20 evaluations.**
  - Random points near the survivor failed hi_class's c_s² test in 4 of 14 cases.
  - In every failing case the most negative c_s² sat on a turning point of the field (φ' changes sign), with values −1.6e-5 … −3e-4. Examples: min at z = 0.31 against a flip at 0.32; 5.17 against 5.22; 3.86 against 3.90; 0.19 against 0.19.
  - Cobaya turns a failed evaluation into −inf, and BOBYQA's interpolation then becomes singular.
- v1 outputs are kept in `Ordered_Simulations/SIM155/fits/` for the record.

### 2026-10-23 — Fixes and SIM155 v2 (method only; criteria unchanged)
- **`cs2_safe_smg` 1e-5 → 1e-2.** That is more than 30× the observed turning-point noise and still catches a real instability. The ghost and tensor tests remain active.
- **New driver `SIM155/fit.py`:**
  - Stage 1: Py-BOBYQA scaled to the prior box (rhobeg 0.05, rhoend 2e-5, at most 2500 evaluations); failures return a finite penalty of 1e6.
  - Stage 2: MIGRAD (iminuit), repeated until two passes agree to 0.01 in −log(posterior).
  - Every evaluation goes through Cobaya's model (theory + likelihoods + priors).
- All models are rerun with v2 so the comparison is uniform.
- Status at the time of writing: a ΛCDM D1 test fit is running, plus an end-to-end CMSTG smoke test (tiny budget, output in `fits_v2_smoke/`, not a result).

---

## Decisions and deviations (with reasons)
1. **Helium in the fits uses F* = ½ + ξ²/(4Λ0)** (the attractor value) instead of the pre-registered F0 (today's value), because F0 is only known after solving the background.
   - At the survivor F0 = 0.511574 and F* = 0.511574.
   - `summarize.py` recomputes C3 with the true F0 at every best fit, so the criterion itself is applied as registered.
2. **v1 queue:** `best_of` was cut from Cobaya's default 2 to 1 for jobs not yet started (three seeds still meet "≥ 3 starts"). v1 was then superseded by v2.
3. **c_s² tolerance and the v2 driver:** see above. Both are numerical-method changes made after seeing failures, documented before any v2 result was examined.

## Open issues
- Measure the GR-limit residual (issue 2 above) as a Δχ² on the full likelihood.
- BBN uses linear fitting formulas (Kneller & Steigman 2004). Confirm the survivors with a BBN code.
- The CMB data are Planck 2018 plik-lite, low-ℓ and lensing. DESI's published w0wa significance uses different CMB likelihoods (Planck PR4 CamSpec and ACT lensing), so numbers will differ somewhat.
- Paper III (UV) and Paper IV (galactic) are not audited.
