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
- **Python venv `/home/aion/cmstg-venv`** (outside the repo): Python 3.10, numpy 2.2.6, scipy 1.15.3, matplotlib, emcee 3.1.6, pytest 9.1.1, cobaya 3.6.2, cython 3.3.0, Py-BOBYQA 1.5.0, iminuit 2.33.0 (used only by the abandoned v2 driver), hiclassy 3.4.0.1 (built from the patched source).
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
| `4421053` | this research log; SIM155 v1 records; driver budget overrides |
| `c0e427f` | SIM155 v1 run configurations (yaml) kept; binary/temporary artefacts ignored |
| (next) | SIM155 v3: tested optimiser `common/fitcore.py`, likelihood diagnostics, input-stage error fix; 27 fits launched |

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
    **Correction (later the same day):** that probe used steps too small to cross hi_class's discrete sampling thresholds.
    With steps of 0.004 proposal widths the noise is 0.07 RMS in χ² (see the v3 entry below). The conclusion still
    holds: disagreements of 2–5 in χ² are 30–70× the noise, so they are early termination.
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

### 2026-10-23 — v2 abandoned (too slow); v3 optimiser built and tested; 27 fits launched (method only; criteria unchanged)
- **Both v2 tests died at their time limits and produced no result:**
  - The CMSTG smoke test (MAXFUN 30, NCALL 30) was killed by its 1200 s `timeout`. MIGRAD's Hessian evaluations are not counted by `ncall`, so its budget does not bound its run time.
  - The ΛCDM D1 production test (`MAXT=1 OMP_NUM_THREADS=4 timeout 3000 python fit.py lcdm D1 1`) was killed after 50 min (exit 143). BOBYQA at rhoend 2e-5 with ≤ 2500 evaluations (1.7 s each) was still in stage 1, and v2 wrote output only at the end.
  - At that rate 27 fits would take days. No v2 output exists, so nothing is lost by replacing it.
- **v3 = `common/fitcore.py` (an optimiser module), driven by `SIM155/fit.py`:**
  1. Coarse: Py-BOBYQA inside the prior box, axes scaled to [0, 1], ≤ 400 evaluations.
  2. Curvature: central-difference Hessian at the current point, in units of the proposal widths. Eigenvalues are folded to |w| and floored at 1e-3, so neither noise nor non-convexity can stretch a direction by more than about 32 widths.
  3. Fine: Py-BOBYQA in whitened coordinates, where −log(post) ≈ ½|u|². Steps 2–3 repeat (at most 3 rounds) until a round gains less than 0.1.
  - The prior box is a continuous wall, f(clip x) + 50·|(x − clip x)/proposal|², so the rotated stage can approach an edge without a jump. This matters for CMSTG's GR edge (ξ = 0).
  - A failed evaluation is retried up to 3 times at a 1e-3-width jitter. If all fail it scores best-so-far + 50; a Hessian stencil that fails is dropped.
  - Output is written after every stage (`fits_v3/<run>.json`), including the inverse-Hessian (Laplace) errors.
- **First v3 launch (19:08 EDT), stopped after about 5 min with no output:**
  - The first draft regularised the Hessian in raw units, where eigenvalues span about 10⁷, so the floor was meaningless.
  - It also used a hard 1e6 penalty at the box edges.
  - Both were fixed as above.
- **Diagnostics on the real D1 likelihood (`SIM155/diagnostics.py`; outputs `diag_noise.json`, `diag_fails_*.json`):**
  - *Noise.* −log(post) was evaluated at 11 points spaced 0.004 proposal widths apart, along H0 and along ω_cdm, around the ΛCDM reference point.
    - Residuals from a quadratic fit are random ±0.04 jumps between neighbouring points.
    - RMS = 0.032 (H0) and 0.037 (ω_cdm) in −log(post), i.e. **0.07 in χ²**.
    - This is hi_class's discrete k/ℓ sampling at default precision. It sets the precision of any best-fit χ² at about ±0.1, far below the pre-registered thresholds (6 and 3).
    - Breakdown by likelihood (rerun with per-likelihood χ²): plik-lite TTTEEE 0.064 / 0.074 (H0 / ω_cdm direction); lensing 0.0003; low-ℓ TT ≤ 0.0002; low-ℓ EE, DESI DR2 and Pantheon+ 0.
    - So the noise is entirely high-ℓ C_ℓ sampling, as expected. A future MCMC could reduce it with tighter hi_class C_ℓ/k-sampling precision.
  - *Failures.* 0 of 40 CMSTG-A points drawn from the fit's `ref` distribution fail with `cs2_safe_smg = 1e-2` (v1 failed 4/14 at 1e-5).
    - −log(post) ranges 1205.5–1719 (median 1293).
    - The best random point is Λ0 = 0.26, ξ = 0.064 (ξ²/4Λ0 = 0.004).
- **Synthetic validation (`tests/test_fitcore.py`, in CI):** Gaussian −log(post) with correlations up to 0.995 and widths spanning five decades, at the measured noise.
  - Over 36 benchmark problems (n = 7 and 9, interior and on-edge minima, noise 0.035–0.07, 0–5% random failures) the final excess over the true minimum was:
    - typically 0.01–0.08;
    - worst 0.17 at the measured noise and 0.53 at twice the noise plus failures;
    - after 280–670 evaluations.
  - Before the retry rule existed, 5% failures left excesses of 2–12. The retry rule removed that.
- **hiclassy input-stage errors.** hiclassy raises `CosmoSevereError` for any failure in its input stage, which includes the Ω_smg = −1 closure shooting. Cobaya re-raises that error, so one extreme point would abort a whole fit, whereas it assigns zero likelihood to `CosmoComputationError`.
  - `cobaya_hiclass._HiClass` now converts input-stage failures to `CosmoComputationError`.
  - Set-up errors (unread parameters, unknown level) stay fatal.
  - `fitcore` raises if no evaluation in the coarse stage succeeds, so a broken set-up cannot pass silently.
- **Launch: 2026-10-23 19:22 EDT.** `SIM155/run_all_v3.sh` runs 27 fits (3 models × D1–D3 × seeds 1–3), 6 parallel jobs × 2 OpenMP threads.
  - Order: CMSTG D1 first, then ΛCDM D1, w0wa D1, CMSTG D2/D3, then the rest.
  - Expected cost: about 3 s per ΛCDM evaluation and 8 s per CMSTG evaluation; about 3.5 h in total.

### 2026-10-23 — Systematic check: χ² offset of the hi_class model in its exact GR limit (`SIM155/gr_limit_chi2.py`)
- **Question:** in the exact GR limit (ξ = 0, so φ ≡ 0 and F ≡ ½) the hi_class CMSTG model differs from plain ΛCDM by 2e-4 … 9e-4 of the C_ℓ peaks (scalar-sector numerics, hi_class_model/README.md). Does that bias the CMSTG-vs-ΛCDM Δχ²?
- **Method:** full D1 likelihood at the same parameters, CMSTG (ξ = 0) minus plain ΛCDM, at four cosmologies (Planck-like reference, then ±1σ moves in H0/ω_cdm, H0/n_s, ω_b/τ/logA).
- **Result:** Δχ² = +0.139, +0.174, −0.129, −0.114. **Mean +0.02 ± 0.08 (sd 0.16).**
  - At the reference point it is all plik-lite (+0.185) and low-ℓ TT (−0.048); DESI, Pantheon+ and low-ℓ EE are 0.
  - It is identical for Λ0 = 0.15 and 0.5 (Λ0 is irrelevant at ξ = 0) and continuous at ξ = 1e-4 (+0.140).
- **Conclusion:** no systematic offset. The scatter is evaluation noise; the difference of two noisy evaluations is expected to scatter by √2 × 0.07–0.1. CMSTG-vs-ΛCDM Δχ² values need no correction, only the ±0.1-level noise allowance. Output: `SIM155/gr_limit_chi2.json`.

### 2026-10-23 (evening) — SIM155 v3 interim results (D1), the CMSTG valley, and the profile likelihood
**Finished fits, D1 = Planck 2018 + DESI DR2 + Pantheon+**. −log(post) is the minimised quantity; χ²_eff = χ² + ((A_planck − 1)/0.0025)².

| model | −log(post), seeds 1/2/3 | spread | best χ² | best χ²_eff | Δχ² vs ΛCDM | ΔAIC vs ΛCDM | best-fit parameters |
|---|---|---|---|---|---|---|---|
| ΛCDM | 1200.929 / 1200.914 / 1200.939 | 0.025 | 2432.605 | 2432.785 | 0 | 0 | H0 = 68.30 ± 0.28, ω_cdm = 0.1179, ω_b = 0.02251 |
| w0waCDM | 1198.643 / 1198.604 / 1198.566 | 0.077 | 2423.415 | 2423.485 | +9.19 | −5.19 | w0 = −0.850, wa = −0.552, H0 = 67.68 |

- **Optimiser performance.** Each fit took 285–620 evaluations (20–46 min) with no failed evaluations. Round 1 gained 6.5, against a predicted 5.7; the last round gained ≈ 0.01, against a predicted ≈ 0.05.
- **Laplace errors match published posteriors.** The ΛCDM inverse-Hessian errors (σ(H0) = 0.28, σ(ω_cdm) = 6.2e-4, σ(n_s) = 0.0033, σ(τ) = 0.0071) match published DESI DR2 + CMB widths. This validates the curvature step.
- **v3 beats v1.** The v3 ΛCDM optimum (−log(post) 1200.914) is 0.20 better than v1's best (1201.112).
  - v1's lower likelihood-only χ² (2432.49 vs 2432.61) came from pulling A_planck 0.83σ off its prior, which cost 0.69 in prior χ² (χ²_eff 2433.18 vs 2432.79).
  - Hence χ²_eff is reported alongside χ², and each model's best seed is chosen by the minimised quantity.
- **w0waCDM agrees with DESI DR2.** w0 = −0.85, wa = −0.55, Δχ² = 9.2 for CMB + DESI DR2 + Pantheon+. DESI DR2 reports w0 = −0.838 ± 0.055, wa = −0.62^{+0.22}_{−0.19} and Δχ²_MAP ≈ −10.7 with different CMB likelihoods (PR4 CamSpec + ACT lensing).
- **Convergence flag.** `fit.py` flags convergence when the last round gains less than 0.05, but `fitcore` stops at tol = 0.1. So w0wa s1, which stopped on a 0.0999 gain, was labelled unconverged.
  - `summarize.py` recomputes the flag from the stored gains with 0.1.
  - `fit.py` is not edited mid-run, because queued jobs load it when they start.

**CMSTG-A on D1: a long, flat, curved valley.** The three seeds stop at nearly the same χ² at different points:

| seed | stage | −log(post) | χ² | Λ0 | ξ | ΔF = ξ²/(4Λ0) | φ* = ξ/(2Λ0) | H0 |
|---|---|---|---|---|---|---|---|---|
| 1 | round 2 (round 3 running) | 1197.660 | 2428.64 | 0.266 | 0.159 | 0.0237 | 0.30 | 67.17 |
| 2 | round 2 | 1198.044 | 2429.43 | 0.110 | 0.063 | 0.0090 | 0.29 | 68.31 |
| 3 | final ("converged") | 1198.025 | 2429.36 | 0.145 | 0.081 | 0.0114 | 0.28 | 68.14 |

- **Shape of the valley.** Along it the attractor field φ* ≈ 0.28–0.30 stays fixed while ΔF and H0 change; it is the F–H0 degeneracy found in the audit, partly broken by the CMB.
  - χ² falls slowly with ΔF: about 0.7 between ΔF = 0.011 and 0.024.
  - The Hessian has 1–2 negative eigenvalues in proposal units, the signature of a curved valley.
- **Seed 3 stopped early.** It "converged" (round gain 0.03) about 0.37 in −log(post) above seed 1, which was still improving. A local optimiser in a curved valley can stop on a small round gain far from the valley's minimum.
- **BBN.** C3 needs ΔF ≲ 0.012 (SIM154). So the best-fit location matters for criterion (b), and where in the valley a seed stops is set by the optimiser, not by the data.

**Decision (method only; criteria unchanged): profile likelihood in ΔF (`SIM155/profile.py`, `run_profile.sh`).**
- **Why.** This follows decision 5 (a seed spread > 0.5 in χ² triggers a re-fit). It was made after seeing the seed spread above, but before any profile result.
- **How.**
  - At fixed ΔF, Λ0 is sampled with its flat prior and ξ = 2√(Λ0 ΔF) is computed. The pre-registered box ξ ≤ 0.3 is kept by Λ0 ≤ 0.0225/ΔF; along the valley this caps ΔF at about 0.042.
  - Everything else is minimised with the same `fitcore` optimiser.
  - Each point starts from the ΛCDM D1 best fit, with Λ0 at the valley guess ΔF/0.28².
- **Grid (D1).** ΔF = 0.004, 0.008, 0.012, 0.016, 0.020, 0.024, 0.032, 0.040.
- **What it gives.**
  - min over ΔF of the profile, compared with the best free seed, is the CMSTG best fit used for criteria (a)–(b).
  - The profile also shows how much improvement survives where BBN allows (ΔF ≲ 0.012). That reading is reported as exploratory, outside the pre-registered verdict.
- **Scheduling.** Profile jobs fill CPU slots as the main fits finish, matched on `^…/python (fit|profile).py` so the xargs wrappers are not counted.
  - A first version of the runner matched the wrappers as well, so it would have waited for all 27 fits. It was stopped before starting anything.
- D2/D3 profiles will be run only if criterion (c) could matter, i.e. only if CMSTG passes (a) and (b) on D1.
- **Update, same evening (before any profile job had started).** Free seed 1 ended after its 3 rounds without converging (last gain 0.29).
  - Its end point: χ² 2428.10 (Δχ² = +4.51 vs ΛCDM), Λ0 = 0.288, ξ = 0.201, **ΔF = 0.035**, φ* = 0.35, H0 = 66.50, ω_b = 0.02227. Still drifting to larger ΔF and lower H0, and the valley bends (φ* rises with ΔF).
  - The grid now extends to ΔF = 0.048 and 0.056. The box ξ ≤ 0.3 ends the valley near ΔF ≈ 0.05.
  - The Λ0 start guess follows the bend: φ*(ΔF) = 0.28 + 4·max(0, ΔF − 0.018).
  - The runner was restarted with 10 points in the order 0.040, 0.048, 0.012, 0.032, 0.024, 0.056, 0.008, 0.016, 0.004, 0.020.

### 2026-10-24 — SIM155 complete: pre-registered verdict FAIL (Ordered_Simulations/SIM155/RESULT.md)
**27 free fits** (finished 23:58 EDT on 2026-10-23) **plus the 10-point D1 profile** (finished on 2026-10-24).

| data | ΛCDM χ² | w0waCDM Δχ² | CMSTG-A Δχ² (source) | CMSTG − w0wa |
|---|---|---|---|---|
| D1 | 2432.60 | +9.19 | **+4.65** (profile, ΔF = 0.040) | −4.54 |
| D2 | 1055.25 | +15.53 | +4.50 (free seed 1) | −11.02 |
| D3 | 2676.22 | +19.60 | +3.19 (free seed 3) | −16.41 |

- **Criteria.**
  - (a) fails: +4.65 < 6.
  - (b) fails: C2 passes (|γ−1| = 2.9e-8), but C3 fails at the best fit (Y_p +3.25σ, D/H +2.24σ; F0 = F* = 0.540).
  - (c) fails: D2 +4.50, D3 +3.19.
  - Not competitive with w0waCDM: χ² gap 4.5–16.4, against the pre-registered "within 3".
- **Profile (D1).** Δχ²: +2.4 (ΔF = 0.004), +3.0 (0.008), +3.0 (0.012), +3.5 (0.016), +3.8 (0.020), +4.1 (0.024), +4.5 (0.032), +4.65 (0.040), +3.5 (0.048), +3.9 (0.056).
  - The BBN edge (Y_p pull = 2σ) is at ΔF ≈ 0.015.
  - The best BBN-allowed point is free seed 3 (ΔF = 0.0113): **Δχ² = +3.25**, ΔAIC = +0.75. Exploratory, not a criterion.
- **Physics** (`growth.py`; χ² breakdown in RESULT.md).
  - At the best fit G at BBN is 8.0% stronger than today, 0.25% stronger at z = 1090, and equal to today for z ≲ 10. So the effect is early (between BBN and recombination).
  - Gains come from plik-lite (−3.35), low-ℓ EE (−1.13) and DESI (−1.86). Pantheon+ gets worse (+0.60); w0waCDM gains −2.59 there. CMSTG cannot mimic late-time dark energy.
  - S8 = 0.871 (+7.5% vs ΛCDM 0.811) at the best fit, and +2.5% at the BBN-allowed point. Both worsen the weak-lensing S8 tension.
- **SIM154's gain shrank.** Its +17.7 (DESI DR2, background level, fixed Planck ω's, against a fixed ΛCDM) shrank to +4.65 in the joint fit. ΛCDM re-adjusts and its DESI χ² drops from ≈ 30 to 12.7. w0waCDM likewise fell from +23.8 to +9.2.
- **Theory consequences** (THEORY_NOTES 4.5). There is a structural tension for Damour–Nordtvedt-class models of this type:
  - local gravity needs the field at the attractor today, while a late-time dark-energy-like signal needs it still moving;
  - fast convergence makes the effect early and BBN-limited;
  - escaping that needs screening, or a separate dark-energy sector.
- **ROADMAP consequences.** MCMC was pre-registered only for a model that passes, so it is not run. 6.4 ("one focused paper if a model survives") does not apply. 5.2 and 5.3 are done at the full-CMB best fits.
- **Runs not repeated.** D2/D3 CMSTG profiles (they cannot change the verdict); a family B fit (not pre-registered; SIM154 shows it is nearly identical to A near the attractor).

### 2026-10-24 — Phase 8: can curvature memory be screened? (SCREENING_NOTES.md, SIM156)
- **The request.** The user asked to explore screening ("can we explore that idea?").
- **Literature checked against the papers:**
  - Wang–Hui–Khoury 2012 (no-go for thin-shell screening);
  - Babichev–Deffayet–Esposito-Farèse 2011 (Vainshtein does not screen Ġ);
  - Burrage & Dombrowski 2020 (|α_M| ≲ 0.002 from LLR);
  - Creminelli et al. 2020 (|α_B| ≲ 1e-2 from GW-induced instabilities);
  - recent DESI + non-minimal coupling papers (arXiv:2606.14014 "very narrow bands"; arXiv:2609.10133; Wolf et al. 2025).
  - The Planck distance-prior table (Chen, Huang & Wang 2019) was read from the PDF, because an automated summary had the correlation signs wrong.
- **Theory map** (SCREENING_NOTES):
  - The family-A attractor plus a quartic potential is exactly a symmetron in the curvature. Screening the Sun and Milky Way needs Λ0 ≳ 8e4, and then ΔF ≲ 3e-6: no-go.
  - A genuine cubic Galileon screens Cassini, but LLR caps α_M.
  - The project's SIM142 "Galileon" was a rescaled kinetic term with no screening.
- **SIM156** (pre-registered and committed in ed37175 before running): pinned-thawing CMSTG on compressed Planck + DESI DR2 + Pantheon+.
  - Validation passed: V1 (K_R = 0.99999, K_lA = 0.99802); V2 (Pantheon+ χ² 1405.58, identical to SIM155); V3 (w0wa +7.68).
  - Within |α_M(0)| ≤ 0.002: +2.77, identical to thawing quintessence (the coupling switches off). FAIL (a), as predicted.
  - The 72-point grid showed a strong-coupling region the free fits had missed. Three extra starts there give the unconstrained optimum +4.95 at α_M(0) = −0.111 (56× LLR), Λ0 = 0.88.
- **Conclusion.** Universally coupled F(Ψ)R curvature memory cannot produce the DESI-preferred dynamics: early it is BBN-limited (SIM155), late it is limited by LLR (Vainshtein) or WHK (thin-shell). The remaining directions are different theories (ROADMAP 8.4).

---

## Decisions and deviations (with reasons)
1. **Helium in the fits uses F* = ½ + ξ²/(4Λ0)** (the attractor value) instead of the pre-registered F0 (today's value), because F0 is only known after solving the background.
   - At the survivor F0 = 0.511574 and F* = 0.511574.
   - `summarize.py` recomputes C3 with the true F0 at every best fit, so the criterion itself is applied as registered.
2. **v1 queue:** `best_of` was cut from Cobaya's default 2 to 1 for jobs not yet started (three seeds still meet "≥ 3 starts"). v1 was then superseded by v2.
3. **c_s² tolerance and the v2 driver:** see above. Both are numerical-method changes made after seeing failures, documented before any v2 result was examined.
4. **v3 driver replaces v2** (too slow; no v2 result was ever produced). Also a method-only change, documented before any v3 result: same data, models, priors, starts and criteria.
5. **Seed spread as the convergence check.** Three seeds per model; the best is reported, and the spread between seeds is reported next to it.
   - At the measured noise, a spread up to about 0.5 in χ² is expected.
   - A larger spread triggers a re-fit from the best point, logged here, before any criterion is applied.
6. **CMSTG best fit from a profile likelihood in ΔF = ξ²/(4Λ0)** (see the evening entry): the free CMSTG fits stop at different points of a flat, curved valley.
   - The best fit is the lower of (i) the best free seed and (ii) the minimum of the profile.
   - Same data, priors (including the box ξ ≤ 0.3) and criteria.

## Open issues
- ~~Measure the GR-limit residual (issue 2 above) as a Δχ² on the full likelihood.~~ Done 2026-10-23: no systematic offset (+0.02 ± 0.08); see the GR-limit entry in the chronology.
- BBN uses linear fitting formulas (Kneller & Steigman 2004). Confirm the survivors with a BBN code.
- The CMB data are Planck 2018 plik-lite, low-ℓ and lensing. DESI's published w0wa significance uses different CMB likelihoods (Planck PR4 CamSpec and ACT lensing), so numbers will differ somewhat.
- Paper III (UV) and Paper IV (galactic) are not audited.
