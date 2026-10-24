# SIM155 — RESULT (2026-10-24)

**Verdict (pre-registered criteria, SPEC.md / ROADMAP 3.4): FAIL.**

- **The test.** The CMSTG attractor (family A) was fitted to Planck 2018 + DESI DR2 + supernovae with a full Boltzmann code (hi_class) and compared with ΛCDM and w0waCDM.
- **(a) fails.** CMSTG improves on ΛCDM by at most Δχ² = +4.65 on the baseline data set; the pre-registered pass line is 6.
- **(b) fails.** That improvement exists only where the change in G since BBN is about 3× larger than primordial helium allows.
- **(c) fails.** It is smaller on the other two supernova sets.
- **Not competitive with w0waCDM**, which improves on ΛCDM by +9.2 … +19.6 on the same data (pre-registered "competitive" = within 3).

## Best fits (χ²; Δ > 0 means better than ΛCDM)

| data | ΛCDM χ² | w0waCDM Δχ² (ΔAIC) | CMSTG-A Δχ² (ΔAIC) | CMSTG − w0waCDM |
|---|---|---|---|---|
| D1: Planck 2018 + DESI DR2 + Pantheon+ | 2432.60 | +9.19 (−5.19) | **+4.65** (−0.65) | −4.54 |
| D2: … + Union3 | 1055.25 | +15.53 (−11.53) | +4.50 (−0.50) | −11.02 |
| D3: … + DES-Y5 | 2676.22 | +19.60 (−15.60) | +3.19 (+0.81) | −16.41 |

- **Fits.** Every model and data set was fitted from 3 random starts with the tested optimiser `common/fitcore.py`. For CMSTG on D1 a 10-point profile likelihood in ΔF was added (see below).
- **Agreement between starts.** ΛCDM starts agree to ≤ 0.05 in χ²_eff; w0waCDM D1 to 0.15.
- **Precision.** Evaluation noise is ±0.07 in χ²; the GR-limit offset is +0.02 ± 0.08. Neither matters at these thresholds (see Systematics).
- **Full tables** (χ²_eff, seed spreads, every profile point): `results.md`; machine-readable: `output.json`.

## Pre-registered criteria (CMSTG on D1)
- **(a) Δχ² ≥ 6 vs ΛCDM: FAIL.** The best of 3 free fits plus the 10-point profile is Δχ² = +4.65 (Δχ²_eff = +4.82), at ΔF = 0.040 (Λ0 = 0.297, ξ = 0.218, H0 = 66.27).
- **(b) C2 and C3 at the best fit: FAIL.**
  - C2 (Cassini, LLR) holds: |γ − 1| = 2.9e-8, because the attractor is reached (F0 = F* = 0.540).
  - C3 (BBN) fails: Y_p pull +3.25σ, D/H pull +2.24σ.
- **(c) Δχ² ≥ 6 on D2 or D3: FAIL.** +4.50 (Union3) and +3.19 (DES-Y5).
- **Competitive with w0waCDM (within 3): no.** CMSTG trails it by 4.5, 11.0 and 16.4.

## Exploratory (not pre-registered): what survives where BBN allows it
The best CMSTG point that passes C2 and C3 on D1 is a free fit at ΔF = 0.0113 (Λ0 = 0.145, ξ = 0.081, H0 = 68.14; Y_p pull +1.80σ).
- It gains **Δχ² = +3.25** (Δχ²_eff +3.33) over ΛCDM with two extra parameters. That is ΔAIC = +0.75, a 1.3σ-level preference for two extra parameters, i.e. none.
- The profile points inside the BBN edge (ΔF = 0.004, 0.008, 0.012) give +2.4 … +3.0, consistent with that.
- Figure: `figures/sim155_profile.png`.

## Why: where the improvement comes from (D1, χ² change vs the ΛCDM best fit)

| likelihood | CMSTG best (ΔF = 0.040) | CMSTG, BBN-allowed (ΔF = 0.011) | w0waCDM |
|---|---|---|---|
| Planck high-ℓ TTTEEE (plik-lite) | −3.35 | −1.55 | −2.60 |
| Planck low-ℓ EE | −1.13 | −0.84 | −0.98 |
| Planck low-ℓ TT | +0.41 | +0.25 | +0.31 |
| Planck lensing | +0.69 | +0.10 | −0.44 |
| DESI DR2 BAO | −1.86 | −1.85 | −2.88 |
| Pantheon+ | +0.60 | +0.63 | −2.59 |
| A_planck prior | −0.16 | −0.08 | −0.11 |
| **total (χ²_eff)** | **−4.82** | **−3.33** | **−9.30** |

(Entries are rounded; totals are from `summarize.py`.)

- **The field settles before recombination.** At the best fit, G at BBN is 8.0% stronger than today, 0.25% stronger at recombination, and equal to today for z ≲ 10 (`growth.json`). At the BBN-allowed point the figures are 2.3% and 0.5%.
  - So CMSTG acts as an early-universe modification of gravity (between BBN and recombination). After recombination it expands like ΛCDM with a bare Λ, at shifted parameters.
- **Where it gains.** The CMB damping tail and the sound horizon (plik, low-ℓ EE, DESI). It also lowers ω_b (0.02226 vs 0.02251) and H0 (66.3 vs 68.3).
- **Where it cannot gain.** It cannot produce the late-time evolution that the supernovae and DESI favour. Pantheon+ gets worse (+0.6), while w0waCDM gains −2.6 there. That is why w0waCDM's advantage grows with the supernova sets that prefer late-time evolution most (Union3, DES-Y5) while CMSTG's does not.
- **Why BBN binds.** The same early change of G is what BBN measures: Y_p and D/H rise with S = √(G_BBN/G0). The profile shows Δχ² growing with ΔF up to ≈ 0.04, but BBN cuts at ΔF ≈ 0.015.
- **Structure growth gets worse.**
  - At the CMSTG best fit S8 = 0.871, against 0.811 for ΛCDM (+7.5%); σ8 is +3.9% and fσ8 is +3.3 … +3.9% at z = 0–1.
  - At the BBN-allowed point S8 is +2.5%.
  - Either way the tension with weak-lensing S8 (≈ 0.76–0.79) gets worse.
  - w0waCDM: S8 +1.5%.

## The valley, and why a profile was needed
- **The degeneracy.** The CMSTG posterior has a long, curved, nearly flat valley in (Λ0, ξ, H0): larger ΔF = ξ²/(4Λ0) goes with lower H0, at an attractor field φ* = ξ/(2Λ0) ≈ 0.28–0.37. This is the F–H0 degeneracy found in the audit, partly broken by the full CMB.
- **Free starts stop at different points.** Two of the three free D1 starts "converged" near ΔF ≈ 0.01. The third drifted to ΔF = 0.035 without converging. The spread was 1.3 in χ²_eff, above the 0.5 re-fit trigger (RESEARCH_LOG decision 5).
- **The profile settles it.** The profile likelihood in ΔF (decision 6; method only, criteria unchanged) fixes ΔF and minimises everything else.
  - It shows a smooth rise from Δχ² ≈ +2.4 (ΔF = 0.004) through ≈ +3.0 … +3.5 at the BBN edge (ΔF ≈ 0.015, where the Y_p pull reaches 2σ) to a broad maximum ≈ +4.5 … +4.65 at ΔF ≈ 0.03–0.04.
  - It falls again beyond: +3.5 at 0.048 and +3.9 at 0.056, where H0 is down to 65.5.
  - At fixed ΔF a second, nearly flat direction (Λ0 at fixed ΔF) leaves ≈ 0.2–0.3 of point-to-point scatter; for this reason both the profile and the free fits are used.
- **Same issue on D2/D3.** The CMSTG D2 and D3 free fits show the same valley (seed spreads 1.2 and 0.6). Their profiles were not run: criterion (c) cannot change the verdict once (a) and (b) fail on D1, as decided in advance in RESEARCH_LOG.

## Why SIM154's large DESI gain shrank
- **SIM154 was not a joint fit.** It found Δχ² ≈ +17.7 (DESI DR2) for this model in a background-level test with Planck's ω_b and ω_m fixed, against a ΛCDM with the same fixed parameters.
- **ΛCDM re-adjusts in a joint fit.** ΛCDM itself moves (H0 68.3, ω_cdm 0.1179) and its DESI χ² drops from ≈ 30 to 12.7.
- **Same for w0waCDM.** It fell from +23.8 (SIM154's setup) to +9.2. SIM154's caveat ("the DESI preference is overstated") is confirmed and quantified.

## Systematics (none changes the verdict)
- **Evaluation noise.** hi_class at default precision jitters χ² by 0.07 RMS, all from plik-lite (`diagnostics.py`). Best-fit χ² values are good to ±0.1–0.3.
- **GR-limit offset.** The hi_class CMSTG model at ξ = 0 vs plain ΛCDM gives Δχ² = +0.02 ± 0.08 over 4 cosmologies (`gr_limit_chi2.py`). No bias.
- **w0waCDM on D2 has its own weak w0–wa degeneracy** (seed spread 1.85). Its true best fit can only be lower, which widens CMSTG's deficit.
- **Data.** Planck plik-lite rather than full plik/CamSpec; DESI's own w0waCDM numbers use PR4 CamSpec + ACT lensing. Our w0waCDM (w0 = −0.85, wa = −0.55, Δχ² = 9.2 on D1) agrees with DESI DR2 (w0 = −0.838 ± 0.055, wa = −0.62^{+0.22}_{−0.19}).
- **BBN.** Linear fitting formulas (Kneller & Steigman 2004). A 3σ Y_p pull at the best fit is far beyond their accuracy; a BBN code would refine the edge (ΔF ≈ 0.015), not move it by 3×.

## Method record
- **v1** (Cobaya `minimize`) stopped early and crashed on turning-point c_s² noise.
- **v2** (BOBYQA + MIGRAD) was too slow.
- **v3.** `fitcore.py` (coarse BOBYQA, then Hessian-whitened BOBYQA rounds, failure retries, continuous prior wall), unit-tested in CI.
- All changes were method-only and documented before the corresponding results; see RESEARCH_LOG.md.

## Files
- `fit.py`: v3 driver; `run_all_v3.sh`: 27 fits; `fits_v3/`: one JSON per fit, with the full χ² breakdown and Laplace errors.
- `profile.py`, `run_profile.sh`; `profile/`: the ΔF profile.
- `summarize.py`: criteria, writes `results.md` and `output.json`; `plot_profile.py`: `figures/sim155_profile.png`.
- `growth.py`: σ8, S8, fσ8(z) and G(z) at the best fits; writes `growth.json`.
- `diagnostics.py` and `gr_limit_chi2.py`: noise, failures, GR-limit systematics.
- `fits/`: v1 runs, kept for the record.
