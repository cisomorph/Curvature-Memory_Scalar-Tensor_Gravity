# DESI Year 3 Trigger Checklist
## Paper VII Activation Protocol

**Purpose:** Evaluate whether DESI Y3 data release triggers Paper VII execution.  
**Expected timing:** DESI Y3 data release ~2027 (exact date TBD).  
**Decision deadline:** Within 4 weeks of the official DESI Y3 cosmological analysis paper.

---

## Trigger Conditions (all five must be satisfied)

### Condition 1: Official DESI Y3 data release
- [ ] DESI collaboration has released the Y3 cosmological constraints paper
      (full galaxy and quasar spectroscopic survey analysis, not preliminary)
- [ ] The paper is publicly available (arXiv or DESI publication)
- [ ] The paper includes w₀wₐ constraints from the combined BAO + CMB analysis

**Do not activate on:** Conference preliminary results, DESI data challenge papers,
or Y3 BAO-only papers without CMB combination.

---

### Condition 2: Tension ≥ 3σ with ΛCDM
- [ ] The combined DESI Y3 + Planck (or equivalent CMB) analysis gives w₀wₐ
      tension with ΛCDM at ≥ 3σ significance
- [ ] Record the exact σ value: _________σ  
- [ ] Record the best-fit w₀: _________  
- [ ] Record the best-fit wₐ: _________  
- [ ] Record the paper's stated significance level and methodology

**Marginal case (2σ ≤ tension < 3σ):** Do not activate Paper VII. Write a 2-page
note updating the locked CMSTG prediction to Y3 data. The CMSTG structural
prediction (2.77σ from SIM144) was derived against Y1 data; update it for Y3 and
compare.

**Enhancement case (tension > 4σ):** Activate Paper VII with elevated urgency.
Move pre-data preparation timeline up if not complete.

---

### Condition 3: Multi-tracer consistency
- [ ] The w₀wₐ tension is present in ≥ 2 independent DESI tracer samples when
      analyzed separately (LRGs, ELGs, QSOs, BGS, Lyα — check whichever are
      available in Y3)
- [ ] The tracer samples give w₀ and wₐ constraints that are mutually consistent
      at ≤ 2σ (i.e., the tension is coherent, not driven by one anomalous tracer)

**If tension is driven by one tracer:** Flag for human review. Do not automatically
activate. The single-tracer anomaly may be a systematic effect.

Record:
- LRG tension: _________σ
- ELG tension: _________σ  
- QSO tension: _________σ
- BGS tension: _________σ
- Lyα tension: _________σ

---

### Condition 4: No identified systematic explanation
- [ ] The DESI Y3 paper has performed systematic checks (fiber assignment, spectral
      pipeline, photometric calibration, stellar contamination, angular
      clustering) and found no uncorrected systematic that would reduce tension
      to < 2σ
- [ ] Independent analysis teams (if any) have reproduced the tension using
      different analysis pipelines
- [ ] The tension is not reduced to < 2σ by any alternative BAO analysis
      approach (e.g., alternative template or fitting range) described in the
      DESI Y3 paper itself

**If systematic concerns are raised but unresolved:** Flag. Delay activation by
4 weeks pending independent analysis. Activate only if systematic is refuted or
shown to be insufficient to reduce tension below 3σ.

---

### Condition 5: No contradicting evidence from independent experiments
- [ ] No other Stage-IV survey (Euclid, DESI+CMB-S4, DES Y6, HST/JWST
      calibrations) has published results contradicting the DESI Y3 DE tension
      at ≥ 2σ significance within 12 months of DESI Y3 release
- [ ] The Planck CMB constraints used in the combined analysis are the currently
      accepted Planck 2018 results (or updated CMB results if available); no
      major Planck re-analysis has changed θ* by more than 0.1%

**If Euclid contradicts DESI Y3 tension:** Flag. Do not activate. The tension
may be an instrumental systematic. Wait for reconciliation (6–12 months).

---

## Decision Outcome

### All 5 conditions satisfied → ACTIVATE Paper VII
Immediately begin:
1. Phase VII-1 (if not complete from pre-data work): start candidate comparison document
2. Or, if Phase VII-1 is pre-completed: immediately proceed to Phase VII-2

### Any condition not satisfied → DO NOT ACTIVATE Paper VII
Write instead: a 2–5 page note titled
"CMSTG Structural Prediction vs. DESI Year 3: [Confirmed/Reduced] Tension"

Content of the short note:
- Report current CMSTG tension level against Y3 data (recompute Δχ² using Y3 BAO values)
- If tension reduced: note that the locked action is now within 2σ of DESI Y3 and
  Paper VII is not warranted
- If tension increased but < 3σ: note that Paper VII remains contingent on future data

---

## Data to Record at Decision Time

Regardless of activation outcome, record the following for the CMSTG archive:

| Quantity | DESI Y1 value | DESI Y3 value |
|----------|--------------|--------------|
| Best-fit w₀ | −0.55 (Y1) | _________ |
| Best-fit wₐ | −1.32 (Y1) | _________ |
| Tension with ΛCDM | 2.5–3σ (Y1) | _________σ |
| r_s (Mpc, from BAO) | 147.09 ± 0.26 | _________ |
| H(z=0.51)/r_s | [DESI Y1] | _________ |
| H(z=0.71)/r_s | [DESI Y1] | _________ |
| H(z=0.93)/r_s | [DESI Y1] | _________ |
| CMSTG DESI tension | 2.77σ (SIM144) | Recompute |

Compute updated CMSTG Δχ²_DESI against Y3 data using the SIM144 methodology
(locked action background integrator + DESI BAO likelihood). Record result as
SIM159 in Ordered_Simulations/SIM_MAP.md.

---

## Notes

- The trigger evaluation should be completed by one person, documented in a
  ~1 page "DESI Y3 Assessment" note added to the CMSTG repository.
- If in doubt about any condition, err toward NOT activating Paper VII prematurely.
  A false positive wastes months. A false negative can be corrected 4 weeks later
  once systematics are resolved.
- If DESI Y3 and DESI Y4 data arrive in quick succession, base the trigger
  evaluation on the Y4 data if Y4 becomes available within 6 months of Y3.
