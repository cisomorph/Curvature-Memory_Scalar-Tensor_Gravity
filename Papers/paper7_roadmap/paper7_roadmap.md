# Paper VII Contingency Roadmap
## CMSTG Response to DESI Year 3 Dark Energy Tension

**Status:** Contingency — activates only if DESI Y3 confirms w₀wₐ tension ≥ 3σ  
**Drafted:** 2026-05-16  
**Trigger evaluation:** See `trigger_checklist.md`  
**Pre-data work:** See `predata_workplan.md`  
**References:** See `bibliography.bib`

---

## Executive Summary

Paper VII extends the CMSTG action with a pre-recombination physics sector
to occupy the structural loophole identified at the end of Paper I (Section 6,
"The Remaining Loophole"). The loophole: if r_s is simultaneously increased
by ~1.5% through pre-recombination modification, late-time H(z) boosts can
satisfy both Planck θ* and DESI w₀wₐ, evading Theorem 2 of Paper I. All
routes to this loophole via the existing Ψ field are closed (SIM147–148). A
new sector is required.

If DESI Y3 confirms tension at ≥ 3σ, occupying this loophole becomes
scientifically urgent. Paper VII is the execution plan. If DESI Y3 reduces
tension to < 2σ, a 2-page note suffices and this roadmap is shelved.

The roadmap runs in four gated phases over ~20–30 weeks from DESI Y3 data
release. Pre-data preparation (see `predata_workplan.md`) can reduce this to
~10–15 weeks by completing Phase VII-1 and partial Phase VII-2 ahead of time.

---

## CMSTG Programme Context

The relevant results from the completed programme:

| Result | Location | Key number |
|--------|----------|-----------|
| Locked action established | Paper II, SIM88–113 | Λ₀ = 0.003 M_Pl⁻², Ψ₀ = 2.62 M_Pl, F_eff = 0.521 |
| DESI tension structural | Paper I, SIM144 | 2.77σ irreducible within action class |
| Theorem 1 (H(z) suppression) | Paper I, SIM131–136 | Curvature-sourced Ψ suppresses H(z); 6 mechanism classes all FAIL |
| Theorem 2 (θ* violation) | Paper I, SIM141–143 | Late-time H(z) boost shifts θ* by ≥ 26σ at fixed r_s |
| Theorem 2 loophole | Paper I, Section 6 | Simultaneous Δr_s/r_s ~ +1.5% could compensate |
| Ψ_pre loophole closed | SIM147–148 | τ_mem = 205,000 yr; pre-bang Ψ memory decays before recombination |
| UV finiteness (Paper III) | SIM115–120 | k_m ≥ 0.01 Mpc⁻¹ required; k_m < 2×10⁻³ Mpc⁻¹ ruled out |

**Consequence:** The loophole cannot be occupied from within the locked Ψ
sector. A new physical sector active near z ~ 1000–3000 is the only viable
route.

The required r_s shift: Paper I quotes Δr_s/r_s ~ +1.5% as sufficient to
compensate a 5% late-time H(z) boost. The precise target depends on the DESI
Y3 tension level (to be determined from the data). For planning:

| DESI Y3 tension | Required Δr_s/r_s | Required Δw₀ or ΔH₀ |
|-----------------|--------------------|----------------------|
| 3.0σ | ~+1.5% | Δw₀ ~ 0.2 or ΔH₀ ~ 1.5 km/s/Mpc |
| 3.5σ | ~+2.0% | Δw₀ ~ 0.3 or ΔH₀ ~ 2.0 km/s/Mpc |
| 4.0σ | ~+2.5% | Δw₀ ~ 0.4 or ΔH₀ ~ 2.5 km/s/Mpc |

These estimates follow from the CMSTG H(z) structure established in Paper II.
Recalibrate against DESI Y3 best-fit parameters once available.

---

## Trigger Conditions

**See `trigger_checklist.md` for the full evaluation protocol.**

Abbreviated: Paper VII proceeds if and only if
(a) DESI Y3 public data release with official cosmological analysis,
(b) w₀wₐ tension ≥ 3σ with ΛCDM in the combined DESI + CMB fit,
(c) consistent signal across ≥ 2 independent tracer samples,
(d) no identified systematic explanation reducing tension to < 2σ, and
(e) no independent experiment (DES, Euclid) contradicting the tension.

---

## Phase VII-1: Candidate Identification (Weeks 1–4)

**Goal:** Rank the three candidate mechanisms by viability within CMSTG.
**Deliverable:** 5–10 page comparison document (see format below).

### The three candidates

#### Candidate A: Early Dark Energy (EDE)

**Mechanism:** A separate scalar field φ with potential V(φ) undergoes a
brief phase of cosmological energy dominance near matter-radiation equality
(z_c ~ 3000–5000), increasing H(z) in this range and thereby decreasing r_s.
After the phase transition, φ rapidly dilutes and is cosmologically invisible.

**Standard EDE reference:** Poulin, Smith, Karwal & Kamionkowski (2019)
[arXiv:1811.04083] — cited in Paper I. Also Smith, Poulin & Amin (2020)
[arXiv:1908.06995]; Murgia, Murgia & Freese (2021) [arXiv:2007.04999].

**Embedding in CMSTG:**
The CMSTG action gains a new term:
```
S_VII = S_CMSTG + ∫d⁴x√(-g)[ -½(∂φ)² - U(φ) ]
```
where S_CMSTG is the locked Phase 1 action (Λ₀, Ψ₀ unchanged). The φ field
is minimally coupled to gravity (no φ²R term); it contributes only to T_μν.

Candidate potential forms (in order of simplicity):
1. **Axion EDE:** U(φ) = M⁴[1 − cos(φ/f)]^n, n ∈ {2,3}
2. **Power-law EDE:** U(φ) = M⁴(φ/M_Pl)^2n with dynamical attractor
3. **Rock'n'roll EDE:** U(φ) = M⁴[1 − cos(φ/f)]² with n=2 kinetic boost
4. **NEDE:** New Early Dark Energy with first-order phase transition [Niedermann & Slothus 2020]

**Parameters:** For axion EDE: {M (energy scale), f (decay constant), n, z_c (phase transition redshift)}.
Approximate viable ranges from Planck + pre-DESI-Y3 constraints:
- f_EDE (peak EDE fraction) ~ 5–15%
- z_c ~ 3000–5000
- f (axion decay constant) ~ 0.1–0.5 M_Pl

**CMSTG-specific constraint checklist for Candidate A:**
1. φ must not develop a non-minimal coupling to R (i.e., no φ²R term generated radiatively). Check Paper III results — if the CMSTG UV structure generates φ²R at loop level, this must be bounded.
2. φ must not couple to Ψ at tree level. If φ couples to Ψ through gravity, verify the coupling is suppressed by M_Pl².
3. The SIM147 τ_mem constraint does not directly apply to φ (φ is a new field, not Ψ). Confirm independence.
4. GW170817: φ must not contribute to tensor propagation speed modification. For minimal gravitational coupling of φ, c_T = c is automatic.
5. BBN constraint: φ energy density at z ~ 10⁸ must be < 10% of radiation energy. Verify U(φ) dilutes faster than radiation after the phase transition (requires w_φ > 1/3 during dilution).

**Structural problems with Candidate A (known from literature):**
- Worsens S8 tension (increases matter power spectrum via enhanced early growth)
- Requires f_EDE ~ 10% which adds Δχ²_CMB ~ 5–10 from Planck high-ℓ TT
- May conflict with Lyman-α forest constraints if z_c too early
- Naturalness: requires two independent hierarchically separated scales (M << f << M_Pl)

**Assessment:** Strongest phenomenological fit to H0 tension literature; most
studied; clear CMSTG embedding; structural problems are known and quantified.

---

#### Candidate B: Extra Relativistic Species (ΔNeff)

**Mechanism:** Additional relativistic energy density at recombination from a
dark sector (dark photons, sterile neutrinos, majorons, etc.) increases H(z)
at z ~ 1000–10⁶, decreasing r_s. Unlike EDE, the effect persists throughout
radiation domination rather than being concentrated near one redshift.

**Standard references:** Schöneberg et al. (2022) [arXiv:2207.07045];
Bernal et al. (2016) [arXiv:1607.05617]; Abazajian et al. (2019) [arXiv:1907.04473].

**Embedding in CMSTG:**
```
S_VII = S_CMSTG + S_dark_radiation
```
where S_dark_radiation contains a dark sector decoupled from the Standard Model
and the Ψ field, minimally coupled to gravity. The dark sector contributes:
ρ_dark(z) = ΔNeff × (7/8) × (4/11)^(4/3) × ρ_γ(z)    [fermion dark radiation]
ρ_dark(z) = ΔNeff × (1/2) × (...)                       [boson dark radiation]

**Parameter:** ΔNeff ∈ [0.1, 0.5] (current Planck 95% CL upper bound: ΔNeff < 0.34).

**r_s shift from ΔNeff:**
Δr_s/r_s ≈ −(1/2) × (ΔNeff / 3.046) × [∂ ln r_s / ∂ ln N_eff]
For ΔNeff = 0.3: Δr_s/r_s ≈ −1.5% (note: negative sign — ΔNeff decreases r_s)

**Critical issue:** ΔNeff > 0 DECREASES r_s (more radiation → smaller sound horizon), which worsens the H0 tension rather than helping it. For DESI w₀wₐ tension specifically (which differs from H0 tension), the direction of the r_s shift must be re-evaluated.

**CMSTG-specific analysis required:**
Check whether the DESI w₀wₐ tension calls for Δr_s/r_s > 0 (loophole direction) or < 0. If the DESI tension is in w₀wₐ space (not H0 space), the sign of the required r_s shift may differ from the H0 tension. This is the key pre-data question for Candidate B.

**Assessment:** Simplest mechanistically; weakest in terms of resolving the DESI tension in the correct direction; most testable by CMB-S4 and future Neff measurements.

---

#### Candidate C: Modified Recombination

**Mechanism:** Changes to baryon-photon recombination physics shift z* (effective decoupling redshift) and modify r_s. Approaches include:
- Modified helium fraction Y_He at recombination
- Varying fundamental constants near z ~ 1000 (e.g., α, m_e)
- Dark matter-baryon scattering modifying photon mean free path
- Early reionization modifying optical depth

**Standard references:** Hart & Chluba (2018) [arXiv:1808.04094];
Lee & Chluba (2022) [arXiv:2204.11832]; Hart & Chluba (2017) [arXiv:1705.03928].

**Embedding in CMSTG:**
Candidate C requires modifying either the baryon sector (Standard Model)
or introducing a baryon–Ψ coupling active at recombination. The latter is
severely constrained by the SIM151 cosmological screening bound:
β_∞ ≲ 2.4×10⁻¹⁰ at z = 0 → β(z ~ 1000) ≲ β_∞ × (1+z)^k for some k
For any reasonable k, β(z ~ 1000) ≪ 1, making the Ψ-baryon coupling
mechanism ineffective at recombination.

Alternative within CMSTG: modify Y_He or m_e via a separate hidden sector
weakly coupled to the baryon sector, not through Ψ. This is effectively
Candidate B with a baryon-coupling added.

**Assessment:** Most difficult to embed naturally within CMSTG without
violating SIM151 screening constraints. Produces sharpest spectral
predictions (measurable with CMB-S4 at 5–10σ). Lower priority than A or B.

---

### Phase VII-1 Week-by-Week Plan

**Week 1 (days 1–7): Literature survey and DESI Y3 data ingest**

Tasks:
1. Download and read DESI Y3 official cosmological analysis paper.
   Record: best-fit w₀, wₐ; tension with ΛCDM in σ; tension with Planck θ*; breakdown by tracer sample.
2. Record DESI Y3 best-fit r_s (from BAO standard ruler) and D_H(z)/r_s, D_M(z)/r_s at each tracer redshift.
3. Compute the required Δr_s/r_s to resolve tension at the DESI Y3 central value:
   - Start from Paper I Section 4 derivation of θ* shift vs Δr_s
   - Target Δr_s = r_s(Y3 best-fit) − r_s(locked CMSTG from SIM88)
4. Check Planck 2018 Neff posterior and CMB high-ℓ χ² for updated constraints.
5. Read Murgia et al. 2021, Schöneberg et al. 2022, and Smith et al. 2020 for
   current EDE constraints with Y1-equivalent data. Note which constraints will
   change with Y3 data.

**Week 2 (days 8–14): Candidate A analysis**

Tasks:
1. Compute r_s shift for axion EDE at n = 2, 3 for f_EDE ∈ [0.05, 0.15]:
   Using r_s = ∫₀^{t*} cs(t)/a dt, estimate Δr_s numerically under
   ρ_EDE(z) = f_EDE × ρ_tot(z_c) × [1 − tanh²((z−z_c)/Δz)]
   (top-hat approximation around z_c).
2. Verify GW170817 compatibility: φ minimal coupling → c_T = c automatically. Flag only.
3. Verify BBN compatibility: ρ_φ(z_BBN) / ρ_rad(z_BBN) < 0.1. Compute for
   axion EDE dilution rate (w_φ ~ 1 after phase transition for n = 2).
4. Check UV finiteness compatibility: does adding φ with V(φ) modify the
   CMSTG one-loop self-energy Σ(k) computed in Paper III? If φ has no coupling
   to Ψ and no R coupling, Paper III's two-loop result is unchanged. Verify.
5. Compute Planck TT/EE penalty (approximate): ΔNeff_eff from EDE injection.
   Use Schöneberg et al. Fig. 3 or equivalent to read off approximate Δχ²_Planck
   for the required f_EDE.

**Week 3 (days 15–21): Candidates B and C analysis**

Tasks:
1. Candidate B: Determine sign and magnitude of Δr_s/r_s from ΔNeff = 0.1, 0.2, 0.3.
   Compare direction with the sign required to resolve DESI Y3 tension.
   If sign is wrong: Candidate B is ruled out as the primary mechanism.
2. Candidate B: Compute Planck constraint on ΔNeff using current posterior (Aghanim+2020).
   Determine whether the ΔNeff required to produce the needed Δr_s is within 2σ.
3. Candidate C: Check β_∞ ≲ 2.4×10⁻¹⁰ bound from SIM151 against any Ψ-baryon
   mechanism at recombination. Document that this rules out Ψ-mediated modified
   recombination.
4. Candidate C: Identify whether a non-Ψ modified recombination (e.g., dark sector
   interacting with electrons) can be embedded without new gravitational coupling.
   Estimate order-of-magnitude r_s shift achievable.

**Week 4 (days 22–28): Comparison document**

Produce the Phase VII-1 comparison document (5–10 pages) with:

1. Summary table (see format below)
2. Recommendation: which candidate to pursue first in Phase VII-2
3. Ranking with explicit go/no-go reasoning per candidate
4. Outstanding uncertainties that Phase VII-2 will resolve

**Comparison document table format:**

| Criterion | Candidate A (EDE) | Candidate B (ΔNeff) | Candidate C (Mod. Recomb.) |
|-----------|-------------------|---------------------|---------------------------|
| Required Δr_s achievable? | f_EDE ~ 7–12% | Depends on sign | Limited by β_∞ bound |
| Planck penalty (Δχ²) | +5 to +15 | < +3 | Varies |
| BBN compatible? | Yes (w_φ > 1/3 after transition) | Yes (DR decoupled) | Yes |
| GW170817 compatible? | Yes (minimal coupling) | Yes | Yes |
| UV-finite in CMSTG context? | Yes (φ not coupled to R,Ψ) | Yes | TBD |
| New parameters | {M, f, n, z_c} | {ΔNeff} | Varies |
| Naturalness | Poor (two scales) | Moderate | Poor |
| Literature status | Active debate | Constrained | Less studied |
| **Recommended?** | **Primary** | **Secondary** | **Tertiary** |

If Candidate A fails Phase VII-2 consistency checks, default to the next
ranked candidate without returning to Phase VII-1.

---

## Phase VII-2: Action Derivation and Consistency Checks (Weeks 4–12)

**Goal:** Establish that the leading candidate action is self-consistent and
passes all seven constraints.
**Deliverable:** Derivation document with pass/fail per constraint.

### Assumed leading candidate: Candidate A (axion EDE)

If Phase VII-1 ranks a different candidate first, adapt the following to
that candidate. The structure is the same.

### Extended action

```
S_VII = ∫d⁴x√(-g) { [(½ + Λ₀Ψ²)R − ½(∂Ψ)² − ½m₀²Ψ²]
                    + [−½(∂φ)² − M⁴(1−cos(φ/f))^n] }
        + S_SM
```

New parameters: {M, f, n}. The phase transition redshift z_c is determined
by M and f through the Klein-Gordon equation:
  φ̈ + 3Hφ̇ + V'(φ) = 0
The field starts at φ_i ~ πf (bottom of potential), rolls when H ~ m_φ,
oscillates, and dilutes. z_c corresponds to H(z_c) ~ m_φ(φ).

### Consistency checks (sequential; failure at any check stops the candidate)

**Check 1: Recombination-era r_s shift**
- Solve the background Friedmann equation with ρ_φ(z) included.
- Compute r_s = ∫₀^{a*} cs(a)/(a²H(a)) da.
- Target: Δr_s/r_s = [r_s computed from DESI Y3 data] − [locked CMSTG r_s].
- Acceptance: |Δr_s/r_s − target| < 0.3% (i.e., within 20% of the required shift).
- Numerical tool: Extend the SIM120 background integrator with φ evolution.

**Check 2: CMB acoustic angle θ***
- Compute θ* = r_s(z*)/D_A(z*) for the extended action.
- Acceptance: |θ* − Planck best fit| < 0.1% (1σ Planck precision is 0.04%).
- This check verifies the EDE-increased r_s is compensated by an appropriate D_A.

**Check 3: CMB TT/EE/TE power spectrum**
- Compute the CMB power spectra using CAMB or CLASS with the extended action.
- Acceptance: Δχ²_CMB ≤ 9 (3σ) relative to locked Phase 2 best fit.
  (Note: f_EDE ~ 10% typically costs Δχ² ~ 5–15; this must be checked.)
- If Δχ² > 9 for any f_EDE that produces the required Δr_s: Candidate A FAILS
  this check. Move to Candidate B.

**Check 4: BAO standard ruler**
- Compute D_H(z)/r_s and D_M(z)/r_s at DESI Y3 tracer redshifts.
- Acceptance: Within 2σ of DESI Y3 measured values at all tracer redshifts.

**Check 5: Structure growth (f_σ8)**
- Compute the linear growth rate f(z) and σ8(z) for the extended action.
- Verify φ does not enhance early structure growth beyond BOSS/DESI RSD constraints.
- Key concern: EDE increases early H(z) → smaller Jeans scale → enhanced early growth → larger σ8.
  The S8 tension (σ8 vs CMB) must not worsen by > 2σ.
- Acceptance: S8 = σ8(Ω_m/0.3)^0.5 within 2σ of KiDS/DES constraints.

**Check 6: GW propagation**
- φ has no ∂φ∂φ coupling to the metric perturbation h_μν at tree level (minimal coupling).
- Verify no loop-level coupling of φ to gravitons modifies c_T. (For M << M_Pl and minimal coupling, this is automatic at 1-loop order.)
- Acceptance: |c_T/c − 1| < 10⁻¹⁵ (GW170817 bound).

**Check 7: BBN and early universe**
- φ must have w_φ > 1/3 after phase transition to dilute faster than radiation.
- For axion EDE: w_φ ≈ n/(n+1) for n = 2 → w_φ = 2/3 > 1/3. ✓
- For n = 3: w_φ = 3/4 > 1/3. ✓
- Compute ρ_φ(z_BBN)/ρ_rad(z_BBN) < 10⁻²: verify dilution rate is sufficient.

**Outcome of Phase VII-2:**
- All 7 checks pass → proceed to Phase VII-3.
- Check 1 or 2 fails → parameter tuning within Candidate A (adjust M, f, n).
- Check 3 fails (CMB penalty too large) → Candidate A is dead; move to Candidate B or C.
- Check 5 fails (S8 tension) → marginal case; document and continue with caveat.
- Check 6 or 7 fails → structural problem; document, reassess candidate.

---

## Phase VII-3: Numerical Pipeline and Parameter Scan (Weeks 12–20)

**Goal:** Build the simulation pipeline and perform joint cosmological fits.
**Deliverable:** Best-fit parameters, posterior distributions, Δχ² improvements.

### Simulation architecture

The Phase VII pipeline extends the existing CMSTG simulation infrastructure:

| Component | Base | Extension |
|-----------|------|-----------|
| Background integrator | SIM120 | Add φ(z) evolution via RK4 KG equation |
| CMB power spectra | CAMB/CLASS (used in SIM113) | Enable EDE module in CAMB (`AxiCLASS` or CLASS with EDE patch) |
| BAO likelihood | SIM111 DESI Y1 pipeline | Update to DESI Y3 data release files |
| Planck likelihood | SIM113 TT/EE/TE | No change needed |
| RSD likelihood | SIM113 f_σ8 | Update to DESI Y3 if new measurements released |

**Recommended tool:** `AxiCLASS` (CLASS extension for axion EDE; Poulin+2019)
or `EDE-CAMB` patch from Smith+2020. Both are public. Choose based on whichever
was used in Phase 2 CMB pipeline (check SIM113 code).

### Parameter scan structure

**Free parameters in Phase VII scan:**
- CMSTG locked: Λ₀ = 0.003 (fixed, not scanned — locked action)
- CMSTG locked: Ψ₀ = 2.62 M_Pl (fixed)
- CMSTG locked: k_m (fixed at Paper III value)
- EDE new: M (energy scale), f (decay constant), n ∈ {2, 3} (discrete)
- Standard cosmological: H₀, ω_b, ω_cdm, A_s, n_s (scanned jointly)

**Decision on locking:** Do NOT scan over Λ₀ or Ψ₀. They are locked by
Papers I–III. Scanning them would reopen questions settled by Phases 1–3.

**Likelihood:** Joint Planck TT/TE/EE (2018) + DESI Y3 (all tracers) + BBN (Ω_b h²).
Do NOT include late-universe weak lensing (KiDS, DES) in the primary fit — include as
a post-hoc consistency check only (S8 tension is a separate problem).

**SIM numbering:**
- SIM154: Background integrator with φ (Candidate A validation)
- SIM155: CMB power spectra with EDE sector
- SIM156: Joint parameter scan (MCMC or nested sampling)
- SIM157: DESI Y3 BAO fit with extended action
- SIM158: Posterior distributions and Δχ² vs. locked CMSTG and ΛCDM

**Decision criteria at Phase VII-3:**
- Δχ²(extended CMSTG vs locked CMSTG) > 9: extended action preferred at > 3σ. → PROCEED to VII-4.
- Δχ²(extended CMSTG vs ΛCDM) > 0: extended CMSTG preferred. → PROCEED.
- Δχ²(extended CMSTG vs locked CMSTG) < 4: preference below 2σ. → Report as "mild improvement"; downgrade Paper VII scope.
- Best-fit EDE parameters require f_EDE > 20%: conflict with current CMB constraints; re-examine candidate selection.

---

## Phase VII-4: Paper VII Drafting (Weeks 20–30)

**Journal target:** JCAP (primary), PRD (backup). Same target as Papers I–IV.

**Page target:** 20–25 pages (longer than companion note; comparable to Paper II).

### Paper structure

**Title (placeholder):**
"Curvature-Memory Scalar-Tensor Gravity with Early Dark Energy:
Resolving the DESI Year 3 Tension via Pre-Recombination Sound Horizon Modification"

**Abstract (~250 words):**
- State the DESI Y3 tension (confirmed σ-level).
- State Theorem 2 loophole from Paper I.
- State that the loophole requires pre-recombination physics; Ψ routes closed (SIM147–148).
- State the EDE solution: new scalar φ minimally coupled to gravity.
- State key numbers: f_EDE, z_c, Δr_s/r_s, DESI Y3 fit improvement (Δχ²).
- State predictions for future surveys.

**Section 1: Introduction**
- Recap of CMSTG programme (Papers I–IV) in 2 paragraphs.
- DESI Y3 result and tension level.
- Statement of the Theorem 2 loophole from Paper I.
- Statement that SIM147–148 closed the Ψ_pre route to the loophole.
- This paper's response: extend the action with φ.
- Paper structure outline.

**Section 2: The Theorem 2 Loophole and Why a New Sector is Required**
- Reproduce Theorem 2 statement from Paper I (cite).
- Reproduce the loophole identification from Paper I Section 6.
- Derive the required Δr_s/r_s quantitatively from DESI Y3 data.
- Show that SIM147–148 closed the Ψ_pre route (summarize, cite).
- Conclude: new sector required.

**Section 3: Extended Action**
- State S_VII explicitly.
- Derive field equations for Ψ (unchanged) and φ (new KG equation).
- Show φ does not modify Ψ field equation at tree level (decoupling argument).
- Show the locked CMSTG cosmological solution (Papers I–III) is preserved in the low-φ regime.

**Section 4: Pre-Recombination Cosmology and r_s Shift**
- Background Friedmann equations with φ.
- φ evolution through z_c; phase transition; dilution.
- Compute Δr_s/r_s as a function of {M, f, n}.
- Show match to DESI Y3 target at best-fit {M, f, n}.

**Section 5: Consistency with Locked Phase 1–3 Results**
- Pass each of Checks 1–7 from Phase VII-2 using SIM154–158 results.
- Tabulate: check name, target, result, status (pass/marginal/fail).
- Discuss any marginally passing checks (especially CMB TT penalty and S8).

**Section 6: DESI Y3 Fit**
- Joint posterior distributions (figures from SIM158).
- Best-fit parameter values with uncertainties.
- Δχ²(extended CMSTG vs locked CMSTG) and interpretation.
- Δχ²(extended CMSTG vs ΛCDM) and interpretation.
- BAO distance ratios D_H/r_s, D_M/r_s at each tracer redshift (figures).

**Section 7: Predictions for Future Surveys**
- DESI Y5/Y6: Precision on w₀wₐ expected to improve by ~√(5/3). Predicted
  tension level if EDE sector is correct: quote expected σ improvement.
- CMB-S4: f_EDE ~ 10% is detectable at ~3σ through high-ℓ TT and Neff constraint.
  Compute predicted ΔNeff_eff from the EDE sector.
- Euclid: Structure growth predictions from f_σ8 posterior in Phase VII-3.
- Simons Observatory: CMB-TE constraint on EDE through damping tail suppression.

**Section 8: Conclusions**
- Summary: DESI Y3 confirmed at Nσ; Theorem 2 loophole occupied by EDE extension;
  5 checks pass; 2 marginal.
- Prediction: CMB-S4 will detect or rule out this mechanism at ~3σ.
- CMSTG programme status: Papers I–IV establish boundaries; Paper VII proposes
  the minimal extension needed to resolve DESI tension.

**Supplementary material (appendix or online):**
- Full SIM154–158 code and output data.
- Posterior MCMC chains.

---

## Fallback Protocol

If Candidate A fails Phase VII-2 or VII-3:

| Failure | Response |
|---------|----------|
| Check 3 (CMB penalty > 3σ) | Move to Candidate B (ΔNeff) |
| Check 5 (S8 tension) | Continue with explicit caveat section in Paper VII |
| Δχ² improvement < 2σ in Phase VII-3 | Downgrade: no Paper VII; write a short note |
| Candidate B also fails CMB penalty | Move to Candidate C (modified recombination) |
| All three candidates fail | Paper VII is not written; write a 5-page "closed loophole" note |

The "closed loophole" note (if all candidates fail) would state:
- The Theorem 2 loophole cannot be occupied by EDE, ΔNeff, or modified recombination
  within the CMSTG framework as currently formulated.
- The DESI Y3 tension is structural within the entire CMSTG action class
  (extending Papers I–IV's conclusion from Y1 to Y3).
- Future paper required to consider qualitatively different architecture changes.

---

## Constraint Reference Table

For use during Phases VII-1 and VII-2. All numbers at time of roadmap drafting
(2026-05-16); update from DESI Y3 release.

| Constraint | Current value | Source | Sensitivity to EDE |
|------------|--------------|--------|-------------------|
| r_s (Planck+BAO) | 147.09 ± 0.26 Mpc | Planck 2018 | Modified by φ |
| θ* (Planck CMB) | (1.04101 ± 0.00029)° | Planck 2018 | Must be preserved |
| H₀ (CMB-indirect) | 67.4 ± 0.5 km/s/Mpc | Planck 2018 | Modified by φ |
| N_eff (Planck) | 2.99 ± 0.17 | Planck 2018 | ΔNeff_eff from EDE |
| σ8 (Planck+lensing) | 0.811 ± 0.006 | Planck 2018 | May worsen |
| S8 = σ8(Ω_m/0.3)^0.5 | 0.766 ± 0.020 | KiDS-1000 | Key tension |
| f_σ8(z=0.38) | 0.497 ± 0.045 | BOSS DR12 | Modified by φ growth |
| |γ_PPN − 1| (Cassini) | < 2.3×10⁻⁵ | Bertotti+2003 | Not affected (φ minimal) |
| c_T/c (GW170817) | 1 ± 10⁻¹⁵ | Abbott+2017 | Not affected (φ minimal) |
| ΔNeff (BBN, 2σ) | < 0.5 | PDG 2022 | Must satisfy |
| DESI Y1 tension | 2.77σ (CMSTG) | SIM144 | Target to resolve |

---

## Notes on Scope and Conservatism

1. **Locked action:** Paper VII does not modify Λ₀, Ψ₀, m₀, k_m, or any other
   parameter of the Phase 1 locked action. If any of these need modification to
   accommodate φ, that constitutes a Phase 6 programme, not a Paper VII extension.
   Paper VII is an additive extension (S_VII = S_CMSTG + S_EDE), not a revision.

2. **Minimality:** The strongly preferred route is a single new parameter (if
   possible). Axion EDE with fixed n = 2 leaves {M, f} as free parameters;
   the ratio M/f is the physically meaningful combination (sets z_c and f_EDE
   jointly). If a one-parameter family within Candidate A works, prefer it.

3. **Pre-data work:** Phase VII-1 and partial Phase VII-2 (derivation, checks 1-2, 6-7)
   can be completed before DESI Y3 release. See `predata_workplan.md`.

4. **If theory landscape changes:** By 2027, new EDE variants or new observational
   constraints may alter the candidate ranking. Revisit Phase VII-1 if more than
   12 months have elapsed since this roadmap was drafted.

5. **Authorship note:** All simulation code and outputs will be archived at the
   CMSTG GitHub repository with SIM154+ numbering consistent with the existing
   SIM_MAP.md.
