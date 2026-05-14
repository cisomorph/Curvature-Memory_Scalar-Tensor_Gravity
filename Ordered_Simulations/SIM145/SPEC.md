# SIM145 — UV Consistency Recheck on Winning Mechanism

**Tier:** 3 (Validation)
**Depends on:** SIM144 in progress or complete
**Prerequisite reads:** RESEARCH_RULES.md, PHASE4_ROADMAP.md, winner's SPEC, SIM102–106 UV programme

---

## Question

Does the winning Phase 4 mechanism preserve the UV finiteness and Ward identity structure established for Phase 1 through SIM102–106?

## Motivation

SIM104 established Ward identity Π_hh(0) = 0 at full one-loop. SIM105 established RG flow with Λ₀ = 0 UV fixed point (negative beta function, GR limit). SIM106 closed the two-loop graviton sector. These are CMSTG's strongest theoretical results — the theory is UV-complete at loop level.

Any new term in the Lagrangian can silently break these results. A winning mechanism that solves DESI but introduces UV divergences is not a physical theory — it is an effective description that needs to be embedded in something better.

This sim's role is to verify that the mechanism survives the same tests Phase 1 survived.

## Scope per mechanism

The tests required depend on which sim won:

### If SIM140 (step potential V(Ψ)) won
- **Low urgency** for UV. Polynomial / bounded potential does not introduce new quartic divergences.
- Verify: V(Ψ) vertex contributions to Π_hh(0) vanish (no new couplings to h_μν).
- Verify: resummed Ψ propagator still regulates Ψ-bubble insertions.
- Deliverable: formal argument + one-loop diagrammatic check. Small.

### If SIM141 (running Λ₀(a)) won
- **High urgency.** The Ward identity derivation in SIM104 assumed constant Λ₀. If Λ₀ varies, the derivation must be redone.
- Verify: Π_hh(0) = 0 when Λ₀ is treated as a background field.
- Check: does the a-dependence of Λ₀ require promoting it to a full field? If yes, SIM141 has silently introduced a new degree of freedom and the whole framework needs rederivation.
- Deliverable: full SIM104-analog calculation with Λ₀(a). Large.

### If SIM142 (Galileon G₃) won
- **Highest urgency.** G₃ introduces new vertices. Even if the theory is classically safe (α_i stable), loop structure may generate new divergences.
- Verify: Ward identity at tree and one-loop for the G₃-extended theory.
- Verify: RG flow of c₃ — is there an IR attractor analogous to Λ₀ ≈ 0.003?
- Check: Horndeski theories are generically non-renormalizable but have calculable effective-theory cutoffs (Λ₃ ~ (M_Pl H₀²)^(1/3) ≈ 10³ km⁻¹). Identify the cutoff and verify it is above the scales of interest.
- Deliverable: full SIM104 + SIM105 analog with G₃ vertex. Very large.

### If SIM143 (decoupled φ) won
- **Medium urgency.** Standard quintessence is perturbatively fine at the cosmological-background level but has EFT cutoff questions at higher energies.
- Verify: U(φ) quartic self-coupling (if exponential, expand in series) is non-negative and stable.
- Verify: Ψ and φ sectors remain decoupled at loop level — no induced Ψ-φ coupling from matter loops.
- Deliverable: formal argument for sector independence + EFT cutoff statement. Medium.

## Procedure (general)

1. Read winner's RESULT.md. Identify the new Lagrangian term(s).
2. Determine new Feynman vertices arising from the term.
3. For each Phase 1 UV result (Ward identity, Λ₀ fixed point, two-loop graviton), check whether the new vertex contributes to the relevant diagrams.
4. For each diagram that does receive new contributions:
   - Compute one-loop correction with memory-regulated propagators
   - Verify finiteness (either exact cancellation, symmetry argument, or explicit regulation by k_m)
5. If any divergence appears: is it renormalizable (shifts existing coupling) or new (requires new counterterm)?

## Success criteria

- (a) Ward identity Π_hh(0) = 0 preserved at one-loop in the new theory
- (b) No new quartic UV divergences
- (c) RG flow of Λ₀ (and c₃ if SIM142) well-behaved — no Landau pole in the range of interest
- (d) Two-loop graviton sector stable (δZ_h finite)
- (e) Memory damping k_m unchanged or consistently redefined

All (a)-(e): SIM145 PASS → winning mechanism is UV-consistent with Phase 1.
Partial: document which results are preserved and which require new work. A partial UV pass can still allow Paper VIII to proceed IF the broken result is acknowledged as a Phase 5 task.

## Failure modes to watch

- **Silent breaks:** the Ward identity looks "obviously preserved" but isn't because the new vertex couples differently. Always work out the actual diagram.
- **Regulator dependence:** memory damping is a specific regulation scheme. Results derived there may not hold in dimensional regularization. Make sure the Phase 1 scheme (memory-regulated) is used consistently.
- **EFT breakdown:** for Galileon especially, the theory is valid only below a cutoff. If the cutoff is at or below Hubble scale today, the theory isn't cosmologically applicable at the classical level either — the classical fit is accidental.
- **Assuming previous results:** do not write "Ward identity preserved by analogy with SIM104." Derive or explicitly cite the argument.

## Deliverables

- `output.json` with per-check verdicts and any numerical loop integrals
- `RESULT.md` with full UV analysis, diagrams (can be Feynman diagrams rendered with TikZ), analytic cancellations
- If PASS: confirmation section for Paper VIII
- If FAIL: "UV obstruction identified" section explaining what breaks and motivating Phase 5

## Estimated time

Small (SIM140 winner) to Very Large (SIM142 winner). Calibrate after SIM144 winner is identified.
