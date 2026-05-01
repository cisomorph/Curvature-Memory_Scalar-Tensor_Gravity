# SIM139 — SIM128 RSD Shape Diagnostic

**Tier:** 1 (Diagnostic)
**Depends on:** SIM128 outputs (S₈ = 0.753 at p=1, ν=0.550; χ²/N = 2.24 on RSD shape)
**Prerequisite reads:** RESEARCH_RULES.md, PHASE4_ROADMAP.md, Paper VI §4

---

## Question

Where does the SIM128 Ĝ_eff/G = 1 − νM̂ᵖ coupling fail the fσ₈(z) shape test? Is the disagreement localized to high-z, low-z, or distributed?

## Motivation

SIM128 is CMSTG's closest approach to resolving the S₈ tension: S₈ = 0.753 lands −0.24σ from KiDS-1000. The only reason SIM128 was classified PARTIAL and not PASS is the RSD shape χ²/N = 2.24 — the integrated amplitude is right but the redshift evolution is wrong.

If the shape breaks only at low-z or only at high-z, the M̂(a) functional form can likely be modified (e.g., a different power p or a late-time cutoff) to fix the shape while preserving S₈. That would turn a PARTIAL into a PASS without Phase 4 work.

## Action spec

No new physics. Decompose existing χ² against the 9 fσ₈ measurements.

## Inputs

- SIM128 best-fit: ν = 0.550, p = 1, Ĝ_eff/G = 1 − 0.550·M̂(a)
- 9 fσ₈ measurements at z ∈ {0.067, 0.15, 0.38, 0.51, 0.60, 0.70, 0.85, 1.05, 1.48} (6dFGRS, SDSS MGS, BOSS, VIPERS, eBOSS LRG/ELG/QSO — see master paper bibliography)
- CMSTG structure growth solver from SIM92 and SIM128
- M̂(a) memory field from SIM125

## Procedure

1. Reproduce SIM128 χ²/N = 2.24 as sanity check.
2. Compute model fσ₈(z) at each of the 9 data points using SIM128 best-fit parameters.
3. Compute per-z residuals (fσ₈_model − fσ₈_obs) / σ_obs.
4. Identify shape failure mode:
   - Low-z excess / high-z deficit (or vice versa)?
   - Single outlier bin?
   - Systematic slope?
5. Test three alternative M̂(a) modifications (phenomenological — label clearly):
   - (a) Late-time cutoff: Ĝ_eff unchanged for z > z_c, flat for z < z_c; scan z_c ∈ {0.3, 0.5, 0.7}
   - (b) Different power: p ∈ {0.5, 1.5, 2.0} with ν re-fit
   - (c) Two-parameter: Ĝ_eff/G = 1 − νM̂ + ηM̂² with both fit
6. For each, compute both S₈ and RSD χ²/N. Look for a modification that preserves S₈ < 0.78 while bringing RSD χ²/N below ~1.2.

## Success criteria

- **SHAPE_FIXABLE:** at least one modification achieves S₈ ∈ [0.75, 0.78] AND RSD χ²/N < 1.2. → SIM128-like mechanism upgraded to PASS; possible Paper VIII result without Phase 4.
- **SHAPE_STRUCTURAL:** no modification in the tested space achieves both. The shape failure is a property of the coupling form, not its parametrization. → continue to Tier 2 as planned.
- **LOCALIZED_OUTLIER:** one measurement dominates the χ² (>60% contribution); excluding it brings χ²/N ≤ 1. → flag for discussion; do NOT silently drop data, but note this as a possible data-systematic.

## Failure modes to watch

- Survey-by-survey systematics: the 9 fσ₈ measurements come from 6 different surveys with different selection and modelling. Some level of inter-survey scatter is expected. Check whether the failure is within one survey (→ modelling) or across all surveys (→ theory).
- p is a phenomenological index (labelled as such in SIM128). Any modification explored here is also phenomenological. Do not claim a derived result.
- S₈ computation depends on σ₈(0) normalization. Use the same normalization convention as SIM128 (check code, do not assume).

## Deliverables

- `output.json` with per-z residuals, modification scan results, and verdict
- `RESULT.md` with:
  - Figure: fσ₈ model vs data for SIM128 best-fit and for best alternative modification
  - Table: (ν, p) or (z_c) scan with S₈ and χ²/N
  - Verdict and recommendation
- If SHAPE_FIXABLE: a draft section for Paper VIII describing the upgrade

## Estimated time

Small. Uses existing SIM128 infrastructure.
