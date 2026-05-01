# SIM142 — P4-A: Galileon G₃(Ψ)□Ψ Sector

**Tier:** 2 (Mechanism)
**Depends on:** Tier 1 complete; SIM140, SIM141 results reviewed
**Prerequisite reads:** RESEARCH_RULES.md, PHASE4_ROADMAP.md, Paper VII direction P4-A, SIM108 GW constraints

---

## Question

Can the simplest kinetic-gravity-braiding Horndeski term G₃(Ψ)□Ψ, added to the Phase 1 canonical action, produce tracker behaviour with effective w(z) crossing −1 at late times while satisfying all GW, CMB, and BAO constraints?

## Motivation

Galileon / kinetic-gravity-braiding is the most theoretically sophisticated of the Paper VII directions and the most powerful — it can naturally produce w_0 < −1, w_a < 0 without fine-tuning. If SIM140 and SIM141 have failed or only partial-passed, this is the remaining mainline mechanism.

The trade-off: G₃ introduces non-trivial dynamics that can re-break the graviton sector constraints that Phase 1 spent SIM102–106 closing.

**Run this only after SIM140 and SIM141 have been documented.** If either produces a clean PASS, SIM142 may be deferred.

## Action spec

```
L = (1 + 2Λ₀Ψ²)/2 · R  −  ½(∂Ψ)²  −  ½m₀²Ψ²  +  G₃(Ψ) □Ψ  +  L_matter
```

with Λ₀ = 0.003, m₀ as in Phase 1, and G₃(Ψ) parametrized as:

- Simplest: G₃(Ψ) = c₃ · Ψ (linear)
- Alternative: G₃(Ψ) = c₃ · Ψ² (quadratic)

The coefficient c₃ is scanned across 4 decades: c₃ ∈ {10⁻⁴, 10⁻³, 10⁻², 10⁻¹} in appropriate natural units.

Note: the G₃ term gives kinetic braiding — mixes Ψ and metric kinetic sectors. Requires careful modification of both Friedmann equations and the scalar EOM.

## Theoretical checks (required before PASS)

1. **GR recovery:** G₃(0) = 0 for linear and quadratic forms. Ψ → 0 limit gives back standard action. OK by construction.

2. **c_T = c (GW170817):** G₃ alone does NOT modify c_T (unlike G₄, G₅ which do). c_T = c is preserved. **Verify this explicitly** — the literature result relies on Horndeski classification; make sure the specific Lagrangian used respects it.

3. **No ghost / no gradient instability:** the kinetic matrix of the scalar-tensor system must be positive-definite. This is the critical check. For kinetic braiding, the conditions are:
   - α_T = 0 (tensor speed)
   - α_M related to F'(Ψ) — inherited from Phase 1
   - α_B (braiding) proportional to G₃'(Ψ)Ψ̇H
   - α_K + (3/2)α_B² / (effective mass terms) > 0 — scalar kinetic positivity
   
   Evaluate α_i at the best-fit background and verify all positive.

4. **No tachyon at Phase 1 attractor:** m² of perturbations around Ψ̄ must be positive.

5. **Ward identity:** G₃ preserves diffeomorphism invariance. Ward identity Π_hh(0) = 0 should hold at tree level. One-loop corrections unknown — flag for SIM145.

6. **UV finiteness:** G₃(Ψ) introduces new vertices. The memory-regulated propagator from Phase 1 should still suppress loops if k_m is unchanged, but G₃ vertex structure may introduce new divergences. Flag for SIM145 full check.

## Inputs

- Horndeski Friedmann and scalar EOM (standard references: Kobayashi 2019 "Horndeski theory and beyond" review; Bellini & Sawicki 2014 for α_i definitions)
- Phase 1 solver as starting point
- CLASS with Horndeski support (hi_class branch) — may require new installation
- Planck plikHM, DESI Y1, RSD, GW170817 constraints

## Procedure

1. Install and validate hi_class or equivalent Horndeski-capable Boltzmann solver. Test that c₃ = 0 reproduces Phase 1 SIM90.
2. Derive the modified Friedmann and scalar EOM symbolically (or use established Horndeski library). Document derivation in RESULT.md.
3. For G₃(Ψ) = c₃·Ψ and G₃(Ψ) = c₃·Ψ², scan c₃ over 4 decades.
4. For each c₃:
   - Background evolution: compute H(z), w_eff(z), F_eff(z)
   - Evaluate α_K, α_B, α_M, α_T and check stability
   - If stable: compute CMB, BAO, RSD χ²
   - If unstable: record and move on
5. Identify best point minimizing joint χ².
6. For best candidate, explicitly verify GW170817 c_T bound: |c_T − c|/c < 6 × 10⁻¹⁵.

## Success criteria (all must hold)

- (a) Background stable: all α_i stability conditions satisfied
- (b) c_T = c verified to GW170817 precision
- (c) CMB preserved: 100θ* ∈ [1.0408, 1.0414], Δ(−2 ln L_plik) < +5
- (d) DESI tension < 2σ
- (e) RSD fσ₈ shape χ²/N < 1.5 (carryover concern from SIM128; check here)
- (f) No tachyon, Ward identity preserved at tree level
- (g) Phase 1 passes preserved (BBN, Solar System — verify γ with screening)

If (a)-(g) all hold: PASS (conditional on SIM145 full UV recheck, which is more urgent here than for SIM140).

## Failure modes to watch

- **Horndeski complexity:** easy to implement incorrectly. Validate against known Galileon cosmology results (e.g., De Felice & Tsujikawa 2010) before running on CMSTG.
- **Stability vs fit tradeoff:** Horndeski parameter space often has narrow "stable AND fits data" regions. A scan that only looks at χ² will miss that nominally-best points are unstable. Check stability FIRST, fit SECOND.
- **c_T bound is extremely tight:** 6 × 10⁻¹⁵. Any mechanism relying on G₅ is dead. G₃ alone is safe but verify numerically; sanity check against SIM108 which passed this bound for Phase 1.
- **Solar System screening:** Vainshtein screening from Galileon can interact strangely with chameleon screening from Phase 1 (SIM107). Verify γ = 1 ± 2.3 × 10⁻⁵ at Cassini sensitivity.
- **Running this sim on a laptop:** hi_class with MCMC is expensive. Consider cloud compute or narrow the c₃ grid based on a coarse-grain pass first.

## Deliverables

- `output.json` with best c₃, G₃ form, full stability report, all χ²
- `RESULT.md` with EOM derivation, scan results, stability analysis, GW check
- Figures: α_i(z), w_eff(z), F_eff(z), H(z) vs data
- If PASS: draft Paper VIII section (this is the highest-impact PASS outcome possible in Phase 4)

## Estimated time

Large. Most complex of the Tier 2 sims. New Lagrangian sector, new solver requirements, harder stability analysis.

## Deferral rule

If SIM140 or SIM141 has already produced a clean PASS (all criteria met, SIM144/145/146 in progress), **defer SIM142 to a future cycle**. Document the deferral in PHASE4_ROADMAP.md. The mechanism-hunt goal is a working mechanism, not all possible mechanisms.
