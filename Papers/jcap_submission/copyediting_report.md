# Copyediting Report — CMSTG Four-Paper Series
## Phase A JCAP Submission Preparation
**Date:** 2026-05-16

---

## Paper I — `paper1_nogo/cmstg_paper1_nogo.tex`

### Applied Fixes

| # | Location | Fix |
|---|----------|-----|
| 1 | `\thanks`, line 27 | `Email:` → `E-mail:` (consistency with Papers II–IV) |
| 2 | `\bibitem{Nicolis2009}` | `arXiv:0811.2399` → `arXiv:0811.2197` (correct Nicolis-Rattazzi-Trincherini 2009 Galileon arXiv ID) |

### Human-Review Flags

| # | Location | Issue |
|---|----------|-------|
| F1 | `tab:p3summary`, lines 229–241 | Overfull hbox ~94 pt overwidth. Table column widths need adjustment (`tabcolsep` reduction or column restructuring). |
| F2 | Abstract, paragraph 3 | "Fourteen simulation-based mechanism tests" — verify count against current SIM_MAP (Phase 5 adds SIM147–SIM153; if Paper I's scope is Paper IV companion, count may be 15–21). |
| F3 | Section headings containing math | `$\sigma$` in headings not wrapped in `\texorpdfstring{$\sigma$}{\unichar{"03C3}}` — may generate PDF bookmark warnings. Low priority. |

---

## Paper II — `paper2_framework/cmstg_paper2_framework.tex`

### Applied Fixes

*None.* Paper II is not affected by Phase 5 results; no copyediting issues reached the fix threshold.

### Human-Review Flags

| # | Location | Issue |
|---|----------|-------|
| F1 | `tab:all-passes`, lines 900–921 | Overfull hbox ~51 pt overwidth. Same remedy as Paper I F1. |
| F2 | Section 5 / SIM88 result | `G_eff/G = 0.999985` (15 ppm deviation). With Λ₀ = 0.003 M_Pl⁻² and Ψ̄ = 2.62 M_Pl, F_eff = ½ + Λ₀Ψ̄² = 0.5 + 0.003×6.8644 ≈ 0.521, giving G_eff/G = 1/(2×0.521) ≈ 0.960 — roughly 4% below G, not 15 ppm. **Potential internal inconsistency.** Either the SIM88 result is for a different parameter set, or the 15 ppm figure refers to something other than the global G_eff/G ratio. Author clarification required before submission. |
| F3 | `\bibitem{Nicolis2009}` | arXiv:0811.2197 — **correct** in this paper (no fix needed; note for completeness). |

---

## Paper III — `paper3_uv/cmstg_paper3_uv.tex`

### Applied Fixes

*None.*

### Human-Review Flags

| # | Location | Issue |
|---|----------|-------|
| F1 | Preamble | `\newcommand{\km}{k_m}` — style inconsistency: Papers I, II, IV write `k_m` directly in math mode; only Paper III defines a macro. Cosmetic; no functional issue. |

---

## Paper IV — `paper4_galactic/cmstg_paper4_galactic.tex`

### Applied Fixes (Copyediting)

| # | Location | Fix |
|---|----------|-----|
| 1 | Line 86 | `confirmed exhaustive across` → `confirmed exhaustively across` |

### Phase 5 Integration Changes (Deliverable 2)

| # | Change | Detail |
|---|--------|--------|
| P1 | Preamble | Added `\theoremstyle{plain}` + `\newtheorem{theorem}{Theorem}` (amsthm environment declarations) |
| P2 | Abstract | Trimmed SIM150 sweep paragraph (76→55 words); replaced paragraph 4 with Phase 5 summary: Option B shape inversion (r=−0.81), Option C mass gap 466×, closure of all three action-extension channels |
| P3 | Section 5 (Minimum New Physics) | Added forward pointer: "Phase 5 simulations (SIM147–SIM153) subsequently test all three options; Section 8 reports the results." |
| P4 | New Section 8 | "Phase 5 Structural Obstructions to Dark Matter Mechanisms" inserted before Conclusions; five subsections (see below) |
| P5 | Section 9 (Conclusions) | Opening paragraph updated; two new bullet items added for SIM151–152 and SIM153; closing paragraph updated to "all mechanism channels exhausted" |
| P6 | Bibliography | Added `\bibitem{Bertotti2003}` (Cassini PPN bound, Nature 425, 2003) and `\bibitem{WilsonCompNote}` (companion short note, in preparation) |
| P7 | Data/code sentence | Updated SIM list to include SIM147–SIM153 |

**New Section 8 subsections:**
- 8.1 Temporal Channel Obstruction (SIM147–SIM148): kernel decay τ_mem ≈ 205,000 yr; Ψ_pre loophole closure (75 scan points, Δχ²_DESI = 0.000)
- 8.2 Coupling-Strength Obstruction: Third Structural No-Go (SIM149–SIM150): Theorem (thm:geff-ceiling); Cassini → Λ₀ < 1.29×10⁻³ M_Pl⁻² vs. Λ₀_req ≥ 0.154 M_Pl⁻² (factor ~119 incompatibility); G_eff/G ceiling 1.041 vs. required 3.115
- 8.3 Shape-Inversion Obstruction (SIM152): Theorem (thm:shapeinversion); monotone β(ρ_b) → E_model↘ while flat curves need E_req↗; Pearson r = −0.81, χ²/dof = 91.9
- 8.4 Condensate Mass-Gap Obstruction (SIM153): three obstructions — 466× mass gap, gap-closing needs r_c ≥ 232 kpc (solid-body profile), soliton requires sign flip in m₀² or λ
- 8.5 Synthesis: Three Channels Exhausted: Theorem (thm:dm-nogo) Structural Irreducibility; DM must come from structurally distinct sector

### Compilation Status

Paper IV compiles cleanly to 17 pages (expanded from 14 pp). All cross-references resolve on second pass. Remaining warnings:
- `Font OT1/cmr/m/scit` substitution — harmless; `\textsc` inside italic theorem environment

### Human-Review Flags

| # | Location | Issue |
|---|----------|-------|
| F1 | `\bibitem{WilsonCompNote}` | Marked "(in preparation, 2026)." Must be updated with journal/arXiv reference before final submission. |
| F2 | Abstract word count | ~235 words — within JCAP 250-word limit, but check after any further edits. |
| F3 | Theorem numbering | Theorems 1, 2 (Section 8.2, 8.3) and Theorem 3 (Section 8.5) number sequentially via amsthm. If Papers I or III theorems are cross-referenced in Paper IV, numbering is local only — no conflict, but worth noting. |

---

## Submission Gap List

The following items must be resolved before uploading to JCAP Editorial Manager.

### Blockers (must fix before submission)

| # | Issue | Papers affected | Action required |
|---|-------|----------------|-----------------|
| G1 | **JCAP document class** | I, II, III, IV | All papers use `\documentclass[12pt,a4paper]{article}` + `geometry`. JCAP requires the `jcap` class (iopart-based). This is a substantial reformatting task. |
| G2 | **ORCID** | I, II, III, IV | All papers have `ORCID: [ORCID-ID]` placeholder. Author must supply actual ORCID. |
| G3 | **Institutional address** | I, II, III, IV | All papers list only "Independent researcher." JCAP typically requires a full mailing address. Author must supply. |
| G4 | **Editor name in cover letter** | cover_letter.tex | `[Editor Name]`, `[ADDRESS LINE 1/2]`, `[City, Country, Postcode]`, and three reviewer slots all contain placeholders. |
| G5 | **WilsonCompNote reference** | IV | Must be updated from "(in preparation)" to actual arXiv or journal reference before/at submission. If the companion note is not yet submitted, a preprint arXiv ID suffices. |

### Non-blocking (should fix before submission)

| # | Issue | Papers affected | Action |
|---|-------|----------------|--------|
| G6 | Overfull hboxes in tables | I (94 pt), II (51 pt) | Adjust `tabcolsep` or restructure columns. Ugly in print but not a rejection reason. |
| G7 | G_eff/G = 15 ppm in SIM88 | II | Author should verify and correct or add explanatory footnote. |
| G8 | BibTeX `.bib` files | I, II, III, IV | Papers use inline `thebibliography`. JCAP accepts this, but a `.bib` + BibTeX workflow is preferred by the editorial system for cross-reference checking. |
| G9 | `\texorpdfstring` in section headings with math | I | Low priority; affects PDF bookmarks only. |
| G10 | Abstract word counts (final check) | I, II, III, IV | JCAP limit is 250 words. Paper IV is ~235 words; others should be spot-checked after any further edits. |

### Submission package inventory

| File | Status |
|------|--------|
| `cmstg_paper1_nogo.pdf` | Ready |
| `cmstg_paper2_framework.pdf` | Ready |
| `cmstg_paper3_uv.pdf` | Ready |
| `cmstg_paper4_galactic.pdf` | Ready (17 pp, Phase 5 integrated) |
| `cmstg_paper1_nogo.tex` | Ready |
| `cmstg_paper2_framework.tex` | Ready |
| `cmstg_paper3_uv.tex` | Ready |
| `cmstg_paper4_galactic.tex` | Ready |
| `paper1_figures/` | Present (verify completeness) |
| `paper2_figures/` | Present (verify completeness) |
| `paper3_figures/` | Present (verify completeness) |
| `paper4_figures/` | Present (verify completeness) |
| `cover_letter.tex` | Template only — placeholders must be filled |
| `cover_letter.pdf` | Not yet compiled |

---

*Report generated during Phase A JCAP submission preparation, 2026-05-16.*
