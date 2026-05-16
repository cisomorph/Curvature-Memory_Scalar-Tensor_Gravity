# Pre-Data Preparation Work Plan
## Paper VII — Parallel with Phase D, Before DESI Y3

**Purpose:** Reduce the time from DESI Y3 data release to Paper VII submission
from ~30 weeks (if starting cold) to ~10–15 weeks (if Phase VII-1 and
partial VII-2 are pre-completed).

**When to execute:** During Phase D (JCAP editorial process for Papers I–IV),
estimated mid-2026 through early 2027.

**Resources required:** ~4–8 weeks of focused work, parallelizable with
editorial correspondence for Papers I–IV.

---

## Pre-Work Item 1: Complete Phase VII-1 (Candidate Ranking)
**Estimated effort:** 3–4 weeks  
**Target completion:** Before end of 2026

Execute the full Phase VII-1 week-by-week plan from `paper7_roadmap.md`
but using DESI Y1 data (arXiv:2404.03002) as the placeholder:

**Placeholder targets (from SIM144/DESI Y1):**
- Required Δr_s/r_s: +1.5% (Paper I loophole estimate)
- DESI tension level: 2.77σ (SIM144)
- H₀ discrepancy: ~5% vs Planck indirect

**Outputs:**
1. Candidate comparison document (5–10 pages) with DESI Y1 placeholder numbers.
   This document needs only the numbers replaced once Y3 data arrives; structure
   and conclusions are likely stable.
2. Ranking: Candidate A (EDE) > Candidate B (ΔNeff) > Candidate C (Mod. Recomb.)
   UNLESS literature developments by end of 2026 change this assessment.

**Key papers to read and summarize during Pre-Work 1:**
- Poulin, Smith, Karwal & Kamionkowski (2019) [arXiv:1811.04083] — foundational EDE paper
- Smith, Poulin & Amin (2020) [arXiv:1908.06995] — EDE oscillations and CMB
- Murgia, Murgia & Freese (2021) [arXiv:2007.04999] — EDE constraints
- Schöneberg et al. (2022) [arXiv:2207.07045] — H0 tension landscape
- Kamionkowski & Riess (2022) [arXiv:2211.04492] — H0 tension review
- Hill et al. (2020) [arXiv:2003.07355] — EDE constraints from CMB alone
- Niedermann & Slothus (2020) [arXiv:1910.10739] — NEDE alternative
- Ye & Zhang (2020) [arXiv:2001.02451] — Cold EDE

For each: note (a) the mechanism, (b) the r_s shift achieved, (c) the Planck
χ² penalty, (d) the S8 impact, (e) whether it embeds naturally in a
scalar-tensor framework with minimal new coupling.

---

## Pre-Work Item 2: Partial Phase VII-2 — Action Derivation
**Estimated effort:** 2–3 weeks  
**Target completion:** Early 2027

Complete the action derivation and consistency checks 6–7 before data arrives:

**Check 6 (GW170817):** Trivially passes for minimal φ coupling. Write 1-paragraph
proof: "φ has no ∂φ∂φ coupling to metric perturbation h_μν at tree level;
c_T = c is automatic." Document and file.

**Check 7 (BBN):** Compute ρ_φ(z_BBN)/ρ_rad(z_BBN) as a function of {M, f}
using the axion EDE dilution rate. For n = 2: w_φ = 2/3 after phase transition;
ρ_φ ∝ a^{-4(1+w_φ)} = a^{-20/3}; ρ_rad ∝ a^{-4}; ratio ∝ a^{-8/3}, falls
rapidly. This check will pass for any z_c > 1000. Document.

**Checks 1–2 partial:** Implement the modified Friedmann integrator with ρ_φ(z)
using a top-hat approximation for the EDE injection. Compute r_s(M, f) for a
grid of {M, f} at n = 2. This is the SIM154 pre-run. Cache the results; once
DESI Y3 Δr_s target is known, read off the {M, f} values immediately.

Output: A grid of r_s(M, f) values at n = 2, 3 for:
- M/H₀ ∈ [10, 1000] (log-spaced, 20 points)
- f/M_Pl ∈ [0.05, 0.5] (log-spaced, 15 points)
- n ∈ {2, 3}

This is ~600 background integrations, each requiring ~1 second. Total: ~10 minutes.
Implement and run before DESI Y3.

**Check 3 approximation:** Using Schöneberg et al. (2022) Fig. 3 (or equivalent),
tabulate approximate Δχ²_Planck as a function of f_EDE (which maps to M, f from
the r_s grid above). This approximation tells us whether a viable {M, f} exists
that simultaneously achieves the required Δr_s and stays within Planck bounds.
If no such point exists at DESI Y1 tension level (2.77σ), Candidate A is
immediately flagged as marginal for Paper VII.

---

## Pre-Work Item 3: Simulation Infrastructure
**Estimated effort:** 1–2 weeks  
**Target completion:** Before DESI Y3 release

**3a: DESI Y3 likelihood placeholder**  
Set up the DESI BAO likelihood code framework using the Y1 data format. When Y3
data arrives, replace the Y1 data files with Y3 data files. No code changes
should be needed if the data format is consistent (expected, given DESI's
published data formats).

Code location: Same pipeline as SIM111/SIM144.  
New files needed: `desi_y3_bao_likelihood.py` (wrapper stub, data files to be
replaced when Y3 releases).

**3b: EDE background integrator (SIM154 pre-run)**  
Extend the SIM120 background integrator with:
- Klein-Gordon equation: φ'' + 3Hφ' + a² V'(φ) = 0 (primes = d/dN, N = ln a)
- Axion EDE potential: V(φ) = M⁴(1 − cos(φ/f))^n
- ρ_φ = ½φ'²/(a²) + V(φ); p_φ = ½φ'²/(a²) − V(φ)
- Updated Friedmann: H² = (8πG/3)(ρ_r + ρ_m + ρ_Λ + ρ_φ)

Initial conditions: φ_i = π f − δφ (just displaced from top of potential),
φ'_i = 0, initialized at z_i = 10⁶. Solve forward through recombination.

Test against the Poulin+2019 Table 1 benchmark values at f_EDE = 0.123, z_c = 3388.

**3c: AxiCLASS or EDE-CAMB installation**  
Install the public EDE extension to CAMB or CLASS. Run a test at known EDE
parameters (Poulin+2019 best fit) and verify CMB power spectra match published
figures. This ensures the CMB pipeline is operational before DESI Y3 data arrives.

---

## Pre-Work Item 4: Literature Bibliography
**Estimated effort:** 2–3 hours  
**Target completion:** Concurrent with Pre-Work 1

The `bibliography.bib` file in this directory is the starting point. As Pre-Work 1
proceeds, add entries for every paper read. When DESI Y3 arrives, add the Y3
paper and any new relevant papers that appeared in 2026–2027.

Fields to maintain per entry: key, author, title, journal, year, arXiv ID.

---

## Pre-Work Completion Checklist

| Item | Description | Estimated completion | Done? |
|------|-------------|---------------------|-------|
| PW1 | Candidate comparison document | End of 2026 | [ ] |
| PW2a | GW170817 check document | 2026-Q3 | [ ] |
| PW2b | BBN check document | 2026-Q3 | [ ] |
| PW2c | r_s(M,f) grid at n=2,3 (SIM154 pre-run) | 2026-Q4 | [ ] |
| PW2d | Δχ²_Planck vs f_EDE table | 2026-Q4 | [ ] |
| PW3a | DESI Y3 likelihood stub | 2026-Q4 | [ ] |
| PW3b | EDE background integrator (SIM154) | 2026-Q4 | [ ] |
| PW3c | AxiCLASS/EDE-CAMB installed and tested | 2027-Q1 | [ ] |
| PW4 | bibliography.bib complete for pre-data literature | Ongoing | [ ] |

---

## Time Savings Summary

| Phase | Cold start | With pre-data work |
|-------|-----------|-------------------|
| VII-1 (Candidate ranking) | 4 weeks | 0 (done) |
| VII-2 (Derivation + checks 1–2) | 4 weeks | 1–2 weeks (grid look-up) |
| VII-2 (Checks 3–7) | 4 weeks | 2–3 weeks |
| VII-3 (Pipeline + scan) | 8 weeks | 4–6 weeks (pipeline ready) |
| VII-4 (Drafting) | 10 weeks | 8 weeks |
| **Total** | **~30 weeks** | **~15–17 weeks** |

The pre-data work investment (~6–8 weeks) buys back ~12–15 weeks from the
post-data execution timeline.
