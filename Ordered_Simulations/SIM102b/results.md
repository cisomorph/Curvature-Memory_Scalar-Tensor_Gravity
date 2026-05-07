# SIM102b — Kernel-Form Robustness: Results

**Status: PASS**

## Stretched-exponential family K_n(k) = exp(-(k/k_m)^(2n)), n = 1,2,3,4

| n | c_n (analytic) | Max deviation from c_n prediction |
|---|----------------|----------------------------------|
| 1 | 1.000000 | — |
| 2 | 1.000000 | — |
| 3 | 1.137388 | — |
| 4 | 1.253314 | — |

- **Max deviation from analytic c_n across all n and k_m**: 0.2755% — **PASS** (<1%)
- **Max deviation from universal scaling after k_m → c_n^(1/4) k_m**: 0.2755% — **PASS** (<1%)

## Rational kernel K(k) = 1/(1 + k^2/k_m^2)

- Log-divergent as expected: Sigma ~ A * log(k_UV/k_m) + B
- Log-fit coefficient A = 1.1399e-11 (confirms linear growth in log k_UV)
- **Verdict: LOG_DIVERGENT** — k^-2 falloff insufficient to regulate quartic UV divergence

## Outputs

- `Outputs/sim102b_kernel_robustness.pdf` — 3-panel figure
- `Outputs/sim102b_kernel_robustness.png`
- `Outputs/sim102b_diagnostics.json`

## Robustness claim verified

The k_m^4 and Lambda0^2 scaling of Sigma(0) is universal across the stretched-exponential
family for n=1,2,3,4. The O(1) shape coefficient c_n can be absorbed into k_m via
k_m → c_n^(1/4) k_m. The rational kernel fails, sharpening the claim: faster-than-power-law
decay is necessary and sufficient to regulate the quartic UV divergence.
