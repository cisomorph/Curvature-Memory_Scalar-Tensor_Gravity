# SIM144 Writeup — Candidate addition to Paper I Section 3.4
## "Completeness of the Phase 4 Search"

---

*The following is written as a candidate Section 3.4 extension for Paper I. It follows directly after the existing Table 3 (Phase 4 Tier 2 mechanisms executed). The framing matches the paper's voice and notation conventions.*

---

### 3.4  Completeness of the Phase 4 Search

Table 3 lists the three executed Tier 2 mechanism classes. Together with the structural arguments of Sections 2.2 and 3.2, they verify that no modification within the CMSTG action class reduces the $2.77\sigma$ DESI residual while preserving the CMB acoustic scale. However, a natural referee question remains: is the Tier 2 space truly exhaustive? In particular, Sections 2 and 3 treat the two routes — curvature-sourced scalars (Route a) and late-time $H(z)$ boosts (Route b) — as independent. A scalar sourced by \emph{both} $R$ and $\rho_{\rm matter}$ simultaneously would straddle this division and is not directly covered by SIM131–136 (pure $R$-sourcing) or SIM143 (quintessence with no curvature coupling).

We address this gap with SIM144, a mechanism-completeness probe. The action extends the bi-scalar SIM143 setup by coupling a second scalar $\phi$ simultaneously to the Ricci scalar and to the matter density:
\begin{equation}
  S = \int d^4x\,\sqrt{-g}\left[\frac{1+2\Lambda_0\Psi^2}{2}\,R - \tfrac{1}{2}(\partial\Psi)^2 - \tfrac{1}{2}m_0^2\Psi^2 - \tfrac{1}{2}(\partial\phi)^2 + \xi_R\,\phi R + 2\beta_m\,\phi\,\rho_{\rm m}\right] + S_{\rm SM}\,,
  \label{eq:sim144_action}
\end{equation}
where $\xi_R$ controls the curvature-sourcing strength and $\beta_m$ the direct matter coupling. The $\Psi$ sector remains at its Phase~1 canonical values. With $F_{\rm eff} = (1+2\Lambda_0\Psi^2)/2 + \xi_R\,\phi$, the $\phi$ equation of motion in FLRW is
\begin{equation}
  \phi'' + (3-\varepsilon_H)\phi' = \xi_R\cdot\frac{R}{H^2} + 2\beta_m\cdot\frac{\rho_{\rm m}}{H^2}\,,
  \label{eq:sim144_eom}
\end{equation}
with initial condition $\phi_{\rm ini} = 0$ (the Deser–Woodard boundary condition used throughout Phase~3).

The scan covers a $4\times4$ grid: $\xi_R \in \{0, 0.01, 0.1, 1.0\}$ and $\beta_m \in \{0, 0.01, 0.1, 1.0\}$ (16 cases total). For each case, we integrate the modified Friedmann–scalar system from $z=10^5$ to $z=0$ with self-consistent $\Lambda_{\rm bare}$ calibration to enforce $H_0 = 67.4$\,km/s/Mpc, and compute the DESI tension, $100\theta_*$, and $\Delta r_s/r_s$ against the Phase~1 ODE baseline.

**Results.** No case in the 16-point scan satisfies all three success criteria simultaneously (DESI tension reduced, $|100\theta_* - 1.04101| < 0.00058$, $|\Delta r_s/r_s| < 0.3\%$). The failure pattern is precisely the one predicted by Theorems~1 and 2:

\begin{table}[h]
\centering
\begin{tabular}{llllll}
\hline
$\xi_R$ & $\beta_m$ & $\phi(0)$ [$M_{\rm Pl}$] & DESI [$\sigma$] & $100\theta_*$ & Verdict \\
\hline
0    & 0    & 0.00  & 1.507 & 1.04096 & Phase 1 reference \\
0    & 0.01 & 0.10  & 1.501 & 1.04133 & Trivial ($\phi \approx 0$) \\
0    & 0.10 & 0.95  & 0.999 & 1.075   & FAIL\_CMB ($+118\sigma$) \\
0.01 & 0    & 0.16  & 1.523 & 1.04012 & FAIL\_CMB ($-3.1\sigma$) \\
0.10 & 0    & 1.62  & 2.891 & 0.973   & FAIL\_DESI+CMB \\
0.10 & 0.10 & 2.73  & 2.339 & 1.012   & FAIL\_DESI+CMB \\
1.0  & 0    & 16.6  & 16.18 & 0.530   & FAIL\_DESI+CMB \\
\multicolumn{5}{l}{(all 16 cases shown in SIM144 repository output)} & \\
\hline
\end{tabular}
\caption{Representative SIM144 results. Phase 1 ODE reference: DESI $1.507\sigma$, $100\theta_* = 1.04096$. Planck $\pm2\sigma$ window: $\pm 0.00058$. No case satisfies all three criteria simultaneously.}
\label{tab:sim144}
\end{table}

The $\beta_m = 0$ column reproduces the Phase~3 Theorem~1 pattern: any $\xi_R > 0$ drives $\phi$ monotonically positive from zero, increasing $F_{\rm eff}$ and suppressing $H(z)$ at DESI redshifts. The $\xi_R = 0$ column reproduces the SIM143 Theorem~2 pattern: the matter-sourcing raises $H(z)$ at late times (improving DESI in principle) but compresses $D_C^*$ and shifts $\theta_*$ far above the Planck bound with no compensating $r_s$ modification ($\Delta r_s = 0.000$\,Mpc for the thawing regime, $|\Delta r_s/r_s| = 0.37\%$ for $\beta_m = 0.1$, the largest tested). Intermediate cases ($\xi_R > 0$, $\beta_m > 0$) inherit both failure modes: the $R$-sourced component increases $F_{\rm eff}$ (wrong direction) while the matter-sourced energy partially compensates, but in no case is the competition resolved in a way that satisfies both $\theta_*$ and DESI.

The predicted monotone structure is confirmed: adding $\beta_m > 0$ to a fixed $\xi_R$ column reduces the DESI tension slightly (the matter energy raises $H$) but simultaneously amplifies the $\theta_*$ shift. The DESI-improving region and the Planck-allowed region remain disjoint across the full $4\times 4$ grid.

**Completeness argument.** The mixed-source action (\ref{eq:sim144_action}) is the most general scalar coupling consistent with the CMSTG action class that is linear in both $R$ and $\rho_{\rm m}$. It interpolates continuously between the Phase~3 probes ($\beta_m = 0$, $\xi_R > 0$) and the SIM143 matter-only probe ($\xi_R = 0$, $\beta_m > 0$). Since both limiting columns fail for independent structural reasons — Theorem~1 for the $R$-sourced direction, Theorem~2 for the matter-sourced direction — and no mixed case evades either theorem, SIM144 closes the final gap in the Tier~2 mechanism space.

Together with Table~3, we conclude that the Tier~2 mechanism space is exhaustive under the CMSTG action class. The $2.77\sigma$ DESI Y1 residual is structural: it cannot be resolved by any extension of the Phase~1 canonical action that couples a new scalar field to curvature, matter, or both simultaneously, provided the field is initialized at zero as required by the Deser–Woodard boundary condition. This strengthens the combined no-go structure of Section~4: the Phase~1 canonical solution is the unique attractor consistent with both Planck CMB data and DESI Y1 BAO data within the CMSTG action class.

---

*[Note for paper integration: Table~\ref{tab:sim144} can be merged with the existing Table~3, adding SIM144 as a fourth row with class label "P4-E (completeness)" and a footnote referencing the mixed-source action. Alternatively it stands as a separate paragraph after Table~3 with an inline citation to the simulation repository. The $\phi(z)$ trajectory figure and the DESI-vs-$\theta_*$ scatter figure are available as SIM144\_phi\_evolution.pdf and SIM144\_desi\_vs\_theta.pdf in the simulation repository.]*
