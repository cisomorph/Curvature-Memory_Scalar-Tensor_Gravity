import numpy as np, re
from scipy import stats
from scipy.optimize import minimize
src = open('cmstg/Ordered_Simulations/SIM121C/Published_Scripts/sim121c_mcmc.py').read()
# load only the model/likelihood definitions (stop before the scan/MCMC)
exec(src[:src.index('# PRE-SCAN')])

print("== 1. Repo's 'DESI tension' metric vs a real significance")
for chi2,n,k in [(41.4917,6,3)]:
    print(f"  repo: sqrt(chi2/N) = {np.sqrt(chi2/n):.2f} 'sigma'")
    for dof in (n, n-k):
        p = stats.chi2.sf(chi2, dof); print(f"  chi2={chi2:.1f}, dof={dof}: p={p:.1e} -> {stats.norm.isf(p/2):.2f} sigma (two-sided)")

print("\n== 2. Same pipeline, plain LCDM (F0=0.5, w0=-1, wa=0)")
def run(F0,w0,wa):
    H0 = H0_from_theta(F0,w0,wa); Hm = np.array([H_CMSTG(z,F0,w0,wa,H0) for z in DESI_z])
    c2 = np.sum(((Hm-DESI_H)/DESI_s)**2); return H0, c2, Hm
for lab,p in [("LCDM",(0.5,-1,0)),("CMSTG Phase-1 (F0=0.5206, Lambda)",(0.5206,-1,0)),("SIM121C MAP",(0.55998,-0.60026,0.48621))]:
    H0,c2,Hm = run(*p)
    print(f"  {lab:34s} H0_param={H0:6.2f}  H(z=0)={H_CMSTG(0,*p,H0):6.2f}  chi2_DESI={c2:6.2f}  repo-metric={np.sqrt(c2/6):.2f}")
print("  note E2(z=0)=0.5/F0, so H(z=0)!=H0 whenever F0!=0.5 -> F0 is just a rescaling of H")

print("\n== 3. Repo's DESI H(z) table vs DESI Y1 D_H/r_d (values from memory of DESI 2024 III, Table 1 -- please verify)")
desi = {0.510:(20.98,0.61),0.706:(20.08,0.60),0.930:(17.88,0.35),1.317:(13.82,0.42),2.330:(8.52,0.17)}
rd = 147.09
for z,H,s in zip(DESI_z,DESI_H,DESI_s):
    if z in desi:
        DH,e = desi[z]; Hd = c_kms/(DH*rd); sd = Hd*e/DH
        print(f"  z={z:5.3f}: repo H={H:6.1f}+-{s:4.1f}   DESI(D_H, r_d={rd}) H={Hd:6.1f}+-{sd:4.1f}   diff={(H-Hd)/sd:+5.1f} sigma")
    else:
        print(f"  z={z:5.3f}: repo H={H:6.1f}+-{s:4.1f}   DESI Y1 has no D_H at this z (BGS is D_V-only)")

print("\n== 4. Cassini PPN at the canonical field value")
L0, Psi = 0.003, 2.62
F = 0.5 + L0*Psi**2; Fp = 2*L0*Psi
g_std = Fp**2/(F+2*Fp**2)                    # standard massless scalar-tensor: |gamma-1| = F'^2/(F+2F'^2)
a = L0*Psi/(1+2*L0*Psi**2); g_repo = 2*a**2/(1+a**2)   # SIM107's own formula
print(f"  F={F:.4f} (G_eff/G_N today = {0.5/F:.4f}, i.e. {100*(1-0.5/F):.1f}% off, vs README '1 - 15 ppm')")
print(f"  |gamma-1| standard formula = {g_std:.2e}  -> {g_std/2.3e-5:.0f}x Cassini bound 2.3e-5")
print(f"  |gamma-1| SIM107's formula = {g_repo:.2e}  -> {g_repo/2.3e-5:.1f}x Cassini bound")
print(f"  SIM107 itself reports Psi0_critical(Cassini)=1.139 and only evaluates Psi0=1.0; canonical Psi=2.62")
