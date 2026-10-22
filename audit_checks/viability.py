import sys, numpy as np
sys.path.insert(0, '/home/aion/CMSTG-fixes/Ordered_Simulations/common')
import cmstg_obs as obs
from st_solver import ScalarTensorBG
XI = 1/6
def F134(p,L0): return 0.5*np.exp(-2*L0*p*p)+XI*p
def Fp134(p,L0): return -2*L0*p*np.exp(-2*L0*p*p)+XI
cases = {
 'SIM130 Psi_ini=2.62': (ScalarTensorBG(F=lambda p,N: 0.5+0.003*p*p, F_psi=lambda p,N: 0.006*p), 2.62),
 'SIM130 Psi_ini=3.5':  (ScalarTensorBG(F=lambda p,N: 0.5+0.003*p*p, F_psi=lambda p,N: 0.006*p), 3.5),
 'SIM132 L0=0.2':       (ScalarTensorBG(F=lambda p,N: 0.5*(1-0.4*p*p)+XI*p, F_psi=lambda p,N: -0.4*p+XI), 0.0),
 'SIM134 L0=0.2':       (ScalarTensorBG(F=lambda p,N: F134(p,0.2), F_psi=lambda p,N: Fp134(p,0.2)), 0.0),
}
H0_per_yr = lambda H0: H0/3.0857e19*3.156e7
print(f"{'case':22} {'H0':>6} {'Psi0':>7} {'Psi0p':>8} {'F0':>7} {'F_psi0':>9} {'|g-1|':>9} {'xCassini':>8} {'Gdot/G [1/yr]':>14} {'F(z*)':>7} {'dG/G BBN-now':>12}")
for k,(m,pi) in cases.items():
    bg = m.theta_matched(pi, 0.0, 0.674)
    p0,y0,F0 = bg['psi0'], bg['y0'], bg['F0']
    Fp0 = m.F_psi(p0,0.0)
    Ft, Ftp = 2*F0, 2*Fp0
    g1 = Ftp**2/(Ft+2*Ftp**2)
    gdot = -(Fp0*y0/F0)*H0_per_yr(bg['H0'])
    Fstar = np.interp(1089.92, bg['z'], bg['F']); Fearly = bg['F'][-1]
    print(f"{k:22} {bg['H0']:6.2f} {p0:7.4f} {y0:8.4f} {F0:7.4f} {Fp0:9.5f} {g1:9.2e} {g1/2.3e-5:8.1f} {gdot:14.2e} {Fstar:7.4f} {F0/Fearly-1:+12.4f}")
print("Bounds: Cassini |gamma-1| < 2.3e-5; LLR Gdot/G = (7.1 +- 7.6)e-14 /yr (Hofmann & Mueller 2018); BBN |dG/G| <~ 0.1")
