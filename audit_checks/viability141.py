import sys, numpy as np
sys.path.insert(0, '/home/aion/CMSTG-fixes/Ordered_Simulations/SIM141'); sys.path.insert(0, '/home/aion/CMSTG-fixes/Ordered_Simulations/common')
import importlib.util
spec = importlib.util.spec_from_file_location('s141', '/home/aion/CMSTG-fixes/Ordered_Simulations/SIM141/run.py')
src = open('/home/aion/CMSTG-fixes/Ordered_Simulations/SIM141/run.py').read()
g = {'__file__': '/home/aion/CMSTG-fixes/Ordered_Simulations/SIM141/run.py', '__name__': 'x'}
exec(src[:src.index('print("=" * 72)')], g)
H0yr = lambda H0: H0/3.0857e19*3.156e7
for form,p,lab in [('constant',{},'Phase 1 const'),('linear',{'beta':0.10},'linear b=0.10'),('tanh',{'a_trans':0.67,'sigma':0.2,'gamma':0.5},'tanh 0.67/0.5'),('tanh',{'a_trans':0.67,'sigma':0.2,'gamma':1.0},'tanh 0.67/1.0')]:
    m = g['model'](form,p); bg = m.theta_matched([g['PSI_BAR']],[0.0],0.6759)
    f0,y0 = bg['f0'],bg['y0']; F0 = bg['F0']
    Fp = m.F_i(f0,0.0)[0]; FN = m.F_N(f0,0.0)
    Ft, Ftp = 2*F0, 2*Fp
    g1 = Ftp**2/(Ft+2*Ftp**2)
    gdot = -((Fp*y0[0]+FN)/F0)*H0yr(bg['H0'])
    print(f"{lab:16} H0={bg['H0']:.2f} Psi0={f0[0]:.3f} F_psi={Fp:.2e} |g-1|={g1:.2e} ({g1/2.3e-5:.1f}x Cassini)  Gdot/G={gdot:+.2e}/yr ({(gdot-7.1e-14)/7.6e-14:+.1f} sigma vs LLR)")
