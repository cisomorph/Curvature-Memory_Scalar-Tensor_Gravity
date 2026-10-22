import numpy as np
src = open('cmstg/Ordered_Simulations/SIM121C/Published_Scripts/sim121c_mcmc.py').read()
exec(src[:src.index('# PRE-SCAN')])
for H0 in (50,60,67.36,75,90):
    print(f"LCDM H0={H0}: 100theta*={theta_star_CMSTG(0.5,-1,0,H0):.5f}  r_s(z_drag)={r_s_CMSTG(0.5,-1,0,H0):.2f} Mpc")
print("Planck 2018 LCDM reference: r_drag = 147.09 Mpc, 100theta* = 1.04110")
