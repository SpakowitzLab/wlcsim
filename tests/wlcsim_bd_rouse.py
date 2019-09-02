import matplotlib.pyplot as plt
import numpy as np

# parameters of Burgess simulation
N = 101; L = 17475; b = 15; D = 16667; dt = 2e-4; Nt = 1e6; Nt_save = 1e4;
# runs in a few seconds on my laptop
t = np.arange(0, Nt*dt, dt); t_save = t[::int(Nt/Nt_save)]

import wlcsim
# perform simulation of discretized Rouse polymer (beads on thermal springs)
x = wlcsim.bd.rouse.jit_rouse(N, L, b, D, t, t_save)
ta_msd, count = wlcsim.bd.rouse._get_bead_msd(x, k=51)
# plot results vs naive diffusion comparisons
plt.plot(t_save[1:], ta_msd[1:]/count[1:], 'g', label='Simulation')
plt.plot(t_save[1:], 6*D*t_save[1:], 'k-.' label='Free Monomer, Theory')
plt.plot(t_save[1:], 6*D*t_save[1:]/N, 'k-.' label='0th Rouse Mode, Theory')

# compare to full analytical theory
tmsd = np.logspace(-1, 2, 101)
msd_rouse = wlcsim.analytical.rouse.rouse_mid_msd(tmsd, b, N, D)
plt.plot(tmsd, msd_rouse, 'k', label='Full Rouse Theory')



