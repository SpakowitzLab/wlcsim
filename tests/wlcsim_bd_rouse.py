import matplotlib.pyplot as plt
import numpy as np


def test_burgess_sim():
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

    plt.yscale('log'); plt.xscale('log')
    plt.legend()
    plt.xlabel('time (s)')
    plt.ylabel('MSD (nm^2/s)')

    plt.title('Simulation matches analytical Rouse theory')



def test_rouse_msd_line_approx():
    """Figure out what Dapp is exactly using our analytical result.

    if we use msd_approx(t) = 3*bhat*np.sqrt(Dhat*t)/np.sqrt(3)*1.1283791(6)
    then ``|msd(t) - msd_approx(t)|/msd(t) = np.sqrt(2)/N*np.power(t, -1/2)``.
    ``msd_approx(t) > msd(t)`` in this range, so that means that
    ``msd_approx(t)/(np.sqrt(2)/N*np.power(t, -1/2) + 1) = msd(t)``.

    if we redefine msd_approx with this correction, the new relative error is
    about ``np.sqrt(2)/N/100*t**(-1/2)``? gotta see if this carries over to
    other chain parameters though...it does not...  this time msd(t) is bigger,
    so ``msd_approx = msd_approx(t)/(1 - np.sqrt(2)/N/100*t**(-1/2))``.

    okay actually looks like extra factor is additive?

    for N = 1e8+1; L = 25; b = 2; D = 1;, we have
    msd_approx(t) - msd(t) = 1.01321183e-07

    for N = 1e8+1; L = 174; b = 150; D = 166;, we have
    msd_approx(t) - msd(t) = 5.2889657(3)e-05

    for N = 1e8+1; L = 174; b = 15; D = 166;, we have
    msd_approx(t) - msd(t) = 5.2889657(3)e-06

    for N = 1e8+1; L = 17.4; b = 15; D = 166;, we have
    msd_approx(t) - msd(t) = 5.2889657(3)e-07

    for N = 1e8+1; L = 17.4; b = 15; D = 16.6;, we have
    msd_approx(t) - msd(t) = 5.2889657(3)e-07 (no change)

    for N = 1e7+1; L = 17.4; b = 15; D = 16.6;, we have
    msd_approx(t) - msd(t) = 5.2889657(3)e-06 (no change)

    so the answer is like
    3*bhat*np.sqrt(Dhat*t)/np.sqrt(3)*1.1283791615 - 0.202642385398*b*L/N
    """
    N = 1e8+1; L = 174; b = 150; D = 166; dt = 1e-2; Nt = 1e6; Nt_save = 1e4;
    t = np.arange(0, Nt*dt, dt); t_save = t[::int(Nt/Nt_save)];
    Nhat = L/b; L0 = L/(N-1); Dhat = D*(N)/Nhat; bhat = np.sqrt(L0*b)
    plt.plot(tmsd, np.abs(msd_rouse - 3*bhat*np.sqrt(Dhat*tmsd)/np.sqrt(3)*1.12837916)/msd_rouse, 'k')
    plt.yscale('log')
    plt.xscale('log')
