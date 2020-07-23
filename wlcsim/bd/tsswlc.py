r"""
Simulation twistable, stretchable-shearable wormlike chains.

Notes
-----
The ssWLC is implemented as described in `Koslover et. al., Soft Matter, 2013,
9, 7016
<https://pdfs.semanticscholar.org/87e5/64cd9b84db3ddf38a869c14de09fa2e05de2.pdf>`_.
This defines the locations (:math:`r^{(i)}(t_j)`) and tangent vectors
(:math:`u^{(i)}(t_j)`) of each (:math:`i`\th) bead at each time
point (:math:`t_j`).

The twist energy is (loosely speaking) quadratic in the angle that one must
rotate the material normals of one bead to get them to align with the next bead
(after the two beads have been rotated once to make the tangent vectors align).

There are other strategies one might take for this, but we will leave them as
"#TODO" in the code for now.
"""
from ..tabulation import dsswlc_from_del

# #TODO Alternate twist versions
# 1) Based on Hermite polynomials.
# The other two material normals, (the "twistable" part of the chain) are
# implemented independently. Since the ssWLC specifies the locations and
# tangent vectors at each bead, the chain is approximated by a cubic Hermite
# polynomial between each pair of beads. This cubic Hermite polynomial has a
# natural (Levi-Civita) connection relating the material normals at each bead
# to a pair of material normals at the next bead corresponding to "zero twist".
# The twist energy of each segment is quadratic in the angle about the tangent
# vector (:math:`\theta = \int \omega_3 ds`) by which the curve deviates from
# this natural frame.
# 2) Include bend/twist and shear/twist coupling.

# # Fortran parameters from old code
# EB=PARA(1)
# EPAR=PARA(2)
# EPERP=PARA(3)
# GAM=PARA(4)
# ETA=PARA(5)
# XIR=PARA(6)
# XIU=PARA(7)
# LBOX=PARA(8)
# LHC=PARA(9)
# VHC=PARA(10)
# MAGR=sqrt(XIR*2.0/DT)
# MAGU=sqrt(XIU*2.0/DT)
# DRPAR = (R[i+1] - R[i]).U[i]
# DRPERP = (R[i+1] - R[i]) - DRPAR*U[i]
# FRAND = MAGR*rnorm()
# TRAND = MAGU*rnorm(), projected to be perpendicular to u
# GR = (U[i+1] - (U[i+1].U[i])U[i] - ETA*DRPERP
# F = -ETA*EB*GR _ EPAR*(DRPAR - GAM)*U[i] + EPERP*DRPERP
# FELAS[i] += F, FELAS[i+1] -= F
# GU = U[i+1] - U[i] - ETA*DRPERP
# T1 = EB*GU
# T2 = -ETA*EB*DRPAR*DU + ETA*EB*(1 - U[i+1].U[i])*DR
#      - EPAR*(DRPAR - GAM)*DR + EPERP*DRPAR*DRPERP
# TELAS[i] += T1 + T2, TELAS[i+1] -= T1
# actual force/torque is /= XI[U/R]

def sswlc(N, L, lp, t, t_save):
    """
    Simulate a stretchable-shearable WLC.

    Parameters
    ----------
    N : int
        The number of beads to use.
    L : float
        The length of the chain.
    lp : float
        The persistence length of the underlying WLC being approximated. (Must
        be same units as *L*.
    t : array_like
        The times for which to simulate the polymer's motion.
    t_save : (M, ) array_like
        The subset of *t* for which we should save the output.

    Returns
    -------
    r : (M, N, 3) array of float
        The positions of each bead at each save time.
    u : (M, N, 3) array of float
        The tangent vectors to the underlying wormlike chain at each bead.
    """
    L0 = L / (N - 1)
    delta = L0/lp
    e_b, gam, e_par, e_perp, eta, xi_u, _ = dsswlc_from_del(delta)
    # same normalization done in the fortran code for now
    e_b = lp*e_b/(delta*lp)
    e_par = e_par/(delta*lp*lp)
    e_perp = e_perp/(delta*lp*lp)
    gam = delta*lp*gam
    eta = eta/lp
    xi_u = xi_u*delta
    xi_r = delta
    # Lena's estimate of minimum dt possible
    max_dt = (1/2)*xi_u/(e_perp*gam**2)
    # e_twist = lt/(delta/lp)
    if np.any(np.diff(t) > max_dt):
        raise ValueError(f"Maximum recommended time step is {max_dt}")
    return _jit_sswlc_lena(N, L, lp, t, t_save, e_b, gam, e_par, e_perp,
                           eta, xi_u, xi_r)


# @jit(nopython=True)
def _jit_sswlc_lena(N, L, lp, t, t_save, e_b, gam, e_par, e_perp,
                    eta, xi_u, xi_r):
    """
    TODO: test
    """
    rtol = 1e-5
    # initial position
    r0, u0 = _init_circle(N, L)
    # pre-alloc output
    if t_save is None:
        t_save = t
    r = np.zeros(t_save.shape + r0.shape)
    u = np.zeros(t_save.shape + r0.shape)
    # setup for saving only requested time points
    save_i = 0
    if t[0] == t_save[save_i]:
        r[0] = r0
        u[0] = u0
        save_i += 1
    # at each step i, we use data (r,t)[i-1] to create (r,t)[i]
    # in order to make it easy to pull into a new functin later, we'll call
    # t[i-1] "t0", old r (r[i-1]) "r0", and t[i]-t[i-1] "h".
    for i in range(1, len(t)):
        h = t[i] - t[i-1]
        F_brown = np.sqrt(2*xi_r/h)*np.random.randn(*r0.shape)
        T_brown = np.sqrt(2*xi_u/h)*np.random.randn(*u0.shape)

        F_elas, T_elas = ft_elas_sswlc(r0, u0, e_b, gam, e_par, e_perp, eta,
                                       xi_u, xi_r)
        dr1 = (F_elas + F_brown)/xi_r  # slope at beginning of time step
        du1 = (T_elas + T_brown)/xi_u
        for j in range(N):
            du1[j] -= (du1[j]@u0[j])*u0[j]
        r_est = r0 + dr1*(h/2)  # first estimate at midpoint
        u_est = u0 + du1*(h/2)
        u_est = u_est/np.linalg.norm(u_est, axis=1)[:, None]

        F_elas, T_elas = ft_elas_sswlc(r_est, u_est, e_b, gam, e_par, e_perp,
                                       eta, xi_u, xi_r)
        dr2 = (F_elas + F_brown)/xi_r  # first estimate slope at midpoint
        du2 = (T_elas + T_brown)/xi_u
        for j in range(N):
            du2[j] -= (du2[j]@u_est[j])*u_est[j]
        r_est = r0 + dr2*(h/2)  # second estimate at midpoint
        u_est = u0 + du2*(h/2)
        u_est = u_est/np.linalg.norm(u_est, axis=1)[:, None]

        F_elas, T_elas = ft_elas_sswlc(r_est, u_est, e_b, gam, e_par, e_perp,
                                       eta, xi_u, xi_r)
        dr3 = (F_elas + F_brown)/xi_r  # second estimate slope at midpoint
        du3 = (T_elas + T_brown)/xi_u
        for j in range(N):
            du3[j] -= (du3[j]@u_est[j])*u_est[j]
        r_est = r0 + dr3*h  # estimate at endpoint
        u_est = u0 + du3*h
        u_est = u_est/np.linalg.norm(u_est, axis=1)[:, None]

        F_elas, T_elas = ft_elas_sswlc(r_est, u_est, e_b, gam, e_par, e_perp,
                                       eta, xi_u, xi_r)
        dr4 = (F_elas + F_brown)/xi_r  # estimate slope at next time point
        du4 = (T_elas + T_brown)/xi_u
        for j in range(N):
            du4[j] -= (du4[j]@u_est[j])*u_est[j]

        # average the slope estimates
        r0 = r0 + h * (dr1 + 2*dr2 + 2*dr3 + dr4)/6
        u0 = u0 + h * (du1 + 2*du2 + 2*du3 + du4)/6
        u0 = u0/np.linalg.norm(u0, axis=1)[:, None]

        # save if at a time in t_save
        if np.abs(t[i] - t_save[save_i]) < rtol*np.abs(t_save[save_i]):
            r[save_i] = r0
            u[save_i] = u0
            save_i += 1
            if save_i >= len(t_save):
                break
    return r, u


# @jit(nopython=True)
def _jit_sswlc_clean(N, L, lp, t, t_save, e_b, gam, e_par, e_perp,
                     eta, xi_u, xi_r):
    """
    WARNING: known to be buggy...
    """
    rtol = 1e-5
    # initial position
    r0 = np.zeros((N, 3))
    r0[:, 0] = np.linspace(0, L, N)
    u0 = np.zeros((N, 3))
    u0[:, 0] = 1
    # pre-alloc output
    if t_save is None:
        t_save = t
    r = np.zeros(t_save.shape + r0.shape)
    u = np.zeros(t_save.shape + r0.shape)
    # setup for saving only requested time points
    save_i = 0
    if t[0] == t_save[save_i]:
        r[0] = r0
        u[0] = u0
        save_i += 1
    # at each step i, we use data (r,t)[i-1] to create (r,t)[i]
    # in order to make it easy to pull into a new functin later, we'll call
    # t[i-1] "t0", old r (r[i-1]) "r0", and t[i]-t[i-1] "h".
    for i in range(1, len(t)):
        h = t[i] - t[i-1]
        F_elas, T_elas = ft_elas_sswlc(r0, u0, e_b, gam, e_par, e_perp, eta,
                                       xi_u, xi_r)
        dW_r = np.random.randn(*r0.shape)
        dW_u = np.random.randn(*u0.shape)
        # -1 or 1, p=1/2
        S_r = 2*(np.random.rand() < 0.5) - 1
        S_u = 2*(np.random.rand() < 0.5) - 1
        # D = sigma^2/2 ==> sigma = np.sqrt(2*D)
        F_brown = np.sqrt(2*xi_r/h)*(dW_r - S_r)
        T_brown = np.sqrt(2*xi_u/h)*(dW_u - S_u)
        for j in range(N):
            T_brown[j] -= (T_brown[j]@u0[j])*u0[j]
        # estimate for slope at interval start
        K1_r = (F_elas + F_brown)/xi_r
        K1_u = (T_elas + T_brown)/xi_u
        r1 = r0 + h*K1_r
        u1 = u0 + h*K1_u
        for j in range(N):
            u1[j] = u1[j]/np.linalg.norm(u1[j])
        # estimate for slope at interval end
        F_brown = np.sqrt(2*xi_r/h)*(dW_r + S_r)
        T_brown = np.sqrt(2*xi_u/h)*(dW_u + S_u)
        for j in range(N):
            T_brown[j] -= (T_brown[j]@u0[j])*u0[j]
        F_elas, T_elas = ft_elas_sswlc(r1, u1, e_b, gam, e_par, e_perp, eta,
                                       xi_u, xi_r)
        K2_r = (F_elas + F_brown)/xi_r
        K2_u = (T_elas + T_brown)/xi_u
        # average the slope estimates
        r0 = r0 + h * (K1_r + K2_r)/2
        u0 = u0 + h * (K1_u + K2_u)/2
        for j in range(N):
            u0[j] = u0[j]/np.linalg.norm(u0[j])
        if np.abs(t[i] - t_save[save_i]) < rtol*np.abs(t_save[save_i]):
            r[save_i] = r0
            u[save_i] = u0
            save_i += 1
            if save_i >= len(t_save):
                break
    return r, u
