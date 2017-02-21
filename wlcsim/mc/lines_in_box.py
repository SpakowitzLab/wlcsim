import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import pandas as pd
from ..utils import path as wpath
import os
from numba import jit
import multiprocessing
import pscan

SMALL_NUM = np.power(10.0, -8)

@jit
def dsegment(xinit0, xfinal0, xinit1, xfinal1):
    """Find distance between two segments defined by initial
    and ending position vectors."""
    u = xfinal0 - xinit0
    v = xfinal1 - xinit1
    w = xinit0 - xinit1
    a = np.dot(u,u)
    b = np.dot(u,v)
    c = np.dot(v,v)
    d = np.dot(u,w)
    e = np.dot(v,w)
    D = a*c - b*b
    sc = D; sN = D; sD = D;
    tc = D; tN = D; tD = D;

    if D < SMALL_NUM:
        sN = 0.0         # force using point P0 on segment S1
        sD = 1.0         # to prevent possible division by 0.0 later
        tN = e
        tD = c
    else:                 # get the closest points on the infinite lines
        sN = (b*e - c*d)
        tN = (a*e - b*d)
        if sN < 0.0:        # sc < 0 => the s=0 edge is visible
            sN = 0.0
            tN = e
            tD = c
        elif sN > sD:  # sc > 1  => the s=1 edge is visible
            sN = sD
            tN = e + b
            tD = c
    if tN < 0.0:            # tc < 0 => the t=0 edge is visible
        tN = 0.0
        # recompute sc for this edge
        if -d < 0.0:
            sN = 0.0
        elif -d > a:
            sN = sD
        else:
            sN = -d
            sD = a
    elif tN > tD:      # tc > 1  => the t=1 edge is visible
        tN = tD
        # recompute sc for this edge
        if (-d + b) < 0.0:
            sN = 0
        elif (-d + b) > a:
            sN = sD
        else:
            sN = (-d +  b)
            sD = a
    # finally do the division to get sc and tc
    sc = sN / sD if np.abs(sN) > SMALL_NUM else 0.0
    tc = tN / tD if np.abs(tN) > SMALL_NUM else 0.0

    # get the difference of the two closest points
    dP = w + (sc * u) - (tc * v)  # =  S1(sc) - S2(tc)

    return np.linalg.norm(dP)

@jit
def uniform_segment_from_unit_cube():
    """Draws uniformly random segments from (unit) cube (in first quadrant)"""
    # sample a unit vector perpendicular to the xy-plane
    # using the "three normals" trick to generate a point on a sphere, then
    # ignoring the sign of the third component
    face_dx = np.random.randn(3)
    face_dx = face_dx/np.linalg.norm(face_dx)
    face_dx[2] = np.abs(face_dx[2]) # ignore sign
    # draw one of the six faces of the cube
    face = np.floor(np.random.rand(1)*6)
    # sample a point uniformly on that face
    face_x = np.random.rand(2)
    xinit = np.zeros(3)
    dxinit = np.zeros(3)
    # transform xinit and face_dx into face's coordinates
    if face == 0: # in xz-plane
        xinit[0] = face_x[0]
        xinit[2] = face_x[1]
        dxinit[1] = face_dx[2]
        dxinit[0] = face_dx[0]
        dxinit[2] = face_dx[1]
    elif face == 1: # in yz-plane
        xinit[1] = face_x[0]
        xinit[2] = face_x[1]
        dxinit[0] = face_dx[2]
        dxinit[1] = face_dx[0]
        dxinit[2] = face_dx[1]
    elif face == 2: # in xy-plane
        xinit[0:2] = face_x
        dxinit = face_dx
    elif face == 3: # opposite xz-plane
        xinit[0] = face_x[0]
        xinit[1] = 1.0
        xinit[2] = face_x[1]
        dxinit[1] = -face_dx[2]
        dxinit[0] = face_dx[0]
        dxinit[2] = face_dx[1]
    elif face == 4: # opposite yz-plane
        xinit[0] = 1.0
        xinit[1] = face_x[0]
        xinit[2] = face_x[1]
        dxinit[0] = -face_dx[2]
        dxinit[1] = face_dx[0]
        dxinit[2] = face_dx[1]
    elif face == 5: # opposite xy-plane
        xinit[0:2] = face_x
        xinit[2] = 1.0
        dxinit = face_dx
        dxinit[2] = -dxinit[2]
    else:
        raise ValueError('internal error: face \\notin {0,1,2,3,4,5}')
    # now get collision point with other side of cube, by pretending our
    # segment is a path and asking at what time the path intersects the cube
    # if dx_i < 0, then x_i + dx_i*t exits the cube when x_i == 0
    # if dx_i > 0, then at x_i == 1
    tout = -xinit/dxinit # waste a couple cycles
    is_dx_pos = dxinit > 0
    if np.any(is_dx_pos):
        tout[is_dx_pos] = np.min((1.0 - xinit[is_dx_pos])/dxinit[is_dx_pos])
    xfinal = xinit + np.min(tout)*dxinit
    if np.any(xfinal > 1.0):
        raise ValueError('internal error: xfinal should be inside unit cube')
    # get the segment initial and final locations
    return xinit, xfinal

@jit
def uniform_segment_from_cube(a):
    """Draws uniformly random segments from cube (in first quadrant) defined by
    its side length."""
    xi,xf = uniform_segment_from_unit_cube()
    return a*xi, a*xf

def volume_of_cylinder(xi, xf, a=1):
    """Calculat the volume of a cylinder with center running from xi to xf
    inside a cube of side length a in the 1st quadrant."""
    #TODO: implement, then replace calculation fo phi in monte_carlo with call
    # to this guy
    # correct method:
    # cases:
    #     1) hitting no edges - naive calculation works, same "in" as out
    #     2) hitting one edge - small modification
    #     3) hitting two+edges - ignored for now, hope these events are "rare"
    # case 1) d^2*pi*L is the correct formula. draw the "extra" volume you'd be
    # adding and subtracting from the "corners" that "jut out" of the box to
    # convince yourselves that they cancel out and you get back the naive
    # formula.
    # case 2) if you drew the first case out, you saw that the intersection of
    # the cylinder with the box edge is an ellipse, with major and minor axes
    # determined by the angle of intersection as (r|sec(\theta)|, r),
    # respectfully, where theta is angle from the normal to the wall to the
    # normal of the cylinder's natural cap. thus, you can consider case two as
    # the same as case 1, but subtracting the part ommited by the "collision"
    # with the wall not intersected by the center of the cylinder. since
    # finding the equation for this volume to be omitted seems hard in general,
    # since the ellipse's axes will in general not line up with the box edge,
    # we approximate by an intersection that is instead a circle with the
    # average of the major and minor axes as
    pass

def calculate_sphericity(lines):
    vectors = [line[1] - line[0] for line in lines]
    # vectors = [vector/np.norm(vector) for vector in vectors]
    outers = [np.outer(a, a) for a in vectors]
    inners = [np.inner(a, a) for a in vectors]
    S = np.zeros((3,3))
    Z = np.sum(inners)
    for s in outers:
        S += s
    S = S/Z
    eig = np.linalg.eig(S)
    evals = eig[0]
    evals.sort()
    return evals, (evals[0]+evals[1])*3/2


def plot_segment(xinits, xfinals):
    """Plot a list of uniformly generated segments."""
    # xfinals = np.stack(xfinals, axis=-1)
    # # 2-by-3-by-num_segments paths arrays
    # paths = np.stack(xinits, xfinals, axis=0)
    # easier than above, just plot one by one
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    for i in range(len(xinits)):
        xinit = xinits[i]
        xfinal = xfinals[i]
        ax.plot([xinit[0], xfinal[0]], [xinit[1], xfinal[1]], [xinit[2], xfinal[2]])
    ax.set_xlim(0,1)
    ax.set_ylim(0,1)
    ax.set_zlim(0,1)
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    ax.set_zlabel('z')

def detect_collision(xi, xf, xi_new, xf_new, d, a=1, periodic=False):
    """Are cylinders defined by start and end centerpoint sof xi,xf and
    xi_new,xf_new, each of diameter d, overlapping? use periodic boundary
    conditions in the box of side length a centered in the first quadrant if
    requested."""
    if not periodic:
        return dsegment(xi, xf, xi_new, xf_new) < d
    # one for each of the 26 cubes adjacent to you that can affect you
    multipliers = [(i-1,j-1,k-1) for i in range(3) for j in range(3) for k in range(3)]

# alternate collision detection, more accurate:
# first, if the point of closest approach of the cylinders' centers is interior
# to the cube, then the #TODO

def monte_carlo_with_removal(nsteps, k_remove, d, a=1, lines=None, phi=None,
                             periodic=False):
    """Perform nsteps MC steps. Half steps try to remove a uniformly random
    line with probability exp(-k_remove). Other half attempts (nsteps times)
    to insert a random line into a box of side length
    a. If the line overlaps with any other line (in the sense that they
    represent the centers of cylinders of radius d), reject the new line.

    Returns table of line insertion attempts and stats about the box at that
    time with a column specifying whether or not the move was accepted.

    TODO: calculate actual volume being added to box instead of just using the
    length*d^2*pi approximation, which only works for thin d."""
    if lines is None:
        lines = []
    num_lines = len(lines)
    if phi is None:
        phi = 0 # volume fraction of cylinders in box
    move_history = pd.DataFrame(index=np.arange(0, nsteps, 1),
                                columns=['xi1', 'xi2', 'xi3', 'xf1', 'xf2',
                                         'xf3', 'phi', 'num_lines',
                                         'did_succeed', 'is_removal'])
    for i in range(nsteps):
        if np.random.rand(1) < 0.5:
            num_lines = len(lines)
            if num_lines == 0:
                # record an unsuccessful removal move
                move_history.loc[i] = [np.nan, np.nan, np.nan,
                                       np.nan, np.nan, np.nan,
                                       phi, num_lines, False, True]

                continue
            # choose a random line in the box
            remove_ind = int(np.random.rand(1)*num_lines)
            xi_rem, xf_rem = lines[remove_ind]
            L = np.linalg.norm(xf_rem - xi_rem)
            did_succeed = False
            if np.random.rand(1) < np.exp(-k_remove*L):
                did_succeed = True
                # delete that line
                line = lines.pop(remove_ind)
                num_lines -= 1
                phi -= np.pi*d*d*np.linalg.norm(xf_rem - xi_rem)
            move_history.loc[i] = [xi_rem[0], xi_rem[1], xi_rem[2],
                                   xf_rem[0], xf_rem[1], xf_rem[2],
                                   phi, num_lines, did_succeed, True]
        else:
            xi_new, xf_new = uniform_segment_from_cube(a)
            did_succeed = True
            for xi, xf in lines:
                if detect_collision(xi, xf, xi_new, xf_new, d, a, periodic):
                    did_succeed = False
                    break
            move_history.loc[i] = [xi_new[0], xi_new[1], xi_new[2],
                                xf_new[0], xf_new[1], xf_new[2],
                                phi, num_lines, did_succeed, False]
            if did_succeed:
                num_lines += 1
                phi += np.pi*d*d*np.linalg.norm(xf_new - xi_new)
                lines += [(xi_new, xf_new)]
    numeric_columns = ['xi1', 'xi2', 'xi3', 'xf1', 'xf2', 'xf3', 'phi',
                       'num_lines']
    move_history[numeric_columns] = move_history[numeric_columns].apply(pd.to_numeric)
    move_history['did_succeed'] = move_history.did_succeed.astype(bool)
    move_history['is_removal'] = move_history.is_removal.astype(bool)
    return lines, phi, move_history

#DEPRECATED!
# def monte_carlo(nsteps, d, a=1, lines=None, phi=None):
#     """Attempt (nsteps times) to put a random lines into a box of side length
#     a. If the line overlaps with any other line (in the sense that they
#     represent the centers of cylinders of radius d), reject the new line.

#     Returns table of line insertion attempts and stats about the box at that
#     time with a column specifying whether or not the move was accepted.

#     TODO: calculate actual volume being added to box instead of just using the
#     length*d^2*pi approximation, which only works for thin d."""
#     if lines is None:
#         lines = []
#     num_lines = len(lines)
#     if phi is None:
#         phi = 0 # volume fraction of cylinders in box
#     move_history = pd.DataFrame(index=np.arange(0, nsteps, 1),
#                                 columns=['xi1', 'xi2', 'xi3', 'xf1', 'xf2',
#                                          'xf3', 'phi', 'num_lines',
#                                          'did_succeed'])
#     for i in range(nsteps):
#         xi_new, xf_new = uniform_segment_from_cube(a)
#         did_succeed = True
#         for xi, xf in lines:
#             if dsegment(xi, xf, xi_new, xf_new) < d:
#                 did_succeed = False
#                 break
#         move_history.loc[i] = [xi_new[0], xi_new[1], xi_new[2],
#                                xf_new[0], xf_new[1], xf_new[2],
#                                phi, num_lines, did_succeed]
#         if did_succeed:
#             num_lines += 1
#             phi += np.pi*d*d*np.linalg.norm(xf_new - xi_new)
#             lines += [(xi_new, xf_new)]
#     move_history = move_history.apply(pd.to_numeric)

#     return lines, phi, move_history

# this kind of refactoring is not nice because we can't save all the phi's in
# move_history this way
# def try_add_uniform_segment(lines, a=1):
#     xi_new, xf_new = uniform_segment_from_cube(a)
#     did_succeed = True
#     for xi, xf in lines:
#         if dsegment(xi, xf, xi_new, xf_new) < d:
#             did_succeed = False
#             break
#     move_history.loc[i] = [xi_new[0], xi_new[1], xi_new[2],
#                             xf_new[0], xf_new[1], xf_new[2],
#                             phi, num_lines, did_succeed]
#     if did_succeed:
#         line = (xi_new, xf_new)
#         lines += [line]
#         return line
#     return None
def get_all_line_lengths(moves_df):
    return np.sqrt(
        np.power(moves_df['xf3'] - moves_df['xi3'], 2)
      + np.power(moves_df['xf2'] - moves_df['xi2'], 2)
      + np.power(moves_df['xf1'] - moves_df['xi1'], 2)
    )

# original function to run a bunch of sims to get the dependence of the
# acceptance probabilities on phi. we will instead be using the equilibrium
# values from an insertion/removal process to get the acceptance probabilities
# at equilibrium. if you want to come back and still do this, you should modify
# this function to work with the new monte_carlo functions
def get_accept_prob(outdir, nPhi, nL, steps_per_batch, k_remove, d, a=1, num_sims=None):
    """Bin acceptance ratio dependence on phi into nP bins. Restart the
    MC sim every time that the rate that phi is increasing per nS steps is less
    than restart_thresh."""
    if steps_per_batch < 100:
        Warning('nS < 100 will result in bad averaging.')
    if nPhi < 2:
        raise ValueError('nPhi < 2 corresponds to no binning!')
    if nL < 2:
        raise ValueError('nL < 2 corresponds to no binning!')
    # get directory to save to
    outdir = wpath.get_unique_folder(outdir)
    # save index, currently not used
    num_restarts = 0
    save_ind = 0
    # file names to save to
    prob_file = "prob.txt"
    prob_file = os.path.join(outdir, prob_file)
    acc_file = "acc.txt"
    acc_file = os.path.join(outdir, acc_file)
    tot_file = "tot.txt"
    tot_file = os.path.join(outdir, tot_file)
    # initialize bin stuffs
    phiBins = np.linspace(0, 1, nPhi+1)
    phibins_file = os.path.join(outdir, "phibins.txt")
    np.savetxt(phibins_file, phiBins)
    LBins = np.linspace(0, 1*np.sqrt(3), nL+1)
    lbins_file = os.path.join(outdir, "lbins.txt")
    np.savetxt(lbins_file, LBins)
    accepted_moves = np.zeros((nPhi, nL))
    total_moves = np.zeros((nPhi, nL))
    # initialize inter-mc save state
    lines = None
    phi = None
    # initialize restart checker so that we always do a second simulation
    prev_phi_ave = -float('inf')
    while True:
        # perform MC
        lines, phi, new_moves = monte_carlo_with_removal(steps_per_batch,
                                                         k_remove, d, a, lines, phi)
        new_moves['L'] = get_all_line_lengths(new_moves)
        # update ratios
        #TODO: check if np.digitize is faster
        for i in range(nPhi):
            for j in range(nL):
                tot_mask = ((phiBins[i] < new_moves['phi'])
                        & (new_moves['phi'] < phiBins[i+1])
                        & (LBins[j] < new_moves['L'])
                        & (new_moves['L'] < LBins[j+1])
                        & (~new_moves['is_removal']))
                acc_mask = tot_mask & (new_moves['did_succeed'])
                accepted_moves[i,j] += np.sum(acc_mask)
                total_moves[i,j] += np.sum(tot_mask)
        # check if should restart
        phi_ave = new_moves['phi'].mean()
        if phi_ave < prev_phi_ave:
            # reset mc state
            lines = None
            phi = None
            # reset to base case
            prev_phi_ave = -float('inf')
            num_restarts += 1
        else:
            # keep track fo previous phi for convergence criterion
            prev_phi_ave = phi_ave
        # write out current values for rate
        np.savetxt(prob_file, np.divide(accepted_moves, total_moves))
        np.savetxt(acc_file, accepted_moves)
        np.savetxt(tot_file, total_moves)
        save_ind += 1
        if num_sims is not None and num_restarts >= num_sims:
            break
    return accepted_moves, total_moves

def get_equilibrium(steps_per_batch, k_remove, d, a=1, lines=None, phi=None):
    """For the values of d and k_remove provide, run until equilibrium. Here,
    equilibrium means run in batches of steps_per_batch MC steps, then if the
    average phi (packing fraction) in a batch drops below the average phi of
    the previous batch. One more batch is then run to calculate the equilibrium
    value.

    This works since the approach to equlibrium is an exponential (monotone)
    plus noise, so as long as steps_per_batch is chosen large enough that the
    reandom fluctuations during the rise to equilibrium get averaged out.

    For cylinder packign of radius 0.01 in a box of side length 1, a reasonable
    number for steps_per_batch is 1000. For all values tested that return
    reasonable answers, a reasonable number for steps_per_batch is 10000.
    """
    prev_phi_ave = -float('inf')
    lines = []; phi = 0;
    while True:
        lines, phi, moves = monte_carlo_with_removal(steps_per_batch, k_remove,
                                                     d, a, lines, phi)
        moves = moves[~moves.is_removal]
        phi_ave = moves['phi'].mean()
        # if termination condition reached, run again so taht we decorrelate
        # from the artificial lower value of phi
        if phi_ave < prev_phi_ave:
            _, _, moves = monte_carlo_with_removal(steps_per_batch, k_remove,
                                                   d, a, lines, phi)
            moves = moves[~moves.is_removal]
            return moves['phi'].mean()
        prev_phi_ave = phi_ave


def equilibrium_mapper(p):
    return (p['k'], p['d'],
            get_equilibrium(p['s'], p['k'], p['d']))

def scan_equilibriums(steps_per_batch, k_removes, ds, num_cores=8):
    """Use pscan to get a list of equilibrium phi levels as a function of the
    removal probability (i.e. k_remove a.k.a. chemical potential of the bath of
    "sticks" feeding our box filling process) and the diameter/box width ratio."""
    script_name = os.path.basename(__file__)
    #print(script_name + ': Running scan_equilibriums!')
    p = multiprocessing.Pool(num_cores)
    scan = pscan.Scan({'s': [steps_per_batch], 'k': k_removes, 'd': ds})
    print("k_remove\twidth_ratio\tequilibrium_phi")
    for k,d,phi in p.imap_unordered(
            equilibrium_mapper, scan.params(), chunksize=10):
        print(str(k) + '\t' + str(d) + '\t' + str(phi))

def plot_equilibrium_scan(outfile):
    df = pd.read_table(outfile)
    hf = plt.figure()
    ha = hf.add_subplot(111, projection='3d')
    ha.scatter(df.width_ratio, df.equilibrium_phi, df.k_remove)
    ha.set_xlabel('$d$: Cylinder diameter/Box side length')
    ha.set_ylabel('$\phi$: Equilibrium packing density')
    ha.set_zlabel('$\mu$: Chemical potential/length')
    ha.set_title('Equilibrium Packing Densities from MC')
    return ha

def plot_accept_prob(outdir):
    files = os.listdir(outdir)
    prob = None; acc = None; tot = None
    phiBins = None; Lbins = None
    for file in files:
        if file.startswith('prob'):
            prob = np.loadtxt(os.path.join(outdir, file))
        elif file.startswith('acc'):
            acc = np.loadtxt(os.path.join(outdir, file))
        elif file.startswith('tot'):
            tot = np.loadtxt(os.path.join(outdir, file))
        elif file.startswith('phibins'):
            phiBins = np.loadtxt(os.path.join(outdir, file))
        elif file.startswith('lbins'):
            Lbins = np.loadtxt(os.path.join(outdir, file))

    if prob is None or acc is None or tot is None \
            or phiBins is None or Lbins is None:
        raise FileNotFoundError('one of [phi/l]bins, prob, acc, tot not found in ' + outdir)
    # due to a bug, accepted rates not seem to be saved correctly
    acc = tot*prob
    # get bin centers for plotting
    dphiBin = phiBins[1] - phiBins[0]
    phiBinC = phiBins[:-1] + dphiBin/2
    dLbin = Lbins[1] - Lbins[0]
    LbinC = Lbins[:-1] + dphiBin/2
    phimesh = np.zeros((len(phiBinC), len(LbinC)))
    Lmesh = np.zeros((len(phiBinC), len(LbinC)))
    for i in range(len(phiBinC)):
        for j in range(len(LbinC)):
            phimesh[i,j] = phiBinC[i]
            Lmesh[i,j] = LbinC[j]
    X = phimesh; Y = Lmesh; Z = -np.log(prob)/Lmesh;
    Zave = -np.log(prob) # *Lmesh/Lmesh
    is_missing = np.isnan(Zave) | np.isinf(Zave)
    Zave[is_missing] = 0.0
    Lnorm = Lmesh;
    Lnorm[is_missing] = 0.0
    Zave = np.sum(Zave, axis=1)/np.sum(Lnorm, axis=1);
    hf = plt.figure()
    ha = hf.add_subplot(111, projection='3d')
    ha.scatter(X,Y,Z)
    ha.set_xlabel('$\phi$: Density Occupied')
    ha.set_ylabel('$\Delta{}L$: length of inserted stick')
    ha.set_zlabel('-ln[Prob(Accept)]/$\Delta{}L$')
    hf2 = plt.figure()
    ha2 = hf2.add_subplot(111)
    plt.scatter(X,Z,s=1,c=Y)
    plt.colorbar(label='$\Delta{}L')
    ha2.set_xlabel('$\phi$: Density Occupied')
    ha2.set_ylabel('-ln[Prob(Accept)]/$\Delta{}L$')
    plt.plot(X,Zave)
    return ha,ha2





if __name__ == '__main__':


    # k_removes = np.linspace(0.01, 5, 51)
    # ds = np.linspace(0.01, 0.1, 11)
    # scan_equilibriums(20000, k_removes, ds, num_cores=32)


    get_accept_prob('tmp/accept_probs_d_0.01_k_5', nPhi=100, nL=50,
                    steps_per_batch=1000, k_remove=5, d=0.01,
                    num_sims=1000)


