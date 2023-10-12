#/usr/bin/python 
# npagane | python utilities file

from pathlib import Path
import numpy as np
#Modules added by ABC
from pylab import *
from scipy.optimize import curve_fit
from scipy.spatial.distance import squareform
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
# set parameters
helix_height = 3.4
bp_per_turn = 10.5
length_per_bp = helix_height/bp_per_turn # length in nm of 1 bp # changed from 0.332 (presumed typo)
default_omega = 2*np.pi/bp_per_turn # default twist of DNA
max_bp_wrapped = 160
nucleosome_height = 5.5 # nm
nucleosome_radius = 5.2 # nm
#nucleosome_center = np.array([4.8455, -2.4445, 0.6694])
default_dir = str(Path(__file__).parent / Path('..'))

# read in translation and rotation files
nucleosome_tran = np.loadtxt('%s/input/nucleosomeT' %(default_dir))
nucleosome_rot = np.loadtxt('%s/input/nucleosomeR' %(default_dir))

# Get nucleosome center
default_wrap = 127 # to do: make dynamic
nucleosome_tran_sub = nucleosome_tran[-1*default_wrap:]
xmid = (np.min(nucleosome_tran_sub[:,0]) + np.max(nucleosome_tran_sub[:,0]))/2
ymid = (np.min(nucleosome_tran_sub[:,1]) + np.max(nucleosome_tran_sub[:,1]))/2
zmid = (np.min(nucleosome_tran_sub[:,2]) + np.max(nucleosome_tran_sub[:,2]))/2
nucleosome_center = np.array([xmid, ymid, zmid])

def rotate_bead(Uin, Vin, Rin, link_bp, wrap_bp):
    """
    Rotate and translate computational bead to account for DNA wrapping nucleosome 
    and expected entry-exit DNA angle
    If the bead is a DNA bead, this function does NOT do any rotations/translations.
    This is essentially just translated from Fortran (src/mc/nucleosome.f90) to Python.
    
    Parameters
    ----------
    Uin : numpy.array 
        orientation vector U of bead, tangent to polymer contour
    Vin : numpy.array 
        orientation vector V of bead, constructed to define twist and perpenducalr to U
    Rin : numpy.array 
        euclidian position of bead
    link_bp : float or int
        discretization of computational bead in bp
    wrap_bp : int
        number of basepairs wrapping bead, i.e. 127 for nucleosomes and 1 for DNA beads

    Returns
    ----------
    Uout, Vout, Rout : list
        list of the resultant vectors from rotating and translating Uin, Vin, and Rin
    """

    # Make a rotational matrix
    mat = np.matrix([Vin, np.cross(Uin, Vin), Uin]).T
    #interpolate for wrap_bp
    ind_down = int(np.floor(wrap_bp))
    ind_up = int(np.ceil(wrap_bp))
    if (ind_up == 0): 
        ratio = 0
    else:
        ratio = wrap_bp/ind_up
    offratio = 1 - ratio
    #max_bp_wrapped - ind_up is equiva
    inter_tran = ratio*nucleosome_tran[max_bp_wrapped - ind_up] + offratio*nucleosome_tran[max_bp_wrapped - ind_down] 
    inter_rot = ratio*nucleosome_rot[(3*(max_bp_wrapped - ind_up)):(3*(max_bp_wrapped - ind_up) + 3),:] \
                    + offratio*nucleosome_rot[(3*(max_bp_wrapped - ind_down)):(3*(max_bp_wrapped - ind_down) + 3),:]
    Rout = Rin + np.matmul(mat, inter_tran)

    # Rotation matrix to propogate the double helix twist
    link_rot = np.matrix([[np.cos(default_omega*link_bp), -np.sin(default_omega*link_bp), 0.],
                        [np.sin(default_omega*link_bp), np.cos(default_omega*link_bp), 0.],
                        [0., 0., 1.]])
    mat = np.matmul(np.matmul(mat, inter_rot), link_rot)
    Uout = mat[: ,2]/np.linalg.norm(mat[:, 2])
    Vout = mat[:, 0]/np.linalg.norm(mat[:, 0])
    return np.squeeze(np.asarray(Uout)), np.squeeze(np.asarray(Vout)), np.squeeze(np.asarray(Rout))
    

def DNAhelix(bp, omega=default_omega, v=length_per_bp):
    """
    Parametrize the double helix and return the phosphate (strand 1), base, and 
    phosphate (strand 2) positions. This function helps to interpolate the "actual" 
    positions of DNA basepairs between coarse-grained computational bead locations.
    Inspired from https://www.slideshare.net/dvidby0/the-dna-double-helix
    
    Parameters
    ----------
    bp : int or float
        number of basepairs to construct a double helix for 
    omega : float, optional
        default : 2 :math:`\pi` / 10.5 bp
        twist of DNA in radians per bp
    v : float, optional
        default : 0.332 nm
        length of DNA per bp in nm

    Returns
    ----------
    phos1, base, phos2 : list
        list of positions of the two phosphate backbone strands and the 
        nitrogenous base positions for the number of basepairs requested

    """     
    r = 1.0 # radius of DNA helix
    alpha = 0 # angular offset for DNA strand 1 
    beta = 2.4 # angular offset for DNA strand 2
    z = v*bp
    phos1 = np.array([r*np.cos(omega*bp + alpha), r*np.sin(omega*bp + alpha), z])
    phos2 = np.array([r*np.cos(omega*bp + beta), r*np.sin(omega*bp + beta), z])
    base = np.array([0, 0, z])
    return phos1, base, phos2

def get_uv_angle(u1, u2):
    """
    Get the UV angle, i.e. the "bond" angle between two consecutive beads.
    For two DNA beads, the UV angle should fluctuate around 0 or multiples of :math:`\pi`
    since DNA should be relatively straight. The UV bond angle for a nucleosome bead should 
    fluctuate around the entry/exit angle of a crystallized nucleosome. 
    
    Parameters
    ----------
    u1 : numpy.array
        orientation vector U of bead 1
    u2 : numpy.array
        orientation vector U of bead 2

    Returns
    ----------
    bond_angle : float
        bond angle between beads 1 and 2 in radians

    """     
    bond_angle = np.dot(u1,u2)
    bond_angle = bond_angle/(np.linalg.norm(u1)*np.linalg.norm(u2))
    # round avoids numerical errors that will result in bond_angle > 1 
    bond_angle = np.round(bond_angle, 6)
    # since bond_angle <=1, we can now take arccos
    return np.arccos(bond_angle)

## Added by Ariana:

### Autocorrelation Analysis
#Input: simulation object with all trials and snapshots (e.g. data=mc.Simulation(...)

#Autocorrelation
# exponential function
def func(x, a, c, d):
    return a*np.exp(-c*x)+d

def calc_autocorr(data, trial_labels, num_snaps=100, max_lag=30, eq = 10):
    num_trials = len(data.returnTrials())
    e2es = np.zeros(shape=(num_trials,num_snaps-eq+1))
#     for i in range(len(trial_labels)):
    for i in range(0,num_trials):
        print(i)
        for snapshot in range(eq,num_snaps+1):
#             print(snapshot-eq)
            data.trials[trial_labels[i]].snapshots[snapshot].pairwiseNucleosomeDistance()
            # Length of pair_dist is ((Num_Nucs)^2)/2 - Num_Nucs 
            pre = data.trials[trial_labels[i]].snapshots[snapshot].pair_dist
#             print(pre)
            # Put in square form to see as a Num_Nucs by Num_Nucs
            # Grab the last value in the first row of squareform array
            e2e = squareform(pre)[0][-1]
#             print(e2e)
            e2es[i,snapshot-eq] = e2e

    autocorrs_arr = np.zeros(shape=(num_trials,max_lag))
    for trial in range(0,num_trials):
        row = int(trial)
#         print(row)
        snapshots = np.arange(num_snaps)
        autocorrs = []
        for lag in range(0,max_lag):
            dists_no_lag = e2es[row,:][lag:]
            if lag == 0:
                dists_with_lag = dists_no_lag
            else:
                dists_with_lag = e2es[row,:][:-1*lag]
            corr = np.corrcoef(dists_no_lag,dists_with_lag)[0, 1]
            autocorrs.append(corr)
        autocorrs_arr[trial, :] = np.array(autocorrs)

        mean_autocorrelations = autocorrs_arr.mean(axis=0)
    return mean_autocorrelations

def calc_autocorr_all(data, trial_labels, num_snaps=100, max_lag=30, eq = 10):
    """
    Modified version of calc_autocorr that returns pearson correlation values for all trials
    not just the average between trials
    """
    num_trials = len(data.returnTrials())
    e2es = np.zeros(shape=(num_trials,num_snaps-eq+1))
#     for i in range(len(trial_labels)):
    for i in range(0,num_trials):
        print(i)
        for snapshot in range(eq,num_snaps+1):
#             print(snapshot-eq)
            data.trials[trial_labels[i]].snapshots[snapshot].pairwiseNucleosomeDistance()
            # Length of pair_dist is ((Num_Nucs)^2)/2 - Num_Nucs 
            pre = data.trials[trial_labels[i]].snapshots[snapshot].pair_dist
            print(pre)
            # Put in square form to see as a Num_Nucs by Num_Nucs
            # Grab the last value in the first row of squareform array
            e2e = squareform(pre)[0][-1]
#             print(e2e)
            e2es[i,snapshot-eq] = e2e

    autocorrs_arr = np.zeros(shape=(num_trials,max_lag))
    for trial in range(0,num_trials):
        row = int(trial)
#         print(row)
        snapshots = np.arange(num_snaps)
        autocorrs = []
        for lag in range(0,max_lag):
            dists_no_lag = e2es[row,:][lag:]
            if lag == 0:
                dists_with_lag = dists_no_lag
            else:
                dists_with_lag = e2es[row,:][:-1*lag]
            corr = np.corrcoef(dists_no_lag,dists_with_lag)[0, 1]
            autocorrs.append(corr)
        autocorrs_arr[trial, :] = np.array(autocorrs)

        mean_autocorrelations = autocorrs_arr.mean(axis=0)
    return autocorrs_arr, mean_autocorrelations

def fit_exponential(mean_autocorrelations, max_lag = 30):
    ''' Return a, b, and tau such that tau = 1/lambda
    And the following function is fit to the data: f(x) = ae^(-lambda * x) + b
    The x series is the lag intervals in terms of snapshots
    The y series is the Pearson Correlation
    '''
    x = np.array(range(0,max_lag))
    y = mean_autocorrelations
    popt, pcov = curve_fit(func, x, y)
    y_vals_of_fit = func(x, *popt)
    lam = popt[1]
    tau = 1/lam
    a = popt[0]
    b = popt[2]
    return tau, a, b, y_vals_of_fit

def plot_auto_corr(filepath, data, trial_labels, num_snaps=100, max_lag=30, eq = 10):
   mean_autos = calc_autocorr(data, trial_labels, eq=eq, num_snaps=num_snaps)
   tau, a, b, y_vals_of_fit = fit_exponential(mean_autos)
   fig, ax = plt.subplots()
   ax.scatter(range(0,len(mean_autos)),mean_autos, alpha = 0.5)
   ax.plot(range(0,len(mean_autos)),y_vals_of_fit, color='r', linestyle = '-', linewidth = 3)
   ax.set_xlabel('Lag Interval [Snapshots]', size=14)
   ax.set_ylabel('Pearson Correlation', size=14)
   ax.set_title('End-to-End Distance Autocorrelation', size=16)
   ax.text(0.75, 0.7,'a$e^{-\lambda x}+b$ \n a='+str(round(a,2))+'\n \N{greek small letter tau}= '+str(round(tau,3))+'\n b= '+str(round(b, 3)), transform = ax.transAxes, size = 12)
   plt.savefig(filepath+'.pdf')
   plt.show()
   return mean_autos, tau, a, b, y_vals_of_fit

def categorizeLocs(locations, nuc_pos):
    '''
    Categorize positions of centered beads, given the index of the nucleosome position.
   
    Parameters
    ----------
    locations (nd.array): centered bead locations. Shape: num beads x 3
    nuc_pos (int): index of nucleosome bead 


    To Do:
    Generalize for multinucleosomes, only for mononucleosomes as is. 
    '''
    linker1 = locations[:nuc_pos]
    nuc_center = locations[nuc_pos]
    linker2 = locations[(nuc_pos+1):]
    return linker1, linker2, nuc_center

def vectToSeries(vect, spread = 8, j = 100):
    """ Take a 3d vector with an x, y, and z coordinate and returns 2 points where inpute are scaled by the spread"""
    linepts = vect * np.mgrid[-1*spread:spread:2j][:, np.newaxis]
    return linepts

def getPrincipalComp(vector, print_info=0):
    datamean = vector.mean(axis=0)
    uu, dd, vv = np.linalg.svd(vector - datamean)
    pc1 = vv[0]  
    # We could stop here, but the PC could be pointing in the negative direction of what is intended
    # I think this may depend on the directions of the next two PCs, I am not sure
    # Adding test to see if PC1 is pointing in intended direction (same direction as vector)
    # If the vector is pointing in the correct direction, the dist1 should be less than dist 2
    # Make a vector of length spread*2 in the dir of pc1 and center it at the mean of the input vector
    # If PC1 is already in the correct direction, the distance between 
    
    spread = vectToSeries(pc1, spread = 1)
    start_point = spread[0] + datamean
    end_point = spread[1] + datamean

    dist1 = np.linalg.norm(start_point-vector[0])
    dist2 = np.linalg.norm(start_point-vector[-1])
    if dist1 < dist2:
        pass
    elif dist1 > dist2:
        pc1 = -1*pc1
    else:
        raise ValueError('dist 1 vs dist 2 logic issue') 

    return pc1

def getOrthoNormal(vin, uin):
    """
    Compute the orthonormal vector for the given nucleosome

    Parameters
    ----------
    vin : np.ndarrary() of len 3 
        v vector of the last linker bead entering the nuclesome core particle (x, y, z directions)
    
    uin : np.ndarrary() of len 3
        u vector of the last linker bead entering the nuclesome core particle (x, y, z directions)
    
    Returns
    -------
    ortho_norm : np.ndarrary() of len 3 (x, y, z directions)
        Orthonormal vector in the direction of the bottom face of the nucleosome (closest to DNA entry) 
        to the top face of the nuceolsome (closest to DNA exit). In other words, if you use the left hand rule to
        wrap the fingers of your left hand in the direction of the DNA wrapping around the nucleosome core particle,
        your thumb will be pointing in the direction of the ortho_norm
    
    """
    # nucleosome_height = 5.5
    # Take the cross product v and u to define the third axis
    # which will be used in the rotation matrix 
    cross = np.cross(vin, uin)
    
    # mat is a rotation matrix. 
    # It's determinant = 1 and its inverse equals its transpose
    mat = np.matrix([vin, cross, uin]).reshape([3,3]).T
    
    # Nucleosome center computed from wlcsim/input/nucleosomeT
    # TO DO: Change this to new computed values
    center = np.asarray([4.1899999999999995, -2.8882080789053335,  0.22996943061202169]) 
    
    # The y-direction is the superhelical axis  
    # nuc_top represents a point: the center of the top face of the nucleosome 
    nuc_top = center + np.asarray([0, 0.5 * nucleosome_height, 0]) 
    # nuc_bottom represents a point: the center of the bottom face of the nucleosome 
    nuc_bottom = center - np.asarray([0, 0.5 * nucleosome_height, 0]) 
    # Apply rotation matrix to nucleosome center of the top face, and center of the bottom face
    # from the origin. No need to translate in the correct place in 3D space because we just need vector direction
    rot_top = np.matmul(mat, nuc_top)
    rot_bottom = np.matmul(mat, nuc_bottom)
    # Define vector between these two points 
#     ortho_norm = np.reshape(np.array(rot_top - rot_bottom), (-1,1))
    ortho_norm = np.array(rot_top - rot_bottom).flatten()
    return ortho_norm 

def dotproduct(v1, v2):
    return sum((a*b) for a, b in zip(v1, v2))

def length(v):
    return math.sqrt(dotproduct(v, v))

def angle(v1, v2):
    try: 
        arccos_input = dotproduct(v1, v2) / (length(v1) * length(v2))
        # Rounding arcos input because decimal imprecisions cause the attempt to 
        # take arcos of impossible value (e.g. -1.000000002)
        ang = math.acos(np.round(arccos_input,2))*180/np.pi
    except ValueError:
        print("Angle computation is throwing an error", flush = True)
        ang = np.nan
    return ang