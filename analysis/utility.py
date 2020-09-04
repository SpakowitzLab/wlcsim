#/usr/bin/python 
# npagane | python utilities file

from pathlib import Path
import numpy as np

# set parameters
length_per_bp = 0.332 # length in nm of 1 bp
default_omega = 2*np.pi/10.5 # default twist of DNA
max_bp_wrapped = 147
nucleosome_height = 5.5 # nm
nucleosome_radius = 5.2 # nm
nucleosome_center = np.array([4.8455, -2.4445, 0.6694])
default_dir = str(Path(__file__).parent / Path('..'))

# read in translation and rotation files
nucleosome_tran = np.loadtxt('%s/input/nucleosomeT' %(default_dir))
nucleosome_rot = np.loadtxt('%s/input/nucleosomeR' %(default_dir))

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
        number of basepairs wrapping bead, i.e. 147 for nucleosomes and 1 for DNA beads

    Returns
    ----------
    Uout, Vout, Rout : list
        list of the resultant vectors from rotating and translating Uin, Vin, and Rin
    """
    mat = np.matrix([Vin, np.cross(Uin, Vin), Uin]).T
    Rout = Rin + np.matmul(mat, nucleosome_tran[max_bp_wrapped - wrap_bp])
    linkRot = np.matrix([[np.cos(default_omega*link_bp), -np.sin(default_omega*link_bp), 0.],
                        [np.sin(default_omega*link_bp), np.cos(default_omega*link_bp), 0.],
                        [0., 0., 1.]])
    mat = np.matmul(np.matmul(mat, nucleosome_rot[(3*(max_bp_wrapped - wrap_bp)):(3*(max_bp_wrapped - wrap_bp) + 3),:]), linkRot)
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
