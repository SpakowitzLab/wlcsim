# NP | risca lab | 01/2020 | visualize chromatin sims

import numpy as np
import sys
import r2pdb
import os

input_folder = '../../data'
output_folder = 'pdb/'
timePts = 110

# gather discretization data and translation/rotation matrix data
disc = np.loadtxt('%s/discretization' %(input_folder))
wrap = disc[0]; bps = disc[1]
angle = 2*np.pi/10.5
nucleosomeTRAN = np.loadtxt('../../input/nucleosomeT') 
nucleosomeROT = np.loadtxt('../../input/nucleosomeR') 

def linkROT(linkBP):
    return np.matrix([np.cos(angle*linkBP), -np.sin(angle*linkBP), 0, np.sin(angle*linkBP), np.cos(angle*linkBP), 0, 0, 0 , 1]).reshape([3,3])

#define topology ('linear' or 'circular')
topo = 'linear'
#define index range for coordinate files
file_inds = range(0,timePts)
os.system('rm -r %s/*' %output_folder)
for ind in file_inds:
    #load file
    r = np.loadtxt('%s/r%sv0' %(input_folder,ind))
    u = np.loadtxt('%s/u%sv0' %(input_folder,ind))
    # find centers
    nNuc = np.sum(wrap>1)
    for i in range(len(r)):
        if (wrap[i] > 1):
            uin = np.asarray(u[i,0:3]); vin = np.asarray(u[i,3:6]); cross = np.cross(uin, vin)
            mat = np.matrix([vin, cross, uin]).reshape([3,3]).T
            tranMat = np.asarray(nucleosomeTRAN[int(147-wrap[i])])
            rtemp = r[i,:] + np.matmul(mat, tranMat)
            r[i,:] = (rtemp + r[i,:]) / 2.0 # reset nuc position
    #get pdb lines
    dna = r2pdb.mkpdb(r, topology = topo)
    #save pdb file
    r2pdb.save_pdb('%s/snap%0.3d.pdb' %(output_folder,ind),dna)
