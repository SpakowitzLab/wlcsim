# NP | risca lab | 01/2020 | visualize chromatin sims

import numpy as np
import sys
import r2pdb
import os
import sys

if ('nicolepagane' not in os.getcwd()):
    print('hi, sorry nicole is lazy and this is configured to run on her computer, for chromatin simulations, and for the code version after date jan ~28. if you want to run on your machine, you will need to change/comment out a few things in this script and movie.py.\nx0x0 np')

if (len(sys.argv) <= 1):
    timePts = 110+1
    channel='0'
else:
    timePts = int(sys.argv[1])+1
    channel = sys.argv[2]
    if (channel == 'PT'):
        nodes = np.loadtxt('../../data/nodeNumber')
        nodes = np.vstack((np.linspace(0, np.shape(nodes)[1]-1, np.shape(nodes)[1]), nodes))
        channel = np.asarray(nodes[:,-1], 'int')
    else:
        channel = [channel]*timePts
input_folder = '../../data'
output_folder = 'pdb/'

# gather translation/rotation matrix data
center = np.asarray([4.492472962875028, -2.223066978638025, 0.3943949777108554])
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
    #load file r and u data
    r = np.loadtxt('%s/r%sv%s' %(input_folder,ind,channel[ind]))
    u = np.loadtxt('%s/u%sv%s' %(input_folder,ind,channel[ind]))
    # load discretization data
    disc = np.loadtxt('%s/d%sv%s' %(input_folder,ind,channel[ind]))
    wrap = disc[0]; bps = disc[1]
    # find centers for nucleosomes
    nNuc = np.sum(wrap>1)
    for i in range(len(r)):
        if (wrap[i] > 1):
            uin = np.asarray(u[i,0:3]); vin = np.asarray(u[i,3:6]); cross = np.cross(uin, vin)
            mat = np.matrix([vin, cross, uin]).reshape([3,3]).T
            tranMat = np.asarray(nucleosomeTRAN[int(147-wrap[i])])
            rtemp = r[i,:] + np.matmul(mat, tranMat)
            r[i,:] = (rtemp + r[i,:]) / 2.0 # reset nuc position
    #get pdb lines
    #dna = r2pdb.mkpdb(r, topology = topo)
    dna = r2pdb.mkpdb(r, wrap, topology = topo)
    #save pdb file
    r2pdb.save_pdb('%s/snap%0.3d.pdb' %(output_folder,ind),dna)

# run pymol
#os.system("~/Applications/pymol/pymol -r movie.py -- " + str(timePts))
os.system("/Applications/PyMOL.app/Contents/MacOS/MacPyMOL -r movie.py -- " + str(timePts))
