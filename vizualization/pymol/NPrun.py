# NP | risca lab | 01/2020 | visualize chromatin sims

import numpy as np
import sys
import r2pdb
import os
import sys
#from pymol import cmd

if ('nicolepagane' not in os.getcwd()):
    print('hi, sorry nicole is lazy and this is configured to run on her computer, for chromatin simulations, and for the code version after date jan ~28. if you want to run on your machine, you will need to change/comment out a few things in this script and movie.py.\nx0x0 np')

if (len(sys.argv) <= 1):
    timePts = 110+1
    channel='0'
    channelOG = channel
else:
    timePts = int(sys.argv[1])+1
    channel = sys.argv[2]
    channelOG = channel
    if (channel == 'PT'):
        nodes = np.loadtxt('../../data/nodeNumber')
        nodes = np.vstack((np.linspace(0, np.shape(nodes)[1]-1, np.shape(nodes)[1]), nodes))
        channel = np.asarray(nodes[:,-1], 'int')
    else:
        channel = [channel]*timePts
input_folder = '../../data'
output_folder = 'pdb/'

#define topology ('linear' or 'circular')
topo = 'linear'
#define index range for coordinate files
file_inds = range(0,timePts)
os.system('rm -r %s/*' %output_folder)
for ind in file_inds:
    #load file r
    r = np.loadtxt('%s/r%sv%s' %(input_folder,ind,channel[ind]))
    dna = r2pdb.mkpdb(r, topology = topo)
    #save pdb file
    r2pdb.save_pdb('%s/snap%0.3d.pdb' %(output_folder,ind),dna)

# run pymol
#os.system("~/Applications/pymol/pymol -r movie.py -- " + str(timePts))
os.system("/Applications/PyMOL.app/Contents/MacOS/MacPyMOL -r movie.py -- " + str(timePts) + " " + channelOG)
