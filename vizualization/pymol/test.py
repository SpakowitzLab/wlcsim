#Sample script for generating a set of pdb files

import numpy as np
import sys
import r2pdb
import os

input_folder = '../../data'
output_folder = 'pdb/test'

#define topology ('linear' or 'circular')
topo = 'linear'
#define index range for coordinate files
file_inds = range(1,101)
os.system('rm -r %s/*' %output_folder)
for ind in file_inds:
    #load file
    r = np.loadtxt('%s/r%s' %(input_folder,ind))
    #get pdb lines
    lines = r2pdb.mkpdb(r, topology = topo)
    #save pdb file
    r2pdb.save_pdb('%s/snap%0.3d.pdb' %(output_folder,ind),lines)
