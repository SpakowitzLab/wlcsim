#Sample script for generating a set of pdb files

import numpy as np
import sys
import r2pdb
import os

input_folder = '../../data'
output_folder = 'pdb'

#define topology ('linear' or 'circular')
topo = 'circular'
#define index range for coordinate files
file_inds = range(1,101,1)
#Define range of linking numbers
LKs = range(0,4)

#loop over files and save to pdb
for lk in LKs:
    #Create output folder if it does not exists
    os.system('mkdir %s/LK_%s/' %(output_folder,lk))
    #clear output folder
    os.system('rm -r %s/LK_%s/*' %(output_folder,lk))
    for ind in file_inds:
        #load file
        r = np.loadtxt('%s/LK_%s/r%s' %(input_folder,lk,ind))
        #get pdb lines
        lines = r2pdb.mkpdb(r, topology = topo)
        #save pdb file
        r2pdb.save_pdb('%s/LK_%s/snap%0.3d.pdb' %(output_folder,lk,ind),lines)
