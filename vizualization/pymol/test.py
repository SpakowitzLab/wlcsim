#Sample script for generating a set of pdb files

import numpy as np
import sys
import r2pdb
import os

input_folder = '../../data'
output_folder = 'pdb/'


LL = float(sys.argv[1])
NB = float(sys.argv[2])
nNuc = float(sys.argv[3])

linkerBead = int(np.round(LL/(LL/((NB-nNuc)/(nNuc+1))))) # for nucs on both ends
#linkerBead = int(np.round(LL/(LL/((NB-nNuc)/(nNuc))))) # for nuc on one end
linkerBead=9

angle = 2*np.pi/10.5

print(linkerBead)

nucleosomeTRAN = np.asarray([2.2877798060913035, -4.828870632041149, -3.297352970952407])
nucleosomeROT = np.matrix([-0.4342885913498671, 0.15622238982734812, -0.8871234324151178, -0.24474180849679034, -0.9682619137169162, -0.05069826051232343, -0.866888035790706, 0.19509851706737052, 0.4587392527798819]).reshape([3,3])
linkBP = 5#linkerBead
linkROT = np.matrix([np.cos(angle*linkBP), -np.sin(angle*linkBP), 0, np.sin(angle*linkBP), np.cos(angle*linkBP), 0, 0, 0 , 1]).reshape([3,3])

#define topology ('linear' or 'circular')
topo = 'linear'
#define index range for coordinate files
file_inds = range(0,110)
os.system('rm -r %s/*' %output_folder)
for ind in file_inds:
    #load file
    r = np.loadtxt('%s/r%sv0' %(input_folder,ind))
    u = np.loadtxt('%s/u%sv0' %(input_folder,ind))
    # find centers
    centers = np.zeros(int(nNuc)*3).reshape([int(nNuc), 3])
    nIter = 0
    for i in range(len(r)):
        if ((i % (linkerBead+1) == 0) and (i != 0) and (i != len(r)-1)):
            uin = np.asarray(u[i,0:3]); vin = np.asarray(u[i,3:6]); cross = np.cross(uin, vin)
            mat = np.matrix([vin, cross, uin]).reshape([3,3]).T
            rtemp = r[i,:] + np.matmul(mat, nucleosomeTRAN)
            #mat = np.matmul(np.matmul(mat, nucleosomeROT), linkROT)
            #uoutNormed = mat[:,2] - mat[:,2].mean(axis=1)
            #uoutNormed = mat[:,2] / np.abs(mat[:,2]).max(axis=1)
            #uout = mat[:,2]/uoutNormed
            centers[nIter,:] = (rtemp + r[i,:]) / 2.0#(1.0*r[i,:] + 1.0*r[i+1,:]) / 2.0
            nIter += 1
    #get pdb lines
    r = r + np.asarray([40,0,-40])
    #centers = centers + np.asarray([40,0,-40])
    lines = r2pdb.mkpdb(r, linkerBead, topology = topo)
    lines2 = r2pdb.mkpdb(centers, 0, topology=topo)
    #save pdb file
    r2pdb.save_pdb('%s/snap%0.3d.pdb' %(output_folder,ind),lines)
    r2pdb.save_pdb('%s/center%0.3d.pdb' %(output_folder,ind),lines2)
    #get pdb lines
    #lines = r2pdb.mkpdb(r, linkerBead, topology = topo)
    #save pdb file
    #r2pdb.save_pdb('%s/snap%0.3d.pdb' %(output_folder,ind),lines)
