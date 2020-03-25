# NP | risca lab | 01/2020 | visualize chromatin sims

import numpy as np
import sys
import r2pdb
import os
import sys
import pandas as pd
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
# define nuc geometry
nucT = pd.read_csv('../../input/nucleosomeT', delim_whitespace=True,header=None)
nucR = pd.read_csv('../../input/nucleosomeR', delim_whitespace=True,header=None)

# define nuc prop function
lengthPerBP = .33
defaultOmega = 34.3*np.pi/180

nucleosomeTran = {}
nucleosomeTran[1] = np.asarray([0., 0., 0.])
nucleosomeTran[147] = np.asarray([2.2877798060913035, -4.828870632041149, -3.297352970952407])

nucleosomeRot = {}
nucleosomeRot[1] = np.eye(3)
nucleosomeRot[147] = np.matrix([[-0.4342885913498671, 0.15622238982734812, -0.8871234324151178], 
                                [-0.24474180849679034, -0.9682619137169162, -0.05069826051232343], 
                                [-0.866888035790706, 0.19509851706737052, 0.4587392527798819]])
def nucleosomeProp(Uin, Vin, Rin, linkBP, wrapBP):
    angle = 2*np.pi/10.5
    mat = np.matrix([Vin, np.cross(Uin, Vin), Uin]).T
    Rout = Rin + np.matmul(mat, nucleosomeTran[wrapBP])
    linkRot = np.matrix([[np.cos(angle*linkBP), -np.sin(angle*linkBP), 0.],
                         [np.sin(angle*linkBP), np.cos(angle*linkBP), 0.],
                         [0., 0., 1.]])
    mat = np.matmul(np.matmul(mat, nucleosomeRot[wrapBP]), linkRot)
    Uout = mat[:,2]/np.linalg.norm(mat[:,2])
    Vout = mat[:,0]/np.linalg.norm(mat[:,0])
    return np.squeeze(np.asarray(Uout)), np.squeeze(np.asarray(Vout)), np.squeeze(np.asarray(Rout))

# inspired from https://www.slideshare.net/dvidby0/the-dna-double-helix
def DNAhelix(bp, omega=defaultOmega, v=0.332): # distances=nm, angles=radians
    r = 1.0 
    alpha = 1 
    beta=2.4 
    if (type(bp)==int or type(bp)==float):
        bp = np.asarray([bp])
    strand1 = np.zeros(len(bp)*3).reshape([len(bp),3])
    strand2 = np.zeros(len(bp)*3).reshape([len(bp),3])
    strand1[:,0] = r*np.cos(omega*bp+alpha); strand2[:,0] = r*np.cos(omega*bp+beta)
    strand1[:,1] = r*np.sin(omega*bp+alpha); strand2[:,1] = r*np.sin(omega*bp+beta)
    z = v*bp; strand1[:,2] = z; strand2[:,2] = z
    return strand1, strand2

def get_UV_angle(r1,r2,u1,u2):
    t = np.linspace(0,1,2)
    pt1 = r1+u1*t[0]
    pt2 = r1+u1*t[1]
    pt3 = r2+u2*t[0]
    pt4 = r2+u2*t[1]
    b1 = pt2 - pt1
    b2 = pt4 - pt3
    return np.arccos(np.round(np.dot(b1,b2)/(np.linalg.norm(b1)*np.linalg.norm(b2)),6))

#define index range for coordinate files
file_inds = range(0,timePts)
os.system('rm -r %s/*' %output_folder)
for ind in file_inds:
    #load file r
    r = np.loadtxt('%s/r%sv%s' %(input_folder,ind,channel[ind]))
    u = np.loadtxt('%s/u%sv%s' %(input_folder,ind,channel[ind]))
    # load discretization data
    disc = np.loadtxt('%s/d%sv%s' %(input_folder,ind,channel[ind]))
    wrap = disc[0]; bps = disc[1]
    dna = r2pdb.mkpdb(r, topology = topo)
    #save pdb file
    #r2pdb.save_pdb('%s/snap%0.3d.pdb' %(output_folder,ind),dna)
    temp = np.sum(np.max(np.asarray([wrap,bps]),0))-1
    bpRes = np.zeros(int(temp)*3*3).reshape([int(temp*3),3])
    indR = 0
    for i in range(len(r)):
        maxBp = bps[i]
        if (wrap[i] != 1):
            Uin = u[i,0:3]; Vin = u[i,3:6]
            for j in range(len(nucT)):
                Rin = np.asarray(nucT.iloc[j,:])
                row = np.zeros(3*3).reshape([3,3])
                Uout, Vout, Rout = nucleosomeProp(Uin, Vin, Rin, 1, 1)
                strand1, strand2 = DNAhelix(j,v=0)
                mat = np.matrix([Vin, np.cross(Uin, Vin), Uin]).T
                strand1 = np.asarray(np.matmul(mat, strand1[0])[0])
                strand2 = np.asarray(np.matmul(mat, strand2[0])[0])
                Vin = Vout; Uin = Uout
                # strand 1 backbone
                row[0,:] = r[i,:] + np.matmul(mat,Rin+strand1[0])
                # strand 2 backbone
                row[2,:] = r[i,:] + np.matmul(mat,Rin+strand2[0])
                # base
                row[1,:] = r[i,:] + np.matmul(mat,Rin+(strand1[0]+strand2[0])/2)
                # plot
                bpRes[indR:indR+3,:] = row
                indR = indR + 3
                #ax.plot3D(row[:,0], row[:,1], row[:,2], '--o', c=cmap[i],alpha=0.5)
        else:
            if (i < len(r)-1):
                omega = get_UV_angle(r[i,:], r[i+1,:], u[i,3:6], u[i+1,3:6])/maxBp
                v = omega/defaultOmega*lengthPerBP
            Uin = u[i,0:3]; Vin = u[i,3:6]; Rin = r[i,:]
            for j in range(int(maxBp)):
                row = np.zeros(3*3).reshape([3,3])
                strand1, strand2 = DNAhelix(j,omega=omega,v=v)
                mat = np.matrix([Vin, np.cross(Uin, Vin), Uin]).T
                strand1 = np.matmul(mat, strand1[0])
                strand2 = np.matmul(mat, strand2[0])
                Uout, Vout, Rout = nucleosomeProp(Uin, Vin, Rin, 1, 1)
                Vin = Vout
                # strand 1 backbone
                row[0,:] = r[i,:]+strand1
                # strand 2 backbone
                row[2,:] = r[i,:]+strand2
                # base
                row[1,:] = r[i,:]+(strand1+strand2)/2
                # plot
                bpRes[indR:indR+3,:] = row
                indR = indR + 3
                #ax.plot3D(row[:,0], row[:,1], row[:,2], '--o', c='grey',alpha=0.75)
    dna = r2pdb.mkpdb(bpRes, topology = topo)
    r2pdb.save_pdb('%s/fine%0.3d.pdb' %(output_folder,ind),dna)
            
# run pymol
#os.system("~/Applications/pymol/pymol -r movie.py -- " + str(timePts) + " " + channelOG)
#os.system("/Applications/PyMOL.app/Contents/MacOS/MacPyMOL -r movie.py -- " + str(timePts) + " " + channelOG)
