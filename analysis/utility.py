#/usr/bin/python 
# npagane | python utilities file

import pathlib
import numpy as np

# set parameters
lengthPerBP = 0.332
defaultOmega = 34.3*np.pi/180
defaultDirectory = str(pathlib.Path(__file__).parent.absolute()).split(__file__)[0]+'/../'

# read in translation and rotation files
nucleosomeTran = np.loadtxt(defaultDirectory+'input/nucleosomeT')
nucleosomeRot = np.loadtxt(defaultDirectory+'input/nucleosomeR')

# define bead translation + rotation
def rotateBead(Uin,Vin,Rin,linkBP,wrapBP):
    angle = 2*np.pi/10.5
    mat = np.matrix([Vin, np.cross(Uin, Vin), Uin]).T
    Rout = Rin + np.matmul(mat, nucleosomeTran[147-wrapBP])
    linkRot = np.matrix([[np.cos(angle*linkBP), -np.sin(angle*linkBP), 0.],
                        [np.sin(angle*linkBP), np.cos(angle*linkBP), 0.],
                        [0., 0., 1.]])
    mat = np.matmul(np.matmul(mat, nucleosomeRot[(3*(147-wrapBP)):(3*(147-wrapBP)+3),:]), linkRot)
    Uout = mat[:,2]/np.linalg.norm(mat[:,2])
    Vout = mat[:,0]/np.linalg.norm(mat[:,0])
    return np.squeeze(np.asarray(Uout)), np.squeeze(np.asarray(Vout)), np.squeeze(np.asarray(Rout))
    
# parametrize the double helix and return the backbone, base, backbone positions
# inspired from https://www.slideshare.net/dvidby0/the-dna-double-helix
def DNAhelix(bp,omega=defaultOmega, v=lengthPerBP): # distances=nm, angles=radians
    r = 1.0 
    alpha = 0 
    beta=2.4
    if (type(bp)==int or type(bp)==float):
        bp = np.asarray([bp])
    phos1 = np.zeros(len(bp)*3).reshape([len(bp),3])
    phos2 = np.zeros(len(bp)*3).reshape([len(bp),3])
    base = np.zeros(len(bp)*3).reshape([len(bp),3])
    phos1[:,0] = r*np.cos(omega*bp+alpha); phos2[:,0] = r*np.cos(omega*bp+beta); base[:,0] = 0
    phos1[:,1] = r*np.sin(omega*bp+alpha); phos2[:,1] = r*np.sin(omega*bp+beta); base[:,1] = 0
    z = v*bp; phos1[:,2] = z; phos2[:,2] = z; base[:,2] = z
    return phos1, base, phos2

# get the change in the UV angle
def getUVangle(r1,r2,u1,u2):
    t = np.linspace(0,1,2)
    pt1 = r1+u1*t[0]
    pt2 = r1+u1*t[1]
    pt3 = r2+u2*t[0]
    pt4 = r2+u2*t[1]
    b1 = pt2 - pt1
    b2 = pt4 - pt3
    return np.arccos(np.round(np.dot(b1,b2)/(np.linalg.norm(b1)*np.linalg.norm(b2)),6))
