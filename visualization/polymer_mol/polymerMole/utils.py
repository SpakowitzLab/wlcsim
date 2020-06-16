import math
import numpy as np
from polymerMole.glob import *

def rotateZ(r,theta,center=np.array([32.0,32.0,32.0])):
    R = np.array([[np.cos(theta), -np.sin(theta), 0.0],
                  [np.sin(theta),  np.cos(theta), 0.0],
                  [0.0          ,  0.0          , 1.0]])
    return np.dot(R,r-center)+center
def rotateX(r,theta,center=np.array([32.0,32.0,32.0])):
    R = np.array([[1.0          ,            0.0,            0.0],
                  [0.0          ,  np.cos(theta), -np.sin(theta)],
                  [0.0          ,  np.sin(theta),  np.cos(theta)]])
    return np.dot(R,r-center)+center

def nameSeries(x,letter='A'):
    if x<10:
        return '  '+letter+str(x)
    elif x<100:
        return ' '+letter+str(x)
    elif x<1000:
        return ''+letter+str(x)
    else:
        raise ValueError('Too many colors')

def SetAtomTypeTwoGroups(nbeads, n1):
    atomType = []
    for n in range(nbeads):
        if n<n1:
            atomType.append('  A1')
        else:
            atomType.append('  A2')
    return atomType

def SetAtomTypeSequentially(nbeads,Ncolors):
    atomType=[]
    for n in range(0,nbeads):
        color_n = int(np.floor(n*(Ncolors-1)/nbeads))
        atomType.append(nameSeries(color_n,letter='A'))
    return atomType


def SetAtomTypesByMethylation(METH):
    nbeads=len(METH)
    atomType=[]
    for n in range(0,nbeads):
        if METH[n]==1:
            atomType.append('  A2')
        elif METH[n]==2:
            atomType.append('  A3')
        elif METH[n]==0:
            atomType.append('  A1')
        elif METH[n] == 8:
            atomType.append(' HL8')
        else:
            raise ValueError('Meth type not recognized, '+str(METH[n]))
    return atomType

def SetAtomTypeByPolymer(polymerLengthFile,index):
    if (polymerLengthFile == None):
        raise ValueError("You must supply a polymerLengthFile")
    f = open(polymerLengthFile)
    firstBeads = [1]
    for line in f:
        firstBeads.append(firstBeads[-1]+int(line))
    polymer = 0
    atomType = []
    for ind in index:
        while ind>=firstBeads[polymer+1]:
            polymer=polymer+1
        atomType.append(nameSeries(polymer+1,letter='A'))
    return atomType

def SetAtomTypeVariableMethyaltionLevel(METH):
    atomType = []
    for val in METH:
        atomType.append(nameSeries(val,letter='A'))
    return atomType

import pdb
def ColorCohesin(atomType,leftends,index):
    for ii in range(len(atomType)):
        if leftends[index[ii]] == 1:
            atomType[ii] = '  C1'

def drawConfinement(R, center, file_obj, Ntot, nboundary, nbeads):
    for N in range(0,nboundary):
         atomName='BLCK'
         x=center[0];
         y=center[1]+R*math.sin(2*3.141592*N/nboundary)
         z=center[2]+R*math.cos(2*3.141592*N/nboundary)
         Ntot=Ntot+1;
         print ('HETATM%5d %s %s          %8.3f%8.3f%8.3f  1.00  1.00           C'
               %(Ntot,atomName,resname,x,y,z),file=file_obj)
         if N != 0:
             print('CONECT%5d%5d'%(Ntot,Ntot-1),file=file_obj)

    print('CONECT%5d%5d'%(nbeads+1,Ntot),file=file_obj)
    return Ntot


def drawCube(low, high, file_obj, Ntot):
    atomName='BLCK'
    corners = []
    if hasattr(low, "__iter__"):
        if not hasattr(high, "__iter__"):
            raise ValueError("low and high must be same length")
        for xx in [low[0],high[0]]:
            for yy in [low[1],high[1]]:
                for zz in [low[2],high[2]]:
                    Ntot=Ntot+1
                    print ('HETATM%5d %s %s          %8.3f%8.3f%8.3f  1.00  1.00           C'
                          %(Ntot,atomName,resname,xx,yy,zz),file=file_obj)
                    corners.append([xx,yy,zz,Ntot])
    else:
        for xx in [low,high]:
            for yy in [low,high]:
                for zz in [low,high]:
                    Ntot=Ntot+1
                    print ('HETATM%5d %s %s          %8.3f%8.3f%8.3f  1.00  1.00           C'
                          %(Ntot,atomName,resname,xx,yy,zz),file=file_obj)
                    corners.append([xx,yy,zz,Ntot])

    for aa in range(8):
        for bb in range(8):
            if aa==bb:
                continue
            if np.sum([corners[aa][0]==corners[bb][0], \
                       corners[aa][1]==corners[bb][1], \
                       corners[aa][2] == corners[bb][2]]) != 2:
                if (not (corners[aa][0] == low and corners[bb][0]) == low) or True:
                    continue
            print('CONECT%5d%5d'%(corners[aa][3],corners[bb][3]),file=file_obj)
    return Ntot

def colorFace(low, high, file_obj, Ntot, atomName='BLCK'):
    corners = []
    for xx in [low,high]:
        for yy in [low,high]:
            for zz in [low,high]:
                Ntot=Ntot+1
                print ('HETATM%5d %s %s          %8.3f%8.3f%8.3f  1.00  1.00           C'
                      %(Ntot,atomName,resname,xx,yy,zz),file=file_obj)
                corners.append([xx,yy,zz,Ntot])

    #for aa in range(8):
    #    for bb in range(8):
    #        if aa==bb:
    #            continue
    #        if np.sum([corners[aa][0]==corners[bb][0], \
    #                   corners[aa][1]==corners[bb][1], \
    #                   corners[aa][2] == corners[bb][2]]) != 2 and \
    #            not (corners[aa][0] == low and corners[bb][0]) == low:
    #            continue
    #        print('CONECT%5d%5d'%(corners[aa][3],corners[bb][3]),file=file_obj)
    return Ntot

from numba import jit

@jit(nopython =True)
def fastSquareFilter(invec, halfWidth=0):
    """Applies square filter.

    Args:
        invec (ndarray): Input vector.
        width (int): With of square_filter
    """
    NT=invec.shape[0]
    out = np.zeros(NT)

    total = invec[0]
    lower = 0
    upper = 0
    for center in range(0,NT):
        while upper < min(center+halfWidth, NT-1):
            upper = upper + 1
            total = total + invec[upper]
        while lower < max(0,center-halfWidth):
            total = total - invec[lower]
            lower = lower + 1
        out[center] = total/float(upper-lower+1)
    return out
