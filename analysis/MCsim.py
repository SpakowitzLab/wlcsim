#/usr/bin/python 
# npagane | python object to store MC data and simple analysis

import sys
import pathlib
from .utility import *
sys.path.append(defaultDirectory+'../vizualization/pymol')
import vizualization.pymol.r2pdb as r2pdb
import numpy as np
import pandas as pd
import os
import scipy
import scipy.special 
import scipy.spatial
import pickle

class Simulation:
    def __init__(self,path_to_data=defaultDirectory+'data/',trials=[''],trajectories=[''],
                 time_min=0,time_max=110,channel=0):
        # path to data directory
        self.path_to_data = path_to_data
        # store Trials in dictionary
        self.trials = {}
        for trial in trials:
            self.trials[trial] = Trial(self.path_to_data+trial,
                    time_min,time_max,trajectories,channel=channel)
            channel += 1
        # linearized snapshots for multiprocessing
        self.linearized_snapshots = []
    def returnTrials(self):
        return self.trials.keys()
    def returnTrajectories(self):
        tempDict = {}
        for i in self.trials.keys():
            tempDict[i] = self.trials[i].trajectories.keys()
        return tempDict
    def returnSnapshots(self):
        tempDict = {}
        for i in self.trials.keys():
            tempDict[i] = {}
            for j in self.trials[i].trajectories.keys():
                tempDict[i][j] = self.trials[i].trajectories[j].snapshots.keys()
        return tempDict
    def linearizeSnapshots(self):
        for i in self.trials.keys():
            for j in self.trials[i].trajectories.keys():
                for k in self.trials[i].trajectories[j].snapshots.keys():
                    self.linearized_snapshots.append(self.trials[i].trajectories[j].snapshots[k])
    def getCenterBeads(self):
        #mulitprocessing does not work with class methods GRRR >:|
        if (self.linearized_snapshots==None):
            self.linearizeSnapshots()
        for i in self.linearized_snapshots:
            i.centerBeads()
    def getPairwiseNucleosomeDistance(self):
        if (self.linearized_snapshots==None):
            self.linearizeSnapshots()
        for i in self.linearized_snapshots:
            i.pairwiseNucleosomeDistance()
    def getReducedPairwiseNucleosomeDistance(self):
        if (self.linearized_snapshots==None):
            self.linearizeSnapshots()
        for i in self.linearized_snapshots:
            i.reducedPairwiseNucleosomeDistance()
    def getInterpolate(self):
        if (self.linearized_snapshots==None):
            self.linearizeSnapshots()
        for i in self.linearized_snapshots:
            i.interpolate()
    def getRICCbreak(self):
        if (self.linearized_snapshots==None):
            self.linearizeSnapshots()
        for i in self.linearized_snapshots:
            i.RICCbreak()

class Trial:
    def __init__(self,path_to_data,time_min,time_max,trajectories,channel):
        # path to data directory
        self.path_to_data = path_to_data
        # channels
        if (trajectories==['']):
            self.channels = np.linspace(0,len(trajectories)-1,len(trajectories),dtype='int')
        else:
            self.channels = [channel]*len(trajectories)
        # store Trajectories in dictionary
        self.trajectories = {}
        for i,trajectory in enumerate(trajectories):
            self.trajectories[trajectory] = Trajectory(self.path_to_data+trajectory,
                    time_min,time_max,self.channels[i])

class Trajectory:
    # constructor 
    def __init__(self,path_to_data,time_min,time_max,channel):
        # path to data directory
        self.path_to_data = path_to_data
        # set time range
        self.time_min = time_min
        self.time_max = time_max+1 # add plus one here for iteration
        self.PT = False
        # store Snapshots in dictionary
        self.snapshots = {}
        # snapshots stats
        self.end_to_end = []
        self.energies = []
        # check for PT
        if os.path.exists(self.path_to_data+'nodeNumber'):
            nodes = np.loadtxt(self.path_to_data+'nodeNumber')
            nodes = np.vstack((np.linspace(0, np.shape(nodes)[1]-1, np.shape(nodes)[1]), nodes))
            self.channel = np.asarray(nodes[:,-1], 'int')
            self.PT = True
        else:
            self.channel = [channel]*(self.time_max-self.time_min)
        # load snapshots in dictionary
        for i,time in enumerate(range(self.time_min,self.time_max)):
            self.snapshots[time] = Snapshot(self.path_to_data,time,self.channel[i])
    def getEndToEnd(self):
        for time in range(self.time_min,self.time_max):
            self.end_to_end.append(self.snapshots[time].end_to_end)
    def getEnergies(self):
        for time in range(self.time_min,self.time_max):
            self.energies.append(self.snapshots[time].energies)
    def setEquilibrium(self,time):
        for time in range(self.time_min,time):
            self.snapshots.pop(time)
    def playFineMovie(self,path=defaultDirectory+'vizualization/pymol/pdb/',topo='linear',pymol='~/Applications/pymol/pymol'):
        for time in range(self.time_min,self.time_max):
            self.snapshots[time].saveFineGrainedPDB(path=path,topo=topo)
        os.system(pymol + " -r "+defaultDirectory+'vizualization/pymol/'+"movieFine.py -- " 
                    + str(self.time_max-self.time_min) + " " + path)
    def playCoarseMovie(self,path=defaultDirectory+'vizualization/pymol/pdb/',topo='linear',pymol='~/Applications/pymol/pymol'):
        for i,time in enumerate(range(self.time_min,self.time_max)):
             self.snapshots[time].saveCoarseGrainedPDB(path=path,topo=topo)
        if (self.PT):
            os.system(pymol + " -r "+defaultDirectory+'vizualization/pymol/'+"movieCoarse.py -- " 
                        + str(self.time_max-self.time_min) + " PT " + path + " " + self.path_to_data)
        else:
            os.system(pymol + " -r "+defaultDirectory+'vizualization/pymol/'+"movieCoarse.py -- " 
                        + str(self.time_max-self.time_min) + " " + str(self.channel[i]) + " " 
                        + path + " " + self.path_to_data)

class Snapshot:
    # constructor
    def __init__(self,path_to_data,time,channel):
        # path to data directory
        self.path_to_data = path_to_data
        # determine time of snapshot
        self.time = time
        # determine channel
        self.channel = channel
        # load position and u data of computational beads
        self.r = np.loadtxt('%sr%sv%s' %(self.path_to_data,self.time,self.channel))
        self.u = np.loadtxt('%su%sv%s' %(self.path_to_data,self.time,self.channel))
        if (np.shape(self.u)[1]>3):
            temp = self.u[:,0:3]
            self.v = self.u[:,3:6]
            self.u = temp
        else:
            self.v = None
        # load discretization data
        disc = np.loadtxt('%sd%sv%s' %(self.path_to_data,self.time,self.channel), dtype='int')
        self.wrap = disc[0]; self.basepairs = disc[1]
        # assign constants from postion data
        self.n_beads = len(self.r)
        self.end_to_end = np.linalg.norm(self.r[-1,:]-self.r[0,:])
        self.n_bps = np.sum(np.max(np.asarray([self.wrap,self.basepairs]),0))-1+np.sum(self.basepairs[np.linspace(0,self.n_beads-1,self.n_beads,dtype='int')[np.asarray(self.wrap>1)]])
        self.end_to_end_norm = self.end_to_end/(self.n_bps*lengthPerBP)
        # energies
        with open('%senergiesv%s' %(self.path_to_data,self.channel)) as fp:
            for i, line in enumerate(fp):
                if i == 0:
                    # get column specifiers
                    cols = line.strip().split()
                elif i == self.time+1:
                    # get time data
                    energies = line.strip().split()
                    break
        cols = cols[4:]
        energies = energies[2:]
        self.energies = pd.DataFrame([energies],columns=cols)
        # centered beads 
        self.center_r = None
        # pairwise nucleosomes
        self.n_nucs = np.sum(self.wrap>1)
        self.n_pair_nucs = int(scipy.special.comb(self.n_nucs,2))
        self.pair_nucs = None
        self.reduced_pair_nucs = None
        # interpolation/ricc-seq stuff
        self.bps = None
        self.n_pair_bps = int(scipy.special.comb(self.n_nucs,2))
        self.break_length_s1 = None; self.break_location_s1 = None
        self.break_length_b = None; self.break_location_b = None
        self.break_length_s2 = None; self.break_location_s2 = None
    # determine center of beads from regular polygons
    def centerBeads(self,nside=12,type='regular'):
        if (type!='regular'):
            print('can only find center of beads for regular shapes at the moment.')
            return 0
        self.center_r = np.zeros((self.n_beads-1)*3).reshape([self.n_beads-1,3])
        for i in range(self.n_beads-1):
            if (self.wrap[i]>1): # nucleosome
                # make rotation matrix
                uin = np.asarray(self.u[i,:]); vin = np.asarray(self.v[i,:]); cross = np.cross(uin, vin)
                mat = np.matrix([vin, cross, uin]).reshape([3,3]).T
                # center of nucleosome material frame
                center = np.asarray([4.8455, -2.4445, 0.6694])
                # rotate into center of material frame 
                poly = self.r[i,:] + np.matmul(mat, center)
            else: # dna
                # rotate into center of material frame 
                poly = (self.r[i,:]+self.r[i+1,:])/2.0
            self.center_r[i,:] = poly
    # determine the pairwise distances between nucleosomes
    def pairwiseNucleosomeDistance(self):
        if (self.center_r == None):
            self.centerBeads()    
        nucLocs = np.asarray(np.linspace(0,self.n_beads-1,self.n_beads)[self.wrap>1],dtype='int')
        self.pair_nucs = scipy.spatial.distance.pdist(self.center_r[nucLocs,:])
    # determine the reduced pairwise distances between nucleosomes
    # such that n+x where x -> 1..self.n_nucs-1
    def reducedPairwiseNucleosomeDistance(self):
        if (self.pair_nucs == None):
            self.pairwiseNucleosomeDistance()
        self.reduced_pair_nucs = np.zeros((self.n_nucs-1)).reshape([self.n_nucs-1])
        # sum up distances
        iterTemp = 0
        for i in range(self.n_nucs-1):
            for j in range(i+1,self.n_nucs):
                self.reduced_pair_nucs[j-i-1] += self.pair_nucs[iterTemp]
                iterTemp += 1
        # normalize
        for i in range(self.n_nucs-1):
            self.reduced_pair_nucs[i] /= (self.n_nucs-i-1)
    # interpolate atoms into coarse grained chain
    def interpolate(self):
        # rotate DNA strand into material frame
        offset = 0
        self.bps = np.zeros(int(self.n_bps)*3*3).reshape([int(self.n_bps),3,3])
        indR = 0
        for i in range(self.n_beads-1):
            maxBp = self.basepairs[i]
            if (i < self.n_beads-1 and self.wrap[i] == 1):
                # find the rotational/height offset between beads
                omega = getUVangle(self.r[i,:], self.r[i+1,:], self.v[i,:], self.v[i+1,:])/maxBp
                v = defaultOmega/omega*lengthPerBP
                offset = offset + (np.pi-omega*maxBp) % np.pi
            if (self.wrap[i] > 1): # nucleosome
                Uin = self.u[i,:]; Vin = self.v[i,:]; Rin = self.r[i,:]
                # define mat for rotation to lab frame
                mat = np.matrix([self.v[i,:], np.cross(self.u[i,:], self.v[i,:]), self.u[i,:]]).T
                for j in range(len(nucleosomeTran)):
                    row = np.zeros(3*3).reshape([3,3])
                    Uout, Vout, Rout = rotateBead(Uin,Vin,Rin,1,1)
                    omega = getUVangle(Rin, Rout, Vin, Vout)
                    strand1, base, strand2 = DNAhelix(j,omega=omega,v=0)
                    Vin = Vout; Uin = Uout; Rin = np.asarray(nucleosomeTran[len(nucleosomeTran)-1-j,:])
                    # strand 1 backbone
                    row[0,:] = self.r[i,:] + np.matmul(mat,Rin+strand1[0])
                    # strand 2 backbone
                    row[2,:] = self.r[i,:] + np.matmul(mat,Rin+strand2[0])
                    # base
                    row[1,:] = self.r[i,:] + np.matmul(mat,Rin+base[0])
                    # save atoms
                    self.bps[indR,:,:] = row
                    indR = indR + 1
                # add the extruding linker from the nucleosome
                Uout, Vout, Rout = rotateBead(self.u[i,:],self.v[i,:],self.r[i,:],self.basepairs[i],self.wrap[i])
                RnucEnd = Rout; Rin = Rout; Vin = Vout; Uin = Uout
                for j in range(int(maxBp)):
                    row = np.zeros(3*3).reshape([3,3])
                    strand1, base, strand2 = DNAhelix(j,omega=0,v=v)
                    mat = np.matrix([Vin, np.cross(Uin, Vin), Uin]).T
                    strand1 = np.matmul(mat, strand1[0])
                    base = np.matmul(mat, base[0])
                    strand2 = np.matmul(mat, strand2[0])
                    Uout, Vout, Rout = rotateBead(Uin, Vin, Rin, 1, 1)
                    Vin = Vout
                    # strand 1 backbone
                    row[0,:] = RnucEnd+strand1
                    # strand 2 backbone
                    row[2,:] = RnucEnd+strand2
                    # base
                    row[1,:] = RnucEnd+base
                    # save atoms
                    self.bps[indR,:,:] = row
                    indR = indR + 1
                # reset to 0 since nucleosome should have reset offset
                offset = 0
            else: # dna bead
                Uin = self.u[i,:]; Vin = self.v[i,:]; Rin = self.r[i,:]
                for j in range(int(maxBp)):
                    row = np.zeros(3*3).reshape([3,3])
                    strand1, base, strand2 = DNAhelix(float(j+offset),omega=0,v=v)
                    mat = np.matrix([Vin, np.cross(Uin, Vin), Uin]).T
                    strand1 = np.matmul(mat, strand1[0])
                    base =  np.matmul(mat, base[0])
                    strand2 = np.matmul(mat, strand2[0])
                    Uout, Vout, Rout = rotateBead(Uin, Vin, Rin, 1, 1)
                    Vin = Vout
                    # strand 1 backbone
                    row[0,:] = self.r[i,:]+strand1
                    # strand 2 backbone
                    row[2,:] = self.r[i,:]+strand2
                    # base
                    row[1,:] = self.r[i,:]+base
                    # save atoms
                    self.bps[indR,:,:] = row
                    indR = indR + 1
    def RICCbreak(self,cutoff=3.5,noise=50.0): # units in nm or bp
        if (self.bps == None):
            self.interpolate()
        # figure out combinatorial correlated cleaves
        nPair = int(scipy.special.comb(self.n_bps,2))
        indPair = np.zeros(nPair*2).reshape([nPair,2])
        ind = 0
        for i in range(self.n_bps-1):
            extension = self.n_bps-i-1
            indPair[ind:ind+extension,0] = int(i)
            indPair[ind:ind+extension,1] = np.linspace(i+1,i+extension, extension, dtype='int')
            ind += extension
        # find 3d distance on both strands and base
        pairS1 = scipy.spatial.distance.pdist(self.bps[:,0,:])
        pairB = scipy.spatial.distance.pdist(self.bps[:,1,:])
        pairS2 = scipy.spatial.distance.pdist(self.bps[:,2,:])
        cutIndS1 = pairS1 <= cutoff; cutIndB = pairB <= cutoff; cutIndS2 = pairS2 <= cutoff
        indBreakS1 = np.linspace(0,nPair-1,nPair,dtype='int')[cutIndS1]
        indBreakB = np.linspace(0,nPair-1,nPair,dtype='int')[cutIndB]
        indBreakS2 = np.linspace(0,nPair-1,nPair,dtype='int')[cutIndS2]
        fragBreakS1 = indPair[indBreakS1,1]-indPair[indBreakS1,0]
        fragBreakB = indPair[indBreakB,1]-indPair[indBreakB,0]
        fragBreakS2 = indPair[indBreakS2,1]-indPair[indBreakS2,0]
        noiseIndS1 = fragBreakS1 >= noise; noiseIndB = fragBreakB >= noise; noiseIndS2 = fragBreakS2 >= noise
        self.break_length_s1 = fragBreakS1[noiseIndS1]; self.break_location_s1 = pairS1[cutIndS1][noiseIndS1]
        self.break_length_b = fragBreakB[noiseIndB]; self.break_location_b = pairB[cutIndB][noiseIndB]
        self.break_length_s2 = fragBreakS2[noiseIndS2]; self.break_location_s2 = pairS2[cutIndS2][noiseIndS2]
    def saveCoarseGrainedPDB(self,path=defaultDirectory+'vizualization/pymol/pdb/',topo='linear'):
        dna = r2pdb.mkpdb(self.r, topology = topo)
        r2pdb.save_pdb('%s/coarse%0.3d.pdb' %(path,self.time),dna)
    def saveFineGrainedPDB(self,path=defaultDirectory+'vizualization/pymol/pdb/',topo='linear'):
        if (self.bps == None):
            self.interpolate()
        dna = r2pdb.mkpdb(np.asarray(self.bps).reshape([3*self.n_bps,3]), topology = topo)
        r2pdb.save_pdb('%s/fine%0.3d.pdb' %(path,self.time),dna)
    def saveRICCbreak(self,path=defaultDirectory+'analysis/data'):
        tempDict = {'strand1': [self.break_length_s1, self.break_location_s1],
                    'base': [self.break_length_b, self.break_location_b],
                    'strand2': [self.break_length_s2, self.break_location_s2]}
        f = open('%s/ricc%0.3d.pkl' %(path,self.time),"wb")
        pickle.dump(tempDict,f)
        f.close()

