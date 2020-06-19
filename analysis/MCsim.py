#/usr/bin/python 
# npagane | python object to store MC data with tools to perform simple analyses

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
from scipy.ndimage import gaussian_filter
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
                    time_min,time_max,trajectories,channel)
            print('read in ' + str(trial))
        self.linearizeSnapshots()
        # remove intermediary data structures if they do not exist
        self.unnest()
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
        self.linearized_snapshots = []
        snapshots=self.returnSnapshots()
        for i in snapshots.keys():
            for j in snapshots[i]:
                for k in snapshots[i][j]:
                    self.linearized_snapshots.append(self.trials[i].trajectories[j].snapshots[k])
    def unnest(self):
        trials=self.returnTrials()
        trajectories=self.returnTrajectories()
        snapshots=self.returnSnapshots()
        if '' in trials and '' not in trajectories['']:
            self.trajectories = {} 
            for i in trajectories['']:
                self.trajectories[i] = self.trials[''].trajectories[i]
        elif '' in trials and '' in trajectories['']:
            self.snapshots = {} 
            for i in snapshots['']['']:
                self.snapshots[i] = self.trials[''].trajectories[''].snapshots[i]
        elif '' not in trials:
            for i in trials:
                if '' in trajectories[i]:
                    self.trials[i] = self.trials[i].trajectories['']
    def getCenterBeads(self):
        #mulitprocessing does not work with class methods GRRR >:|
        for i in self.linearized_snapshots:
            i.centerBeads()
    def getPairwiseNucleosomeDistance(self):
        for i in self.linearized_snapshots:
            i.pairwiseNucleosomeDistance()
    def getReducedPairwiseNucleosomeDistance(self):
        for i in self.linearized_snapshots:
            i.reducedPairwiseNucleosomeDistance()
    def getInterpolate(self):
        for i in self.linearized_snapshots:
            i.interpolate()
    def getRICCbreak(self):
        for i in self.linearized_snapshots:
            i.RICCbreak()

class Trial:
    def __init__(self,path_to_data,time_min,time_max,trajectories,channel):
        # path to data directory
        self.path_to_data = path_to_data
        if trajectories != ['']:
            self.channels = trajectories
        else:
            self.channels = [channel]*len(trajectories)
        if type(time_max) != int:
            self.time_maxs = time_max            
        else:
            self.time_maxs = [time_max]*len(trajectories)
        if type(time_min) != int:
            self.time_mins = time_min
        else:
            self.time_mins = [time_min]*len(trajectories)
        # store Trajectories in dictionary
        self.trajectories = {}
        for i,trajectory in enumerate(trajectories):
            self.trajectories[trajectory] = Trajectory(self.path_to_data,
                    self.time_mins[i],self.time_maxs[i],self.channels[i])
        # check for PT
        if os.path.exists(self.path_to_data+'nodeNumber'):
            for i in range(1,max(self.channels)+1):
                self.trajectories['PT'+str(i)] = Trajectory(self.path_to_data,
                        max(self.time_mins),min(self.time_maxs),1,temperature=i)
            print('PT configured')

class Trajectory:
    # constructor 
    def __init__(self,path_to_data,time_min,time_max,channel,temperature=None):
        # path to data directory
        self.path_to_data = path_to_data
        # set time range
        self.time_min = time_min
        self.equilibrium_time = self.time_min
        self.time_max = time_max+1 # add plus one here for iteration
        # store Snapshots in dictionary
        self.snapshots = {}
        # snapshots stats
        self.end_to_end = []
        self.reduced_pair_nucs = []
        self.temperature = temperature
        if (self.temperature != None):
            nodes = np.loadtxt(self.path_to_data+'nodeNumber')
            nodes = np.vstack((np.linspace(0, np.shape(nodes)[1]-1, np.shape(nodes)[1]), nodes))
            self.channel = np.asarray(nodes[:,temperature], 'int')
        else:
            self.channel = [channel]*(self.time_max-self.time_min)
        # load snapshots in dictionary
        for i,time in enumerate(range(self.time_min,self.time_max)):
            self.snapshots[time] = Snapshot(self.path_to_data,time,self.channel[i])
    def setEquilibriumTime(self,time):
        self.equilibrium_time = time
    def getEndToEnd(self):
        for time in range(self.equilibrium_time,self.time_max):
            self.end_to_end.append(self.snapshots[time].end_to_end)
    def getEnergies(self):
        self.energies = self.snapshots[self.equilibrium_time].energies
        for time in range(self.equilibrium_time+1,self.time_max):
            self.energies = self.energies.append(self.snapshots[time].energies)
    def getReducedPairwiseNucleosomeDistance(self):
        for time in range(self.equilibrium_time,self.time_max):
            self.snapshots[time].reducedPairwiseNucleosomeDistance()
            self.reduced_pair_nucs.append(self.snapshots[time].reduced_pair_nucs)
        nnuc = self.snapshots[time].n_nucs # assume all snapshots have the same number of nucleosomes
        self.reduced_pair_nucs = np.asarray(self.reduced_pair_nucs).reshape([self.time_max-self.equilibrium_time,nnuc-1])
    def playFineMovie(self,path=defaultDirectory+'vizualization/pymol/pdb/',topo='linear',pymol='~/Applications/pymol/pymol'):
        for time in range(self.time_min,self.time_max):
            self.snapshots[time].saveFineGrainedPDB(path=path,topo=topo)
        os.system(pymol + " -r "+defaultDirectory+'vizualization/pymol/'+"movieFine.py -- " 
                    + str(self.time_max-self.time_min) + " " + path)
    def playCoarseMovie(self,path=defaultDirectory+'vizualization/pymol/pdb/',topo='linear',pymol='~/Applications/pymol/pymol'):
        for i,time in enumerate(range(self.time_min,self.time_max)):
            self.snapshots[time].saveCoarseGrainedPDB(path=path,topo=topo)
        if (self.temperature != None):
            os.system(pymol + " -r "+defaultDirectory+'vizualization/pymol/'+"movieCoarse.py -- " 
                        + str(self.time_max-self.time_min) + " PT " + str(self.temperature) + " " + path + " " + self.path_to_data)
        else:
            os.system(pymol + " -r "+defaultDirectory+'vizualization/pymol/'+"movieCoarse.py -- " 
                        + str(self.time_max-self.time_min) + " " + str(self.channel[i]) + " 1 "  
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
        disc = np.loadtxt('%sd%sv%s' %(self.path_to_data,self.time,self.channel), dtype='float')
        self.wrap = np.asarray(disc[0],dtype='int'); self.basepairs = disc[1]
        # assign constants from postion data
        self.n_beads = len(self.r)
        self.end_to_end = np.linalg.norm(self.r[-1,:]-self.r[0,:])
        self.n_bps = int(np.round(np.sum(self.basepairs[self.basepairs!=0])+np.sum(self.wrap[self.wrap>1])))
        self.end_to_end_norm = self.end_to_end/(self.n_bps*lengthPerBP)
        # energies
        with open('%senergiesv%s' %(self.path_to_data,self.channel)) as fp:
            for i, line in enumerate(fp):
                if i == 0:
                    # get column specifiers
                    cols = line.strip().split()
                else:
                    temp = line.strip().split()
                    if int(temp[0]) == self.time:
                        # get time data
                        energies = temp
                        break
        cols = cols[4:]
        try:
            energies = energies[2:]
        except:
            energies = [np.nan]*len(cols)
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
        self.break_length_s1 = None; self.break_location_s1 = None; self.break_distance_s1 = None
        self.break_length_b = None; self.break_location_b = None;  self.break_distance_b = None
        self.break_length_s2 = None; self.break_location_s2 = None;  self.break_distance_s2 = None
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
        try :
            if (self.center_r == None):
                self.centerBeads()    
        except:
            print('warning: this has already been run')
            return
        nucLocs = np.asarray(np.linspace(0,self.n_beads-1,self.n_beads)[self.wrap>1],dtype='int')
        self.pair_nucs = scipy.spatial.distance.pdist(self.center_r[nucLocs,:])
    # determine the reduced pairwise distances between nucleosomes
    # such that n+x where x -> 1..self.n_nucs-1
    def reducedPairwiseNucleosomeDistance(self):
        try:
            if (self.pair_nucs == None):
                self.pairwiseNucleosomeDistance()
        except:
            print('warning: this has already been run')
            return
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
        self.bps = np.zeros(self.n_bps*3*3).reshape([self.n_bps,3,3])
        indR = 0
        connect = []; chain = []
        chainNum = 1
        for i in range(self.n_beads):
            if self.basepairs[i] != 0:
                maxBp = self.basepairs[i]
                omega = defaultOmega + (getUVangle(self.r[i,:], self.v[i+1,:], self.v[i,:], self.v[i+1,:])%(2*np.pi) - (maxBp*defaultOmega)%(2*np.pi))
                v = omega/defaultOmega*lengthPerBP
                Uout, Vout, Rout = rotateBead(self.u[i,:], self.v[i,:], self.r[i,:], self.basepairs[i], self.wrap[i])
                matIn = np.matrix([self.v[i,:], np.cross(self.u[i,:],self.v[i,:]), self.u[i,:]]).T
                mat = np.matrix([Vout, np.cross(Uout,Vout), Uout]).T
                if (self.wrap[i] > 1): # nucleosome
                    for j in range(len(nucleosomeTran)):
                        row = np.zeros(3*3).reshape([3,3])
                        strand1, base, strand2 = DNAhelix(j,v=0)
                        Rin = np.asarray(nucleosomeTran[len(nucleosomeTran)-1-j,:])
                        # strand 1 backbone
                        row[0,:] = self.r[i,:] + np.matmul(matIn,Rin+strand1[0])
                        # strand 2 backbone
                        row[2,:] = self.r[i,:] + np.matmul(matIn,Rin+strand2[0])
                        # base
                        row[1,:] = self.r[i,:] + np.matmul(matIn,Rin+base[0])
                        # save atoms
                        self.bps[indR,:,:] = row
                        indR = indR + 1                    
                        #connect.append((indR,indR+1))
                        chain.extend([str(chainNum)]*3)
                        if (indR == self.n_bps): break
                    # add the extruding linker from the nucleosome
                    for j in range(int(np.round(maxBp))):
                        row = np.zeros(3*3).reshape([3,3])
                        strand1, base, strand2 = DNAhelix(j)#,omega=0,v=v)
                        # strand 1 backbone
                        row[0,:] = Rout + np.matmul(mat, strand1[0])
                        # strand 2 backbone
                        row[2,:] = Rout + np.matmul(mat, strand2[0])
                        # base
                        row[1,:] = Rout + np.matmul(mat, base[0])
                        # save atoms
                        self.bps[indR,:,:] = row
                        indR = indR + 1
                        chain.extend([str(chainNum)]*3)
                        if (indR == self.n_bps): break
                else: # dna bead
                    for j in range(int(np.round(maxBp))):
                        row = np.zeros(3*3).reshape([3,3])
                        strand1, base, strand2 = DNAhelix(j)#,omega=omega,v=v)
                        # strand 1 backbone
                        row[0,:] = Rout + np.matmul(mat, strand1[0])
                        # strand 2 backbone
                        row[2,:] = Rout + np.matmul(mat, strand2[0])
                        # base
                        row[1,:] = Rout + np.matmul(mat, base[0])
                        # save atoms
                        self.bps[indR,:,:] = row
                        indR = indR + 1
                        chain.extend([str(chainNum)]*3)
                        if (indR == self.n_bps): break
            else:
                chainNum +=1
        return chain
    # distance constraint ricc-seq
    def RICCbreak(self,cutoff=3.5,noise=50.0): # units in nm or bp
        try :
            if (self.bps == None):
                self.interpolate()
        except:
            print('warning: this has already been run')
            return
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
        self.break_length_s1 = fragBreakS1[noiseIndS1]; self.break_location_s1 = indPair[indBreakS1][noiseIndS1]; self.break_distance_s1 = pairS1[cutIndS1][noiseIndS1]
        self.break_length_b = fragBreakB[noiseIndB]; self.break_location_b = indPair[indBreakB][noiseIndB]; self.break_distance_b = pairB[cutIndB][noiseIndB]
        self.break_length_s2 = fragBreakS2[noiseIndS2]; self.break_location_s2 = indPair[indBreakS2][noiseIndS2]; self.break_distance_s2 = pairS2[cutIndS2][noiseIndS2]
    def saveCoarseGrainedPDB(self,path=defaultDirectory+'vizualization/pymol/pdb/',topo='linear'):
        chain = []; connect =  []
        chainNum = 1
        for i in range(self.n_beads):
            chain.append(chainNum) 
            if self.basepairs[i]==0:
                chainNum += 1
            else:
                connect.append((i,i+1))
        dna = r2pdb.mkpdb(self.r,topology=topo,chain=chain,connect=connect)
        r2pdb.save_pdb('%scoarse%0.3d.pdb' %(path,self.time),dna)
    def saveFineGrainedPDB(self,path=defaultDirectory+'vizualization/pymol/pdb/',topo='linear'):
        chain = self.interpolate()
        dna = r2pdb.mkpdb(np.asarray(self.bps).reshape([3*self.n_bps,3]),topology=topo,chain=chain)#,connect=connect)
        r2pdb.save_pdb('%sfine%0.3d.pdb' %(path,self.time),dna)
    def saveRICCbreak(self,path=defaultDirectory+'analysis/data'):
        self.RICCbreak()
        tempDict = {'strand1': [self.break_length_s1, self.break_location_s1],
                    'base': [self.break_length_b, self.break_location_b],
                    'strand2': [self.break_length_s2, self.break_location_s2]}
        f = open('%s/riccBreak%0.3d.pkl' %(path,self.time),"wb")
        pickle.dump(tempDict,f)
        f.close()
    def RICCmatrix(self,blur=20,removeLow=True,init=0.5):
        self.RICCbreak()
        mat = init*np.ones(self.n_bps**2).reshape([self.n_bps,self.n_bps])
        for iterN in range(len(self.break_location_s1)):
            j = int(self.break_location_s1[iterN,0])
            k = int(self.break_location_s1[iterN,1])
            mat[j,k] += 1
        for iterN in range(len(self.break_location_s2)):
            j = int(self.break_location_s2[iterN,0])
            k = int(self.break_location_s2[iterN,1])
            mat[j,k] += 1
        mat = gaussian_filter(mat, sigma=blur)
        if (removeLow):
            mat[mat < 1] = np.nan
        return mat
    def saveRICCmat(self,path=defaultDirectory+'analysis/data'):
        mat = self.RICCmatrix(blur=0,removeLow=False,init=0)
        with open('%s/riccMat%0.3d.txt' %(path,self.time), 'w') as f:
            for i in range(np.shape(mat)[0]):
                for j in range(np.shape(mat)[1]):
                    f.write(str(i) + ' ' + str(j) + ' ' + str(mat[i,j])+ '\n')
