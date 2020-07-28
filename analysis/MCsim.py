#/usr/bin/python 
# npagane | simulation object class and subclasses to parse wlcsim output
"""
Simulation objects to help read in wlcsim data and perform simple analyses/visualizations.
This was specifically designed for Risca lab usage, but feel free to use/edit if it helps you
parse your wlcsim output, too. 

Generally, the data structure is assumed to be organized as follows:

MAIN_SIMULATION_FOLDER
    TRIALS_<I>
        wlcsim
            data
                r<K>v<J>

where there are <I> number of trials and <J> number of trajectories for <K> number of snapshots. Differing 
trials represent different experiements (i.e. changing condtitions/parameters) while different trajectories 
are essentially technical replicates--UNLESS you are parallel tempering, which then means that the different 
trajectories are different "temperatures". The snapshots are wlcsim outputs throughout the simulation "time". 

Therefore the hierarchy of these classes are nested as follows:
Simulation
    Trajectory*
        Trial*
            Snapshot*
* can have multiple instances

---------------
EXAMPLE USAGE 1 
---------------
# we are in a Python instance from within the wlcsim directory 
# we ran a simulation on 1 core (i.e. no parallelizaiton/parallel tempering) and generated 10 snapshots
dat = Simulation(time_max=10) # read in data

# since the Trial and Trajectory objects are null, we reference them with the empty string "''"
dat.trials[''].trajectories[''].playCoarseMovie(pymol=PATH_TO_PYMOL) # visualize trajectory in PyMol

# we want to look at the end-to-end distance of the snapshot at "time point" 1
dat.trials[''].trajectories[''].snapshots[1].end_to_end # get end-to-end distance

---------------
EXAMPLE USAGE 2 
---------------
# we are in a Python instance from NOT within the wlcsim directory 
# we ran two simulations labeled "trial1" and "trial2" that live under the PATH_TO_DATA directory.
# each simulation was run on 2 cores (for parallelizaiton) and generated 10 snapshots
dat = Simulation(path_to_data=PATH_DO_DATA, trials=["trial1", "trial2"], trajectories=[1, 2], time_max=10) # read in data

# we want to look at the movie of the "trial1" simulation's first trajectory
dat.trials['trial1'].trajectories[1].playCoarseMovie(pymol=PATH_TO_PYMOL) # visualize trajectory in PyMol

# we want to look at the end-to-end distance of the snapshot at "time point" 9 of the second trajetory in "trial2"
dat.trials['trial2'].trajectories[2].snapshots[9].end_to_end # get end-to-end distance
"""
import sys
from .utility import *
from .r2pdb import *
import numpy as np
import pandas as pd
import os
import scipy
import scipy.special 
import scipy.spatial
from scipy.ndimage import gaussian_filter
import pickle

class Simulation:
    """
    `Simulation` object to store all your data. You should only directly interface with THIS object when
    loading data. It assumes you are in the top level directory of the wlcsim repo; otherwise, you can 
    be anywhere and specify location of data (i.e path_to/wlcsim/data/) in `path_to_data`. 
    
    You must specify the number of "time points" to read in as `time_max`, and if you are running this
    from within the wlcsim directory, you do not need to specify anything else. 
    If you are elsewhere, you need to not only specify `path_to_data` but also `trials`if you ran simulations 
    across different parameter values or conditions.
    If you ran several instances of the same simulation across several cores (i.e. either for PT/replica exchange
    or for general parallelization), then you must specify the `trajectories`.
    """
    def __init__(self, path_to_data = "%s/data/" %(default_dir), trials = [''], trajectories = [''],
                 time_min = 0, time_max = 110, channel = 0):
        """
        `Simulation` constructor to read in simulation data. 

        Upon reading in the data, this class will try to "unnest" the nested classes as much as possible.
        For example, if you have several trajectories but only one trial, then you can reference the `Trajectory`
        object directly from the `Simulation` class without having to pass through a `Trial` class. 
        
        Parameters
        ----------
        path_to_data : string
            path to either the wlcsim/data directory or the top level directory above the nested `trials`
        trials : list of strings
            list of the subdirectories in `path_to_data` that specify the different trials of the simulations 
        trajectories : list 
            list of the "channel"/thread values of the simulation replicas, i.e. [1,2,3,4]
        time_min : int, optional
            default : 0
            minimum "time point" that you want to start reading in the data from
        time_max : int, optional
            default : 110
            maximum "time point" that you want to end reading in the data from. 
        channel : int
            "channel"/thread value of a specific replica you want to read in  
        """
        # path to data directory
        self.path_to_data = path_to_data
        # store Trials in dictionary
        self.trials = {}
        for trial in trials:
            self.trials[trial] = Trial("%s%s" %(self.path_to_data, trial),
                    time_min, time_max, trajectories, channel)
            print('read in %s' %(str(trial)))
        self.linearizeSnapshots()
        # remove intermediary data structures if they do not exist
        self.unnest()

    def returnTrials(self):
        """Return a list of the `Trial` subdirectory names."""
        return self.trials.keys()

    def returnTrajectories(self):
        """Return a dictionary with all of the `Trajectory` values for all the `Trial` objects."""
        tempDict = {}
        for i in self.trials.keys():
            tempDict[i] = self.trials[i].trajectories.keys()
        return tempDict

    def returnSnapshots(self):
        """Return a dictionary with all of the `Snapshot` values for all `Trajectory` classes of all 
        `Trial` classes."""
        tempDict = {}
        for i in self.trials.keys():
            tempDict[i] = {}
            for j in self.trials[i].trajectories.keys():
                tempDict[i][j] = self.trials[i].trajectories[j].snapshots.keys()
        return tempDict

    def linearizeSnapshots(self):
        """Append all `Snapshot` objects from all `Trajectory` and `Trial` classes into the `Simulation`
        object to gather metrics all at once."""
        self.linearized_snapshots = []
        snapshots=self.returnSnapshots()
        for i in snapshots.keys():
            for j in snapshots[i]:
                for k in snapshots[i][j]:
                    self.linearized_snapshots.append(self.trials[i].trajectories[j].snapshots[k])

    def unnest(self):
        """Remove unnecessary nesting of objects, i.e. if there is 1 `Trial` but 2 `Trajectory` classes, 
        have the `Trajectory` objects reachable from the main `Simulation` object."""
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
        """Find the center of each bead positon's volume for each `Snapshot`  in `linearized_snapshots`.
        See the centerBeads method in the `Snapshot` class for more informtion on the calculation."""
        for i in self.linearized_snapshots:
            i.centerBeads()

    def getPairwiseNucleosomeDistance(self):
        """Find the pairwise distances of each center of a nucleosome bead for each `Snapshot`  in `linearized_snapshots`.
        See the pairwiseNucleosomeDistance method in the `Snapshot` class for more informtion on the calculation."""
        for i in self.linearized_snapshots:
            i.pairwiseNucleosomeDistance()

    def getReducedPairwiseNucleosomeDistance(self):
        """Find the reduced pairwise distances of each center of a nucleosome bead for each `Snapshot`  in `linearized_snapshots`.
        See the reducedPairwiseNucleosomeDistance method in the `Snapshot` class for more informtion on the calculation."""
        for i in self.linearized_snapshots:
            i.reducedPairwiseNucleosomeDistance()

    def getInterpolate(self):
        """Interpolate the helical nature of DNA (to single basepair resolution) onto the coarse-grained bead positions for
        each `Snapshot`  in `linearized_snapshots`.
        See the interpolate method in the `Snapshot` class for more informtion on the calculation."""
        for i in self.linearized_snapshots:
            i.interpolate()

    def getRICCbreak(self):
        """Generate the RICC-seq break patterns from the interpolated structures for each `Snapshot`  in `linearized_snapshots`.
        See the RICCbreak method in the `Snapshot` class for more informtion on the calculation."""
        for i in self.linearized_snapshots:
            i.RICCbreak()

class Trial:
    """
    `Trial` object to store all the different parameter/condition settings if you want to analyze several different 
    simulations at once. It is nested within the `Simulation` object and further nests the `Trajectory` class.

    The `Trial` class will automatically search for and detect whether if there is parallel tempering or not and
    then instanitate nested `Trajectory` classes for each "temperature". The lowest relative ranked "temperature", 
    i.e. "PT1" corresponds to the main trajectory with the enhanced sampling. 

    There is no direct usability of this object. 
    """
    def __init__(self, path_to_data, time_min, time_max, trajectories, channel):
        """
        `Trial` constructor to read in trial (differing parameters/conditions) data. 

        This constructor is called through the `Simulation` class, so you should not have to directly
        instantiate any `Trial` objects directly.
        
        Parameters
        ----------
        path_to_data : string
            path to either the wlcsim/data directory or the top level directory above the nested `trials`
        trajectories : list 
            list of the "channel"/thread values of the simulation replicas, i.e. [1,2,3,4]
        time_min : int, optional
            default : 0
            minimum "time point" that you want to start reading in the data from
        time_max : int, optional
            default : 110
            maximum "time point" that you want to end reading in the data from. 
        channel : int
            "channel"/thread value of a specific replica throughout the simulation "timecouse"
        """
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
                self.trajectories['PT%s'%(str(i))] = Trajectory(self.path_to_data,
                        max(self.time_mins), min(self.time_maxs), 1, temperature = i)
            print('PT configured')

class Trajectory:
    """
    `Trajectory` object to store all the replicates of a simulation run across several threads with or without
    parallel tempgering. It is nested within the `Trial` object and further nests the `Snapshot` class.

    There are a few functions that are usable from this class, namely setting up PyMol movies or determining the 
    "evolution" of metrics throughout the simulation.
    """
    def __init__(self, path_to_data, time_min, time_max, channel, temperature = None):
        """
        `Trajectory` constructor to read in trajectory (technical replicates/parallel tempering) data. 

        This constructor is called through the `Trial` class, so you should not have to directly
        instantiate any `Trajectory` objects directly.
        
        Parameters
        ----------
        path_to_data : string
            path to either the wlcsim/data directory or the top level directory above the nested `trials`
        time_min : int, optional
            default : 0
            minimum "time point" that you want to start reading in the data from
        time_max : int, optional
            default : 110
            maximum "time point" that you want to end reading in the data from. 
        channel : int
            "channel"/thread value of a specific replica throughout the simulation "timecouse"
        temperature : int, optional
            default : None
            the relative ranked "temperature" value of the trajectory if parallel tempering
        """
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
        """Set an "equilibrium" time after which "time point" you accept all succeeding snapshots to be equilibriated"""
        self.equilibrium_time = time

    def getEndToEnd(self):
        """Find the end-to-end distance of the polymer for each `Snapshot`  in the `Trajectory`.
        References the end_to_end field in the `Snapshot` class for more informtion on the calculation."""
        for time in range(self.equilibrium_time,self.time_max):
            self.end_to_end.append(self.snapshots[time].end_to_end)

    def getEnergies(self):
        """Determine the energetics of the polymer for each `Snapshot`  in the `Trajectory`.
        References the energies field in the `Snapshot` class for more informtion on the calculation."""
        self.energies = self.snapshots[self.equilibrium_time].energies
        for time in range(self.equilibrium_time+1,self.time_max):
            self.energies = self.energies.append(self.snapshots[time].energies)

    def getReducedPairwiseNucleosomeDistance(self):
        """Find the reduced pairwise distances of each center of a nucleosome bead for each `Snapshot`  in the `Trajectory`.
        See the reducedPairwiseNucleosomeDistance method in the `Snapshot` class for more informtion on the calculation."""
        for time in range(self.equilibrium_time,self.time_max):
            self.snapshots[time].reducedPairwiseNucleosomeDistance()
            self.reduced_pair_nucs.append(self.snapshots[time].reduced_pair_nucs)
        nnuc = self.snapshots[time].n_nucs # assume all snapshots have the same number of nucleosomes
        self.reduced_pair_nucs = np.asarray(self.reduced_pair_nucs).reshape([self.time_max-self.equilibrium_time,nnuc-1])

    def playFineMovie(self,path=default_dir+'/analysis/pdb/',topo='linear',pymol='pymol'):
        """Play PyMol movie of the polymer throughout the simulation "timecourse" after interpolating into a more fine-grained
        structure. See the saveFineGrainedPDB and interpolate methods in the `Snapshot` class for more informtion on the calculations.

        Currently, this function CAN NOT determine the simulation "timecourse" trajectory with parallel tempering. 
        
        Parameters
        ----------
        path : string, optional
            default : "wlcsim/analysis/pdb/"
            path to where you want to store the fine-grained PDB files that PyMol will then read in and visualize
        topo : string, optional
            default : 'linear'
            topology of polymer structure (for the risca lab, it will almost always be 'linear')
        pymol : string, optional
            default : 'pymol'
            exectable command to initiaite PyMol, i.e. "~/Applications/pymol/pymol" 
        """
        for time in range(self.time_min,self.time_max):
            self.snapshots[time].saveFineGrainedPDB(path=path,topo=topo)
        os.system(pymol + " -r "+default_dir+"/analysis/movieFine.py -- " 
                    + str(self.time_max-self.time_min) + " " + path)
                    
    def playCoarseMovie(self,path=default_dir+'/analysis/pdb/',topo='linear',pymol='pymol'):
        """Play PyMol movie of the polymer throughout the simulation "timecourse" visualizing the excluded volume of the chain. 
        See the saveCoarseGrainedPDB method in the `Snapshot` class for more informtion on the calculation.

        This function can determine the simulation "timecourse" trajectory with parallel tempering. 
        
        Parameters
        ----------
        path : string, optional
            default : "wlcsim/analysis/pdb/"
            path to where you want to store the coarse-grained PDB files that PyMol will then read in and visualize
        topo : string, optional
            default : 'linear'
            topology of polymer structure (for the risca lab, it will almost always be 'linear')
        pymol : string, optional
            default : 'pymol'
            exectable command to initiaite PyMol, i.e. "~/Applications/pymol/pymol" 
        """
        for time in range(self.time_min,self.time_max):
            self.snapshots[time].saveCoarseGrainedPDB(path=path,topo=topo)
        if (self.temperature != None):
            os.system(pymol + " -r "+default_dir+"/analysis/movieCoarse.py -- " 
                        + str(self.time_max-self.time_min) + " PT " + str(self.temperature) + " " + path + " " + self.path_to_data)
        else:
            os.system(pymol + " -r "+default_dir+"/analysis/movieCoarse.py -- " 
                        + str(self.time_max-self.time_min) + " " + str(self.channel[-1]) + " 1 "  
                        + path + " " + self.path_to_data)

class Snapshot:
    """
    `Snapshot` object to store all the positions and orientations of the computational beads for a given snapshot. 
    This object also calculates different polymer metrics. It is nested within the `Trajectory` object.

    This class contains most of the fields and functions that are useful for wlcsim analysis. Some useful class 
    fields are:
        `r` : computational bead positions
        `u` : computational bead U orientation vectors
        `v` : computational bead V orientation vectors
        `basepairs` : discreitzation (in bp) of each computational bead
        `wrap` : how much DNA (in bp) is wrapped around a computational bead, i.e. nucleosome=147 and DNA=1
        `n_beads` : number of coarse-grained computational beads in polymer
        `end_to_end` : end-to-end distance of polymer
        `n_bps` : number of DNA basepairs throughout polymer (including DNA wrapped around nucleosomes)
        `end_to_end_norm` : end-to-end distance normalized by its contour length `n_bps`
        `energies` : wlcsim defined energetics of the polymer chain (in kT)
        `center_r` : center positions of computational beads with respect to their excluded volume 
        `n_nucs` : number of nucleosomes on the chain
        `pair_nucs` : nucleosome pairwise distances, i.e. all combinatorics of n choose k
        `reduced_pair_nucs` : reduced nucleosome pairwise distances, i.e. n & n+1, n & n+2, etc. 
        `bps` : location of each interpolated phosphate-sugar-phosphate basepair throughout the chain
        `break_length_s1`, `break_length_b`, `break_length_2` : RICC-seq pairwise fragment lengths using `bps`
        `break_location_s1`, `break_location_b`, `break_location_s2` : RICC-seq pairwise break locations using `bps`
        `break_distance_s1`, `break_distance_b`, `break_distance_s2` : RICC-seq pairwise break distances using `bps`

    """
    def __init__(self,path_to_data,time,channel):
        """
        `Snapshot` constructor to read in snaphot (bead positions and orientations) data. 

        This constructor is called through the `Trajectory` class, so you should not have to directly
        instantiate any `Snapshot` objects directly.
        
        Parameters
        ----------
        path_to_data : string
            path to either the wlcsim/data directory or the top level directory above the nested `trials`
        time : int, optional
            "time point" of the current wlcsim output structure
        channel : int
            "channel"/thread value of the current wlcsim output structure
        """
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
        try: 
            disc = np.loadtxt('%sd%sv%s' %(self.path_to_data,self.time,self.channel), dtype='float')
            self.wrap = np.asarray(disc[0],dtype='int'); self.basepairs = disc[1]
        except:
            self.basepairs = np.array([10.5]*len(self.r))
            self.wrap = np.array([1]*len(self.r))
        # assign constants from postion data
        self.n_beads = len(self.r)
        self.end_to_end = np.linalg.norm(self.r[-1,:]-self.r[0,:])
        self.n_bps = int(np.round(np.sum(self.basepairs[self.basepairs!=0])+np.sum(self.wrap[self.wrap>1])))
        self.end_to_end_norm = self.end_to_end/(self.n_bps*length_per_bp)
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
    def centerBeads(self,nside=16,type='regular'):
        """Get the location of the center of each bead with respect to its excluded volume.
        
        Parameters
        ----------
        nside : int, optional
            default : 16
            number of sides of the modeled excluded volume
        type : string, optional
            default : 'regular'
            whether the excluded geometry is a regular shape or not (currently can only handle regular shapes)

        Generates
        ---------
        center_r : center positions of computational beads with respect to their excluded volume 
        """
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
        """Get the pairwise distance between the center of each nucleosome on the chain, i.e. n choose k
        
        Generates
        ---------
        pair_nucs : nucleosome pairwise distances, i.e. all combinatorics of n choose k
        """
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
        """Get the pairwise distance between the center of each nucleosome on the chain, i.e. rather than 
        n choose k, average over n & n+1, n & n+2, etc.
        
        Generates
        ---------
        reduced_pair_nucs : reduced nucleosome pairwise distances, i.e. n & n+1, n & n+2, etc. 
        """
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
        """Interpolate the phosphate-sugar-phosphate basepair locations between the coarse-grained computational beads. 
        See the rotate_bead() and DNAhelix() functions in the utility module for more information on the calculations.
        
        Generates
        ---------
        bps : location of each phosphate-sugar-phosphate basepair throughout the chain
        """
        # rotate DNA strand into material frame
        self.bps = np.zeros(self.n_bps*3*3).reshape([self.n_bps,3,3])
        indR = 0
        connect = []; chain = []
        chainNum = 1
        for i in range(self.n_beads):
            if self.basepairs[i] != 0:
                maxBp = self.basepairs[i]
                omega = default_omega + (get_uv_angle(self.v[i,:], self.v[i+1,:])%(2*np.pi) - (maxBp*default_omega)%(2*np.pi))
                v = omega/default_omega*length_per_bp
                Uout, Vout, Rout = rotate_bead(self.u[i,:], self.v[i,:], self.r[i,:], self.basepairs[i], self.wrap[i])
                matIn = np.matrix([self.v[i,:], np.cross(self.u[i,:],self.v[i,:]), self.u[i,:]]).T
                mat = np.matrix([Vout, np.cross(Uout,Vout), Uout]).T
                if (self.wrap[i] > 1): # nucleosome
                    for j in range(len(nucleosome_tran)):
                        row = np.zeros(3*3).reshape([3,3])
                        strand1, base, strand2 = DNAhelix(j,v=0)
                        Rin = np.asarray(nucleosome_tran[len(nucleosome_tran)-1-j,:])
                        # strand 1 backbone
                        row[0,:] = self.r[i,:] + np.matmul(matIn,Rin+strand1)
                        # strand 2 backbone
                        row[2,:] = self.r[i,:] + np.matmul(matIn,Rin+strand2)
                        # base
                        row[1,:] = self.r[i,:] + np.matmul(matIn,Rin+base)
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
                        row[0,:] = Rout + np.matmul(mat, strand1)
                        # strand 2 backbone
                        row[2,:] = Rout + np.matmul(mat, strand2)
                        # base
                        row[1,:] = Rout + np.matmul(mat, base)
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
                        row[0,:] = Rout + np.matmul(mat, strand1)
                        # strand 2 backbone
                        row[2,:] = Rout + np.matmul(mat, strand2)
                        # base
                        row[1,:] = Rout + np.matmul(mat, base)
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
        """Determine the RICC-seq break patterns from the interpolated fine-graied polymer structure.
        
        Parameters
        ----------
        cutoff : float, optional
            default : 3.5 nm
            RICC-seq radiolytic cleavage radius under which there will be simulated breaks
        noise : float, optional
            default : 50.0 bp
            fragment length under which any breaks are considered noise and thus thrown out

        Generates
        ---------
        break_length_s1, break_length_b, break_length_2 : RICC-seq pairwise fragment lengths using `bps`
        break_location_s1, break_location_b, break_location_s2 : RICC-seq pairwise break locations using `bps`
        break_distance_s1, break_distance_b, break_distance_s2 : RICC-seq pairwise break distances using `bps`
        """
        try :
            if (self.bps == None):
                self.interpolate()
        except:
            pass
        try :
            if (self.break_length_s1 == None):
                pass
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
    
    def saveCoarseGrainedPDB(self,path=default_dir+'/analysis/pdb/',topo='linear'):
        """Save the coarse-grained wlcsim output in a PDB format.
        
        Parameters
        ----------
        path : string, optional
            default : "wlcsim/analysis/pdb/"
            path to where you want to store the fine-grained PDB files 
        topo : string, optional
            default : 'linear'
            topology of polymer structure (for the risca lab, it will almost always be 'linear')
        """
        chain = []; connect =  []
        chainNum = 1
        for i in range(self.n_beads):
            chain.append(chainNum) 
            if self.basepairs[i]==0:
                chainNum += 1
            else:
                connect.append((i,i+1))
        dna = mkpdb(self.r,topology=topo,chain=chain,connect=connect)
        save_pdb('%scoarse%0.3d.pdb' %(path,self.time),dna)
    
    def saveFineGrainedPDB(self,path=default_dir+'/analysis/pdb/',topo='linear'):
        """Save the interpolated fine-grained wlcsim structure in a PDB format.
        
        Parameters
        ----------
        path : string, optional
            default : "wlcsim/analysis/pdb/"
            path to where you want to store the coarse-grained PDB files 
        topo : string, optional
            default : 'linear'
            topology of polymer structure (for the risca lab, it will almost always be 'linear')
        """
        chain = self.interpolate()
        dna = mkpdb(np.asarray(self.bps).reshape([3*self.n_bps,3]),topology=topo,chain=chain)#,connect=connect)
        save_pdb('%sfine%0.3d.pdb' %(path,self.time),dna)
    
    def saveRICCbreak(self,path=default_dir+'/analysis/data'):
        """Save the RICC-seq break patterns (strand1/strand2/base lengths and locatons) to pickle files.
        
        Parameters
        ----------
        path : string, optional
            default : "wlcsim/analysis/data/"
            path to where you want to store the RICC-seq break pattern pickle files
        """
        self.RICCbreak()
        tempDict = {'strand1': [self.break_length_s1, self.break_location_s1],
                    'base': [self.break_length_b, self.break_location_b],
                    'strand2': [self.break_length_s2, self.break_location_s2]}
        f = open('%s/riccBreak%0.3d.pkl' %(path,self.time),"wb")
        pickle.dump(tempDict,f)
        f.close()
    
    def RICCmatrix(self,blur=20,removeLow=True,init=0.5):
        """Generate the fragment-break frequency matrix from the interpolated chain.
        
        Parameters
        ----------
        blur : int, optional
            default : 20
            Gaussian blur to apply to the frequency matrix to add "noise" and make the matrix less sparse
        removeLow : boolean, optional
            default : True
            whether to remove frequency values less than 1 and replace them as NAN
        init : float, optional
            default : 0.5
            value to initialize the entries of the frequency matrix with if trying to avoid zeros in cells

        Returns
        -------
        mat : matrix
            RICC-seq fragment break frequency matrix 
        """
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
    
    def saveRICCmat(self,path=default_dir+'/analysis/data'):
        """Save the RICC-seq fragment break frequency matrix to a text file.
        
        Parameters
        ----------
        path : string, optional
            default : "wlcsim/analysis/data/"
            path to where you want to store the RICC-seq fragment break frequency matrix file
        """
        mat = self.RICCmatrix(blur=0,removeLow=False,init=0)
        with open('%s/riccMat%0.3d.txt' %(path,self.time), 'w') as f:
            for i in range(np.shape(mat)[0]):
                for j in range(np.shape(mat)[1]):
                    f.write(str(i) + ' ' + str(j) + ' ' + str(mat[i,j])+ '\n')
