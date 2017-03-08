""" This module "understands" wlcsim.exe's output, and will make it into
nice, pretty numpy/pandas arrays for you."""
import numpy as np
import pandas as pd
import os
import logging
import itertools
import glob
import re
from .input import ParsedInput
from .utils.cached_property import cached_property
import shutil

from numba import double
from numba.decorators import jit, autojit

logger = logging.getLogger(__name__)

# keep these arround for cache reasons
# initialize with size 1
i_arr = np.array([[1]])
j_arr = np.array([[1]])
ui = (np.array([0]), np.array([0]))

@autojit
def pairwise_distances3(r):
    """
    Returns np.array of pairwise distances of points in R^3 stored in a
    3-by-num_points numpy array.
    """
    n = r.shape[1]
    dr = np.empty((n, n), dtype=np.float)
    for i in range(n):
        for j in range(n):
            d = 0.0
            for k in range(3):
                tmp = r[k, i] - r[k, j]
                d += tmp*tmp
            dr[i, j] = np.sqrt(d)
    return dr


def df_from_coltimes_method3(coltimes):
    # for caching
    global i_arr
    global j_arr
    global ui
    # method 3: manipulate numpy arrays directly. about same speed as
    # method #2, but insanely more complex to reason about
    N = coltimes.shape[0]
    if not np.all(j_arr.shape == coltimes.shape):
        j_arr = np.linspace(0, N-1, N)
        j_arr = np.tile(j_arr, (N, 1))
        i_arr = j_arr.T
        j_arr = np.reshape(j_arr, (np.product(j_arr.shape), 1))
        i_arr = np.reshape(i_arr, (np.product(i_arr.shape), 1))
        ui = np.triu_indices(N, k=0)
    coltimes[ui] = np.nan
    coltimes_lin = np.reshape(coltimes, (np.product(coltimes.shape), 1))
    to_keep_idx = np.logical_not(np.isnan(coltimes_lin))
    coltimes_lin = coltimes_lin[to_keep_idx]
    i_lin = i_arr[to_keep_idx]
    j_lin = j_arr[to_keep_idx]
    coltimes = np.vstack((coltimes_lin, i_lin, j_lin)).T
    if coltimes.ndim == 1:
        coltimes = coltimes.reshape(1, coltimes.shape[0])
    return pd.DataFrame(coltimes, columns=('coltime', 'i', 'j'))

def df_from_coltimes_method2(coltimes):
    ## # turns out this is slow compared to pandas or textreader
    ##coltimes = np.loadtxt(coltimes_file)
    # method 2: still slow from creating Series?
    coltimes_serieses = []
    # make a Series out of the lower triangular part of coltimes matrix
    for i in range(coltimes.shape[0]):
        for j in range(i-1):
            # method 2: still slow from creating Series?
            coltimes_serieses.append(pd.Series({
                    'i': i, 'j': j,
                    'coltime': coltimes[i][j]}))
            # # method 1: this causes N reallocations, super slow
            # self.coltimes = self.coltimes.append(
            #     pd.Series({'i': i, 'j': j,
            #                'coltime': coltimes[i][j]}),
            #     ignore_index=True
            # )
    # method 2: still slow from creating Series?
    return pd.DataFrame(coltimes_serieses)

def df_from_coltimes_method1(coltimes):
    df = pd.DataFrame()
    # make a Series out of the lower triangular part of coltimes matrix
    for i in range(coltimes.shape[0]):
        for j in range(i-1):
            # method 1: this causes N reallocations, super slow
            self.coltimes = df.append(
                pd.Series({'i': i, 'j': j,
                           'coltime': coltimes[i][j]}),
                ignore_index=True
            )
    # method 2: still slow from creating Series?
    return df

def non_overlapping_pairs(length, spacing, max_pairs=None, min_pairs=0):
    """Return the "most spaced" possible sets of pairs of indices that are a
    specific spacing apart, for getting statistics that aren't biased by the
    overlap between neighboring segments of the polymer."""
    # use the number of possible left endpoints to get the num possible pairs
    num_possible_pairs = np.floor((length - spacing)/spacing).astype(int)
    if num_possible_pairs < min_pairs:
        raise ValueError("Requested too large of a spacing for a given minimum "
                         + "number of overlapping pairs")
    # by default, return all possible non-overlapping pairs
    if max_pairs is None or max_pairs > num_possible_pairs:
        max_pairs = num_possible_pairs
    num_leftover_indices = length - max_pairs*spacing
    # num skips is of course num_pairs-1, fenceposts problem
    skip_per_step = int(num_leftover_indices/(max_pairs-1))
    num_leftover_indices -= skip_per_step*(max_pairs-1)
    # the leftover indices should pad the outside of the things, so we just use
    # half of them to pad on the left
    num_leftover_indices = int(num_leftover_indices/2)
    left, right = (0, spacing)
    left, right = (left + num_leftover_indices, right + num_leftover_indices)
    yield (left, right)
    for i in range(max_pairs-1):
        # if num_leftover_indices > 0:
        #     num_leftover_indices -= 1
        #     left = right+skip_per_step+1
        #     right = left + spacing
        left, right = (right + skip_per_step, right + skip_per_step + spacing)
        yield (left, right)

class Sim:
    """
    Sim provides an interface to the heterogenous data produced by a
    simulation that saves coltimes, positions, and momentums.

    Unless specified during initialization, the coltimes, r, and u are loaded
    on demand upon first access.

    Accessing properties of this class for the first time can return an IOError
    if a file is moved, or doesn't exist. Afterwards, the values will be cached.
    """

    def __init__(self, sim_path, input_file=None, overwrite_coltimes=False, load_method=3):
        # since these are all lazily evaluated, we throw exceptions right away
        # if the files we need do not exist, to prevent long running programs
        # from throwing an exception midway
        if not os.path.isdir(sim_path):
            raise FileNotFoundError('Directory (' + sim_path +
                                    ') provided when constructing'
                                    ' wlcsim.data.Sim does not exist!')
        self.sim_path = sim_path
        self.sim_dir = os.path.basename(sim_path)
        self.input_file = Sim.guess_input_name(sim_path, input_file)
        # without input file, we don't really have a chance of using the
        # simulation data from a simulation anyway, so we can safely error out
        # here
        if not os.path.isfile(self.input_file):
            raise FileNotFoundError(self.input_file)
        self.data_dir = os.path.join(sim_path, 'data')
        self.coltimes_file = os.path.join(self.data_dir, 'coltimes')
        self.overwrite_coltimes = overwrite_coltimes
        #TODO correctly set rfile base for quinn's simulations
        self._rfile_base = os.path.join(self.data_dir, 'r')
        self._ufile_base = os.path.join(self.data_dir, 'u')
        self.load_method = load_method
        self._has_loaded_r = False

    def guess_input_name(sim_path, input_file):
        bruno_input_name = os.path.join(sim_path, 'input', 'input')
        quinn_input_name = os.path.join(sim_path, 'input', 'params')
        if input_file:
            if os.path.isfile(input_file):
                return input_file
            elif os.path.is_file(os.path.join(sim_path, input_file)):
                return os.path.join(sim_path, input_file)
            elif os.path.is_file(os.path.join(sim_path, 'input', input_file)):
                return os.path.join(sim_path, 'input', input_file)
        if os.path.isfile(bruno_input_name):
            return bruno_input_name
        elif os.path.isfile(quinn_input_name):
            return quinn_input_name
        else:
            return None

    @cached_property
    def parsed_input(self):
        return ParsedInput(self.input_file)

    @cached_property
    def coltimes(self):
        coltimes = pd.read_csv(self.coltimes_file, delim_whitespace=True,
                                header=None).values
        coltimes = self.df_from_coltimes(coltimes, method=self.load_method)
        return coltimes

    @cached_property
    def r0(self):
        """ Numpy array of positions of each particle pre-first time point. """
        rfile_name = self._rfile_base + str(0)
        r0 = pd.read_csv(rfile_name, delim_whitespace=True, header=None)
        return r0.values.reshape((3, self.num_polymers, self.num_beads))

    @cached_property
    def rend(self):
        """ Numpy array of positions of each particle pre-first time point. """
        last_saved_ti = np.max(self.rtime_idxs)
        rfile_name = self._rfile_base + str(last_saved_ti)
        rend = pd.read_csv(rfile_name, delim_whitespace=True, header=None)
        return rend

    @property
    def rend2end(self):
        # there's a fast way to do this if there's only one polymer
        #TODO use known linesize to extend this to work no matter how many
        # polymers there are
        if self.num_polymers == 1:
            rend2end = np.empty((self.num_time_points, 3, self.num_polymers))
            rend2end[:] = np.nan
            for ti in self.rtime_idxs:
                rfile_name = self._rfile_base + str(ti)
                with open(rfile_name, 'rb') as f:
                    first_line = f.readline().rstrip()
                    line_len = len(first_line)
                    r1 = np.array([float(rij) for rij in first_line.split()])
                    f.seek(-line_len*2, os.SEEK_END)
                    last_line = f.readlines()[-1].rstrip()
                    rN = np.array([float(rij) for rij in last_line.split()])
                    rend2end[ti,:,0] = rN - r1
        else:
            rend2end = self.r[:,:,:,-1] - self.r[:,:,:,0]
        return rend2end

    @cached_property
    def u0(self):
        """ Numpy array of positions of each particle at each time point. """
        ufile_name = self._ufile_base + str(0)
        u0 = pd.read_csv(rfile_name, delim_whitespace=True, header=None)
        return u0.values.reshape((3, self.num_polymers, self.num_beads))

    @cached_property
    def r(self):
        """ Numpy array of positions of each particle at each time point.
        shape:
        r.shape == (self.num_time_points, 3, self.num_polymers, self.num_beads)
        """
        #TODO make only load saved time points
        r = np.full((self.num_time_points, 3, self.num_polymers,
                     self.num_beads), np.nan)
        for i in range(self.num_time_points):
            rfile_name = self._rfile_base + str(i)
            try:
                ri = pd.read_csv(rfile_name, delim_whitespace=True, header=None)
            except OSError:
                continue
            r[i,:,:,:] = ri.values.T.reshape((1, 3, self.num_polymers,
                                            self.num_beads))
        self._has_loaded_r = True
        return r

    @cached_property
    def u(self):
        """ Numpy array of positions of each particle at each time point. """
        u = np.full((self.num_time_points, 3, self.num_polymers,
                     self.num_beads), 0.0)
        for i in range(self.num_time_points):
            ufile_name = self._ufile_base + str(i)
            try:
                ui = pd.read_csv(ufile_name, delim_whitespace=True, header=None)
            except OSError:
                continue
            u[i,:,:,:] = ui.values.T.reshape((1, 3, self.num_polymers,
                                            self.num_beads))
        return u

    @property
    def num_time_points(self):
        return int(self.params['NUMSAVEPOINTS']) + 1 # account for time 0

    @property
    def num_polymers(self):
        return int(self.params['NP'])

    @property
    def num_beads(self):
        return int(self.params['NB'])

    @cached_property
    def center_of_mass(self):
        # "equal mass" beads, so just average positions
        return np.mean(self.r, axis=3)

    @cached_property
    def rtime_idxs(self):
        rfiles = glob.glob(self._rfile_base + '*')
        return self._time_idxs(self._rfile_base, rfiles)

    @cached_property
    def utime_idxs(self):
        ufiles = glob.glob(self._ufile_base + '*')
        return self._time_idxs(self._ufile_base, ufiles)

    def _time_idxs(self, base, file_names):
        number_from_filename = re.compile(re.escape(base) + '([0-9]+)')
        return np.array([int(re.search(number_from_filename,
                                       file_name).groups()[0])
                         for file_name in file_names])

    @cached_property
    def param_names(self):
        param_names, _ = self.parsed_input
        return param_names

    @property # caching won't help a call that's just a redirection
    def params(self):
        return self.parsed_input.params

    @cached_property
    def rtimes(self):
        return self.rtime_idxs * self.params['DT']

    @cached_property
    def utimes(self):
        return self.utime_idxs * self.params['DT']

    @cached_property
    def polymer_length(self):
        """The sum of interbead distances in the polymer at each time"""
        sqd = np.sum(np.power(self.r[:,:,:,1:] - self.r[:,:,:,0:-1], 2),
                     axis=(1,3))
        return np.sqrt(sqd)


    def dr_mean_by_distance(self, lengths, max_samples_per_len=1000):
        distances = np.full((len(lengths), self.num_time_points, self.num_polymers), np.nan)
        for i,length in enumerate(lengths):
            mean_dist = 0
            num_dists = 0
            for left, right in non_overlapping_pairs(self.polymer_length,
                    length, max_pairs=max_samples_per_len):
                mean_dist += np.linalg.norm(self.r[:,:,:,right] - self.r[:,:,:,left], axis=1)
                num_dists += 1
            distances[i,:,:] = np.squeeze(mean_dist/num_dists)
        return distances

    # the properties of coltimes are all coded as a function of the linear
    # distance between the beads of interest, as opposed to as a function of
    # their absolute location, although that may not be super good?
    @property
    def fraction_collided(self):
        return 1.0 - self.fraction_uncollided

    @property
    def fraction_uncollided(self):
        return np.array([float(self.num_uncollided[i])/float(self.num_possible[i])
                         for i in range(length(self.num_possible))])

    @property
    def num_uncollided(self):
        c = self.coltimes
        distances = c['linear_distance'].unique()
        return np.array([sum(pd.isnull(c.loc[c['linear_distance'] == i, 'coltime']))
                         for i in distances])

    @property
    def num_collided(self):
        return self.num_possible - self.num_uncollided

    @property
    def num_possible(self):
        c = self.coltimes
        distances = c['linear_distance'].unique()
        return np.array([len(c.loc[c['linear_distance'] == i, 'coltime'])
                         for i in distances])

    def position_in_axis(self, axis=0):
        """Project position as a function of N along the chain into one of it's
        dimensions to see how it's oriented.

        Useful for e.g. visualizing if bacterial DNA in confinement aligns
        linearly along the cell's length as in e.g. Caulobacter."""
        return np.squeeze(self.r[:,axis,:,:])

    def df_from_coltimes(self, coltimes, method=3):
        # coltimes_file = os.path.join(self.data_dir, 'coltimes.csv')
        # if os.path.isfile(coltimes_file) and not self.overwrite_coltimes:
        #     coltimes = pd.read_csv(coltimes_file, index_col=0)
        #     return coltimes
        # elif method == 1:
        if method == 1:
            coltimes = df_from_coltimes_method1(coltimes)
        elif method == 2:
            coltimes = df_from_coltimes_method2(coltimes)
        elif method == 3:
            coltimes = df_from_coltimes_method3(coltimes)
        else:
            raise KeyError('load_method=' + str(method) + ' does not exist!')
        nan_fill_coltimes(coltimes)
        self.add_useful_cols_to_coltimes(coltimes)
        # coltimes.to_csv(coltimes_file)
        return coltimes

    def add_useful_cols_to_coltimes(self, coltimes):
        """Only works on coltimes that come from combine_coltimes_csvs."""
        dr = self.params['L']/(self.params['NB'] - 1)
        coltimes['linear_distance'] = dr*np.abs(coltimes['i'] - coltimes['j'])


def nan_fill_coltimes(coltimes):
    """Takes a raw coltimes pandas array from e.g. write_coltimes_csvs and puts
    nan's where the -1.0's are, i.e. where the simulation did not go out far
    enough."""
    coltimes.loc[coltimes['coltime'] == -1.0, 'coltime'] = np.nan



class Scan:
    """Can tell you what dirs in a run_dir are simulation directories, and what
    all the input files said, and more. """

    def __init__(self, run_dir, limit=None, complete=False, *args, **kwargs):
        self.limit = limit
        self.complete = complete
        self.path = run_dir
        self._prepend_path = lambda d: os.path.join(self.path, d)
        self._is_simdir = lambda d: Scan.is_simdir(d, self.path)
        self._is_complete = lambda d: Scan.is_complete(d, self.path)
        # self.wlcsim_data = [Sim(folder, *args, **kwargs) for folder in
        #                     glob.glob(os.path.join(run_dir, '*.*'))
        #                     if os.path.isdir(folder)]
        # self.param_names = self.wlcsim_data[0].param_names
        # for sim in self.wlcsim_data:
        #     sim.coltimes['run_name'] = sim.sim_dir
        #     for param in sim.params:
        #         sim.coltimes[param] = sim.params[param]
        # coltimes = [sim.coltimes for sim in self.wlcsim_data]
        # self.coltimes = pd.concat(coltimes)
        # self.calculate_linear_distance()

    @property
    def sim_dirs(self):
        if self.complete:
            is_sim = lambda d: self._is_simdir(d) and self._is_complete(d)
        else:
            is_sim = self._is_simdir
        sim_dirs = filter(is_sim, os.listdir(self.path))
        if self.limit is None:
            return list(sim_dirs)
        else:
            return list(itertools.islice(sim_dirs, self.limit))

    @property
    def sim_paths(self):
        return list(map(self._prepend_path, self.sim_dirs))

    @staticmethod
    def is_simdir(dirname, run_dir):
        sim_name_re = re.compile('.*\.[0-9]+')
        abs_dirname = os.path.join(run_dir, dirname)
        input_file = os.path.join(abs_dirname, 'input', 'input')
        return sim_name_re.search(dirname) \
                and os.path.isdir(abs_dirname) \
                and os.path.isfile(input_file)

    @staticmethod
    def is_complete(dirname, run_dir):
        abs_dirname = os.path.join(run_dir, dirname)
        complete_file = os.path.join(abs_dirname, 'complete')
        return os.path.isfile(complete_file)

    @cached_property
    def sims(self):
        return [Sim(folder) for folder in self.sim_paths]

    @cached_property
    def inputs(self):
        df = pd.DataFrame([sim.params for sim in self.sims])
        df.index = [sim.sim_dir for sim in self.sims]
        return df

    @cached_property
    def all_param_values(self):
        return {col: self.inputs[col].unique() for col in self.inputs.columns}

    @property
    def run_dir(self):
        return self.path

    @property
    def scan_dir(self):
        return self.path

    # from back when we used this class to combine coltimes's from a scan
    # into one big pandas df
    # def calculate_linear_distance(self):
        # # actual distance along polymer is dx*(number_separating)
        # self.coltimes['dx'] = self.coltimes['L']/(self.coltimes['NB'] - 1)
        # self.coltimes['linear_distance'] = self.coltimes['dx']*np.abs(self.coltimes['i'] - self.coltimes['j'])
    # def calculate_for_each_param(self, f):
        # """unimplemented"""
        # pass

    @property
    def rend2end(self):
        """End bead distances for each simulation as numpy array"""
        example_sim = self.sims[0]
        rend2end = np.empty([len(self.sims), example_sim.num_time_points,
                             3, example_sim.num_polymers])
        rend2end[:] = np.nan
        print('This is a time consuming function, I will print out a number'
              'for each sim I process.')
        for i,sim in enumerate(self.sims):
            print(i)
            rend2end[i,:,:,:] = sim.rend2end
        return rend2end


def write_coltimes_csvs(run_dir, overwrite=False, *args, **kwargs):
    for folder in Scan(run_dir).sim_paths:
        coltimes_file = os.path.join(folder, 'coltimes.csv')
        if not overwrite and os.path.isfile(coltimes_file):
            continue
        try:
            sim = Sim(folder, *args, **kwargs)
        except FileNotFoundError:
            logger.warning('write_coltimes_csvs: Unable to construct sim from' + folder)
            continue
        try:
            sim.coltimes.to_csv(coltimes_file)
        except OSError:
            logger.warning('No coltimes file in: ' + folder)
            continue
        except EmptyDataError:
            logger.warning('coltimes file is empty in: ' + folder)
            continue

def write_initial_dists_csvs(run_dir, overwrite=True, *args, **kwargs):
    for folder in Scan(run_dir).sim_paths:
        dists_file = os.path.join(folder, 'initial_dists.csv')
        if not overwrite and os.path.isfile(dists_file):
            continue
        try:
            sim = Sim(folder, *args, **kwargs)
        except FileNotFoundError:
            logger.warning('write_initial_dists_csvs: Unable to construct sim from' + folder)
            continue
        try:
            np.savetxt(dists_file, pairwise_distances3(sim.r0), delimiter=',')
        except OSError:
            logger.warning('r0 not found for: ' + folder)
            continue
        except EmptyDataError:
            logger.warning('r0 file is empty in: ' + folder)
            continue

def write_coltimes_pkls(run_dir, *args, **kwargs):
    for folder in Scan(run_dir).sim_paths:
        try:
            sim = Sim(folder, *args, **kwargs)
        except FileNotFoundError:
            logger.warning('write_coltimes_pkls: Unable to construct sim from' + folder)
        try:
            sim.coltimes.to_pickle(os.path.join(folder, 'coltimes.pkl'))
        except OSError:
            logger.warning('No coltimes file in: ' + folder)
            continue
        except EmptyDataError:
            logger.warning('coltimes file is empty in: ' + folder)
            continue

# def combine_coltimes_pkls(run_dir, *args, **kwargs):
#     pkl_glob = os.path.join(run_dir, '*', 'coltimes.pkl')
#     for coltimes_pkl in glob.glob(pkl_glob):

def combine_coltimes_csvs(run_dir, conditionals=[],
                          outfile_name='coltimes.csv'):
    """takes a directory that scan_wlcsim.py has output to and an optional list
    of functions that takes a coltimes array and decides whether it may be
    used, and an output location"""
    #TODO: finish converting to following:
    # csv_files = [os.path.join(folder, 'coltimes.csv')
    #              for folder in Scan(run_dir).sim_paths]
    # csv_files = [f for f in csv_files if os.path.isfile(f)]
    csv_glob = os.path.join(run_dir, '*', 'coltimes.csv')
    csv_files = iter(glob.glob(csv_glob))
    first_file = next(csv_files)
    coltimes = pd.read_csv(first_file, index_col=0)
    coltimes['run_name'] = os.path.basename(os.path.dirname(first_file))
    # loop and a half. create file uncoditinoally on pre-first iteration so
    # that we can append from then on
    input_file = os.path.join(os.path.dirname(first_file), 'input', 'input')
    parsed_input = ParsedInput(input_file)
    for param in parsed_input.ordered_param_names:
        coltimes[param] = parsed_input.params[param]
    coltimes.to_csv(os.path.join(run_dir, 'coltimes.csv'), index=False)
    # now append the remaining
    for csv_file in csv_files:
        print(csv_file)
        coltimes = pd.read_csv(csv_file, index_col=0)
        coltimes['run_name'] = os.path.basename(os.path.dirname(csv_file))
        # throw away all coltimes that don't pass our conditional checks
        keep_this_coltimes = True
        for conditional in conditionals:
            if not conditional(coltimes):
                keep_this_coltimes = False
                break
        if not keep_this_coltimes:
            continue
        input_file = os.path.join(os.path.dirname(first_file), 'input', 'input')
        parsed_input = ParsedInput(input_file)
        for param in parsed_input.ordered_param_names:
            coltimes[param] = parsed_input.params[param]
        coltimes.to_csv(os.path.join(run_dir, outfile_name),
                                     index=False, mode='a', header=False)

def combine_complete_coltime_csvs(run_dir):
    combine_coltimes_csvs(run_dir, [], 'complete_coltimes.csv')

def add_useful_cols_to_coltimes(coltimes):
    """Only works on coltimes that come from combine_coltimes_csvs."""
    coltimes['deltaR'] = coltimes['L']/(coltimes['NB'] - 1)
    coltimes['linear_distance'] = coltimes['deltaR']*np.abs(coltimes['i'] - coltimes['j'])

def combine_initial_dists_csvs(run_dir, *args, **kwargs):
    """Not yet implemented."""
    pass

def aggregate_r0_group_NP_L0(run_dir, overwrite=True, *args, **kwargs):
    e2e_folder = os.path.join(run_dir, 'data', 'end_to_end_r0')
    os.makedirs(e2e_folder, mode=0o755, exist_ok=True)
    existing_r0_files = {}
    for folder in glob.glob(os.path.join(run_dir, '*.*')):
        if not os.path.isdir(folder):
            continue
        try:
            sim = Sim(folder, *args, **kwargs)
        except FileNotFoundError:
            logger.warning('aggregate_r0_group_NP_L0: Unable to construct sim from ' + folder)
            continue
        # #beads
        N = sim.params['NB']
        # length in units of simulation
        L = sim.params['L']
        r0_agg_file = os.path.join(e2e_folder, 'r0_' + str(N) + '_' + str(L) + '.csv')
        r0_file = os.path.join(folder, 'data', 'r0')
        if os.path.isfile(r0_agg_file) and not overwrite:
            continue
        if r0_agg_file not in existing_r0_files:
            existing_r0_files[r0_agg_file] = True
            try:
                shutil.copyfile(r0_file, r0_agg_file)
            except OSError:
                logger.warning('r0 not found for: ' + folder)
                open(r0_agg_file, 'a').close() # make sure file exists later in loop
            continue
        else:
            try:
                r0_in = open(r0_file, 'r')
            except:
                logger.warning('r0 not found for: ' + folder)
                continue
            with open(r0_agg_file, 'a') as r0_out:
                for line in r0_in:
                    r0_out.write(line)
            r0_in.close()
