""" This module "understands" wlcsim.exe's output, and will make it into
nice, pretty numpy/pandas arrays for you."""
import numpy as np
import pandas as pd
import os
import logging
import glob
import re
from . import input as winput
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


class Sim:
    """
    Sim provides an interface to the heterogenous data produced by a
    simulation that saves coltimes, positions, and momentums.

    Unless specified during initialization, the coltimes, r, and u are loaded
    on demand upon first access.

    Accessing properties of this class for the first time can return an IOError
    if a file is moved, or doesn't exist. Afterwards, the values will be cached.
    """

    def __init__(self, sim_path, overwrite_coltimes=False, load_method=3):
        # since these are all lazily evaluated, we throw exceptions right away
        # if the files we need do not exist, to prevent long running programs
        # from throwing an exception midway
        if not os.path.isdir(sim_path):
            raise FileNotFoundError('Directory (' + sim_path +
                                    ') provided when constructing'
                                    ' wlcsim.data.Sim does not exist!')
        self.sim_path = sim_path
        self.sim_dir = os.path.basename(sim_path)
        self.input_file = os.path.join(sim_path, 'input', 'input')
        # without input file, we don't really have a chance of using the
        # simulation data from a simulation anyway, so we can safely error out
        # here
        if not os.path.isfile(self.input_file):
            raise FileNotFoundError(self.input_file)
        self.data_dir = os.path.join(sim_path, 'data')
        self.coltimes_file = os.path.join(self.data_dir, 'coltimes')
        self.overwrite_coltimes = overwrite_coltimes
        self._rfile_base = os.path.join(self.data_dir, 'r')
        self._ufile_base = os.path.join(self.data_dir, 'u')
        self.load_method = load_method

    @cached_property
    def coltimes(self):
        coltimes = pd.read_csv(self.coltimes_file, delim_whitespace=True,
                                header=None).values
        coltimes = self.df_from_coltimes(coltimes, method=self.load_method)
        return coltimes

    @cached_property
    def r0(self):
        """ Numpy array of positions of each particle at each time point. """
        rfile_name = self._rfile_base + str(0)
        r0 = pd.read_csv(rfile_name, delim_whitespace=True, header=None)
        return r0.values.T

    @cached_property
    def u0(self):
        """ Numpy array of positions of each particle at each time point. """
        ufile_name = self._ufile_base + str(0)
        u0 = pd.read_csv(rfile_name, delim_whitespace=True, header=None)
        return u0.values.T

    @cached_property
    def r(self):
        """ Numpy array of positions of each particle at each time point. """
        num_time_points = len(self.rtime_idxs)
        r = np.full((num_time_points, 3, self.params['N']), 0.0)
        for i, ti in enumerate(self.rtime_idxs):
            rfile_name = self._rfile_base + str(ti)
            ri = pd.read_csv(rfile_name, delim_whitespace=True, header=None)
            r[i,:,:] = ri.values.T
        return r

    @cached_property
    def u(self):
        """ Numpy array of positions of each particle at each time point. """
        num_time_points = len(self.utime_idxs)
        u = np.full((num_time_points, 3, self.params['N']), 0.0)
        for i, ti in enumerate(self.utime_idxs):
            ufile_name = self._ufile_base + str(ti)
            ui = pd.read_csv(ufile_name, delim_whitespace=True, header=None)
            u[i,:,:] = ui.values.T
        return u

    @cached_property
    def center_of_mass(self):
        # "equal mass" beads, so just average positions
        return np.mean(self.r, axis=2)

    @cached_property
    def rtime_idxs(self):
        rfiles = glob.glob(self._rfile_base + '*')
        return self.time_idxs(self._rfile_base, rfiles)

    @cached_property
    def utime_idxs(self):
        ufiles = glob.glob(self._ufile_base + '*')
        return self.time_idxs(self._ufile_base, ufiles)

    def time_idxs(self, base, file_names):
        number_from_filename = re.compile(re.escape(base) + '([0-9]+)')
        return np.array([re.search(number_from_filename, file_name).groups()[0]
                         for file_name in file_names])

    @cached_property
    def param_names(self):
        param_names, _ = winput.read_file(self.input_file)
        return param_names

    @cached_property
    def params(self):
        _, params = winput.read_file(self.input_file)
        return params

    @cached_property
    def rtimes(self):
        return self.rtime_idxs * self.params['DT']

    @cached_property
    def utimes(self):
        return self.utime_idxs * self.params['DT']

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
        dr = self.params['L']/(self.params['N'] - 1)
        coltimes['linear_distance'] = dr*np.abs(coltimes['i'] - coltimes['j'])


def nan_fill_coltimes(coltimes):
    """Takes a raw coltimes pandas array from e.g. write_coltimes_csvs and puts
    nan's where the -1.0's are, i.e. where the simulation did not go out far
    enough."""
    coltimes.loc[coltimes['coltime'] == -1.0, 'coltime'] = np.nan



class Scan:
    """Can tell you what dirs in a run_dir are simulation directories, and what
    all the input files said, and more. """

    def __init__(self, run_dir, *args, **kwargs):
        is_sim = lambda d: Scan.is_simdir(d, run_dir)
        self.sim_dirs = list(filter(is_sim, os.listdir(run_dir)))
        prepend_path = lambda d: os.path.join(run_dir, d)
        self.sim_paths = list(map(prepend_path, self.sim_dirs))
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

    @staticmethod
    def is_simdir(dirname, run_dir):
        sim_name_re = re.compile('.*\.[0-9]+')
        abs_dirname = os.path.join(run_dir, dirname)
        input_file = os.path.join(abs_dirname, 'input', 'input')
        return sim_name_re.search(dirname) \
                and os.path.isdir(abs_dirname) \
                and os.path.isfile(input_file)

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

    # from back when we used this class to combine coltimes's from a scan
    # into one big pandas df
    # def calculate_linear_distance(self):
        # # actual distance along polymer is dx*(number_separating)
        # self.coltimes['dx'] = self.coltimes['L']/(self.coltimes['N'] - 1)
        # self.coltimes['linear_distance'] = self.coltimes['dx']*np.abs(self.coltimes['i'] - self.coltimes['j'])
    # def calculate_for_each_param(self, f):
        # """unimplemented"""
        # pass

def write_coltimes_csvs(run_dir, overwrite=False, *args, **kwargs):
    for folder in glob.glob(os.path.join(run_dir, '*.*')):
        if not os.path.isdir(folder):
            continue
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
    for folder in glob.glob(os.path.join(run_dir, '*.*')):
        if not os.path.isdir(folder):
            continue
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
    for folder in glob.glob(os.path.join(run_dir, '*.*')):
        if not os.path.isdir(folder):
            continue
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

def combine_coltimes_csvs(run_dir, *args, **kwargs):
    csv_glob = os.path.join(run_dir, '*', 'coltimes.csv')
    csv_files = iter(glob.glob(csv_glob))
    first_file = next(csv_files)
    coltimes = pd.read_csv(first_file, index_col=0)
    coltimes['run_name'] = os.path.basename(os.path.dirname(first_file))
    param_names, params = winput.read_file(os.path.join(os.path.dirname(first_file), 'input', 'input'))
    for param in param_names:
        coltimes[param] = params[param]
    coltimes.to_csv(os.path.join(run_dir, 'coltimes.csv'), index=False)
    # now append the remaining
    for csv_file in csv_files:
        coltimes = pd.read_csv(csv_file, index_col=0)
        coltimes['run_name'] = os.path.basename(os.path.dirname(csv_file))
        param_names, params = winput.read_file(os.path.join(
                os.path.dirname(csv_file), 'input', 'input'))
        for param in param_names:
            coltimes[param] = params[param]
        coltimes.to_csv(os.path.join(run_dir, 'coltimes.csv'),
                                        index=False, mode='a', header=False)

def add_useful_cols_to_coltimes(coltimes):
    """Only works on coltimes that come from combine_coltimes_csvs."""
    coltimes['dr'] = coltimes['L']/(coltimes['N'] - 1)
    coltimes['linear_distance'] = coltimes['dr']*np.abs(coltimes['i'] - coltimes['j'])

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
        N = sim.params['N']
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
