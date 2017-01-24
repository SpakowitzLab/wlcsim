""" This module "understands" wlcsim.exe's output, and will make it into
nice, pretty numpy/pandas arrays for you."""
import numpy as np
import pandas as pd
import os
import logging
import glob
import re
from . import input as winput

logger = logging.getLogger(__name__)

# keep these arround for cache reasons
# initialize with size 1
i_arr = np.array([[1]])
j_arr = np.array([[1]])
ui = (np.array([0]), np.array([0]))

def df_from_coltimes(coltimes, method=3):
    if method == 1:
        return df_from_coltimes_method1(coltimes)
    elif method == 2:
        return df_from_coltimes_method2(coltimes)
    elif method == 3:
        return df_from_coltimes_method3(coltimes)

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
    """Sim provides an interface to the heterogenous data produced by a
    simulation that saves coltimes, positions, and momentums.
    
    Unless specified during initialization, the coltimes, r, and u are loaded
    on demand upon first access."""
    def __init__(self, sim_dir, load_coltimes=False, load_ru=False,
                 load_method=3):
        # unpythonic?
        # if not os.path.isdir(sim_dir):
        #     raise IOError('Directory (' + sim_dir +
        #                   ') provided when constructing'
        #                   ' WlcsimData does not exist!')
        self.sim_dir = sim_dir
        input_file = os.path.join(sim_dir, 'input')
        self.param_names, self.params = winput.read_file(input_file)
        data_dir = os.path.join(sim_dir, 'data')
        coltimes_file = os.path.join(data_dir, 'coltimes')
        rfile_base = os.path.join(data_dir, 'r')
        ufile_base = os.path.join(data_dir, 'u')
        self.coltimes = pd.DataFrame()
        if not os.path.isfile(coltimes_file):
            logger.warning(coltimes_file + ' does not exist!')
        else:
            coltimes = pd.read_csv(coltimes_file, delim_whitespace=True,
                                   header=None).values
            self.coltimes = df_from_coltimes(coltimes, method=load_method)
        if not load_ru:
            self.r = np.zeros((0,0))
            self.u = np.zeros((0,0))
            return
        ti = 0
        rfiles = glob.glob(rfile_base + '*')
        ufiles = glob.glob(ufile_base + '*')
        number_from_filename = re.compile(re.escape(rfile_base) + '([0-9]+)')
        time_points = [re.search(number_from_filename, file_name).groups()[0]
                       for file_name in rfiles]
        try:
            last_ti = np.argmax(time_points)
            r0 = pd.read_csv(rfile_base + '0', delim_whitespace=True,
                             header=None).values
            self.r = np.full(r0.shape + (last_ti+1,), np.nan)
            self.u = np.full(r0.shape + (last_ti+1,), np.nan)
        # if time_points was empty, no r or u was saved: ValueError
        except ValueError:
            last_ti = -1
            self.r = np.full((0,0), np.nan)
            self.u = np.full((0,0), np.nan)
            logger.warning('Initial data not found, do both r0 and u0 exist?')
        # if other time points were saved, but for some reason teh initial one
        # was not...this actually happened.. -.- OSError if reading with pd
        # FileNotFoundError if usign loadtxt
        except (FileNotFoundError, OSError):
            logger.warning('Initial data not found, do both r0 and u0 exist?')
        for ti in range(last_ti+1):
            rfile_name = rfile_base + str(ti)
            ufile_name = ufile_base + str(ti)
            try:
                ri = pd.read_csv(rfile_name, delim_whitespace=True, header=None).values
                self.r[:,:,ti] = ri
            except OSError:
                logger.warning(rfile_name + ' does not appear to exist')
            try:
                ui = pd.read_csv(ufile_name, delim_whitespace=True, header=None).values
                self.u[:,:,ti] = ui
            except OSError:
                logger.warning(ufile_name + ' does not appear to exist')

    def load_r(self, ti, index_type="index"):
        """ loads the ti'th time point. second argument currently does nothing.
            TODO: allow index_type=sim_time to trigger parsing of input file
            and load the time point corresponding to the simulation time
            provided instead of the index of the saved time point.  """
        try
            ri = pd.read_csv(


class Scan:
    """ Combine all the input params and output data from a parameter scan """
    def __init__(self, par_run_dir, load_ru=False, load_method=3):
        self.wlcsim_data = [Sim(folder, load_ru=load_ru) for folder in
                            glob.glob(os.path.join(par_run_dir, '*'))
                            if os.path.isdir(folder)]
        self.param_names = self.wlcsim_data[0].param_names
        for sim in self.wlcsim_data:
            sim.coltimes['run_name'] = os.path.basename(sim.sim_dir)
            for param in sim.params:
                sim.coltimes[param] = sim.params[param]
        coltimes = [sim.coltimes for sim in self.wlcsim_data]
        self.coltimes = pd.concat(coltimes)
        self.calculate_linear_distance()
    def calculate_linear_distance(self):
# actual distance along polymer is dx*(number_separating)
        self.coltimes['dx'] = self.coltimes['L']/(self.coltimes['N'] - 1)
        self.coltimes['linear_distance'] = self.coltimes['dx']*np.abs(self.coltimes['i'] - self.coltimes['j'])
    def calculate_for_each_param(self, f):
        """unimplemented"""
        pass

    @staticmethod
    def write_coltimes_csvs(par_run_dir, overwrite=False, *args, **kwargs):
        for folder in glob.glob(os.path.join(par_run_dir, '*')):
            if not os.path.isdir(folder):
                continue
            coltimes_file = os.path.join(folder, 'coltimes.csv')
            if not overwrite and os.path.isfile(coltimes_file):
                continue
            sim = Sim(folder, *args, load_ru=False, **kwargs)
            sim.coltimes.to_csv(coltimes_file)

    @staticmethod
    def write_initial_dists_csvs(par_run_dir, overwrite=True, *args, **kwargs):
        for folder in glob.glob(os.path.join(par_run_dir, '*')):
            if not os.path.isdir(folder):
                continue
            dists_file = os.path.join(folder, 'initial_dists.csv')
            if not overwrite and os.path.isfile(dists_file):
                continue
            sim = Sim(



    @staticmethod
    def write_coltimes_pkls(par_run_dir, *args, **kwargs):
        for folder in glob.glob(os.path.join(par_run_dir, '*')):
            if not os.path.isdir(folder):
                continue
            sim = Sim(folder, *args, load_ru=False, **kwargs)
            sim.coltimes.to_pickle(os.path.join(folder, 'coltimes.pkl'))

    # @staticmethod
    # def combine_coltimes_pkls(par_run_dir, *args, **kwargs):
    #     pkl_glob = os.path.join(par_run_dir, '*', 'coltimes.pkl')
    #     for coltimes_pkl in glob.glob(pkl_glob):

    @staticmethod
    def combine_coltimes_csvs(par_run_dir, *args, **kwargs):
        csv_glob = os.path.join(par_run_dir, '*', 'coltimes.csv')
        csv_files = iter(glob.glob(csv_glob))
        first_file = next(csv_files)
        coltimes = pd.read_csv(first_file, index_col=0)
        coltimes['run_name'] = os.path.basename(os.path.dirname(first_file))
        param_names, params = winput.read_file(os.path.join(os.path.dirname(first_file), 'input', 'input'))
        for param in param_names:
            coltimes[param] = params[param]
        coltimes.to_csv(os.path.join(par_run_dir, 'coltimes.csv'), index=False)
        # now append the remaining
        for csv_file in csv_files:
            coltimes = pd.read_csv(csv_file, index_col=0)
            coltimes['run_name'] = os.path.basename(os.path.dirname(csv_file))
            param_names, params = winput.read_file(os.path.join(
                    os.path.dirname(csv_file), 'input', 'input'))
            for param in param_names:
                coltimes[param] = params[param]
            coltimes.to_csv(os.path.join(par_run_dir, 'coltimes.csv'),
                                         index=False, mode='a', header=False)

