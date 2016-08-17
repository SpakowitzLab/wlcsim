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



class Sim:
    def __init__(self, sim_dir):
        if not os.path.isdir(sim_dir):
            raise IOError('Directory (' + sim_dir +
                          ') provided when constructing'
                          ' WlcsimData does not exist!')
        self.sim_dir = sim_dir
        input_file = os.path.join(sim_dir, 'input')
        _, self.params = winput.read_file(input_file)
        data_dir = os.path.join(sim_dir, 'data')
        coltimes_file = os.path.join(data_dir, 'coltimes')
        rfile_base = os.path.join(data_dir, 'r')
        ufile_base = os.path.join(data_dir, 'u')
        if not os.path.isfile(coltimes_file):
            logger.warning(coltimes_file + 'does not exist!')
        else:
            coltimes = np.loadtxt(coltimes_file)
            self.coltimes = pd.DataFrame()
            # make a Series out of the lower triangular part of coltimes matrix
            for i in range(coltimes.shape[0]):
                for j in range(i-1):
                    self.coltimes = self.coltimes.append(
                        pd.Series({'i': i, 'j': j,
                                   'coltime': coltimes[i][j]}),
                        ignore_index=True
                    )
        ti = 0
        rfiles = glob.glob(rfile_base + '*')
        ufiles = glob.glob(ufile_base + '*')
        number_from_filename = re.compile(re.escape(rfile_base) + '([0-9]+)')
        time_points = [re.search(number_from_filename, file_name).groups()[0]
                       for file_name in rfiles]
        try:
            last_ti = np.argmax(time_points)
            r0 = np.loadtxt(rfile_base + '0')
            self.r = np.zeros(r0.shape + (last_ti+1,))
            self.u = np.zeros(r0.shape + (last_ti+1,))
        # if time_points was empty, no r or u was saved
        except ValueError:
            last_ti = -1
            self.r = np.zeros((0,0))
            self.u = np.zeros((0,0))
            logger.warning('Initial data not found, do both r0 and u0 exist?')
        for ti in range(last_ti+1):
            rfile_name = rfile_base + str(ti)
            ufile_name = ufile_base + str(ti)
            self.r[:,:,ti] = np.loadtxt(rfile_name)
            self.u[:,:,ti] = np.loadtxt(ufile_name)

class Scan:
    """ Combine all the input params and output data from a parameter scan """
    def __init__(self, par_run_dir):
        self.wlcsim_data = [Sim(folder) for folder in
                            glob.glob(os.path.join(par_run_dir, '*'))]
        for sim in self.wlcsim_data:
            sim.coltimes['run_name'] = os.path.basename(sim.sim_dir)
            for param in sim.params:
                sim.coltimes[param] = sim.params[param]
        coltimes = [sim.coltimes for sim in self.wlcsim_data]
        self.coltimes = pd.concat(coltimes)
# actual distance along polymer is dx*(number_separating)
        self.coltimes['dx'] = self.coltimes['L']/(self.coltimes['N'] - 1)
        self.coltimes['linear_distance'] = self.coltimes['dx']*np.abs(self.coltimes['i'] - self.coltimes['j'])
