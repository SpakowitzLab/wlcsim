#!/usr/bin/env python3
import pscan
import numpy as np
import re
import os
import socket # for getting hostname
import subprocess # for git hash and moving files
import shutil # for copying files
import multiprocessing
import datetime
import wlcsim.input

params = {}
jparams = []
count_funcs = []

## START PARAM SCANNING CONFIG

# specify where simulations will be run and where output will go
# you'll probably want to CHANGE run_name FOR EACH PARAMETER SCAN YOU DO to
# preserve your sanity
run_name = 'increasing-chain-length-0.1_fpt_dist'
output_base = 'par-run-dir'

# the following uses numpy to create arrays of parameters
# (http://docs.scipy.org/doc/numpy/reference/routines.array-creation.html)
# don't forget to specify dtype=np.int32 if you need integers!

# to vary parameters combinatorially, list all the values for all parameters
# you want like this, all combinations will be exected automatically
params['FPT_DIST'] = np.array([0.1])
params['DT'] = np.array([0.5])
# 2 bead chains basically never loop since no bend
# similarly, < 1 bead per persistence length seems to give silly answers
# similarly, < 2 bead per persistence length quantitatively incorrect answers
params['TF'] = np.array([10000000]) # max sim time
params['INDMAX'] = np.array([10000]) # max num saves (determines save freq)
params['SAVE_RU'] = np.array([0])
params['EXIT_WHEN_COLLIDED'] = np.array([1])
# to vary parameters jointly, make dictionaries with values of matching size
# like this. see pscan.py for more details.
jparam = {}
jparam['L'] = np.linspace(10, 100, 10, dtype=np.int32)
jparam['N'] = np.linspace(20, 200, 10, dtype=np.int32) + 1
jparams.append(jparam)

# to change how many times each parameter set is run, change number below
# for more advanced control of how many times to run sims based on param
# values, see the docs of pscan.py
default_repeats_per_param = 10
count_funcs.append(lambda p: default_repeats_per_param)
#count_funcs.append(lambda p: max(1,int(default_repeats_per_param/50)) if p['DT'] < 0.1 else None)

# how many cores to use on this computer
num_cores = multiprocessing.cpu_count()

## END PARAM SCANNING CONFIG

script_dir = os.path.dirname(os.path.realpath(__file__))
os.chdir(script_dir)

# read in parameter names and "default" values
# from files in Andy's input format
input_file = os.path.join(script_dir, 'input', 'input')
ordered_param_names, simulation_params = wlcsim.input.read_file(input_file)
simulation_params.update(params)
scan = pscan.Scan(simulation_params)
for jparam in jparams:
    scan.add_jparam(jparam)
for count_func in count_funcs:
    scan.add_count(count_func)

hostname = socket.gethostname()
run_base = os.path.join(output_base, run_name, hostname)
def run_wlcsim(params):
    run_num = 0
    while True:
# as long as the OS is sane, only one thread can successfully mkdir
        try:
# periods aren't allowed in hostnames, so use as delimiter
            run_dir = run_base + '.' + str(run_num)
            os.makedirs(run_dir, mode=0o755)
            break
        except OSError:
            run_num = run_num + 1
    with open(os.path.join(run_dir, 'commit_hash'), 'w') as commit_hash:
        subprocess.run(['git', 'rev-parse', 'HEAD'], stdout=commit_hash)
    shutil.copyfile(os.path.join(script_dir, 'wlcsim.exe'),
                    os.path.join(run_dir, 'wlcsim.exe'))
    # copy over ourselves to make inspecting an old scan easier
    shutil.copyfile(os.path.join(script_dir, 'scan_wlcsim.py'),
                    os.path.join(run_dir, 'scan_wlcsim.py'))
    os.chmod(os.path.join(run_dir, 'wlcsim.exe'), 0o755)
    # make pre-filled input directory
    shutil.copytree('input', os.path.join(run_dir, 'input'))
    # make empty output directory
    os.mkdir(os.path.join(run_dir, 'data'))
    os.chdir(run_dir)
    # make input file from parameters provided
    wlcsim.input.write_file('input/input', ordered_param_names, params)
    # now we're in the right directory, with input and output ready, go ahead
    # and run the simulation
    with open('./data/wlcsim.log', 'w') as f:
        subprocess.run(['time', './wlcsim.exe'], stdout=f, stderr=subprocess.STDOUT)
    os.chdir(script_dir)
    return params

if __name__ == '__main__':
    p = multiprocessing.Pool(num_cores)
    script_name = os.path.basename(__file__)
    print(script_name + ': Running scan!')
    for params in p.imap_unordered(run_wlcsim, scan.params(), chunksize=1):
        print(script_name + ": " + datetime.datetime.now().isoformat()
              + ": completed run with params: " + str(params))
    print(script_name + ': Completed scan!')
