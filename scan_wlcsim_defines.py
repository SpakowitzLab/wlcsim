#!/usr/bin/env python3
import pscan
import numpy as np
import re
import os
import socket # for getting hostname
import subprocess # for git hash and moving files
import shutil # for copying files
from distutils.dir_util import copy_tree # for copying directories
from pathlib import Path # for "touching" files
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
run_name = 'test-nuc'
output_base = 'par-run-dir'

# the following uses numpy to create arrays of parameters
# (http://docs.scipy.org/doc/numpy/reference/routines.array-creation.html)
# don't forget to specify dtype=np.int32 if you need integers!

# # to vary parameters combinatorially, list all the values for all parameters
# # you want like this, all combinations will be exected automatically
# params['DT'] = np.array([0.1, 1.0, 2.0])
# params['collisionRadius'] = np.array([0.1, 1.0, 2.0])
# # to vary parameters jointly, make dictionaries with values of matching size
# # like this. see pscan.py for more details.
# jparam['L'] = np.linspace(10, 100, 10, dtype=np.int32)
# jparam['N'] = np.linspace(20, 200, 10, dtype=np.int32) + 1
# jparams.append(jparam)

# to change how many times each parameter set is run, change number below
# for more advanced control of how many times to run sims based on param
# values, see the docs of pscan.py
default_repeats_per_param = 100
count_funcs.append(lambda p: default_repeats_per_param)
#count_funcs.append(lambda p: max(1,int(default_repeats_per_param/50)) if p['DT'] < 0.1 else None)

# how many cores to use on this computer
num_cores = multiprocessing.cpu_count() - 10

## END PARAM SCANNING CONFIG

script_dir = os.path.dirname(os.path.realpath(__file__))
os.chdir(script_dir)

# read in parameter names and "default" values
# from files in Andy's input format
params_file = os.path.join('src', 'defines.inc')
input_file = os.path.join(script_dir, params_file)
print(input_file)
parsed_input = wlcsim.input.ParsedInput(input_file)
ordered_param_names = parsed_input.ordered_param_names
# update params dict with the values requested by the user
# remembering that everything must be in allcaps, and that pscan wants a list
# of params
simulation_params = {name.upper(): [value] for name, value in
                     parsed_input.params.items()}
print(simulation_params)
simulation_params.update(params)
scan = pscan.Scan(simulation_params)
for jparam in jparams:
    scan.add_jparam(jparam)
for count_func in count_funcs:
    scan.add_count(count_func)

hostname = socket.gethostname()
run_base = os.path.join(output_base, run_name, hostname)

# now make a "skeleton" directory that we can clone each time we run a new
# simulation. Make sure the skeleton directory is unique to this process, so
# that we can run multiple copies of this process at once outputing to the same
# directory
skel_num = 0
while True:
# as long as the OS is sane, only one thread can successfully mkdir
    try:
# periods aren't allowed in hostnames, so use as delimiter
        skeleton_dir = run_base + '.skeleton' + str(skel_num)
        os.makedirs(skeleton_dir, mode=0o755)
        break
    except OSError:
        skel_num = skel_num + 1
with open(os.path.join(skeleton_dir, 'commit_hash'), 'w') as commit_hash:
    subprocess.run(['git', 'rev-parse', 'HEAD'], stdout=commit_hash)
# some versions of the code have compile-time parameters, in which case we need
# to copy over the whole src directory and makefile
if parsed_input.input_format == wlcsim.input.InputFormat.DEFINES:
    #TODO fill this in
    shutil.copyfile(os.path.join(script_dir, 'Makefile'),
                    os.path.join(skeleton_dir, 'Makefile'))
    copy_tree('src', os.path.join(skeleton_dir, 'src'))
# some have an input file read at runtime, in which case we only need to copy
# over the executable file
else:
    shutil.copyfile(os.path.join(script_dir, 'wlcsim.exe'),
                    os.path.join(skeleton_dir, 'wlcsim.exe'))
    os.chmod(os.path.join(skeleton_dir, 'wlcsim.exe'), 0o755)
# copy over ourselves to make inspecting an old scan easier
shutil.copyfile(os.path.join(script_dir, 'scan_wlcsim.py'),
                os.path.join(skeleton_dir, 'scan_wlcsim.py'))
# make pre-filled input directory
copy_tree('input', os.path.join(skeleton_dir, 'input'))
# make empty output directory
os.mkdir(os.path.join(skeleton_dir, 'data'))

# the function that our process mapper will pass each set of params to
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
    copy_tree(skeleton_dir, run_dir)
    os.chdir(run_dir)
    # save parameter defaults used for the run
    defaults_file = os.path.join('input', 'input.defaults')
    shutil.copyfile(params_file, defaults_file)
    # make new input file from parameters provided
    wlcsim.input.ParsedInput(params=params).write(params_file)
    # now we're in the right directory, with input and output ready, go ahead
    # and run the simulation
    with open('./data/wlcsim.log', 'w') as f:
        if not os.path.exists('./wlcsim.exe'):
            subprocess.run(['make'], stdout=f, stderr=subprocess.STDOUT)
        subprocess.run(['time', './wlcsim.exe'], stdout=f, stderr=subprocess.STDOUT)
    # finally, stick a file flag saying that we've finished our simulation
    # here, for ease of sifting through runs
    Path('complete').touch()
    os.chdir(script_dir)
    return params

if __name__ == '__main__':
    p = multiprocessing.Pool(num_cores)
    script_name = os.path.basename(__file__)
    print(script_name + ': Running scan!')
    # # non-parallel version for testing
    # for params in scan.params():
    #     run_wlcsim(params)
    #     print(script_name + ": " + datetime.datetime.now().isoformat()
    #           + ": completed run with params: " + str(params))
    # parallel version
    for params in p.imap_unordered(run_wlcsim, scan.params(), chunksize=1):
        print(script_name + ": " + datetime.datetime.now().isoformat()
              + ": completed run with params: " + str(params))
    print(script_name + ': Completed scan!')
