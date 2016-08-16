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

params = {}
jparams = []
count_funcs = []

## START PARAM SCANNING CONFIG

# specify where simulations will be run and where output will go
# you'll probably want to CHANGE run_name FOR EACH PARAMETER SCAN YOU DO to
# preserve your sanity
run_name = 'timing-test-small-chains'
output_base = 'par-run-dir'

# the following uses numpy to create arrays of parameters
# (http://docs.scipy.org/doc/numpy/reference/routines.array-creation.html)
# don't forget to specify dtype=np.int32 if you need integers!

# to vary parameters combinatorially, list all the values for all parameters
# you want like this, all combinations will be exected automatically
params['FPT_DIST'] = np.linspace(0.1, 1.5, 15)
params['DT'] = np.array([0.1, 0.01])
# to vary parameters jointly, make dictionaries with values of matching size
# like this. see pscan.py for more details.
jparam = {}
jparam['L'] = np.linspace(2, 50, 49, dtype=np.int32)
jparam['N'] = np.linspace(3, 51, 49, dtype=np.int32)
jparams.append(jparam)
# to change how many times each parameter set is run, change number below
# for more advanced control of how many times to run sims based on param
# values, see the docs of pscan.py
default_repeats_per_param = 10
count_funcs.append(lambda p: default_repeats_per_param)
count_funcs.append(lambda p: max(1,int(default_repeats_per_param/10)) if p['DT'] < 0.1 else None)

# how many cores to use on this computer
num_cores = multiprocessing.cpu_count() - 1

## END PARAM SCANNING CONFIG

# read in parameter names and "default" values
# from files in Andy's input format
record_re = re.compile('!-Record (\d+)')
name_re = re.compile('! *([_A-Za-z0-9]+)')
contains_period_re = re.compile('\.')
next_line_is_val = False
ordered_param_names = []
simulation_params = {}
with open('input/input') as f:
    # first three lines are garbage
    for i in range(3):
        f.readline()
    for line in f:
        if next_line_is_val:
            if contains_period_re.search(line):
                value = float(line.strip())
            else:
                value = int(line.strip())
            simulation_params[name] = value
            ordered_param_names.append(name)
            name = None
            record_number = None
            value = None
            next_line_is_val = False
        record_match = record_re.search(line)
        if record_match:
            record_number = int(record_match.groups()[0])
            continue
        name_match = name_re.search(line)
        if name_match:
            name = name_match.groups()[0]
            next_line_is_val = True
            continue

simulation_params.update(params)
scan = pscan.Scan(simulation_params)
for jparam in jparams:
    scan.add_jparam(jparam)
for count_func in count_funcs:
    scan.add_count(count_func)

script_dir = os.path.dirname(os.path.realpath(__file__))
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
    os.chmod(os.path.join(run_dir, 'wlcsim.exe'), 0o755)
    # make pre-filled input directory
    shutil.copytree('input', os.path.join(run_dir, 'input'))
    # make empty output directory
    os.mkdir(os.path.join(run_dir, 'data'))
    os.chdir(run_dir)
    # make input file from parameters provided
    with open('input/input', 'w') as f:
        # write the three garbage lines
        f.write('!\n!\n\n')
        for i,name in enumerate(ordered_param_names):
            f.write('!-Record ' + str(i) + '\n!  ' + name + '\n')
            f.write(str(params[name]) + '\n\n')
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
