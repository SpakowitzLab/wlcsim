""" This module "understands" the input format of wlcsim.exe """
import re
import os
# def __init__(self, sim_dir):
#     input_dir = os.path.join(sim_dir, 'input')
#     input_file = os.path.join(input_dir, 'input')
#     if not os.path.isdir(sim_dir) \
#             or not os.path.isdir(input_dir) \
#             or not os.path.isfile(input_file):
#         raise IOError('Simulation directory ' + sim_dir + ' does not'
#                       ' contain input/input')

def write_file(input_file, ordered_param_names, params):
    """ writes out a valid input file for wlcsim.exe given parameters
        in the format returned by read_file """
    with open(input_file, 'w') as f:
        # write the three garbage lines
        f.write('!\n!\n\n')
        for i,name in enumerate(ordered_param_names):
            f.write('!-Record ' + str(i) + '\n!  ' + name + '\n')
            f.write(str(params[name]) + '\n\n')

def read_file(input_file):
    """ returns a list: ordered_param_names, and a dict: simulation_params """
    # for a spec of the input file, see how Andy uses it in wlcsim.f90 and
    # similar basically triplets of lines with record number, name, then
    # value
    record_re = re.compile('!-Record (\d+)')
    name_re = re.compile('! *([_A-Za-z0-9]+)')
    contains_period_re = re.compile('\.')
    next_line_is_val = False
    ordered_param_names = []
    simulation_params = {}
    # if we get the sim dir or the input dir, resolve to the actual input file
    if not os.path.isfile(input_file) and os.path.isdir(input_file):
        input_file = os.path.join(input_file, 'input')
    if not os.path.isfile(input_file) and os.path.isdir(input_file):
        input_file = os.path.join(input_file, 'input')
    with open(input_file) as f:
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
    return ordered_param_names, simulation_params
