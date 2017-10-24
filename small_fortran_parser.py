import re
import os


num_to_move = {
    1: 'crank_shaft'
    2: 'slide_move'
    3: 'pivot_move'
    4: 'rotate_move'
    5: 'full_chain_rotation'
    6: 'full_chain_slide'
    7: 'change_binding_state'
    8: 'chain_flip'
    9: 'chain_exchange'
    10: 'reptation'
    11: 'super_reptation'
}
move_to_num = {move: num for num, move in num_to_move.items()}
num_to_dim = {1: 'x', 2: 'y', 3: 'z'}

def variable_from_fortran_line(line):
    var_line = re.compile("""
        \s*                                 # leading whitespace
        (real|integer|character|logical)    # type
        (\(dp|MAXPARAMLEN\))?               # optional "precision"
        \s+                                 # some whitespace
        ([_0-9a-zA-Z]+)                     # variable name
        (\([_a-zA-Z0-9]+\))?                # shape of the variable, opt.
        \s*(!.*)?                           # optional following comment
                          """, re.VERBOSE)
    match = var_line.search(line)
    if match is None:
        return None
    type, _, name, shape, comment = match.groups()
    if shape is not None:
        return None
    return name.upper(), shape, comment

def get_function_body(function_name, file):
    function_lines = []
    with open(file) as f:
        start_re = re.compile('subroutine\s' + function_name)
        end_re = re.compile('end subroutine')
        appending = False
        for line in f:
            if start_re.search(line):
                appending = True
                continue
            if end_re.search(line):
                break
            if appending:
                function_lines.append(line)
    return function_lines

def is_definition(line):
    return parse_definition(line) is not None
def parse_definition(line):
    definition_line = re.compile("""
        \s*                     # leading whitespace
        wlc_p%([_0-9a-zA-Z]+)   # variable name
        (\([_a-zA-Z0-9]\))?     # indexer
        \s*=\s*
        (\S+)                   # value, I made sure whitespace delineates def'n by hand
        \s*(!.*)?               # optional following comment
                                """, re.VERBOSE)
    match = definition_line.search(line)
    if not match:
        return None
    name, indexer, value, comment = match.groups()
    if indexer:
        return None
    return name.upper(), value, comment

def get_default_value(name):
    if name not in params_vals:
        raise ValueError(str(name) + " not in param_vals!")
    return params_vals[name]

def get_type_body(type_name, file):
    type_lines = []
    with open(file) as f:
        start_re = re.compile('type\s' + type_name)
        end_re = re.compile('end type')
        appending = False
        for line in f:
            if start_re.search(line):
                appending = True
                continue
            if end_re.search(line):
                break
            if appending:
                type_lines.append(line)
    return type_lines

param_defaults_body = get_function_body('set_param_defaults', 'src/wlcsim/params.f03')
params_defaults_definitions = [parse_definition(line) for line in
                               param_defaults_body if is_definition(line)]
params_vals = {p[0]: p[1] for p in params_defaults_definitions}
param_lines = get_type_body('wlcsim_params', 'src/wlcsim/params.f03')
all_params = [variable_from_fortran_line(line) for line in param_lines]
all_params = [p for p in all_params if p is not None]
pval_d = {p[0]: params_vals[p[0]] for p in all_params}
pcom_d = {p[0]: p[1] for p in all_params}

# print(all_params)



