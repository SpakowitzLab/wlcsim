import re
import os
import subprocess

import pandas as pd

self_dir, self_file = os.path.split(os.path.realpath(__file__))
src_dir = os.path.join(self_dir, '..', 'src')

mc_moves = pd.Series({
    1: 'crank_shaft',
    2: 'slide_move',
    3: 'pivot_move',
    4: 'rotate_move',
    5: 'full_chain_rotation',
    6: 'full_chain_slide',
    7: 'change_binding_state',
    8: 'chain_flip',
    9: 'chain_exchange',
    10: 'reptation',
    11: 'super_reptation'
})
xyz_names = pd.Series({1: 'x', 2: 'y', 3: 'z'})

#####{{{
# parsing functions

var_re = "[a-zA-Z0-9][_a-zA-Z0-9]*"
float_re = "[-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?"

def extract_from_lines(lines, object_type):
    """Return a DataFrame containing the appropriate columns for the requested
    object to be extracted linewise from a list of lines.

    Objects are parsed only into their text representations, like a
    pre-processor. Ask __dict__.object_types for valid objects to request."""
    object_types = {'declaration': declaration,
                    'wlc_p_default': wlc_p_default}
    if object_type not in object_types:
        raise ValueError(str(object_type) + ' should be one of the strings: '
                         + ', '.join(object_types.keys()))
    objects = [object_types[object_type](line) for line in lines]
    parsed = pd.DataFrame([obj for obj in objects if obj is not None])
    parsed.columns = object_types[object_type].column_names
    return parsed

def declaration(line):
    line_re = re.compile("""
        \s*                                 # leading whitespace
        (real|integer|character|logical)    # type
        (\((""" + var_re + """)\))?         # optional "precision", int or varname
        \s+                                 # some whitespace
        ([_0-9a-zA-Z]+)                     # variable name
        (\((""" + var_re + """)\))?         # shape of the variable, opt.
        \s*(!.*)?                           # optional following comment
                          """, re.VERBOSE)
    match = line_re.search(line)
    if match is None:
        return None
    # some groupings used for optionality, these are ignored
    type, _, type_size, name, _, shape, comment = match.groups()
    return type, type_size, name, shape, comment
declaration.column_names = ['type', 'type_size', 'name', 'shape', 'comment']

def wlc_p_default(line):
    line_re = re.compile("""
        \s*                     # leading whitespace
        wlc_p%(""" + var_re + """)   # variable name
        (\((""" + var_re + """)\))?  # indexer, capture also w/o parens
        \s*=\s*
        ([^\s!]+)               # value, terminated by whitespace or a comment
        \s*(!.*)?               # optional following comment
                                """, re.VERBOSE)
    match = line_re.search(line)
    if not match:
        return None
    name, _, indexer, value, comment = match.groups()
    return name, indexer, value, comment
wlc_p_default.column_names = ['name', 'indexer', 'value', 'comment']

# end parsing functions
#####}}}

#####{{{
# extract parts of the code

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

# end extract parts of the code
#####}}}


#####{{{
# do big code modifications

def wlc_p_to_defines():
    """gfortran is apparently bad at expression folding if there's not constants
    involved. So instead of precomputing all divisions by hand, we make all
    parameters compile time constants to try to offload this job into the
    compiler's constant folding step.

    This function extracts the wlcsim_params type, gets its default values from
    params.set_param_defaults, deletes unused variables from this type and from
    the defaults function, creates a 'defines.h' file, removes all mentions of
    wlc_p from parameter lists, and replaces all uses of 'wlc_p%NAME' with
    'wlc_p__NAME' by brute force regex.

    Special care must be taken with variables like minWindow. These are used as
    arrays, one index for each move. We keep the defaults lines for these
    variables, unrolling the for loops manually to have the defaults be like
    ```
    wlc_p%minWindow(2) = wlc_p__minWindow_slide
    ```

    Some variables are kept regardless, since they are defined by reading in
    from a table that Elena generated or depend in some "turing-complete" way
    on other input parameters.
    """
    vars_to_keep = ['EB', 'EPAR', 'EPERP', 'ESELF', 'GAM', 'ETA', 'XIR',
            'SIGMA', 'XIU', 'DEL', 'LHC', 'VHC', 'REND']
    wlc_p_usage_re = re.compile("[^!]*wlc_p%("+ var_re +")\(("+ var_re +")\)")
    wlc_p_usage_re_grep =      "^[^!]*wlc_p%" + var_re + "\(" + var_re + "\)"
    defaults = extract_from_lines(
            get_function_body('set_param_defaults', 'src/wlcsim/params.f03'),
            'wlc_p_default')
    for col in ['name', 'indexer', 'value']:
        defaults[col] = defaults[col].str.upper()
    variables = extract_from_lines(
            get_type_body('wlcsim_params', 'src/wlcsim/params.f03'),
            'declaration')
    for col in ['type', 'type_size', 'name', 'shape']:
        variables[col] = variables[col].str.upper()
    array_usages = subprocess.check_output(['grep', '-ohER',
            wlc_p_usage_re_grep, src_dir]).splitlines()
    array_usages = [wlc_p_usage_re.search(str(u)) for u in array_usages]
    array_usages = (pd.DataFrame(
            [m.groups() for m in array_usages if m is not None])
            .applymap(str.upper) # "parse" into caps, fortran case-insensitive
            .drop_duplicates())  # get unique name,indexer pairs
    array_usages.columns = ['name', 'indexer']
    is_looped_on = array_usages.indexer.str.contains('[A-Z]')
    looped_names = array_usages.name[is_looped_on]

    defines = []
    defines_columns = ['define_name', 'name', 'indexer', 'value',
                       'declaration comment', 'defaults comment',
                       'type', 'type_size', 'shape', 'is_looped_on',
                       'keep_regardless']
    # not_found = []
    for variable in variables.itertuples():
        if variable.shape is None:
            define_name = 'WLC_P__' + variable.name
            var_ix = defaults.name == variable.name
            if variable.name in vars_to_keep:
                value = None
            elif not var_ix.any():
                raise ValueError('Could not find default for variable ' + variable.name)
                # not_found.append(variable.name)
            else:
                value = defaults.value[var_ix]
            dcom = defaults.comment[var_ix]
            defines.append((define_name, variable.name, None,
                    value, variable.comment, dcom, variable.type,
                    variable.type_size, variable.shape,
                    variable.name in looped_names.values,
                    variable.name in vars_to_keep))
            continue
        if variable.shape == 'NMOVETYPES':
            name_map = mc_moves
        elif variable.shape == 'NDIM':
            name_map = xyz_names
        else:
            name_map = {i: str(i) for i in range(1, int(variable.shape)+1)}
        for idx, name in name_map.items():
            define_name = 'WLC_P__' + variable.name + '_' + name.upper()
            indexer = str(idx)
            var_ix = (defaults.name == variable.name) & (defaults.indexer == indexer)
            # it's possible that the variable is set in a loop to the same
            # value
            if not var_ix.any():
                var_ix = (defaults.name == variable.name) & defaults.indexer.str.contains('[A-Z]')
            # that it's set using vector syntax
            if not var_ix.any():
                var_ix = (defaults.name == variable.name) & defaults.indexer.isnull()
            # but all variables should have a default, even if it's a default
            # that forces an error, unless they're "derived"
            if variable.name in vars_to_keep:
                value = None
            elif not var_ix.any():
                raise ValueError('Could not find default for variable ' + variable.name)
                # not_found.append(variable.name)
            else:
                value = defaults.value[var_ix]
            dcom = defaults.comment[var_ix]
            defines.append((define_name, variable.name, indexer,
                    value, variable.comment, dcom, variable.type,
                    variable.type_size, variable.shape,
                    variable.name in looped_names.values,
                    variable.name in vars_to_keep))
    defines = pd.DataFrame(defines)
    defines.columns = defines_columns

    active_defines = defines.loc[~(defines.is_looped_on | defines.keep_regardless)]
    delete_declaration_lines(active_defines.name, os.path.join(src_dir, 'wlcsim', 'params.f03'))
    delete_definition_lines(active_defines.name, os.path.join(src_dir, 'wlcsim', 'params.f03'))

    defines_file = os.path.join(src_dir, 'defines.h')
    if os.path.exists(defines_file):
        raise OSError(defines_file + ' already exists!')
    with open(defines_file, 'w') as f:
        f.write("\n\n#define PI 4*atan(1.0_dp)\n\n")
        for new_define in active_defines.itertuples():
            f.write("VARIABLE COMMENT: {}\n")
            f.write("DEFAULTS COMMENT: {}\n")
            f.write("#define {} {}")
    #TODO awk command here




def defines_to_wlc_p():
    """The inverse of wlc_p_to_defines. //TODO"""
    #TODO implement
    pass



# end do big code modifications
#####}}}

# # first get all parameter defaults
# param_defaults_body = get_function_body('set_param_defaults', 'src/wlcsim/params.f03')
# param_defaults = pd.Series([parse_definition(line) for line in
#                             param_defaults_body if is_definition(line)])
# param_vals = {p[0]: [] for p in param_defaults_definitions}
# for p in param_defaults_definitions:
#     param_vals[p[0]].append((p[1], p[2]))
# # param_vals = {(p[0], p[1]): p[2] for p in param_defaults_definitions}

# # now get all parameter definitions
# param_lines = get_type_body('wlcsim_params', 'src/wlcsim/params.f03')
# all_params = [variable_from_fortran_line(line) for line in param_lines]
# all_params = [p for p in all_params if p is not None]

# # now use the parameter definitions to unroll into the right number of
# # variables for the array variables
# unrolled_params = {}
# unrolled_defaults = {}
# for name, shape, comment in all_params:
#     if shape is None:
#         unrolled_params[name] = comment
#         unrolled_defaults[name] = param_vals[name][1] # [1] gets just value
#     elif shape == '3':
#         unrolled_defaults[name] = param_vals[name][1] # [1] gets just value
#         for dimension in num_to_dim.values():
#             unrolled_params[name + '_' + dimension.upper()] = comment
#     elif shape == 'NMOVETYPES':
#         for move in mc_moves.values():
#             unrolled_params[name + '_' + move.upper()] = comment
#     else:
#         raise ValueError("Unknown shape '" + shape + "' for variable '" + name + "'.")
# param_name_to_comment = {p[0]: p[1] for p in all_params}

# pval_d = {p[0]: param_vals[p[0]] for p in all_params}


# print(all_params)



