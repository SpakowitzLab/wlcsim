import re
import os
import subprocess

import pandas as pd
import numpy as np

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
    params.set_param_defaults, making a DataFrame with all the information
    needed for make_defines_inc to do it's job.

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
            'SIGMA', 'XIU', 'DEL', 'LHC', 'VHC', 'REND', 'L0', 'EPS',
            'DT', 'NBIN', 'NT', 'SIMTYPE',
                   # parallel tempered variables, quinn
            'CHI', 'MU', 'HA', 'HP1_BIND', 'KAP', 'CHI_L2',
            'RECENTER_ON', 'KAP_ON', 'CHI_ON', 'COUPLE_ON', 'FIELD_INT_ON',
            'BIND_ON', "CHI_L2_ON",
                   # parallel tempered variables, brad
            'LK' ]
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
                       'declaration_comment', 'defaults_comment',
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
                value = defaults.value[var_ix].values[0]
            if value is None:
                dcom = None
            else:
                dcom = defaults.comment[var_ix].values[0]
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
                value = defaults.value[var_ix].values[0]
            dcom = defaults.comment[var_ix].values[0]
            defines.append((define_name, variable.name, indexer,
                    value, variable.comment, dcom, variable.type,
                    variable.type_size, variable.shape,
                    variable.name in looped_names.values,
                    variable.name in vars_to_keep))
    defines = pd.DataFrame(defines)
    defines.columns = defines_columns
    return defines

def make_defines_inc(active_defines):
    """Takes defines from wlc_p_to_defines,

    deletes unused variables from the wlcsim_params type and from the defaults
    function, creates a 'defines.inc' file, and replaces all uses of 'wlc_p%NAME'
    with 'wlc_p__NAME' by brute force regex.
    """
    defines_file = os.path.join(src_dir, 'defines.inc')
    if os.path.exists(defines_file):
        raise OSError(defines_file + ' already exists!')
    with open(defines_file, 'w') as f:
        # f.write("\n\n#define PI 4*atan(1.0_dp)\n\n")
        for new_define in active_defines.itertuples():
            f.write("! VARIABLE COMMENT: {}\n".format(new_define.declaration_comment))
            f.write("! DEFAULTS COMMENT: {}\n".format(new_define.defaults_comment))
            if new_define.type_size:
                f.write("! TYPE: {}({})\n".format(new_define.type, new_define.type_size))
            else:
                f.write("! TYPE: {}\n".format(new_define.type))
            f.write("#define {} {}\n".format(new_define.define_name, new_define.value))

def temp_rename(active_defines, undo=False):
    names_in_replace_order = get_replace_order(active_defines.name.values)
    # now replace all instances of wlc_p variables with their new values
    for name in names_in_replace_order:
        if undo:
            program = ['find', src_dir, '-type', 'f', '-exec', 'sed', '-i',
                    's/TMPTMPTMP' + name + '/wlc_p%' + name+ '/gi',
                    '{}', ';']
        else:
            program = ['find', src_dir, '-type', 'f', '-exec', 'sed', '-i',
                    's/wlc_p%' + name + '/TMPTMPTMP' + name + '/gi',
                    '{}', ';']
        subprocess.check_output(program)

def replace_wlc_p_instances(active_defines):
    # def shape_to_int(shape):
    #     if shape == 'NDIM':
    #         return len(xyz_names)
    #     elif shape == 'NMOVETYPES':
    #         return len(mc_moves)
    #     else:
    #         return int(shape)
    affected_files = subprocess.check_output(['grep', '-iRl', 'wlc_p%', src_dir]).splitlines()
    for file in affected_files:
        # if os.path.basename(file) == b'defines.inc':
        #     continue
        if not re.search("\.f[0-9][0-9]'$", str(file)) and not re.search("\.inc'$", str(file)):
            continue
        # add the include line at the top for defines.inc
        if os.path.basename(file) != b"defines.inc":
            with open(file, 'r') as original: data = original.read()
            with open(file, 'w') as modified: modified.write('#include "../defines.inc"\n' + data)
        names_in_replace_order = get_replace_order(active_defines.name.values)
        # now replace all instances of wlc_p variables with their new values
        for name in names_in_replace_order:
            defines_with_name = active_defines[active_defines.name == name]
            for define in defines_with_name.itertuples():
                if define.shape is None:
                    subprocess.check_output(['sed', '-i',
                            's/wlc_p%' + define.name + '/' + define.define_name + '/gi',
                            file])
                else:
                    subprocess.check_output(['sed', '-i',
                            's/wlc_p%' + define.name + '(' + define.indexer + ')/'
                            + define.define_name + '/gi', file])

def get_replace_order(strings):
    """Return a list of tuples corresponding to the defines that we are going to
    replace, but in order such that all superstrings are replaced first."""
    lengths = np.array([-len(string) for string in strings])
    ix = np.argsort(lengths)
    return np.array(strings)[ix]

def perform_wlc_p_to_define_transform():
    defines = wlc_p_to_defines()
    to_delete = defines['shape'].isnull() & (~defines.keep_regardless)
    declarations_to_delete = defines[to_delete]
    declarations_to_keep = defines[~to_delete]
    delete_declaration_lines(declarations_to_delete.name.values, os.path.join(src_dir, 'wlcsim', 'params.f03'),
                             in_type='wlcsim_params')
    # delete_definition_lines(defines.name.values, os.path.join(src_dir, 'wlcsim', 'params.f03'),
    #                         in_subroutine='set_param_defaults')
    delete_subroutine('set_param_defaults', os.path.join(src_dir, 'wlcsim', 'params.f03'))
    delete_subroutine('read_input_file', os.path.join(src_dir, 'wlcsim', 'params.f03'))
    defines_to_be_written = defines.loc[~defines.keep_regardless]
    make_defines_inc(defines_to_be_written)
    temp_rename(declarations_to_keep)
    replace_wlc_p_instances(declarations_to_delete)
    temp_rename(declarations_to_keep, undo=True)
    add_looped_defaults(defines)

def delete_subroutine(name, file):
    with open(file, 'r') as f:
        lines = f.readlines()
    with open(file, 'w') as f:
        pos = "before subroutine"
        for line in lines:
            if pos == "before subroutine" and re.search("subroutine\s+" + name, line):
                    pos = "inside subroutine"
            elif pos == "inside subroutine" and re.search("end\s+subroutine", line):
                    pos = "after subroutine"
                    continue
            if pos != 'inside subroutine':
                f.write(line)

def add_looped_defaults(defines):
    params_file = os.path.join(src_dir, 'wlcsim', 'params.f03')
    with open(params_file, 'r') as f:
        old_file = f.readlines()
    with open(params_file, 'w') as f:
        pos = 'before module'
        for line in old_file:
            if pos == 'before module' and re.search('module\s+params', line):
                pos = 'before contains'
            if pos == 'before contains' and re.match('^\s*contains\s*$', line):
                pos = 'after contains'
                f.write(line)
                write_looped_defaults(f, defines)
                continue
            f.write(line)

def write_looped_defaults(f, defines):
    f.write("""
    subroutine set_param_defaults(wlc_p)
        implicit none
        ! WARNinG: changing this to intent(out) means that unassigned values
        ! here will become undefined upon return, due to Fortran's weird
        ! intent(out) semantics for records, this would require that a default
        ! value always be given to new parameters in wlc_p, else we would get a
        ! compile time catchable runtime error that is not caught by gcc as of v5.0
        !
        ! this is almost definitely undesireable, since "undefined" means the
        ! behavior will depend on which compiler is used
        type(wlcsim_params), intent(inout) :: wlc_p
""")
    for define in defines.itertuples():
        if define.indexer is None:
            continue
        f.write("        wlc_p%{}({}) = {}\n".format(
                define.name, define.indexer, define.define_name))
    f.write("""
    end subroutine set_param_defaults

""")

def delete_declaration_lines(names, file, in_type=None):
    with open(file, 'r') as f:
        lines = f.readlines()
    with open(file, 'w') as f:
        pos = "before type"
        for line in lines:
            if in_type and pos == "before type":
                if re.search("type\s+" + in_type, line):
                    pos = "inside type"
            elif in_type and pos == "inside type":
                if re.search("end\s+type", line):
                    pos = "after type"
            if not in_type or pos == "inside type":
                match = declaration(line)
                if match and match[2].upper() in names:
                    continue
            f.write(line)

def delete_definition_lines(names, file, in_subroutine=None):
    with open(file, 'r') as f:
        lines = f.readlines()
    with open(file, 'w') as f:
        pos = "before subroutine"
        for line in lines:
            if in_subroutine and pos == "before subroutine":
                if re.search("subroutine\s+" + in_subroutine, line):
                    pos = "inside subroutine"
            elif in_subroutine and pos == "inside subroutine":
                if re.search("end\s+subroutine", line):
                    pos = "after subroutine"
            if not in_subroutine or pos == "inside subroutine":
                match = wlc_p_default(line)
                if match and match[0].upper() in names:
                    continue
            f.write(line)

def defines_to_wlc_p():
    """The inverse of wlc_p_to_defines. //TODO"""
    #TODO implement
    pass



# end do big code modifications
#####}}}
