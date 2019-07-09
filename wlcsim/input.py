""" This module "understands" the input format of wlcsim.exe """
from __future__ import print_function

import re
import os
from pathlib import Path
from enum import Enum

# def __init__(self, sim_dir):
#     input_dir = os.path.join(sim_dir, 'input')
#     input_file = os.path.join(input_dir, 'input')
#     if not os.path.isdir(sim_dir) \
#             or not os.path.isdir(input_dir) \
#             or not os.path.isfile(input_file):
#         raise IOError('Simulation directory ' + sim_dir + ' does not'
#                       ' contain input/input')

class InputFormat(Enum):
    ORIGINAL=1
    LENA=2
    DEFINES=3

renamer = {'COL_TYPE': 'COLLISIONDETECTIONTYPE', 'COLTYPE': 'COLLISIONDETECTIONTYPE', 'FPT_DIST': 'COLLISIONRADIUS',
           'INTON': 'INTERPBEADLENNARDJONES', 'N': 'NB', 'INDMAX': 'NUMSAVEPOINTS'}

def correct_param_name(name):
    """Takes messy param names from different generations of the simulator and
    converts them to their newest forms so that simulations from across
    different years can be tabulated together."""
    name = name.upper()
    if name not in renamer:
        return name
    else:
        return renamer[name]

def brown_to_codename(name, value):
    new_name = 'CODENAME'
    if int(value) == 1:
        return (new_name, 'bruno')
    else:
        return (new_name, 'brad')
# for each old param name that needs to go through correct_param_value, the
# function that will convert its value. recall all "values" are strings at this
# point
revalueer = {'BROWN': brown_to_codename}

def correct_param_value(name, value):
    """Some old param names also have new types. This takes messy param names
    from different generations of the simulator and converts their names and
    values to their newest forms so that simulations from across different
    years can be tabulated together."""
    if name in revalueer:
        return revalueer[name](name, value)
    else:
        return (name, value)

class ParsedInput(object):
    """Knows how to handle various input file types used by wlcsim simulator
    over the years, and transparently converts into new parameter naming
    conventions.

    input = ParsedInput(file_name)
    print(input.ordered_param_names) # see params in order defined
    print(input.ordered_param_values) # to see values
    input.write(outfile_name) # write clone of input file

    """

    def __init__(self, input_file=None, params=None):
        """Can be constructed from an input file or from params directly. If an
        input file is provided, params are just ignored if passed in."""
        self.params = {}
        self.ordered_param_names = []
        self.file_format = InputFormat.DEFINES
        self.input_file = Path(input_file)
        if input_file is not None:
            # if we get the sim dir or the input dir, resolve to the actual input file
            if not os.path.isfile(input_file) and os.path.isdir(input_file):
                input_file = os.path.join(input_file, 'input')
            if not os.path.isfile(input_file) and os.path.isdir(input_file):
                input_file = os.path.join(input_file, 'input')
            self.decide_input_format()
            self.parse_params_file()
        elif params is not None:
            for name, value in params.items():
                self.append_param(name, value)
        else:
            Warning('ParsedInput: no params or input file provided!')

    def __repr__(self):
        return str(self.params)

    @property
    def ordered_param_values(self):
        return [self.params[name] for name in self.ordered_param_names]

    def append_param(self, name, value):
        name = correct_param_name(name)
        name, value = correct_param_value(name, value)
        self.ordered_param_names.append(name)
        self.params.update({name: value})

    def write(self, output_file_name):
        """ writes out a valid input file for wlcsim.exe given parameters
            in the format returned by read_file """
        with open(output_file_name, 'w') as f:
            if self.file_format == InputFormat.ORIGINAL:
                # write the three garbage lines
                f.write('!\n!\n\n')
                for i,name in enumerate(self.ordered_param_names):
                    f.write('!-Record ' + str(i) + '\n!  ' + name + '\n')
                    f.write(str(self.params[name]) + '\n\n')
            elif self.file_format == InputFormat.LENA:
                for name in self.ordered_param_names:
                    f.write(str(name) + ' ' + str(self.params[name]) + '\n')
            elif self.file_format == InputFormat.DEFINES:
                for name in self.ordered_param_names:
                    f.write('#define WLC_P__' + str(name) + ' ' +
                            str(self.params[name]) + '\n')
            else:
                raise ValueError('wlcsim.input: attempt to print a ParsedInput'
                                 ' with unknown file_format.')

    def decide_input_format(self):
        """Decide between the two input formats we know of.
        Not too hard, since one uses Fortran-style comments, which we can look
        out for, and the other uses bash style comments. Further, the former
        specifies param names and values on separate lines, while the latter
        specifies them on the same line."""
        # first see if the file has the expected name for the defines file
        if self.input_file.name == 'defines.inc':
            self.input_format = InputFormat.DEFINES
            return
        # then see if there are any comment lines. if so, we immediately know
        # the file type
        with open(self.input_file) as f:
            for line in f:
                if not line:
                    continue
                elif re.match('\s*!', line):
                    self.input_format = InputFormat.ORIGINAL
                    return
                elif line[0] == '#':
                    self.input_format = InputFormat.LENA
                    return
                else:
                    continue
        # if there are no comments, then for now, it must in fact be Lena's
        # intput file type, otherwise, we would not be able to infer the param
        # names, since these are in comment lines in the original-type input
        # files
        self.input_format = InputFormat.LENA
        return


    def parse_params_file(self):
        """Parse and populate ParsedInput's params, ordered_param_names
        This parser currently understands three file formats:
        1) "ORIGINAL" is the input method understood by Andy's hardcoded
        input reader in the original code used the Spakowitz lab was
        founded.
        2) "LENA" is a spec using slightly more general input reader
        written by Elena Koslover while Andy's student.
        3) "DEFINES" is the format of the src/defines.inc file.
        """
        if self.input_format == InputFormat.ORIGINAL:
            self.parse_original_params_file()
        elif self.input_format == InputFormat.LENA:
            self.parse_lena_params_file()
        elif self.input_format == InputFormat.DEFINES:
            self.parse_defines_params_file()

    def parse_lena_params_file(self):
        """Lena-style input files have comment lines starting with a "#". Any
        other non-blank lines must be of the form
        "[identifier][whitespace][value]",
        where an identifier is of the form "[_a-zA-Z][_a-zA-Z0-9]*", and a
        value can be a boolean, float, int or string. They will always be
        stored as strings in the params dictionary for downstream parsing as
        needed.

        Identifiers, like fortran variables, are interpreted in a
        case-insensitive manner by the wlcsim program, and so will be store in
        all-caps within the ParsedInput to signify this."""
        name_value_re = re.compile('([_A-Za-z][_A-Za-z0-9]*)\s*(.*)\s*')
        with open(self.input_file) as f:
            for line in f:
                if not line or line[0] == '#':
                    continue
                match = name_value_re.match(line)
                if match is None:
                    continue
                name, value = match.groups()
                self.append_param(name, value)


    def parse_defines_params_file(self):
        """Parse file in the format of src/defines.inc. Each line begins with
        #define WLC_P__[PARAM_NAME] [_a-zA-Z0-9]
        where WLC_P__[A-Z]* is the parameter name and the piece after the space is the value of
        the parameter.

        TODO:test
        """
        name_value_re = re.compile('#define WLC_P__([_A-Z]*) ([_a-zA-Z0-9]*)')
        with open(self.input_file) as f:
            for line in f:
                if not line:
                    continue
                match = name_value_re.match(line)
                if match is None:
                    continue
                name, value = match.groups()
                self.append_param(name, value)


    def parse_original_params_file(self):
        """Original-style input files have three garbage lines at the top used
        for commenting then triplets of lines describing one variable each. The
        first line of the triplet is a counter, for ease of hardcoding input
        reader, the second line contains the variable name (which is not used
        by the input reader) in order to make parsing possible outside of the
        ahrdcoded wlcsim input reader, and the third line contains the value of
        the parameter itself. The first two lines of each triplet have a fixed
        form that we use to extract the record numbers and parameter names, but
        these forms are not used by wlcsim itself, which ignores these lnies
        completely. Thus, it is possible that the user specified a particular
        name for a parameter but that name does not match what wlcsim
        interpreted it as, since wlcsim simply uses the order of the parameters
        in this file to determine their identities."""
        record_re = re.compile('!-Record (\d+)')
        name_re = re.compile('! *([_A-Za-z0-9]+)')
        contains_period_re = re.compile('\.')
        next_line_is_val = False
        with open(self.input_file) as f:
            # first three lines are garbage
            for i in range(3):
                f.readline()
            for line in f:
                if next_line_is_val:
                    if contains_period_re.search(line):
                        value = float(line.strip())
                    else:
                        value = int(line.strip())
                    self.append_param(name, value)
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
                    name = name.upper()
                    next_line_is_val = True
                    continue
