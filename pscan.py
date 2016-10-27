"""A module whose sole purpose is to facilitate parameter scans via a generic
framework for specifying parameter values to scan.
Combinatorial parameter sweeps, joint parameter sweeps, and advanced options
for how many times to repeat simulations are available.
TODO: implement "scientific" sweeps, where one parameter is varied at a time
from baseline."""
from functools import reduce # for roll-your-own product()
import operator # for operator.mul in map in roll-your-own product()
import numpy as np # just for unravel_index
import collections

# class P:
#     """A param for PScan"""
#     def __init__(self, name, values, counts=1):
#         self.name = name
#         self.values = values
#         self.counts = counts

# class PList:
#     """Jointly varying parameters"""
#     def __init__(self, params):
#         self.params = params
#     @classmethod
#     def from_dict(cls, dic, keys=None):
#         if keys is None:
#             keys = dic.keys()
#         return cls([P(key, dic[key]) for key in keys])

def get_first(iterable, default=None):
    if iterable:
        for item in iterable:
            return item
    return default

class JointParameterListSizeError(Exception):
    """Raised when two parameters that are meant to vary jointly have a
    different number of values that they are supposed to take."""
    pass

def _check_joint_params(jparam):
    param_names = []
    num_params = []
    for key,val in jparam.items():
        param_names.append(key)
        num_params.append(len(val))
    matches_first = lambda i: num_params[i] != num_params[0]
    bad_inds = filter(matches_first, range(len(num_params)))
    bad_ind = get_first(bad_inds)
    if bad_ind:
        raise JointParameterListSizeError( \
            "Param {} set to vary jointly with {}, but their" \
            "sizes ({},{}) do not match".format(
                param_names[bad_ind], param_names[0],
                num_params[bad_ind], num_params[0]
            )
        )

# this code could replace the two places where we force the comb_params to be
# iterable if we thought this needed to be done by the user to prevent bugs
# for now seems like "just work" is the right choice.

# class ParameterListNotIterableError(Exception):
#     """Raised when a value in the input dictionary (supposed to be a list of
#     parameters to take) is not iterable."""
#     pass

# def _check_comb_param(key, val):
#     if val is None or not isinstance(val, collections.Iterable):
#         raise ParameterListNotIterableError("Param {} asked to take values " \
#                 "{} which should be an iterable list of values (len == 1 is " \
#                 "fine).".format(key, val))

class Scan:
    """Class designed for parameter scanning. Designed under the assumption
    that there are basically three types of "parameter sweeps" that need to be
    done.
    1) run the same parameters many times (e.g. stoch simulation)
    2) vary certain parameters jointly (e.g. (i,j) = (1,2), (2,3), (3,4), ... )
    3) vary parameters combinatorially (e.g. (i,j) = (1,1), (1,2), (2,1), (2,2))
    If a Scan is stopped midway, it will remember its state, allowing you to
    load it back up, change the scan parameters, and continue the run.
    A Scan is *not* useful if your code needs to regularly generate the next
    parameters based on previous simulations.

    If you have several "jointly" varying parameters, each "group" of them will
    interact combinatorially. More abstractly, you should think of each set of
    jointly varying parameters as being the same as a "single" parameter.

    e.g. Suppose we're running a simulation 2 times for each set of (a,b,c),
    except that if 'c' is too big, we run less repeats to save time:

    ```
    >>> import numpy as np
    >>> def f(a,b,c):
    >>>     print('a = ', a, ', b = ', b, ', c = c')
    >>> p['a'] = np.linspace(0, 10, 5)
    >>> p['b'] = np.linspace(0.1, -0.5, 5)
    >>> p['c'] = np.linspace(10, 20, 3)
    >>> default_count = lambda p: 2
    >>> # None specifies to default to previous value
    >>> big_c_count = lambda p: 1 if p['c'] >= 14 else None
    >>>
    >>> s = Scan.from_dict(p, joint_lists=[['a','b']])
    >>> s.add_count(default_count) # these will be called
    >>> s.add_count(big_c_count) # in the order they were added
    >>> s.run_scan(f) # verify scan validity and run

    ```

    will print out

    ```
        0.0 0.1 10.0
        0.0 0.1 10.0
        2.5 -0.05 10.0
        2.5 -0.05 10.0
        5.0 -0.2 10.0
        5.0 -0.2 10.0
        7.5 -0.35 10.0
        7.5 -0.35 10.0
        10.0 -0.5 10.0
        10.0 -0.5 10.0
        0.0 0.1 15.0
        2.5 -0.05 15.0
        5.0 -0.2 15.0
        7.5 -0.35 15.0
        10.0 -0.5 15.0
        0.0 0.1 20.0
        2.5 -0.05 20.0
        5.0 -0.2 20.0
        7.5 -0.35 20.0
        10.0 -0.5 20.0
    ```
    """

    def __init__(self, dic={}, joint_lists=[], default_repeats=1,
                 count_funcs=[]):
# how many times to repeat simulation of a specific parameter by default
        self.default_repeats = default_repeats
# list of functions to iteratively determine how many repeats to actually use
        self.count_funcs = count_funcs
# parameters that need to be combinatorially scanned
        self.comb_params = dict(dic)
# parameters that need to be varied together
        self.joint_params = []
        for jlist in joint_lists:
            # append sub-dictionary with keys from each jlist
            jparam = {key: dic[key] for key in jlist}
            # check self-consistency of joint params before insertion
            _check_joint_params(jparam)
            self.joint_params.append(jparam)
            for key in jlist:
                # remove jointly varying parameters from comb_params as needed
                # as opposed to only adding params to comb_params if they're not
                # in flatten(flatten((joint_lists))
                del self.comb_params[key]
        # make singleton variables iterable
        for key,val in self.comb_params.items():
            if not isinstance(val, collections.Iterable):
                self.comb_params[key] = [self.comb_params[key]]
                #_check_comb_param(key, val)

    def run_scan(self, f):
        """Run f the requested number of times for each set of parameters
        requested."""
        for params in self.params():
            f(**params)

    def params(self):
        """A generator that iterates through all parameters requested the
        correct number of times each. Returns them as a dict with form
        {'param_name': param_value, ... } for use as f(**params)."""
        comb_sizes = []
        comb_keys = []
        for key,val in self.comb_params.items():
            comb_keys.append(key)
            comb_sizes.append(len(val))
        # length of first values array, since lengths should match for sets of
        # jointly varying parameters
        joint_sizes = [len(next(iter(jlist.values()))) for jlist in self.joint_params]
        all_sizes = comb_sizes + joint_sizes
        total_combinations = reduce(operator.mul, all_sizes, 1)
        for i in range(total_combinations):
            sub = np.unravel_index(i, all_sizes)
            params = {}
            # get the comb_params
            for j in range(len(comb_sizes)):
                key = comb_keys[j]
                val_arr = self.comb_params[key]
                params[key] = val_arr[sub[j]]
            # get the joint_params
            for j in range(len(joint_sizes)):
                subj = j + len(comb_sizes)
                for key,val_arr in self.joint_params[j].items():
                    params[key] = val_arr[sub[subj]]
            # how we have one parameter set, check how many times to repeat it
            num_repeats = self.default_repeats
            for func in self.count_funcs:
                c = func(params)
                if c:
                    num_repeats = c
            for i in range(num_repeats):
                yield params

    def add_count(self, func):
        """Add (a) new function(s) to determine how many times to repeat a parameter
        set. Each function should take a dictionary of values (your parameters)
        and return an integer if it wants to determine how many times to run
        those parameters or None if it wants to defer to the default in that
        case.

        The functions will be called in the order that they are added.
        e.g.
        >>> s.Scan(p)
        >>> default_count = lambda p: 5
        >>> # None specifies to default to previous value
        >>> big_c_count = lambda p: 2 if p['c'] >= 14 else None
        >>> s.add_count(default_count)
        >>> s.add_count(big_c_count)

        """
        self.count_funcs.append(func)

    def add_params(self, p):
        """Add new variables or change values to be used for particular
        parameter names that will not be jointly varied."""
        new_params = p.copy()
        for key,val in new_params.items():
            if not isinstance(val, collections.Iterable):
                new_params[key] = [new_params[key]]
            #_check_comb_param(key, val)
        self.comb_params.update(new_params)

    def add_jparam(self, jparam):
        """Add a new set of parameters that should be varied together. Should
        simply be a dict of "param_name": param_values pairs, where
        all param_values are the same size."""
        _check_joint_params(jparam)
        self.joint_params.append(jparam)

    @classmethod
    def from_dict(cls, dic, joint_lists=[]):
        return cls(dic, joint_lists)

class TestScan:
    def test_error_if_param_val_not_iterable():
        pass
    def test_error_if_joint_size_wrong():
        pass
    def test_docs_example():
        ans = """0.0 0.1 10.0
0.0 0.1 10.0
2.5 -0.05 10.0
2.5 -0.05 10.0
5.0 -0.2 10.0
5.0 -0.2 10.0
7.5 -0.35 10.0
7.5 -0.35 10.0
10.0 -0.5 10.0
10.0 -0.5 10.0
0.0 0.1 15.0
2.5 -0.05 15.0
5.0 -0.2 15.0
7.5 -0.35 15.0
10.0 -0.5 15.0
0.0 0.1 20.0
2.5 -0.05 20.0
5.0 -0.2 20.0
7.5 -0.35 20.0
10.0 -0.5 20.0"""
        f = lambda a,b,c: print(a,b,c)
        p['a'] = np.linspace(0, 10, 5)
        p['b'] = np.linspace(0.1, -0.5, 5)
        p['c'] = np.linspace(10, 20, 3)
        default_count = lambda p: 2
        # None specifies to default to previous value
        big_c_count = lambda p: 1 if p['c'] >= 14 else None

        s = Scan.from_dict(p, joint_lists=[['a','b']])
        s.add_count(default_count) # these will be called
        s.add_count(big_c_count) # in the order they were added
        s.run_scan(f) # verify scan validity and run

