.. _output:

About Simulation Output
#######################



Bead positions and identity
===========================

Periodically while running the simulation will take a save point and save the configuration of the polymer(s).  These are saved in the `data/` folder in files such as `data/r100` where the number 100 denotes what savepoint it is.  When using MPI to parallel temper, the data will be saved in the format `data/r110v6` where 100 refers to the save point and 6 refers to the replica.  This data will have three columns for the x, y, and z positions of each bead.  If there are multiple polymers they will be printed one after the next with no delineation.

If `WLC_P__SAVEAB` is specified as true than `data/r100v6` will contain a fourth
column specifying the bead type.  For two tail methylation profiles a value of 0
in this column refers to either tail bound by HP1, a value of 1 or 2 refers to
one or the other tail bound by HP1, while a value of 3 refers to both tails
being bound by HP1.  If `WLC_P__CHANGINGCHEMICALIDENTITY` is true than a fifth
column will specify the number of methylated tails (0, 1 or 2).

For example
::

    21.702    39.628    14.415  0  0
    21.715    39.617    14.892  3  1
    21.205    39.449    14.758  3  2
    20.718    39.472    14.672  3  2
    20.463    39.763    14.803  0  0
    20.209    39.667    15.418  1  0
    19.888    39.643    14.826  2  0
    ...


Other files 
===========

If `WLC_P__SAVEU` is true then the orientation unit vectors for each bead will be stored as in files of a similar format to the positions.  For example `data/u100` or `data/u100v6`.   If `WLC_P__LOCAL_TWIST` is also true then there will be three additional columns in the `data/u100`.

The file `data/energiesv6` contains the various energies and corresponding coefficients.  All energies are in kT’s.  There is one line in this file per savepoint.  The file `data/adaptationsv6` contains the adapted window size, success rate, and so on the various move types.  These files are useful for diagnosing problems and determining if the simulation is equilibrated.

A file `data/repHistory` contains various information about replica coupling if it is in use.

Visualization
=============

We provide visualization code based on pymol.  To run, navigate to ‘visualization/polymer_mole’, enter the savepoint and other settings into `run_visualization.py` and then type `python3 run_visualization.py`.  This will require both python3 and pymol.  The way it works is that the output described above is converted into a pdb file which is then passed to pymol.  Various sets of settings that have been used in the past are included in `run_visualization.py`.  To use one of thes set `image` to the name of the corrisponing set.

If `closePymol` is set to true in `run_visualization` than the output will be
saved to the `data/` directiory. By editing the line `for rep in ...` you can
save many snapshots.
