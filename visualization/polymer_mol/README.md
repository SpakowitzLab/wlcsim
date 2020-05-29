# polymer_mole
Use Pymol to generate many snapshots of polymer configurations, emphasis on chromosomal images.

## Dependancies
Python
Pymol
If you are ssh'ing you may need to have the ability to x forward.

## General Workflow
As seen in example.py, the idea is to programatically do the following for many
snapshots:
1) Use polymerMole.mainDefs.r2bdb to turn a file of x,y,z coordinates into a .pdb file.
2) Create a pml file that specifies viewing options with
polymerMole.mainDefs.makepmlFile
3) Run pymol and save the output to file

## Examples
Run example.py to see example of how the code works.  This example code puts the
output into the floder where the data came from.

Change the path and snapshot number in baseNames to the file you would like to
plot.

Change the plotting details after #Example1/2 to those that you would like.

Most likely you will want to added new options to the code in mainsDefs.r2pdb
and mainDefs.makepmlFile to your likeing.

Have fun!
