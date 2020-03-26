from pymol import cmd
from sys import argv

# this script is run automatically by NPrun.py and will need to be changed for someone else's use

numFrames = int(argv[1])

for idx in range(0,numFrames): cmd.load("pdb/fine%03d.pdb"%idx,"snap")
cmd.intra_fit("snap")
cmd.mset("1 -%d" % numFrames)
cmd.show('spheres', 'resn DNA')
cmd.show('lines', 'resn DNA')
cmd.hide('sticks', 'resn DNA')
cmd.spectrum('count', 'rainbow', 'resn DNA')
cmd.set("sphere_scale", 0.15)
cmd.mplay()
