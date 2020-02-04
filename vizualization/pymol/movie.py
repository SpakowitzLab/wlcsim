from pymol import cmd
from sys import argv

# this script is run automatically by NPrun.py and will need to be changed for someone else's use

numFrames = int(argv[1])

for idx in range(0,numFrames): cmd.load("pdb/snap%03d.pdb"%idx,"snap")
cmd.intra_fit("snap")
cmd.mset("1 -%d" % numFrames)
cmd.show('spheres', 'resn NUC')
cmd.spectrum('count', 'rainbow', 'resn NUC')
cmd.spectrum('count', 'rainbow', 'resn DNA')
cmd.set("sphere_scale", 1.5)
cmd.set("sphere_transparency", 0.2)
cmd.mplay()
