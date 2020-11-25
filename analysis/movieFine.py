from pymol import cmd
from sys import argv

# this script is run automatically by MCsim and will need to be changed for someone else's use

numFrames = int(argv[1])
pathPDB = argv[2]

for idx in range(0,numFrames): cmd.load(pathPDB+"fine%03dv0.pdb"%idx,"snap")
cmd.intra_fit("snap")
cmd.mset("1 -%d" % numFrames)
cmd.show('spheres', 'resn DNA')
#cmd.show('lines', 'resn DNA')
cmd.hide('sticks', 'resn DNA')
cmd.spectrum('count', 'rainbow', 'resn DNA')
cmd.set("sphere_scale", 0.15)
cmd.mplay()
