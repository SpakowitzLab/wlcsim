from pymol.cgo import *
from pymol import cmd
#from sys import argv

#x = int(argv[1])
#y = int(argv[2])
#z = int(argv[3])
#d = int(argv[4])
#print(argv[1])
spherelist = [
   COLOR,    1,   0,    0,
   SPHERE,   50,   50,   50, 250.0/2,
    ]

cmd.load_cgo(spherelist, 'sphere',   1)
cmd.set('cgo_transparency', 0.6)
cmd.set('cgo_sphere_quality', 100)
