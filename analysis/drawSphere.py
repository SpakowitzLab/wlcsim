from pymol.cgo import *
from pymol import cmd


spherelist = [
   COLOR,    1,   0,    0,
   SPHERE,   250/2, 250/2, 250/2, 250/2,
    ]

cmd.load_cgo(spherelist, 'sphere',   1)
cmd.set('cgo_transparency', 0.6)
cmd.set('cgo_sphere_quality', 100)
