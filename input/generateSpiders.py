from spiderMove import *

if len(sys.argv) ==2:
    fname= sys.argv[1]
    outName=None
elif len(sys.argv) ==3:
    fname= sys.argv[1]
    outName=sys.argv[2]

bindpairs = np.loadtxt(fname).astype(int)
bindpairs[bindpairs != -1] = bindpairs[bindpairs != -1] - 1
spiders = makeSpiders(bindpairs,leglength=5,moveable_loop=25)
print_spiders_for_fortran(spiders,filename=outName,offset=1)
