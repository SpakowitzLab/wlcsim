from spiderMove import *

if len(sys.argv) ==2:
    fname= sys.argv[1]
    outName=None
elif len(sys.argv) ==3:
    fname= sys.argv[1]
    outName=sys.argv[2]
elif len(sys.argv) ==4:
    fname= sys.argv[1]
    outName=sys.argv[2]
    leglen=int(sys.argv[3])


bindpairs = np.loadtxt(fname).astype(int)
bindpairs[bindpairs != -1] = bindpairs[bindpairs != -1] - 1
spiders = makeSpiders(bindpairs,leglength=leglen, moveable_loop=leglen*5)
print_spiders_for_fortran(spiders,filename=outName,offset=1)
