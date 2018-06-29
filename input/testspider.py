import numpy as np
import os

length =200
bindpairs = np.zeros(length).astype(int)-1
def makeloop(a,b,bindpairs):
    bindpairs[a]=b
    bindpairs[b]=a


makeloop(30,70,bindpairs)
makeloop(40,60,bindpairs)


#np.savetxt('testPairs',bindpairs,fmt='%d')
#os.system('python spiderMove.py testPairs spiders_test')

import spiderMove as sp

spiders = sp.makeSpiders(bindpairs,leglength=5,moveable_loop=25)

print("+++++++++++++++++++")
for spider in spiders:
    sp.print_sections(spider)
    print("-------------------")
