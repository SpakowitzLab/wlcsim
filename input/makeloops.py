import numpy as np
import bisect as bs
def insertLoop(loops,ends,limits=(0,101)):
    while True:
        # insert a point in a region with probability propto size of region
        pt = np.random.randint(limits[0],limits[1])
        left_n = bs.bisect_left(ends,pt) -1 # left index
        right_n = left_n+1 # right index
        if ends[right_n]-ends[left_n] > 2:
            break

    left_lim = ends[left_n]
    right_lim = ends[right_n]

    # Insert two points uniformally in region
    while True:
        low = np.random.randint(left_lim+1,right_lim)
        high = np.random.randint(left_lim+1,right_lim)
        if low != high: # make sure they aren't the same
            break

    if low>high: # order ends
        temp=high
        high = low
        low = temp

    # add to list of loops and ends
    loops.append((low,high))
    ends.insert(right_n,high)
    ends.insert(right_n,low)

def makeLoopFile(L,nloops,name=None,autoName=False,npoly=1):
    if not isinstance(L,int):
        raise ValueError('L must be an integer')
    if not isinstance(nloops,int):
        raise ValueError('nloops must be an integer')
    loops=[]
    #ends=[0, L+1]
    ends=[]
    for ipoly in range(0,npoly):
        ends.append(1+ipoly*L)
        ends.append(L+ipoly*L)
        for ii in range(0,nloops):
            insertLoop(loops,ends,limits=(1+ipoly*L,L+ipoly*L))
    toprint=np.zeros(L*npoly) - 1
    for ii in range(0,nloops*npoly):
        (a,b)=loops[ii]
        toprint[a-1]=b
        toprint[b-1]=a
    if autoName or name==None or name == '':
        if name==None:
            name = ''
        name="L"+str(L)+"nloops"+str(nloops)+name
    np.savetxt(name,toprint,fmt='%d')
    return toprint

#import sys
#nloops=50000
#L=393216
#npoly=1
#name = None
#if len(sys.argv) == 2:
#    nloops= int(sys.argv[1])
#elif len(sys.argv) == 3:
#    nloops = int(sys.argv[1])
#    L = int(sys.argv[2])
#elif len(sys.argv) == 4:
#    nloops = int(sys.argv[1])
#    L = int(sys.argv[2])
#    name = str(sys.argv[3])
#elif len(sys.argv) == 5:
#    nloops = int(sys.argv[1])
#    L = int(sys.argv[2])
#    npoly = int(sys.argv[3])
#    name = str(sys.argv[4])
#
#makeLoopFile(L,nloops,name=name,npoly=npoly)
