#  ----------------------------------------
#
#    Parses output file and makes pdb
#    Written by Quin MacPherson
#    Started: 11/11/16
#
#   Example usage:
#   for i in `seq 1 100`;
#       do for j in `seq 1 12`;
#           do python expand_r.py r[$i]v[$j].npz > r[$i]v[$j]
#       done
#   done
#-------------------------------------------

import sys
import numpy as np

if len(sys.argv) != 2:
      sys.exit('one argument please')
fname= sys.argv[1]
#fname='full/r109v1'


# ---------------------
#    Read from file
# --------------------

#data = np.loadtxt(fname,np.float)
#R = data[:,0:3].astype(np.float16)
#extra = data[:,3:].astype(np.bool_)
#np.savez_compressed(outname,R=R,extra=extra)

with np.load(sys.argv[1]) as data:
    R=data['a']
    extra=data['b']
N=R.shape[0]


nextra=extra.shape[1]

extra=extra.astype(np.int)

if nextra is 0:
    for ii in range(0,N):
        print('%10.3f%10.3f%10.3f' % (
              R[ii,0],R[ii,1],R[ii,2]))
if nextra is 1:
    for ii in range(0,N):
        print('%10.3f%10.3f%10.3f%d2' % (
              R[ii,0],R[ii,1],R[ii,2],extra[ii,0]))
if nextra is 2:
    for ii in range(0,N):
        print('%10.3f%10.3f%10.3f%2d %2d' % (
              R[ii,0],R[ii,1],R[ii,2],extra[ii,0],extra[ii,1]))
if nextra is 3:
    for ii in range(0,N):
        print('%10.3f%10.3f%10.3f%d2 %d2 %d2' %(
              R[ii,0],R[ii,1],R[ii,2],extra[ii,0],extra[ii,1],extra[ii,2]))


#X=[];
#Y=[];
#Z=[];
#AB=[];
#METH=[];
#index=[];
#n=0;
#count=0;
#skip=1
#with open(fname) as f:
#   for line in f:
#       count=count+1
#       if count%skip != 0:
#           continue
#       temp=line.split()
#       x=float(temp[0])
#       y=float(temp[1]) 
#       z=float(temp[2])
#       
#       X.append(x)  
#       Y.append(y) 
#       Z.append(z)
#       if len(temp)==3:
#           continue  
#       AB1.append(int(temp[3]))
#       if len(temp)==4:
#           continue  
#       AB2.append(int(temp[4]))
#       if len(temp)==5:
#           continue  
#       AB3.append(int(temp[5]))
#
#nbeads=len(X)
#
#
## check to make sure lines are the same length
#
