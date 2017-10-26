#  ----------------------------------------
#
#    Parses output file and makes pdb
#    Written by Quin MacPherson
#    Started: 11/11/16
#
#   Example usage:
#   for i in `seq 1 100`;
#       do for j in `seq 1 12`;
#           do python r2pdb.py r[$i]v[$j] r[$i]v[$j]
#       done
#   done
#-------------------------------------------

import sys
import numpy as np

if len(sys.argv) != 3:
      sys.exit('two argument please')
fname= sys.argv[1]
outname= sys.argv[2]
#fname='full/r109v1'


# ---------------------
#    Read from file
# --------------------

data = np.loadtxt(fname,np.float)

R = data[:,0:3].astype(np.float32)
extra = data[:,3:].astype(np.bool_)

np.savez_compressed(outname,a=R,b=extra)



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
