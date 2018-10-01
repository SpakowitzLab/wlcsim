import numpy as np
import loopSim
import makeloops

lengthSet=[]
with open("polyLengths") as lengthsFile:
    for line in lengthsFile:
        lengthSet.append(int(line))

looptype = "insert"
beadsPerLoop = 120/5

if looptype == "insert":
    for IP in range(0,len(lengthSet)):
        L = lengthSet[IP]
        nloops = int(np.floor(L/beadsPerLoop))
        name = "bindpairs_chrom"+str(IP)
        npoly = 1
        makeloops.makeLoopFile(L,nloops,name=name,npoly=npoly)
elif looptype == "growth":
    for IP in range(0,len(lengthSet)):
        L = lengthSet[IP]
        nloops = int(np.floor(L/beadsPerLoop))
        name = "bindpairs_chrom"+str(IP)
        chain = loopSim.LoopExtrusionChain(L,n_cohesins)

        forward_rate = 1.0
        reverse_rate = 0.0
        falloff_rate = 1.0/(120/5)
        delta_t = 300
        npts = 20
        distances = []
        for ii in range(npts):
            if ii%10==0:
                print("working on",ii," of ",npts)
            chain.run_for_time(delta_t, forward_rate, reverse_rate,
                               falloff_rate)
            distances.append(chain.total_spand_beads())

        np.savetxt("loop_growth_equilibration_"+str(IP),distances)
        chain.print_for_Monte_Carlo("bindpairs_chrom"+str(IP))






