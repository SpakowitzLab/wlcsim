import numpy as np
import loopSim
import makeloops

lengthSet=[]
with open("polyLengths") as lengthsFile:
    for line in lengthsFile:
        lengthSet.append(int(line))

looptype = "growth"
beadsPerLoop = 500
processivity = 500.0

if looptype == "insert":
    for IP in range(0,len(lengthSet)):
        L = lengthSet[IP]
        nloops = int(np.floor(L/beadsPerLoop))
        name = "bindpairs_chrom"+str(IP)
        npoly = 1
        makeloops.makeLoopFile(L,nloops,name=name,npoly=npoly)
elif looptype == "growth":
    for IP in range(0,len(lengthSet)):
        print("working on polymer "+str(IP))
        L = lengthSet[IP]
        nloops = int(np.floor(L/beadsPerLoop))
        name = "bindpairs_chrom"+str(IP)
        chain = loopSim.LoopExtrusionChain(L,nloops)

        forward_rate = 1.0
        reverse_rate = 0.0
        falloff_rate = (forward_rate-reverse_rate)/processivity
        delta_t = 2.0/falloff_rate
        npts = 10
        distances = []
        for ii in range(npts):
            if ii%2==0:
                print("working on",ii," of ",npts)
            chain.run_for_time(delta_t, forward_rate, reverse_rate,
                               falloff_rate)
            distances.append(chain.total_spand_beads())

        np.savetxt("loop_growth_equilibration_"+str(IP),distances)
        chain.print_for_Monte_Carlo("bindpairs_chrom"+str(IP))






