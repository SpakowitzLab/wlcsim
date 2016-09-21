#!/usr/bin/python

# convert an output file of snapshots to a pdb with multiple structures

from pdbutils import *
import sys, os

def lines2Struct(lines,branchscl=0.1,scl=10):
    """Parse lines containing info on chain bead position and orientation. Connect atoms appropriately. Return structure. All coordinates are scaled by scl."""

    beads = []; uvec = []; obstacles=[]
    chainnum = 0; newchain = 1;

    # use occupancies and b factors to track the weight of each bead
    # based on integrated r and u
    # given in the last 2 columns of the snapshot file
    occupancies = []
    bfactors = []
    for line in lines:
        data = line.split()
        if data[0]=='C': # new chain
            chainnum += 1
            newchain = 1
            continue

        coords = [float(x)*scl for x in data[1:7]]

        if (data[0] != 'A'):
            print 'Bad line, skipping:', line
        else:
            
            crd = coords[:3]
            atm = Atom(coords=crd,name=data[0])
            if (len(data)>13):
                tmp = [float(data[13]),float(data[14])]
                if (tmp[0]>0):
                    atm.occupancy = -log(tmp[0])
                if (tmp[1]>0):
                    atm.bfactor = -log(tmp[1])
#                occupancies.append(float(data[13]))
#                bfactors.append(float(data[14]))
                
            beads.append(atm)
            beads[-1].chain = "%d"%chainnum
            if (not newchain):
                beads[-2].conect.append(beads[-1])
                beads[-1].conect.append(beads[-2])

            # orientation vector atoms
            ucrd = list(array(crd) + array(coords[3:])*branchscl)
            uvec.append(Atom(coords=ucrd,name=data[0]+'U'))
            uvec[-1].conect.append(beads[-1])
            beads[-1].conect.append(uvec[-1])
            
            newchain = 0
      
    #if len(occupancies)>0:
        # occ = [x for x in occupancies if x>0]
        # if (len(occ)>0):            
        #     minocc= min(occ)
        #     for c in range(len(beads)):
        #         beads[c].occupancy = occupancies[c]/minocc
        
        # bfact = [x for x in bfactors if x>0]
        # if (len(bfact)>0):
        #      minbfact= min(bfact)
        #      for c in range(len(beads)):
        #          beads[c].bfactor = bfactors[c]/minbfact

    struct = makeBareStruct()
    struct.endlines=[]
#    struct.startlines = ['HEADER    struct\n']
    struct.startlines = []
    struct.atoms.extend(beads)
    struct.atoms.extend(uvec)
    struct.atoms.extend(obstacles)
    struct.renumAtoms()

    return struct

def runcode(argv):
	if len(argv) < 2:
            print """
usage: beadobst2pdb.py infile(s) -o <outfile>
All argument pairs besides the input file are optional. Can also supply a list or glob (eg: file.*.out) of input files. Input files may not start with "-"
-o Gives the output pdb file name. Default output file is: infileroot.pdb
-scl is the scaling factor to convert from .out to .pdb coordinates 
             (default 10, as the .out file is expected to be in nm while the .pdb file is in angstroms)
-branchscl is the factor for scaling branch lengths. Default is 0.1
-skip # skip the first few snapshots
"""
            sys.exit()

        # parse input arguments
        infile = sys.argv[1]                

	if '-o' in argv:
	    ind = argv.index('-o')
	    outfile = argv[ind+1]
	else:
            outfile = infile
            (root,ext) = os.path.splitext(infile)
            outfile = root+'.pdb'

        if '-branchscl' in argv:
            ind = argv.index('-branchscl')
            branchscl = float(argv[ind+1])
        else:
            branchscl = 0.1

        if '-scl' in argv:
            ind = argv.index('-scl')
            scl = float(argv[ind+1])
        else:
            scl = 10

        if '-skip' in argv:
            ind = argv.index('-skip')
            skip = int(argv[ind+1])
            print "Skipping first %d snapshots" %skip
        else:
            skip = 0

        #read in file and get structures
        structs = []
        IF = open(infile)
        lines = IF.readlines()
        infolines = [c for c in range(len(lines)) if lines[c][0]=='I']
        starting = 1;
        for c in range(skip,len(infolines)):
            data = lines[infolines[c]].split()
            nbead = int(data[1])

            start = infolines[c]+1
            if c==len(infolines)-1:
                end = len(lines)
            else:
                end = infolines[c+1]
            
            beadlines = lines[start:end]
#            if len(beadlines) != nbead:
#                print 'Inappropriate number of beads. Skipping struct', c
#            else:
            
            print 'Structure %d, with %d chains. %d total atoms.' %(c,nbead,len(beadlines))
            
            struct = lines2Struct(beadlines,branchscl,scl)   
            append = not (starting)
            starting = 0

            struct.outputPDB(outfile,append,c+1)
            #structs.append(st)

        IF.close()

        # for cf in range(len(infiles)):
        #     infile = infiles[cf]; outfile = outfiles[cf]
        #     print "Input file, output file:", infile, outfile
            
        #     struct = coordFile2Struct(infile,branchscl,scl)
        #     struct.outputPDB(outfile)

if __name__ == '__main__':
	runcode(sys.argv)

