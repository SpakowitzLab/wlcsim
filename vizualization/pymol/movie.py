from pymol import cmd

numFrames = 110

for idx in range(0,numFrames): cmd.load("pdb/nuc%03d.pdb"%idx,"nuc")
cmd.intra_fit("nuc")
cmd.mset("1 -%d" % numFrames)
cmd.mplay()

for idx in range(0,numFrames): cmd.load("pdb/dna%03d.pdb"%idx,"dna")
cmd.intra_fit("dna")
cmd.mset("1 -%d" % numFrames)
cmd.hide('sticks', 'nuc')
cmd.show('spheres', 'nuc')
cmd.mplay()
