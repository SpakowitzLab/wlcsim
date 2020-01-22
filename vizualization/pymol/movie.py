from pymol import cmd

numFrames = 110

for idx in range(0,numFrames): cmd.load("pdb/snap%03d.pdb"%idx,"mov")
cmd.intra_fit("mov")
cmd.mset("1 -%d" % numFrames)
cmd.mplay()
