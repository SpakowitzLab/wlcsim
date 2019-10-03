from pymol import cmd

for idx in range(1,110): cmd.load("pdb/test/snap%03d.pdb"%idx,"mov")
cmd.intra_fit("mov")
cmd.mset("1 -%d" % 110)
cmd.mplay()
