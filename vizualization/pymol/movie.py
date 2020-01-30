from pymol import cmd

numFrames = 110

for idx in range(0,numFrames): cmd.load("pdb/snap%03d.pdb"%idx,"snap")
cmd.intra_fit("snap")
cmd.mset("1 -%d" % numFrames)
cmd.show('spheres', 'resn NUC')
cmd.spectrum('count', 'rainbow', 'resn NUC')
cmd.spectrum('count', 'rainbow', 'resn DNA')
cmd.set("sphere_scale", 1.5)
cmd.set("sphere_transparency", 0.2)
cmd.mplay()
