delete *
reset


movie.load pdb/test/snap*.pdb,snap

hide all

set stick_radius=1.0
show spheres, (name A1)
show sphere, (name A2)
color blue, (name A1)
color red, (name A2)
alter (name A1), vdw=1.0
alter (name A2), vdw = 3.0
show sticks, (name A1)
show sticks, (name A2)
rebuild

#run color_b.py
#color_b snap,mode=hist,gradient=bgr,minimum=0,maximum=1,nbins=40,sat=1.,value=1.

bg_color white
reset

zoom center 200
#cd ../png
#set ray_trace_frames = 1
#mpng snap