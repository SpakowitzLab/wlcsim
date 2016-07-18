delete *
reset

cd pdb
movie.load snap*.pdb,snap
cd ../pdb0
movie.load snap*.pdb,snap0

rotate y,90,state=0
rotate z,90,state=0
rotate x,-90,state=0
rotate z,90,state=0

hide all

set stick_radius=0.0075
show spheres, (name A1)
#color blue, (name A1)
alter (name A1), vdw=0.017
show sticks, (name A1)
hide spheres, snap0
color black, snap0

rebuild

run ../color_b.py
color_b snap,mode=hist,gradient=bgr,minimum=0,maximum=1,nbins=40,sat=1.,value=1.

bg_color white
reset

cd ..

#cd ../png
#set ray_trace_frames = 1
#mpng snap