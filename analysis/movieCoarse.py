from pymol import cmd
from sys import argv
import numpy as np
import colorsys
import sys
from pymol.cgo import *
sys.path.append('../analysis/')
from analysis.utility import *

# this script is run automatically by MCsim and will need to be changed for someone else's use

numFrames = int(argv[1])
channel = argv[2]
temperature = int(argv[3])
pathPDB = argv[4]
pathData = argv[5]
hull = argv[6]
rad = float(argv[7])
if (channel == 'PT'):
    nodes = np.loadtxt(pathData+'nodeNumber')
    nodes = np.vstack((np.linspace(0, np.shape(nodes)[1]-1, np.shape(nodes)[1]), nodes))
    channel = np.asarray(nodes[:,temperature], 'int')
else:
    channel = [channel]*numFrames
input_folder = pathData
side = 16
incr = 2*np.pi/(side/2.0)

for idx in range(0,numFrames): cmd.load(pathPDB+"coarse%03d.pdb"%idx,"snap")
cmd.mset("1 -%d" % numFrames)
cmd.color('gray80', 'snap')

file_inds = range(0,numFrames)
for ind in file_inds:
    #load file r and u data
    r = np.loadtxt('/%sr%sv%s' %(input_folder,ind,channel[ind]))
    u = np.loadtxt('/%su%sv%s' %(input_folder,ind,channel[ind]))
    # load discretization data
    disc = np.loadtxt('/%sd%sv%s' %(input_folder,ind,channel[ind]))
    wrap = disc[0]; bps = disc[1]
    if (hull == 'True'):
        for i in range(len(r)):
            if bps[i] != 0:
                poly = np.zeros(side*3).reshape([side,3])
                uin = np.asarray(u[i,0:3]); vin = np.asarray(u[i,3:6]); cross = np.cross(uin, vin)
                mat = np.matrix([vin, cross, uin]).reshape([3,3]).T
                if (wrap[i]>1):
                    space = 0
                    height = 5.5
                    radius = 5.2
                    center = np.asarray([4.8455, -2.4445, 0.6694])
                    for j in range(int(side/2.0)):
                        # rotate into material frame 
                        vec = np.asarray([radius*np.cos(space), -height/2.0, radius*np.sin(space)])
                        poly[j,:] = r[i,:] + np.matmul(mat, center+vec)
                        space = space + incr
                    space = 0
                    for j in range(int(side/2.0),side):
                        # rotate into material frame 
                        vec = np.asarray([radius*np.cos(space), height/2.0, radius*np.sin(space)])
                        poly[j,:] = r[i,:] + np.matmul(mat, center+vec)
                        space = space + incr
                    # make cylinders
                    x1,y1,z1 = np.mean(poly[:int(side/2),0]), np.mean(poly[:int(side/2),1]), np.mean(poly[:int(side/2),2])
                    (re, g, b) = colorsys.hsv_to_rgb(float(i)/(len(r)-1), 1.0, 1.0)
                    x2,y2,z2 = np.mean(poly[int(side/2):,0]), np.mean(poly[int(side/2):,1]), np.mean(poly[int(side/2):,2])
                    cmd.load_cgo( [ 25.0, 0.25, 9.0, x1, y1, z1, x2, y2, z2, radius, re, g, b, re, g, b ], "seg"+str(i+1)+'nuc')
                    # add extruding linker
                    space = 2*np.pi/side
                    tempU, tempV, tempR = rotate_bead(uin,vin,r[i,:],int(bps[i]),int(wrap[i]))
                    height = np.sqrt(np.dot(r[i+1,:]-tempR, r[i+1,:]-tempR))
                    radius = 1.0
                    center = np.zeros(3)
                    tempMat = np.matrix([tempV, np.cross(tempU, tempV), tempU]).reshape([3,3]).T
                    for j in range(int(side/2.0)):
                        # rotate into material frame 
                        vec = np.asarray([radius*np.sin(space), radius*np.cos(space), -height/2.0])
                        poly[j,:] = (tempR+r[i+1,:])/2.0 + np.matmul(tempMat, center+vec)
                        space = space + incr
                    space = 2*np.pi/side
                    for j in range(int(side/2.0),side):
                        # rotate into material frame 
                        vec = np.asarray([radius*np.sin(space), radius*np.cos(space), height/2.0])
                        poly[j,:] = (tempR+r[i+1,:])/2.0 + np.matmul(tempMat, center+vec)
                        space = space + incr
                else:
                    space = 2*np.pi/side
                    height = np.sqrt(np.dot(r[i+1,:]-r[i,:], r[i+1,:]-r[i,:]))
                    radius = 1.0
                    center = np.zeros(3)
                    for j in range(int(side/2.0)):
                        # rotate into material frame 
                        vec = np.asarray([radius*np.sin(space), radius*np.cos(space), -height/2.0])
                        poly[j,:] = (r[i,:]+r[i+1,:])/2.0 + np.matmul(mat, center+vec)
                        space = space + incr
                    space = 2*np.pi/side
                    for j in range(int(side/2.0),side):
                        # rotate into material frame 
                        vec = np.asarray([radius*np.sin(space), radius*np.cos(space), height/2.0])
                        poly[j,:] = (r[i,:]+r[i+1,:])/2.0 + np.matmul(mat, center+vec)
                        space = space + incr
                # make cylinders
                x1,y1,z1 = np.mean(poly[:int(side/2),0]), np.mean(poly[:int(side/2),1]), np.mean(poly[:int(side/2),2])
                (re, g, b) = colorsys.hsv_to_rgb(float(i)/(len(r)-1), 1.0, 1.0)
                x2,y2,z2 = np.mean(poly[int(side/2):,0]), np.mean(poly[int(side/2):,1]), np.mean(poly[int(side/2):,2])
                cmd.load_cgo( [ 25.0, 0.25, 9.0, x1, y1, z1, x2, y2, z2, radius, re, g, b, re, g, b ], "seg"+str(i+1))
#cmd.mplay()
cmd.orient('snap')

# add sphere volume if set 
if (rad>0):
   spherelist = [
      COLOR,    1,   0,    0,
      SPHERE,   rad, rad, rad, rad,
   ]
   cmd.load_cgo(spherelist, 'sphere',   1)
   cmd.set('cgo_transparency', 0.6)
   cmd.set('cgo_sphere_quality', 100)

