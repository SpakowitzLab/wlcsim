"""For now, just a widget that requires hardcoding everythign about hte
simulation. Eventually, a tool to load in simulation data from pysim or
wlcsim and vizualize it intelligently."""
import numpy as np
import matplotlib
matplotlib.use("Qt4Agg") # This program works with Qt only
from matplotlib.widgets import Slider, Button, RadioButtons
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from PyQt4 import QtGui
from PyQt4 import QtCore
from PyQt4.QtCore import Qt
import os
import re

num_beads = 10
num_time_points = 100 + 1 # + 1 for init cond from MC
ndim = 6

coords = np.zeros((ndim, num_beads, num_time_points))
# define the actual axes layout
fig = plt.figure()
ax3d = plt.axes([0.05, 0.15, 0.9, 0.8], axisbg='w', projection='3d')
ax = plt.axes([0.1, 0.025, 0.8, 0.05], axisbg='lightgoldenrodyellow')

t = 0
line, = ax3d.plot(coords[0,:,t], coords[1,:,t], coords[2,:,t])
def update_drawing():
    x, y, z = coords[0,:,t], coords[1,:,t], coords[2,:,t]
    line.set_data(x, y)
    line.set_3d_properties(z)
    fig.canvas.draw_idle()
update_drawing()

root = fig.canvas.manager.window
panel = QtGui.QWidget()
hbox = QtGui.QHBoxLayout(panel)
textbox = QtGui.QLineEdit(parent = panel)
def update_coords():
    global coords
    global num_time_points
    global num_beads
    num_time_points = 1
    for fname in os.listdir(data_dir):
        nums = re.findall(r'r(\d+)', fname)
        if not nums:
            continue
        ti = int(nums[-1])
        if ti > num_time_points:
            num_time_points = ti
    num_time_points += 1 # account for "zero" data point, i.e. 1-indexing
    # there should always be an initial time point
    r0 = np.loadtxt("{}/r{}".format(data_dir, 0))
    # use the initial time point to infer number of beads and dimension of sim
    num_beads, rdims = np.shape(r0)
    coords = np.zeros((ndim, num_beads, num_time_points))
    for ti in range(num_time_points):
        coords[0:3,:,ti] = np.loadtxt("{}/r{}".format(data_dir, ti)).T
        # center of mass coords
        # for dim in range(3):
        #     coords[dim,:,ti] = coords[dim,:,ti] - np.mean(coords[dim,:,ti])
    ax3d.set_xlim((np.min(coords[0,:,:]), np.max(coords[0,:,:])))
    ax3d.set_ylim((np.min(coords[1,:,:]), np.max(coords[1,:,:])))
    ax3d.set_zlim((np.min(coords[2,:,:]), np.max(coords[2,:,:])))

def update_dir():
    global data_dir
    text = textbox.text()
    if os.path.isdir(text):
        data_dir = textbox.text()
        update_coords()
        update_drawing()
textbox.textChanged.connect(update_dir)
hbox.addWidget(textbox)
panel.setLayout(hbox)

dock = QtGui.QDockWidget("Data Directory", root)
root.addDockWidget(Qt.BottomDockWidgetArea, dock)
dock.setWidget(panel)

slider_t = Slider(ax, 't', 0, num_time_points-1, valinit=0)
def update_t(val):
    global t
    t = int(round(slider_t.val))
    update_drawing()
slider_t.on_changed(update_t)

plt.show()

# ###### EX sliders
# fig, ax = plt.subplots()
# plt.subplots_adjust(left=0.25, bottom=0.25)
# t = np.arange(0.0, 1.0, 0.001)
# a0 = 5
# f0 = 3
# s = a0*np.sin(2*np.pi*f0*t)
# l, = plt.plot(t, s, lw=2, color='red')
# plt.axis([0, 1, -10, 10])

# axcolor = 'lightgoldenrodyellow'
# axfreq = plt.axes([0.25, 0.1, 0.65, 0.03], axisbg=axcolor)
# axamp = plt.axes([0.25, 0.15, 0.65, 0.03], axisbg=axcolor)

# resetax = plt.axes([0.8, 0.025, 0.1, 0.04])
# button = Button(resetax, 'Reset', color=axcolor, hovercolor='0.975')


# def reset(event):
#     sfreq.reset()
#     samp.reset()
# button.on_clicked(reset)

# rax = plt.axes([0.025, 0.5, 0.15, 0.15], axisbg=axcolor)
# radio = RadioButtons(rax, ('red', 'blue', 'green'), active=0)


# def colorfunc(label):
#     l.set_color(label)
#     fig.canvas.draw_idle()
# radio.on_clicked(colorfunc)

# plt.show()

