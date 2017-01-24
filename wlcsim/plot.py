"""For now, just a widget that requires hardcoding everythign about hte
simulation. Eventually, a tool to load in simulation data from pysim or
wlcsim and vizualize it intelligently."""
from utils import math as wlcmath
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
num_time_points = 10
ndim = 3

coords = np.zeros((num_beads, ndim, num_time_points))
# define the actual axes layout
fig = plt.figure()
ax3d = plt.axes([0.05, 0.15, 0.9, 0.8], axisbg='w', projection='3d')
ax = plt.axes([0.1, 0.025, 0.8, 0.05], axisbg='lightgoldenrodyellow')

t = 0
lines = ax3d.plot(coords[:,0,t], coords[:,1,t], coords[:,2,t])
def update_drawing():
    global lines
    lines.pop(0).remove()
    lines = ax3d.plot(coords[:,0,t], coords[:,1,t], coords[:,2,t])
    fig.canvas.draw_idle()
update_drawing()

root = fig.canvas.manager.window
panel = QtGui.QWidget()
hbox = QtGui.QHBoxLayout(panel)
textbox = QtGui.QLineEdit(parent = panel)

def update_coords():
    global num_beads
    global num_time_points
    global ndim
    global coords
    x = np.fromfile(data_file + '.out')
    shape = [int(s) for s in np.loadtxt(data_file + '.shape')]
    num_beads = shape[0]
    ndim = shape[1]
    num_time_points = shape[2]
    coords = x.reshape((num_beads, ndim, num_time_points))
    coords = wlcmath.center_by_mass(coords, particle_axis=0)
    ax3d.set_xlim(np.min(coords[:,0,:]), np.max(coords[:,0,:]))
    ax3d.set_ylim(np.min(coords[:,1,:]), np.max(coords[:,1,:]))
    ax3d.set_zlim(np.min(coords[:,2,:]), np.max(coords[:,2,:]))


def update_file():
    global data_file
    text = textbox.text()
    if os.path.isfile(text + '.out') and os.path.isfile(text + '.shape'):
        data_file = textbox.text()
        update_coords()
        update_drawing()
textbox.textChanged.connect(update_file)
hbox.addWidget(textbox)
panel.setLayout(hbox)

dock = QtGui.QDockWidget("Data Directory", root)
root.addDockWidget(Qt.BottomDockWidgetArea, dock)
dock.setWidget(panel)

slider_t = Slider(ax, 't', 0, 1, valinit=0)
def update_t(val):
    global t
    t = int(round((num_time_points-1)*slider_t.val))
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

