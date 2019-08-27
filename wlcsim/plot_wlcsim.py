"""For now, just a widget that requires hardcoding everythign about hte
simulation. Eventually, a tool to load in simulation data from pysim or
wlcsim and vizualize it intelligently."""
import numpy as np
import matplotlib
matplotlib.use("Qt4Agg") # This program works with Qt only
from matplotlib.widgets import Slider, Button, RadioButtons
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from PyQt5 import QtCore, QtWidgets
from PyQt5.QtCore import Qt
import os
import re
from .data import Sim


class SimulationViewer(object):

    def __init__(self, x=None):
        # will keep track of which *index* in the time axis should be used to
        # draw the current scene
        self.t = 0
        self.sim = None
        # define the actual axes layout
        self.fig = plt.figure()
        # area to plot polymers
        self.ax3d = plt.axes([0.05, 0.15, 0.9, 0.8], facecolor='w', projection='3d')
        # area to plot "slider" for selecting what time to plot
        self.ax = plt.axes([0.1, 0.025, 0.8, 0.05], facecolor='lightgoldenrodyellow')
        # set up the QtGui panel and textbox
        self.root = self.fig.canvas.manager.window
        self.panel = QtWidgets.QWidget()
        self.hbox = QtWidgets.QHBoxLayout(self.panel)
        self.textbox = QtWidgets.QLineEdit(parent=self.panel)
        self.hbox.addWidget(self.textbox)
        self.panel.setLayout(self.hbox)
        self.dock = QtWidgets.QDockWidget("Simulation Directory", self.root)
        self.root.addDockWidget(Qt.BottomDockWidgetArea, self.dock)
        self.dock.setWidget(self.panel)
        # handler to check if new simulation directory has been entered
        self.textbox.textChanged.connect(self.update_dir)
        # slider to control what "time" is plotted, (0,1) is rescaled total
        # simulation time units
        self.slider_t = Slider(self.ax, 't', 0, 1, valinit=0)
        # handler to check if the "time" slider has been moved
        self.slider_t.on_changed(self.update_t)
        plt.show()

    # should never be called before a Sim has been loaded successfully
    def update_drawing(self):
        num_polymers = self.sim.num_polymers
        if num_polymers != len(self.ax3d.lines):
            self.ax3d.lines = []
            for i in range(num_polymers):
                self.ax3d.plot(self.sim.r[self.t,0,i,:],
                               self.sim.r[self.t,1,i,:],
                               self.sim.r[self.t,2,i,:])
        for i,line in enumerate(self.ax3d.lines):
            x, y, z = (self.sim.r[self.t,0,i,:], self.sim.r[self.t,1,i,:],
                    self.sim.r[self.t,2,i,:])
            line.set_data(x, y)
            line.set_3d_properties(z)
        self.fig.canvas.draw_idle()

    def update_ax_limits(self):
        self.ax3d.set_xlim((np.nanmin(self.sim.r[:,0,:,:]), np.nanmax(self.sim.r[:,0,:,:])))
        self.ax3d.set_ylim((np.nanmin(self.sim.r[:,1,:,:]), np.nanmax(self.sim.r[:,1,:,:])))
        self.ax3d.set_zlim((np.nanmin(self.sim.r[:,2,:,:]), np.nanmax(self.sim.r[:,2,:,:])))

    def update_dir(self):
        text = self.textbox.text()
        try:
            self.sim = Sim(text)
        except FileNotFoundError:
            return
        self.update_t(self.slider_t.val)
        self.update_ax_limits()
        self.update_drawing()

    def update_t(self, val):
        if self.sim is None:
            return
        # calculate value on slider
        new_t = int(round(val*self.sim.num_time_points))
        # bound it to 0, num_time_points-1
        new_t = max(0, min(self.sim.num_time_points-1, new_t))
        self.t = new_t
        self.update_drawing()

if __name__ == '__main__':
    sim_viewer = SimulationViewer()

