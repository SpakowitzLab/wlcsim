"""Module to plot polymers, cylinders, and spheres."""
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
import pandas as pd
from .utils.utils import well_behaved_decorator, make_decorator_factory_args_optional
from scipy import stats
from matplotlib.widgets import Slider, Button, RadioButtons
from PyQt5 import QtCore, QtWidgets
from PyQt5.QtCore import Qt


def make_standard_axes(is_3D=False):
    if is_3D:
        return plt.figure().add_subplot(111, projection='3d')
    else:
        return plt.figure().add_subplot(111)

@well_behaved_decorator(has_params=True)
@make_decorator_factory_args_optional
def make_axes_if_blank(is_3D=False):
    """Decorator for plotting function that take an "axes" kwarg, which will
    create an axes for them to use if one isn't passed to them."""
    def wrap(plot_function):
        """make_axes_if_blank's function that wraps plot_function."""
        def wrapped_f(*args, **kwargs):
            """The plot_function wrapped by make_axes_if_blank's wrap."""
            if 'axes' not in kwargs:
                kwargs['axes'] = make_standard_axes(is_3D)
            return plot_function(*args, **kwargs)
        return wrapped_f # return decorated function
    return wrap

@make_axes_if_blank(is_3D=False)
def testplot2(x, y, **kwargs):
    """Testing2 make_axes_if_blank decorator."""
    plt.plot(x, y, axes=kwargs['axes'])

@make_axes_if_blank
def testplot3(x, y, **kwargs):
    """Testing3 make_axes_if_blank decorator."""
    plt.plot(x, y, axes=kwargs['axes'])

@make_axes_if_blank(is_3D=True)
def testplot4(x, y, z, **kwargs):
    """Testing4 make_axes_if_blank decorator."""
    plt.plot(x, y, axes=kwargs['axes'])

@make_axes_if_blank(is_3D=True)
def draw_sphere(x0, r, **kwargs):
    """Draw a 3D sphere with center x0 and radius r."""
    # how densely gridding should happen
    longitudes = kwargs.pop('longitudes', 20)
    longitude_count = longitudes*1j
    latitude_count = longitude_count/2
    u, v = np.mgrid[0:2*np.pi:longitude_count, 0:np.pi:latitude_count]
    x = x0[0] + r*np.cos(u)*np.sin(v)
    y = x0[1] + r*np.sin(u)*np.sin(v)
    z = x0[2] + r*np.cos(v)
    ax = kwargs.pop('axes', None)
    ax.plot_wireframe(x, y, z, **kwargs)
    return ax

@make_axes_if_blank(is_3D=False)
def locally_linear_fit(x, y, window_size=5, **kwargs):
    if window_size % 2 == 0:
        raise ValueError('window_size must be odd')
    if window_size < 3:
        raise ValueError('window_size must be at least 3, so that we have at'
                         ' least two points to fit at endpoints of array.')
    lenx = len(x)
    if lenx < 2:
        raise ValueError('Can\'t fit less than 2 points!')
    inc = int(window_size / 2)
    slope = np.zeros_like(x)
    for i in range(lenx):
        imin = np.max([0, i - inc])
        imax = np.min([lenx, i + inc + 1])
        slope[i], intcpt, r_val, p_val, std_err = stats.linregress(x[imin:imax], y[imin:imax])
    ax = kwargs.pop('axes', None)
    if ax is not None:
        ax.plot(x, slope, **kwargs)
    return ax, slope

class PolymerViewer(object):

    def __init__(self, r=None):
        # will keep track of which *index* in the time axis should be used to
        # draw the current scene
        self.t = 0
        self.r = r
        # define the actual axes layout
        self.fig = plt.figure()
        # area to plot polymers
        self.ax3d = plt.axes([0.05, 0.15, 0.9, 0.8], facecolor='w', projection='3d')
        # area to plot "slider" for selecting what time to plot
        self.ax = plt.axes([0.1, 0.025, 0.8, 0.05], facecolor='lightgoldenrodyellow')
        # # set up the QtGui panel and textbox
        # self.root = self.fig.canvas.manager.window
        # self.panel = QtWidgets.QWidget()
        # self.hbox = QtWidgets.QHBoxLayout(self.panel)
        # self.textbox = QtWidgets.QLineEdit(parent=self.panel)
        # self.hbox.addWidget(self.textbox)
        # self.panel.setLayout(self.hbox)
        # self.dock = QtWidgets.QDockWidget("Simulation Directory", self.root)
        # self.root.addDockWidget(Qt.BottomDockWidgetArea, self.dock)
        # self.dock.setWidget(self.panel)
        self.update_ax_limits()
        # slider to control what "time" is plotted, (0,1) is rescaled total
        # simulation time units
        self.slider_t = Slider(self.ax, 't', 0, 1, valinit=0)
        # handler to check if the "time" slider has been moved
        self.slider_t.on_changed(self.update_t)
        plt.show()

    @property
    def r(self):
        return self._r

    @r.setter
    def r(self, new_r):
        self._r = new_r
        self._num_time_points, self._num_beads, self._d = self._r.shape
        if hasattr(self, 'ax3d'):
            self.update_ax_limits()


    # should never be called before a Sim has been loaded successfully
    def update_drawing(self):
        num_polymers = 1
        if num_polymers != len(self.ax3d.lines):
            self.ax3d.lines = []
            for i in range(num_polymers):
                self.ax3d.plot(self.r[self.t,:,0],
                               self.r[self.t,:,1],
                               self.r[self.t,:,2])
        for i,line in enumerate(self.ax3d.lines):
            x, y, z = (self.r[self.t,:,0], self.r[self.t,:,1],
                       self.r[self.t,:,2])
            line.set_data(x, y)
            line.set_3d_properties(z)
        self.fig.canvas.draw_idle()

    def update_ax_limits(self):
        self.ax3d.set_xlim((np.nanmin(self.r[:,:,0]), np.nanmax(self.r[:,:,0])))
        self.ax3d.set_ylim((np.nanmin(self.r[:,:,1]), np.nanmax(self.r[:,:,1])))
        self.ax3d.set_zlim((np.nanmin(self.r[:,:,2]), np.nanmax(self.r[:,:,2])))

    def update_t(self, val):
        if self.r is None:
            return
        # calculate value on slider
        new_t = int(round(val*self._num_time_points))
        # bound it to 0, num_time_points-1
        new_t = max(0, min(self._num_time_points-1, new_t))
        self.t = new_t
        self.update_drawing()
