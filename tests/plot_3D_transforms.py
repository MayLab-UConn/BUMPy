# -*- coding: utf-8 -*-

from matplotlib import pyplot
import pylab
from mpl_toolkits.mplot3d import Axes3D

def check_3D_structure(coords):
    fig = pylab.figure()
    ax = Axes3D(fig)
    ax.scatter(coords[:,0],coords[:,1],coords[:,2])
    pyplot.show()
    return fig