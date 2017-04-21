# -*- coding: utf-8 -*-

from matplotlib import pyplot
import pylab
from mpl_toolkits.mplot3d import Axes3D
import numpy as np

def check_3D_structure(coords, show=True, axeq=True):
    ''' automated function for checking 3D structures
        Can suppress image with show=false,
        axes should be equal sized unless axeq=False
    '''

    fig = pylab.figure()
    ax = Axes3D(fig)
    ax.scatter(coords[:,0],coords[:,1],coords[:,2])
    if axeq:
        plotmax = np.max(coords)
        plotmin = np.min(coords)
        ax.set_xlim(plotmin,plotmax)
        ax.set_ylim(plotmin,plotmax)
        ax.set_zlim(plotmin,plotmax)
    if show:
        pyplot.show()
    return fig,ax
