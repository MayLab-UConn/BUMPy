''' Run this script to test every aspect of utils.rigid_body_transformations
    Run it from memcurv main folder.

    This is for IPYTHON!!!
'''
# Allows for on the fly modification of code
%load_ext autoreload
%autoreload 2

import numpy as np
import utils.rigid_body_transforms as rbtrans
import testing.gen_random_coordinates as gen
import testing.plot_3D_transforms as pt

# make a flat bilayer
flat_coords = gen.gen_random_coords(10000,100,100,5)
#(flat_base_fig,flat_base_ax) = pt.check_3D_structure(flat_coords)

# center flat bilayer
flat_centered = rbtrans.center_coordinates_3D(flat_coords)
rot_90x       = rbtrans.rotate_coordinates(flat_centered,[90,0,0])
rot_neg90x    = rbtrans.rotate_coordinates(flat_centered,[-90,0,0])
rot_45y       = rbtrans.rotate_coordinates(flat_centered,[0,45,0])
rot_yandz     = rbtrans.rotate_coordinates(flat_centered,[0,45,90])
