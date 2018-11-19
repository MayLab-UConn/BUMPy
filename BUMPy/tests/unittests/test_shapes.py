import unittest
import numpy as np
import sys
import os

sys.path.insert(0, os.path.abspath("../.."))   # hacky way to get access to bumpy.py

from bumpy import Molecules, shapes


# -----------------------------------------------------
# interface for all other shapes
# -----------------------------------------------------
class test_base_shape(unittest.TestCase):
    pass


# -----------------------------------------------------
# basic building blocks
# -----------------------------------------------------
class test_flat_bilayer(unittest.TestCase):
    pass


class test_semisphere(unittest.TestCase):
    pass


class test_cylinder(unittest.TestCase):
    pass


class test_inner_quarter_torus(unittest.TestCase):
    pass


class test_outer_quarter_torus(unittest.TestCase):
    pass


# -----------------------------------------------------
# shapes built from blocks
# -----------------------------------------------------
class test_sphere(unittest.TestCase):
    pass


class test_torus(unittest.TestCase):
    pass


class test_semicylinder_plane(unittest.TestCase):
    pass


class test_semisphere_plane(unittest.TestCase):
    pass


class test_double_bilayer_cylinder(unittest.TestCase):
    pass


class test_capped_cylinder(unittest.TestCase):
    pass


class test_sphere_cylinder(unittest.TestCase):
    pass


class test_buckle(unittest.TestCase):
    pass
