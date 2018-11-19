import unittest
import numpy as np
import sys
import os

sys.path.insert(0, os.path.abspath("../.."))   # hacky way to get access to bumpy.py

from bumpy import Molecules


# TODO : some of these operations will be free-standing function in the future
class test_center_xy(unittest.TestCase):
    pass


class test_center_on_zero(unittest.TestCase):
    pass


class test_translate(unittest.TestCase):
    pass


class test_rotate(unittest.TestCase):
    pass


class test_scale_coordinates_rectangular(unittest.TestCase):
    pass


class test_cylindrical_transform(unittest.TestCase):
    pass


class test_spherical_transform(unittest.TestCase):
    pass


class test_inner_toroidal_transform(unittest.TestCase):
    pass


class test_outer_toroidal_transform(unittest.TestCase):
    pass


class test_scale_flat_to_spherical(unittest.TestCase):
    pass


class test_scale_flat_to_inner_partial_toroid(unittest.TestCase):
    pass


class scale_flat_to_outer_partial_toroid(unittest.TestCase):
    pass
