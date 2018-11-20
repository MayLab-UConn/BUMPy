import unittest
import numpy as np
import sys
import os

sys.path.insert(0, os.path.abspath("../.."))   # hacky way to get access to bumpy.py

from bumpy import Molecules


class test_assign_leaflets(unittest.TestCase):
    pass


class test_get_bilayer_center(unittest.TestCase):
    pass


class test_reorder_by_leaflet(unittest.TestCase):
    pass


class test_reorder_within_leaflet(unittest.TestCase):
    pass


class test_reorder_with_dummies_in_back(unittest.TestCase):
    pass
