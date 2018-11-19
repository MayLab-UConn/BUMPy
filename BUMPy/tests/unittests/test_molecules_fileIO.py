import unittest
import numpy as np
import sys
import os

sys.path.insert(0, os.path.abspath("../.."))   # hacky way to get access to bumpy.py

from bumpy import Molecules


class test_read_input(unittest.TestCase):

    def test_argument_sanity_checks(self):
        pass

    def test_read_pdb(self):
        pass

    def test_read_gro(self):
        pass

    def test_ignore_resnames(self):
        pass


class test_write_coordinates(unittest.TestCase):

    def test_output_file_availability(self):
        # do once functionality in place
        pass

    def test_write_pdb(self):
        pass

    def test_write_gro(self):
        pass

    def test_write_header(self):
        pass

    def test_position_shifter(self):
        pass


class test_write_topology(unittest.TestCase):

    def test_basic_functionality(self):
        pass

    def test_multiple_molecule_types(self):
        pass


class test_write_index(unittest.TestCase):
    pass
