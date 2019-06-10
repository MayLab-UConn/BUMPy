import unittest
import numpy as np
from bumpy import gen_dummy_grid


class test_dummy_grid_generation(unittest.TestCase):

    def test_basic_functionality(self):
        grid = gen_dummy_grid()
        # check size
        self.assertEqual(grid.coords.shape, (200, 3))
        # check naming
        self.assertEqual(grid.metadata.atomname[0], "DUMY")
        self.assertEqual(grid.metadata.resname[199], "DUMY")
        # check leaflet assignment
        self.assertTrue(np.all(grid.metadata.leaflets[:100] == 0))
        self.assertTrue(np.all(grid.metadata.leaflets[100:] == 1))
        # check coordinates
        self.assertAlmostEqual(np.sqrt( ((grid.coords[1, 0:2] - grid.coords[0, 0:2]) ** 2).sum()), 5)
        self.assertTrue(np.all(grid.coords[:100, 2] == -25))
        self.assertTrue(np.all(grid.coords[100:, 2] ==  25))

    # TODO : populate once written in main code
    def test_argument_sanity_check(self):
        pass

    def test_lateral_distance(self):
        grid = gen_dummy_grid(lateral_distance=7.5)
        self.assertAlmostEqual(np.sqrt( ((grid.coords[1, 0:2] - grid.coords[0, 0:2]) ** 2).sum()), 7.5)

    def test_thickness(self):
        grid = gen_dummy_grid(thickness=100)
        self.assertTrue(np.all(grid.coords[:100, 2] == -50))
        self.assertTrue(np.all(grid.coords[100:, 2] ==  50))

    def test_atomname(self):
        grid = gen_dummy_grid(atomname="TEST")
        self.assertEqual(grid.metadata.atomname[0], "TEST")

    def test_resname(self):
        grid = gen_dummy_grid(resname="TEST")
        self.assertEqual(grid.metadata.resname[0], "TEST")
