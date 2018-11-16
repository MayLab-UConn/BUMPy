import unittest
import numpy as np
import sys, os

sys.path.insert(0, os.path.abspath("../.."))   # hacky way to get access to bumpy.py

from bumpy import cart2pol, pol2cart, inner_toroid_angle_from_area, outer_toroid_angle_from_area

class test_cart2pol(unittest.TestCase):

    def test_output_shape(self):
        nparts = 20
        cart_coords = np.zeros((nparts, 3))
        theta, rho, z = cart2pol(cart_coords)
        self.assertEqual(theta.size, nparts)
        self.assertEqual(theta.shape[0], nparts)
        self.assertEqual(rho.size, nparts)
        self.assertEqual(z.size, nparts)

    def test_output_values(self):
        coords = np.array([-4, 3, 2])[np.newaxis, :]
        theta, rho, z = cart2pol(coords)
        self.assertEqual(z[0], 2)
        self.assertEqual(rho[0], 5.0)  # 3 4 5 triangle

        thetacoords = np.array([[0, 0, 0], [1, 1, 0], [-1, 0, 0], [0, -1, 0]])
        theta, rho, z = cart2pol(thetacoords)
        self.assertEqual(theta[0], 0)
        self.assertAlmostEqual(theta[1], np.pi / 4)
        self.assertAlmostEqual(theta[2], np.pi)
        self.assertAlmostEqual(theta[3], - np.pi / 2)  # goes negative after pi

class test_pol2cart(unittest.TestCase):
    def test_output_shape(self):
        pass

    def test_output_values(self):
        pass

class test_inner_toroid_angle_from_area(unittest.TestCase):
    pass

class test_inner_toroid_angle_from_area(unittest.TestCase):
    pass
