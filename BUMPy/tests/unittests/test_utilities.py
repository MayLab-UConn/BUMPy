import unittest
import numpy as np
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
        nparts = 20
        theta = rho = z = np.zeros(nparts)
        cartcoords = pol2cart(theta, rho, z)
        self.assertEqual(cartcoords.shape, (nparts, 3))

    def test_output_values(self):
        rho = np.array((0, 1, 10))
        theta = np.array((0, np.pi / 2, -np.pi))
        z = np.array((-5, 0, 10))
        cartcoords = pol2cart(theta, rho, z)
        # test simple cases
        self.assertEqual(      cartcoords[0, 0],   0)
        self.assertEqual(      cartcoords[0, 1],   0)
        self.assertEqual(      cartcoords[0, 2],  -5)
        self.assertAlmostEqual(cartcoords[1, 0],   0)
        self.assertAlmostEqual(cartcoords[1, 1],   1)
        self.assertAlmostEqual(cartcoords[1, 2],   0)
        self.assertAlmostEqual(cartcoords[2, 0], -10)
        self.assertAlmostEqual(cartcoords[2, 1],   0)
        self.assertAlmostEqual(cartcoords[2, 2],  10)
        # test intermediate cases - split xy, and gt 2pi
        rho   = np.array([10, 10])
        theta = np.array([np.pi / 4, 3 * np.pi])
        z     = np.zeros(2)
        cartcoords = pol2cart(theta, rho, z)
        self.assertAlmostEqual(cartcoords[0, 0], 10 / np.sqrt(2))
        self.assertAlmostEqual(cartcoords[0, 1], 10 / np.sqrt(2))
        self.assertAlmostEqual(cartcoords[1, 0], -10)
        self.assertAlmostEqual(cartcoords[1, 1], 0)


class test_inner_toroid_angle_from_area(unittest.TestCase):
    def test_edge_cases(self):
        r_torus = 20
        r_tube = 10
        A_quarter_torus_inner = (np.pi ** 2) * r_torus * r_tube - 2 * np.pi * (r_tube ** 2)
        self.assertEqual(0, inner_toroid_angle_from_area(r_torus, r_tube, 0))
        self.assertAlmostEqual(np.pi / 2, inner_toroid_angle_from_area(r_torus, r_tube, A_quarter_torus_inner)[0])


class test_outer_toroid_angle_from_area(unittest.TestCase):
    def test_edge_cases(self):
        r_torus = 20
        r_tube = 10
        A_quarter_torus_outer = (np.pi ** 2) * r_torus * r_tube + 2 * np.pi * (r_tube ** 2)
        self.assertEqual(0, outer_toroid_angle_from_area(r_torus, r_tube, 0))
        self.assertAlmostEqual(np.pi / 2, outer_toroid_angle_from_area(r_torus, r_tube, A_quarter_torus_outer)[0])


if __name__ == "__main__":
    unittest.main()
