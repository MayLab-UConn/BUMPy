import unittest
import numpy as np
from bumpy import Molecules, shapes
from tests.testutils import get_relative_path


# -----------------------------------------------------
# interface for all other shapes
# -----------------------------------------------------
class test_base_shape(unittest.TestCase):
    def test_raises_execeptions(self):
        flat_reference = Molecules(get_relative_path("reference_files/test_classes/reference.pdb"))
        baseshape = shapes.shape()
        with self.assertRaises(NotImplementedError):
            baseshape.gen_shape(flat_reference, 0)
        with self.assertRaises(NotImplementedError):
            baseshape.dimension_requirements()
        with self.assertRaises(NotImplementedError):
            baseshape.final_dimensions()


# -----------------------------------------------------
# basic building blocks
# -----------------------------------------------------
class test_shapes_run(unittest.TestCase):

    def setUp(self):
        self.flat_reference = Molecules(get_relative_path("reference_files/test_classes/reference.pdb"))
        self.flat_reference.duplicate_laterally(3, 3)  # should be large enough for all shapes

    def test_shapes_run(self):
        shapes.buckle.gen_shape(                 self.flat_reference, (10, 10), r_buckle=50, l_buckle=50)
        shapes.capped_cylinder.gen_shape(        self.flat_reference, (10, 10), r_cylinder=50, l_cylinder=50)
        shapes.cylinder.gen_shape(               self.flat_reference, (10, 10), r_cylinder=50, l_cylinder=50)
        shapes.double_bilayer_cylinder.gen_shape(self.flat_reference, (10, 10), r_cylinder=50, l_cylinder=50, r_junction=50, l_flat=200)
        shapes.inner_quarter_torus.gen_shape(    self.flat_reference, (10, 10), r_tube=50, r_torus=75)
        shapes.outer_quarter_torus.gen_shape(    self.flat_reference, (10, 10), r_tube=50, r_torus=75)
        shapes.semicylinder_plane.gen_shape(     self.flat_reference, (10, 10), r_cylinder=50, l_cylinder=50, r_junction=50, l_flat=100)
        shapes.semisphere.gen_shape(             self.flat_reference, (10, 10), r_sphere=50)
        shapes.semisphere_plane.gen_shape(       self.flat_reference, (10, 10), r_sphere=50, r_junction=50, l_flat=200)
        shapes.sphere.gen_shape(                 self.flat_reference, (10, 10), r_sphere=50)
        shapes.torus.gen_shape(                  self.flat_reference, (10, 10), r_tube=50, r_torus=150)


class test_shapes_final_dims(unittest.TestCase):
    def test_buckle_final_dims(self):
        # buckle
        # x = l_buckle
        # y = 4 * r_buckle
        # z = 2 * (r_buckle + buffer)
        self.assertTrue(np.allclose(shapes.buckle.final_dimensions(100, 200), np.array((200, 400, 300))))

    def test_capped_cylinder_final_dims(self):
        # capped cylinder
        # x = l_cylinder + 2 * (r_cylinder + buff)
        # y = 2 * (r_cylinder + buff)
        # z = 2 * (r_cylinder + buff)
        self.assertTrue(np.allclose(shapes.capped_cylinder.final_dimensions(100, 200), np.array((500, 300, 300))))

    def test_cylinder_final_dims(self):
        # cylinder
        # x = l_cylinder
        # y = 2 * (buff + r_cylinder)
        # z = 2 * (buff * r_cylinder)
        self.assertTrue(np.allclose(shapes.cylinder.final_dimensions(100, 200), np.array((200, 300, 300))))

    def test_double_bilayer_cylinder_final_dims(self):
        # double_bilayer_cylinder
        # x = l_flat
        # y = l_flat
        # z = r_cylinder + (2 * r_junction) + (2 * buff)
        self.assertTrue(np.allclose(shapes.double_bilayer_cylinder.final_dimensions(100, 200, 60, 600), np.array((600, 600, 420))))

    def test_inner_quarter_torus_final_dims(self):
        # inner_quarter_torus
        # x = 2 * (r_torus + buff)
        # y = 2 * (r_torus + buff)
        # z = (r_tube + buff)
        self.assertTrue(np.allclose(shapes.inner_quarter_torus.final_dimensions(100, 60), np.array((300, 300, 110))))

    def test_outer_quarter_torus_final_dims(self):
        # outer_quarter_torus
        # x = 2 * (r_torus + r_tube + buff)
        # y = 2 * (r_torus + r_tube + buff)
        # z = (r_torus + buff)
        self.assertTrue(np.allclose(shapes.outer_quarter_torus.final_dimensions(100, 60), np.array((420, 420, 110))))

    def test_semicylinder_plane_final_dims(self):
        # semicylinder_plane
        # x = l_cylinder
        # y = l_flat + (2 * r_cylinder) + (2 * r_junction)
        # z = r_cylinder + r_junction +  buff
        self.assertTrue(np.allclose(shapes.semicylinder_plane.final_dimensions(100, 200, 60, 300), np.array((200, 620, 210))))

    def test_semisphere_plane_final_dims(self):
        # semisphere_plane
        # x = l_flat
        # y = l_flat
        # z = r_sphere + r_junction + buff
        self.assertTrue(np.allclose(shapes.semisphere_plane.final_dimensions(100, 60, 400), np.array((400, 400, 210))))

    def test_sphere_final_dims(self):
        # sphere
        # x = 2 * (r_sphere + buff)
        # y = 2 * (r_sphere + buff)
        # z = 2 * (r_sphere + buff)
        self.assertTrue(np.allclose(shapes.sphere.final_dimensions(100), np.array((300, 300, 300))))

    def test_torus_final_dims(self):
        # torus
        # x = 2 * (r_torus + r_tube + buff)
        # y = 2 * (r_torus + r_tube + buff)
        # z = 2 * (r_tube + buff)
        self.assertTrue(np.allclose(shapes.torus.final_dimensions(100, 60), np.array((420, 420, 220))))
