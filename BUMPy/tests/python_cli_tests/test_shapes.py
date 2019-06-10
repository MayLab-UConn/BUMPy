import unittest
import sys
import os
import bumpy
from tests.testutils import PDBComp, std_checker, get_relative_path


class test_shapes(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        print("\nTesting command line : shapes")
        cls.tempstdout = sys.stdout
        sys.stdout = std_checker()

    @classmethod
    def tearDownClass(cls):
        sys.stdout = cls.tempstdout

    def setUp(self):
        self.geometries = { 'buckle' : ['r_buckle:50', 'l_buckle:60'],
                            'capped_cylinder' : ['r_cylinder:50', 'l_cylinder:60'],
                            'cylinder' : ['r_cylinder:50', 'l_cylinder:60'],
                            'double_bilayer_cylinder' : ['r_cylinder:50', 'l_cylinder:50', 'r_junction:40', 'l_flat:250'],
                            'flat_bilayer' : ['x_dimension:50', 'y_dimension:40'],
                            'inner_quarter_torus' : ['r_torus:75', 'r_tube:50'],
                            'outer_quarter_torus' : ['r_torus:75', 'r_tube:50'],
                            'semicylinder_plane' : ['r_cylinder:50', 'l_cylinder:30', 'r_junction:40', 'l_flat:50'],
                            'semisphere' : ['r_sphere:50'],
                            'semisphere_plane' : ['r_sphere:50', 'r_junction:30', 'l_flat:200'],
                            'sphere' : ['r_sphere:50'],
                            'torus' : ['r_torus:75', 'r_tube:50']
                            }
        self.storeArgv = sys.argv
        sys.argv = ["bumpy.py", "-z", "10",
                    "-f", get_relative_path("reference_files/input/input_asymm.gro")]

    def tearDown(self):
        sys.argv = self.storeArgv
        written_test_files = [get_relative_path("test_{:s}.pdb".format(shape)) for shape in self.geometries.keys()]
        for file in written_test_files:
            if os.path.exists(file):
                os.remove(file)

    def test_shapes(self):
        base_args = len(sys.argv)
        for shape, geometry in self.geometries.items():
            sys.argv.extend(["-s", shape, "-o", get_relative_path("test_{:s}.pdb".format(shape)), "-g"])
            sys.argv.extend(geometry)
            bumpy.main()

            success = PDBComp.compareAtomFields(get_relative_path("test_{:s}.pdb".format(shape)),
                                                get_relative_path("reference_files/test_shapes/test_{:s}.pdb".format(shape)))
            if not success:
                os.rename("test_{:s}.pdb".format(shape), "failed_test_{:s}.pdb".format(shape))
            self.assertTrue(success, "Failure in shape {:s}".format(shape))
            sys.argv = sys.argv[:base_args]
