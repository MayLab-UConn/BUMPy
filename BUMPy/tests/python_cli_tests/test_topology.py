import unittest
import sys
import os
import bumpy
from tests.testutils import FileComp, std_checker, get_relative_path


class test_topology_written(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        print("\nTesting command line : topology")
        cls.tempstdout = sys.stdout
        sys.stdout = std_checker()

    @classmethod
    def tearDownClass(cls):
        sys.stdout = cls.tempstdout

    def setUp(self):
        self.storeArgv = sys.argv
        sys.argv = ["bumpy.py", "-s", "sphere", "-o", "temp.gro", "-z", "10", "-g", "r_sphere:50"]

    def tearDown(self):
        sys.argv = self.storeArgv
        written_test_files = ["temp.gro", "topol_complex.top", "topol_with_dummy.top"]
        for file in written_test_files:
            if os.path.exists(file):
                os.remove(file)

    def test_complex_topology(self):
        sys.argv.extend(["-f", get_relative_path("reference_files/input/input_asymm.gro"),
                         "-p", "topol_complex.top"])
        bumpy.main()
        self.assertTrue(FileComp.filesMatchExactly("topol_complex.top",
                                                   get_relative_path("reference_files/test_topology/topol_complex.top")))

    def test_dummy_topology(self):
        sys.argv.extend(["-f", get_relative_path("reference_files/input/input_asymm.gro"),
                         "-p", "topol_with_dummy.top", "--gen_dummy_particles", "--dummy_grid_thickness", "50"])
        bumpy.main()
        self.assertTrue(FileComp.filesMatchExactly("topol_with_dummy.top",
                        get_relative_path("reference_files/test_topology/topol_with_dummy.top")))
