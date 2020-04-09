import unittest
import sys
import os
import bumpy.bumpy as bumpy
from bumpy.tests.testutils import FileComp, std_checker, get_relative_path


class test_index_written(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        print("\nTesting command line : index")
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
        written_test_files = ["temp.gro", "index_complex.ndx", "index_with_dummy.ndx"]
        for file in written_test_files:
            if os.path.exists(file):
                os.remove(file)

    def test_complex_index(self):
        sys.argv.extend(["-f", get_relative_path("reference_files/input/input_asymm.gro"),
                         "-n", "index_complex.ndx"])
        bumpy.main()
        self.assertTrue(FileComp.filesMatchExactly("index_complex.ndx", get_relative_path("reference_files/test_index/index_complex.ndx")))

    def test_dummy_index(self):
        sys.argv.extend(["-f", get_relative_path("reference_files/input/input_asymm.gro"),
                         "-n", "index_with_dummy.ndx", "--gen_dummy_particles", "--dummy_grid_thickness", "50"])
        bumpy.main()
        self.assertTrue(FileComp.filesMatchExactly("index_with_dummy.ndx", get_relative_path("reference_files/test_index/index_with_dummy.ndx")))
