import unittest
import numpy as np
import sys
import os

sys.path.insert(0, os.path.abspath("../.."))   # hacky way to get access to bumpy.py

from bumpy import check_argument_sanity
from testutils import stdout_checker


class test_argument_sanity_checker(unittest.TestCase):

    def setUp(self):
        class args:
            pass
        self.validArgs = args()
        self.validArgs.s = "sphere"
        self.validArgs.f = "reference_files/test_user_interface/reference.pdb"
        self.validArgs.z = 10
        self.validArgs.g = "r_sphere:10"
        self.validArgs.o = "argument_check_output.pdb"
        self.validArgs.p = None
        self.validArgs.n = None

        # redirect stdout for text checking
        self.stdout = stdout_checker()
        self.stored_stdout = sys.stdout
        sys.stdout = self.stdout

    def tearDown(self):
        # restore stdout
        sys.stdout = self.stored_stdout
        written_test_files = ["argument_check_output.pdb", "argument_check_output.top", "argument_check_output.ndx"]
        for file in written_test_files:
            if os.path.exists(file):
                os.remove(file)
    def test_succeeds_with_good_parameters(self):
        check_argument_sanity(self.validArgs)

    def test_bad_shape_throws(self):
        invalidArgs = self.validArgs
        # spelling error
        invalidArgs.s = "spere"
        with self.assertRaises(SystemExit):
            check_argument_sanity(invalidArgs)
        self.assertEqual('Invalid shape selected with argument -s. "spere" is not a valid shape\n', str(self.stdout))

    def test_input_suffix_check(self):
        invalidArgs = self.validArgs
        # spelling error
        invalidArgs.f = "nonexistent"
        with self.assertRaises(SystemExit):
            check_argument_sanity(invalidArgs)
        self.assertEqual('User input for option -f was "nonexistent", does not end in .gro or .pdb\nPlease use a valid suffix\n', str(self.stdout))

    def test_input_is_not_a_file(self):
        invalidArgs = self.validArgs
        # spelling error
        invalidArgs.f = "nonexistent.pdb"
        with self.assertRaises(SystemExit):
            check_argument_sanity(invalidArgs)
        self.assertEqual('Error: input file does not exist\nUser input for -f : {:s}\n'.format(invalidArgs.f), str(self.stdout))

    def test_input_permissions_error(self):

        # set up a file with no read permissions
        permissionless_file = "nopermissions.pdb"
        open(permissionless_file, 'a').close()
        os.chmod(permissionless_file, 000)
        invalidArgs = self.validArgs
        # spelling error
        invalidArgs.f = permissionless_file
        with self.assertRaises(SystemExit):
            check_argument_sanity(invalidArgs)
        os.remove(permissionless_file)
        self.assertEqual('I/O error while trying to read from {:s}\n'.format(invalidArgs.f), str(self.stdout))


class test_display_parameters(unittest.TestCase):
    pass
