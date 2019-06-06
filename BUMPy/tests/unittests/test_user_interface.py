import unittest
import sys
import os

sys.path.insert(0, os.path.abspath("../.."))   # hacky way to get access to bumpy.py
sys.path.insert(0, os.path.abspath(".."))
from bumpy import check_argument_sanity, check_pdb_dimension_overflow, shapes
from testutils import stderr_checker


class emptyArgs:
    pass


class test_argument_sanity_checker_input_output(unittest.TestCase):

    def setUp(self):
        class args:
            pass
        self.validArgs = args()
        self.validArgs.s = "sphere"
        self.validArgs.f = "reference_files/test_user_interface/reference.pdb"
        self.validArgs.z = "10"
        self.validArgs.g = ["r_sphere:100"]
        self.validArgs.o = "argument_check_output.pdb"
        self.validArgs.p = None
        self.validArgs.n = None
        self.validArgs.gen_dummy_particles = False
        # redirect stderr for text checking
        self.stderr = stderr_checker()
        self.stored_stderr = sys.stderr
        sys.stderr = self.stderr

    def tearDown(self):
        # restore stderr
        sys.stderr = self.stored_stderr
        written_test_files = ["argument_check_output.pdb", "argument_check_output.top", "argument_check_output.ndx",
                              "nodir/output.pdb"]
        for file in written_test_files:
            if os.path.exists(file):
                os.remove(file)

    def test_succeeds_with_good_parameters(self):
        check_argument_sanity(self.validArgs)

    def test_no_shape_throws(self):
        noargs = emptyArgs()
        with self.assertRaises(SystemExit):
            check_argument_sanity(noargs)
        self.assertEqual("No shape was selected. Pick a shape to build using the -s flag\n", str(self.stderr))

    def test_bad_shape_throws(self):
        invalidArgs = self.validArgs
        invalidArgs.s = "spere"
        with self.assertRaises(SystemExit):
            check_argument_sanity(invalidArgs)
        self.assertEqual('Invalid shape selected with argument -s. "spere" is not a valid shape\n', str(self.stderr))

    def test_input_suffix_check(self):
        invalidArgs = self.validArgs
        invalidArgs.f = "nonexistent"
        with self.assertRaises(SystemExit):
            check_argument_sanity(invalidArgs)
        self.assertEqual('User input for option -f was "nonexistent", does not end in .gro or .pdb\nPlease use a valid suffix\n', str(self.stderr))

    def test_input_is_not_a_file(self):
        invalidArgs = self.validArgs
        invalidArgs.f = "nonexistent.pdb"
        with self.assertRaises(SystemExit):
            check_argument_sanity(invalidArgs)
        self.assertEqual('Error: input file does not exist\nUser input for -f : {:s}\n'.format(invalidArgs.f), str(self.stderr))

    def test_input_permissions_error(self):
        # set up a file with no read permissions
        permissionless_file = "nopermissions.pdb"
        open(permissionless_file, 'a').close()
        os.chmod(permissionless_file, 000)
        invalidArgs = self.validArgs
        invalidArgs.f = permissionless_file
        with self.assertRaises(SystemExit):
            check_argument_sanity(invalidArgs)
        os.remove(permissionless_file)
        self.assertEqual('I/O error while trying to read from {:s}\n'.format(invalidArgs.f), str(self.stderr))

    def test_output_error(self):
        # write to a directory that doesn't exist
        nonexistent_dir_path = 'nodir/output.pdb'
        invalidArgs = self.validArgs
        invalidArgs.o = nonexistent_dir_path
        with self.assertRaises(SystemExit):
            check_argument_sanity(invalidArgs)
        self.assertEqual('I/O error while trying to write to {:s}\n'.format(invalidArgs.o), str(self.stderr))


class test_argument_sanity_checker_zo(unittest.TestCase):

    def setUp(self):
        self.invalidArgs = emptyArgs()
        self.invalidArgs.s = "sphere"
        self.invalidArgs.f = "reference_files/test_user_interface/reference.pdb"
        self.invalidArgs.z = None
        self.invalidArgs.o = "output.pdb"
        self.invalidArgs.p = None
        self.invalidArgs.n = None
        self.invalidArgs.g = ["r_sphere:100"]
        self.invalidArgs.gen_dummy_particles = False
        # redirect stderr for text checking
        self.stderr = stderr_checker()
        self.stored_stderr = sys.stderr
        sys.stderr = self.stderr

    def tearDown(self):
        # restore stderr
        sys.stderr = self.stored_stderr
        written_test_files = ["output.pdb"]
        for file in written_test_files:
            if os.path.exists(file):
                os.remove(file)

    def test_zo_warning(self):
        # first case with no zo
        check_argument_sanity(self.invalidArgs)
        self.assertEqual("WARNING : zo was not set with the -z flag. Will use a default value of 10 angstroms. This should lead to " +
                         "sufficient accuracy of lipid areas for most purposes, but you should refer to the BUMPy publication to " +
                         "ensure that the default value is sufficient for your simulation purposes\n", str(self.stderr))

    def test_zo_negative_error(self):
        self.invalidArgs.z = "-1"
        with self.assertRaises(SystemExit):
            check_argument_sanity(self.invalidArgs)
        self.assertEqual("zo cannot be negative\n", str(self.stderr))

    def test_zo_multiple_negatives(self):
        self.invalidArgs.z = "1:-1"
        with self.assertRaises(SystemExit):
            check_argument_sanity(self.invalidArgs)
        self.assertEqual("zo cannot be negative\n", str(self.stderr))

    def test_zo_too_many_values(self):
        # now a zo too long
        self.invalidArgs.z = "1:2:3"
        with self.assertRaises(SystemExit):
            check_argument_sanity(self.invalidArgs)
        self.assertEqual("Too many zo values selected\n", str(self.stderr))

    def test_zo_inconvertible_to_float(self):
        self.invalidArgs.z = "a"
        with self.assertRaises(SystemExit):
            check_argument_sanity(self.invalidArgs)
        self.assertEqual('your input of "a" for zo could not be converted to a floating point number\n', str(self.stderr))


class test_argument_sanity_geometry(unittest.TestCase):

    def setUp(self):
        self.invalidArgs = emptyArgs()
        self.invalidArgs.s = "sphere"
        self.invalidArgs.f = "reference_files/test_user_interface/reference.pdb"
        self.invalidArgs.z = "10"
        self.invalidArgs.o = "output.pdb"
        self.invalidArgs.p = None
        self.invalidArgs.n = None

        # redirect stderr for text checking
        self.stderr = stderr_checker()
        self.stored_stderr = sys.stderr
        sys.stderr = self.stderr

    def tearDown(self):
        sys.stderr = self.stored_stderr
        written_test_files = ["output.pdb"]
        for file in written_test_files:
            if os.path.exists(file):
                os.remove(file)

    def test_no_g_option(self):
        with self.assertRaises(SystemExit):
            check_argument_sanity(self.invalidArgs)
        self.assertEqual("Geometry field -g was not specified. Please specify geometry. For example, ./bumpy -s sphere -g r_sphere:100\n", str(self.stderr))

    def test_g_value_error(self):
        self.invalidArgs.g = ["r_sphere:five"]
        with self.assertRaises(SystemExit):
            check_argument_sanity(self.invalidArgs)
        self.assertEqual("Could not convert one or more -g arguments to float - check you inputs\n", str(self.stderr))

    def test_g_missing_arguments(self):
        self.invalidArgs.s = "cylinder"
        self.invalidArgs.g = ["r_cylinder:100"]
        with self.assertRaises(SystemExit):
            check_argument_sanity(self.invalidArgs)
        self.assertEqual('The following required geometry argument(s) for the shape "{:s}" are missing :\n{:s}\n'.format("cylinder", "l_cylinder"), str(self.stderr))

    def test_g_unknown_args(self):
        self.invalidArgs.g = ["r_sphere:100", "l_cylinder:100"]
        with self.assertRaises(SystemExit):
            check_argument_sanity(self.invalidArgs)
        self.assertEqual('You submitted the following geometric argument(s) using the -g flag, but shape "{:s}" does not use those parameters\n{:s}\n'.format("sphere", "l_cylinder"), str(self.stderr))


class test_argument_sanity_dummy(unittest.TestCase):

    def setUp(self):
        class args:
            pass
        self.validArgs = args()
        self.validArgs.s = "sphere"
        self.validArgs.f = "reference_files/test_user_interface/reference.pdb"
        self.validArgs.z = "10"
        self.validArgs.g = ["r_sphere:100"]
        self.validArgs.o = "argument_check_output.pdb"
        self.validArgs.p = None
        self.validArgs.n = None

        # redirect stderr for text checking
        self.stderr = stderr_checker()
        self.stored_stderr = sys.stderr
        sys.stderr = self.stderr

    def tearDown(self):
        # restore stderr
        sys.stderr = self.stored_stderr
        written_test_files = ["argument_check_output.pdb"]
        for file in written_test_files:
            if os.path.exists(file):
                os.remove(file)

    def test_runs_successfully(self):
        self.validArgs.gen_dummy_particles = True
        self.validArgs.dummy_grid_thickness = 20
        check_argument_sanity(self.validArgs)

    def test_fails_with_no_thickness(self):
        invalidArgs = self.validArgs
        invalidArgs.gen_dummy_particles = True
        with self.assertRaises(SystemExit):
            check_argument_sanity(invalidArgs)
        self.assertEqual('Error: you requested dummy particles with the --gen_dummy_particles option, but did not specify the thickness with --dummy_grid_thickness\n', str(self.stderr))

    def test_fails_dummy_unit_conversion(self):
        invalidArgs = self.validArgs
        invalidArgs.gen_dummy_particles = True
        invalidArgs.dummy_grid_thickness = "twenty"
        with self.assertRaises(SystemExit):
            check_argument_sanity(invalidArgs)
        self.assertEqual('Error: your input of "{:s}" for dummy_grid_thickness could not be converted to a floating point number\n'.format(invalidArgs.dummy_grid_thickness), str(self.stderr))


class test_pdb_overflow_checker(unittest.TestCase):

    def setUp(self):
        self.shape = shapes.flat_bilayer
        self.geometric_arguments = {"x_dimension" : 1000, "y_dimension" : 1000}
        self.output_format = "output.pdb"

        self.stderr = stderr_checker()
        self.stored_stderr = sys.stderr
        sys.stderr = self.stderr

    def test_does_not_trigger_correct_dimensions(self):
        check_pdb_dimension_overflow(self.shape, self.geometric_arguments, self.output_format)

    def test_triggers_with_bad_dimensions(self):
        self.geometric_arguments["x_dimension"] = 11000
        with self.assertRaises(SystemExit):
            check_pdb_dimension_overflow(self.shape, self.geometric_arguments, self.output_format)
        self.assertEqual("Error: You have requested a shape with a long dimension of {:d} angstroms, but the PDB file format ".format(self.geometric_arguments["x_dimension"]) +
                         "is limited to 10k angstroms. We suggest using the .gro format, which allows for up to 100k angstroms.\n",
                         str(self.stderr))

    def test_does_not_trigger_gro_file(self):
        self.geometric_arguments["x_dimension"] = 11000
        self.output_format = "output.gro"
        check_pdb_dimension_overflow(self.shape, self.geometric_arguments, self.output_format)


class test_display_parameters(unittest.TestCase):
    pass
