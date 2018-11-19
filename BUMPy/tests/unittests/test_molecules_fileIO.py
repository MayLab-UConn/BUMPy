import unittest
import numpy as np
import sys
import os

sys.path.insert(0, os.path.abspath("../.."))   # hacky way to get access to bumpy.py

from bumpy import Molecules
referenceFilePath = "reference_files/test_molecules_fileIO/"


class test_read_input(unittest.TestCase):

    # TODO : add checks for ignoring HETATOMS, etc
    def test_argument_sanity_checks(self):
        pass

    def test_read_pdb(self):
        inputPDBPath = referenceFilePath + "read_pdb.pdb"
        testMolecule = Molecules(infile=inputPDBPath)

        # test coordinates
        self.assertEqual(testMolecule.coords.shape, (3, 3))
        self.assertAlmostEqual(testMolecule.coords[0, 0], 0)
        self.assertAlmostEqual(testMolecule.coords[1, 1], -10)
        self.assertAlmostEqual(testMolecule.coords[2, 2], 0.7)

        # test atomnames
        self.assertEqual(testMolecule.metadata.atomname[0].strip(), "A")
        self.assertEqual(testMolecule.metadata.atomname[1].strip(), "B")
        self.assertEqual(testMolecule.metadata.atomname[2].strip(), "A")

        # test resnames
        self.assertEqual(testMolecule.metadata.resname[0], "MOLA")
        self.assertEqual(testMolecule.metadata.resname[1], "MOLA")
        self.assertEqual(testMolecule.metadata.resname[2], "MOLB")

    def test_read_gro(self):
        inputGROPath = referenceFilePath + "read_gro.gro"
        testMolecule = Molecules(infile=inputGROPath)

        # test coordinates
        self.assertEqual(testMolecule.coords.shape, (3, 3))
        self.assertAlmostEqual(testMolecule.coords[0, 0], 0)
        self.assertAlmostEqual(testMolecule.coords[1, 1], -10)
        self.assertAlmostEqual(testMolecule.coords[2, 2], 0.7)

        # test atomnames
        self.assertEqual(testMolecule.metadata.atomname[0].strip(), "A")
        self.assertEqual(testMolecule.metadata.atomname[1].strip(), "B")
        self.assertEqual(testMolecule.metadata.atomname[2].strip(), "A")

        # test resnames
        self.assertEqual(testMolecule.metadata.resname[0], "MOLA")
        self.assertEqual(testMolecule.metadata.resname[1], "MOLA")
        self.assertEqual(testMolecule.metadata.resname[2], "MOLB")

    def test_ignore_resnames(self):
        inputPDBPath = referenceFilePath + "read_pdb.pdb"
        testMolecule = Molecules(infile=inputPDBPath, ignore="MOLB")
        self.assertEqual(testMolecule.coords.shape, (2, 3))
        self.assertAlmostEqual(testMolecule.coords[0, 0], 0)
        self.assertAlmostEqual(testMolecule.coords[1, 1], -10)

        # test atomnames
        self.assertEqual(testMolecule.metadata.atomname[0].strip(), "A")
        self.assertEqual(testMolecule.metadata.atomname[1].strip(), "B")

        # test resnames
        self.assertEqual(testMolecule.metadata.resname[0], "MOLA")
        self.assertEqual(testMolecule.metadata.resname[1], "MOLA")


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
