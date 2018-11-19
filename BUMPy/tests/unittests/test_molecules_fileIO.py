import unittest
import numpy as np
import sys
import os

sys.path.insert(0, os.path.abspath("../.."))   # hacky way to get access to bumpy.py

from bumpy import Molecules, Metadata
from testutils import PDBComp, GROComp
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

        # TODO: test resnums
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

    def setUp(self):
        ''' Sets up a molecule with equivalent fields to the reference files used for reading tests. Not just
            loading a molecule directly from file because if the reading gets broken then the writing tests would
            likely break too
        '''
        metadata = Metadata(atomname=np.array(("A", "B", "A"), dtype="<U4"),
                            resname=np.array(("MOLA", "MOLA", "MOLB"), dtype="<U4"),
                            ressize=np.array((2, 0, 1)))
        coords = np.array([[0.000, -5.000, 0.500],
                           [5.000, -10.000, 0.600],
                           [10.000, -20.000, 0.700]])
        boxdims = np.array([100, 100, 80])
        self.simple_system = Molecules(metadata=metadata, coords=coords, boxdims=boxdims)

    def tearDown(self):
        ''' Removes temporary files generated during output tests '''

        written_test_files = ["test_write_simple.pdb", "test_write_simple.gro"]
        for file in written_test_files:
            if os.path.exists(file):
                os.remove(file)

    def test_output_file_availability(self):
        # do once functionality in place
        pass

    def test_write_pdb(self):
        ''' Want exact match with reference PDB, turn off position changing'''
        temporaryOutputPath = "test_write_simple2.pdb"
        self.simple_system.write_coordinates(temporaryOutputPath, position=False, reorder=False)
        self.assertTrue(PDBComp.compareAtomFields(temporaryOutputPath, referenceFilePath + "read_pdb.pdb"))

    def test_write_gro(self):
        ''' Want exact match with reference GRO, turn off position changing'''
        temporaryOutputPath = "test_write_simple.gro"
        self.simple_system.write_coordinates("test_write_simple.gro", position=False, reorder=False)
        self.assertTrue(GROComp.compareAtomFields(temporaryOutputPath, referenceFilePath + "read_gro.gro"))

    def test_write_header(self):
        ''' Tests that the input commands are written to pdb and gro files '''
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


if __name__ == "__main__":
    unittest.main()
