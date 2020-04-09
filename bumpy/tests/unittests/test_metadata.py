import unittest
import numpy as np
from bumpy.bumpy import Metadata


# TODO checks for numpy in metadata
class test_Metadata(unittest.TestCase):
    def setUp(self):
        self.metadata = Metadata(atomname=np.array(("A1", "A2", "A3", "A4"), dtype="<U4"),
                                 resname=np.array(("R1", "R2", "R2", "R3"),  dtype="<U4"),
                                 leaflets=np.array((0, 1, 1, 1)),
                                 ressize=np.array((1, 2, 0, 1)))

    def test_append(self):
        mdata2 = Metadata(atomname=np.array(("B1"), dtype="<U4"), resname=np.array(("B1"), dtype="<U4"),
                          leaflets=np.array((0)), ressize=np.array((1)))
        self.metadata.append(mdata2)
        self.assertTrue(np.alltrue(self.metadata.atomname == np.array(("A1", "A2", "A3", "A4", "B1"), dtype="<U4")))
        self.assertTrue(np.alltrue(self.metadata.resname == np.array(("R1", "R2", "R2", "R3", "B1"), dtype="<U4")))
        self.assertTrue(np.alltrue(self.metadata.leaflets == np.array((0, 1, 1, 1, 0))))
        self.assertTrue(np.alltrue(self.metadata.ressize == np.array((1, 2, 0, 1, 1))))

    def test_reorder(self):
        self.metadata.reorder(np.array((3, 1, 2, 0)))
        self.assertTrue(np.alltrue(self.metadata.atomname == np.array(("A4", "A2", "A3", "A1"), dtype="<U4")))
        self.assertTrue(np.alltrue(self.metadata.resname  == np.array(("R3", "R2", "R2", "R1"), dtype="<U4")))
        self.assertTrue(np.alltrue(self.metadata.leaflets == np.array((1, 1, 1, 0))))
        self.assertTrue(np.alltrue(self.metadata.ressize ==  np.array((1, 2, 0, 1))))

    def test_duplicate(self):
        self.metadata.duplicate(2)
        self.assertTrue(np.alltrue(self.metadata.atomname == np.array(("A1", "A2", "A3", "A4", "A1", "A2", "A3", "A4"), dtype="<U4")))
        self.assertTrue(np.alltrue(self.metadata.resname  == np.array(("R1", "R2", "R2", "R3", "R1", "R2", "R2", "R3"), dtype="<U4")))
        self.assertTrue(np.alltrue(self.metadata.leaflets == np.array((0, 1, 1, 1, 0, 1, 1, 1))))
        self.assertTrue(np.alltrue(self.metadata.ressize ==  np.array((1, 2, 0, 1, 1, 2, 0, 1))))

    def test_slice(self):
        # TODO warn on empty slice - happens a lot with code errors
        newMeta = self.metadata.slice((np.array([1, 3])))
        self.assertTrue(np.alltrue(newMeta.atomname == np.array(("A2", "A4"), dtype="<U4")))
        self.assertTrue(np.alltrue(newMeta.resname  == np.array(("R2", "R3"), dtype="<U4")))
        self.assertTrue(np.alltrue(newMeta.leaflets == np.array((1, 1))))
        self.assertTrue(np.alltrue(newMeta.ressize  == np.array((2, 1))))


if __name__ == "__main__":
    unittest.main()
