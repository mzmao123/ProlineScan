import unittest
from Bio.PDB import *
from DSSPscan import prolineSuperimpose
import numpy as np
parser = PDBParser(PERMISSIVE=1)
class TestProlineSuperimpose(unittest.TestCase):
    def test_prolineSuperimpose(self):
        structure = parser.get_structure("7dwyTest","7dwyTest.pdb")
        targetResidue = structure[0]["A"][807]
        movingResidue = structure[0]["A"][809]
        result = prolineSuperimpose(targetResidue,movingResidue)
        self.assertAlmostEqual(result, 0, delta=0.01)