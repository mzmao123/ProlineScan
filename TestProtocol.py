import unittest
from Bio.PDB import *
from DSSPscan import backboneSuperimpose
from DSSPscan import backboneCompatibility as BC
from DSSPscan import prolineConformation as PC
import numpy as np
parser = PDBParser(PERMISSIVE=1)
class TestProlineSuperimpose(unittest.TestCase):
    def test_prolineSuperimpose(self):
        structure = parser.get_structure("7dwyTest","test/7dwyTest.pdb")
        targetResidue = structure[0]["A"][807]
        movingResidue = structure[0]["A"][809]
        result = backboneSuperimpose(targetResidue,movingResidue)
        self.assertAlmostEqual(result, 0, delta=0.01)

class TestBackboneCompatibility(unittest.TestCase):
    def test_probability(self):
        result = BC("test/rama8000-transpro.data")
        phiPsi = (-99.0,-19.0)
        result = BC.returnProb(phiPsi,result)
        self.assertEqual(result,0.001341088541980592)
    def test_compatiblePositions(self):
        structure = parser.get_structure("test", "test/test.pdb")
        result = BC("test/rama8000-transpro.data")
        result = BC.compatiblePositions(result,structure,"test/DSSPtest.dssp")
        expected = []
        self.assertEqual(expected,result)
class TestProlineConformation(unittest.TestCase):
    def test_similarPair(self):
        result = PC("7dwy", "test/DSSP7dwy.dssp", "test/7dwy.pdb")
        result = PC.similarPair([-76.7,167.9],result)
        self.assertEqual(result,[-76.7,167.9])
