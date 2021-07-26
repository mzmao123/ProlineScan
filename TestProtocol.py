import unittest
from Bio.PDB import *
from DSSPscan import backboneSuperimpose
from DSSPscan import backboneCompatibility as BC
from DSSPscan import prolineConformation as PC
from DSSPscan import mutateSite
from DSSPscan import distanceBetweenResidues
import numpy as np
parser = PDBParser(PERMISSIVE=1)
class TestProlineSuperimpose(unittest.TestCase):
    def test_prolineSuperimpose(self):
        structure = parser.get_structure("7dwyTest","test/7dwyTestSuperimpose.pdb")
        targetResidue = structure[0]["A"][807]
        movingResidue = structure[0]["A"][809]
        result = backboneSuperimpose(targetResidue,movingResidue)
        self.assertAlmostEqual(result, 0, delta=0.01)

class TestBackboneCompatibility(unittest.TestCase):
    def test_probability(self):
        result = BC("test/rama8000-transpro.data")
        phiPsi = (round(-84.9),round(-7.1))
        result = BC.returnProb(phiPsi,result)
        self.assertEqual(result,0.09410222409740221)
    def test_compatiblePositions(self):
        structure = parser.get_structure("test", "test/test.pdb")
        result = BC("test/rama8000-transpro.data")
        result = BC.compatiblePositions(result,structure,"test/DSSPtest.dssp")
        expected = [('A', (' ', 18, ' ')), ('A', (' ', 22, ' ')), ('A', (' ', 23, ' ')), ('A', (' ', 25, ' ')), ('A', (' ', 26, ' ')), ('A', (' ', 32, ' ')), ('A', (' ', 38, ' ')), ('A', (' ', 39, ' ')), ('A', (' ', 43, ' ')), ('A', (' ', 54, ' '))]
        self.assertEqual(expected,result)
class TestProlineConformation(unittest.TestCase):
    def test_similarPair(self):
        result = PC("7dwy", "test/DSSP7dwy.dssp", "test/7dwy.pdb")
        result = PC.similarPair([-76.7,167.9],result)
        self.assertEqual(result,[-76.7,167.9])
class TestMutateSite(unittest.TestCase):
    def test_mutateSite(self):
        testStructure = parser.get_structure("7dwy","test/7dwyTestMutateSite.pdb")
        structure1 = parser.get_structure("7dwy","test/7dwy.pdb")
        modifiedGly = testStructure[0]["A"][809]
        result = mutateSite(testStructure,modifiedGly.get_full_id())
        expected = structure1[0]["A"][809]
        self.assertEqual(result,expected)
class TestDistanceBetween(unittest.TestCase):
    def test_distanceBetween(self):
        structure = parser.get_structure("7dwy", "test/7dwy.pdb")
        proline = structure[0]["A"][809]
        result = distanceBetweenResidues(proline)
        print(result)