import unittest
from Bio.PDB import *
from DSSPscan import prolineSuperimpose
import numpy as np
parser = PDBParser(PERMISSIVE=1)
class TestProlineSuperimpose(unittest.TestCase):
    def test_prolineSuperimpose(self):
        structure = parser.get_structure("7dwy","7dwy.pdb")
        targetResidue = structure[0]["A"][807]
        movingResidue = structure[0]["A"][809]
        result = prolineSuperimpose(targetResidue,movingResidue)
        expectedCoords = [[128.554, 126.011,170.703]
                         [128.179 ,125.745 ,169.309]
                         [126.79 , 125.16  ,169.124]
                         [126.301 ,124.402 ,169.964]
                         [129.252, 124.756, 168.853]
                         [130.408, 125.068 ,169.687]
                         [129.878 ,125.46 , 171.022]]