import Bio.PDB.DSSP
from PDBtoDSSP import pdb_to_dssp
from Bio.PDB import *
from Bio.PDB.DSSP import *
import numpy as np
import json
parser = PDBParser(PERMISSIVE=1)
def dsspList(dsspFile,id):
    def tupleToList(dict):
        returnDict = {}
        for key in dict:
            returnDict[key] = list(dict[key])
        return returnDict
    def dsspDictionarty(dsspFile):
        returnDict = {}
        returnDict = make_dssp_dict(dsspFile)
        return returnDict
    def listRemove(dsspFile):
        dict = dsspFile
        for key in dict:
            del dict[key][5:]
            dict[key].pop(0)
        return dict
    dssp_tuple = dsspDictionarty(dsspFile)
    temp_dict = dssp_tuple[0]
    dssp_dict = tupleToList(temp_dict)
    dssp_dict = listRemove(dssp_dict)
    return dssp_dict[id]
def prolineSuperimpose(targetResidue, movingResidue):
    sup = Superimposer()
    proline1 = targetResidue
    proline2 = movingResidue
    fixedAtoms = [atom for atom in proline1 if atom.name in ["C", "N", "O", "CA"]]
    movingAtoms = [atom for atom in proline2 if atom.name in ["C", "N", "O", "CA"]]
    print("Fixed Atoms ")
    for atom in fixedAtoms:
        print(str(atom) + " " + str(atom.get_coord()))
    print("Moving Atoms ")
    for atom in movingAtoms:
        print(str(atom) + " " + str(atom.get_coord()))
    sup.set_atoms(fixedAtoms,movingAtoms)
    movingModel = proline1.get_parent()
    print(sup.rotran)
    print(sup.rms)
    sup.apply(movingModel.get_atoms())
    return proline2
