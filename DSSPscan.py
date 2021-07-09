import Bio.PDB.DSSP
from PDBtoDSSP import pdb_to_dssp
from Bio.PDB import *
from Bio.PDB.DSSP import *
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
def prolineSuperimpose(targetResidue):
    sup = Superimposer()
    structure = parser.get_structure("7dwy","7dwy.pdb")
    proline1 = structure[0]["A"][807]
    proline2 = targetResidue
    fixedAtoms = [atom for atom in proline1 if atom.name in ["C", "N", "O", "CA"]]
    movingAtoms = [atom for atom in proline2 if atom.name in ["C", "N", "O", "CA"]]
    sup.set_atoms(fixedAtoms,movingAtoms)
    movingModel = proline1.get_parent()
    print(sup.rotran)
    print(sup.rms)
    sup.apply(movingModel.get_atoms())


"""
def prolineList():
    returnDict = {}
    structure = parser.get_structure("7dwy", "7dwy.pdb")
    model = structure[0]
    for chain in model:
        resID = []
        for residue in chain:
            if residue.get_resname() == "PRO":
                resID.append(residue.get_id)
            else:
                pass
        returnDict[chain] = resID
    return returnDict
"""