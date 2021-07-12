from PDBtoDSSP import pdb_to_dssp
from Bio.PDB import *
from Bio.PDB.DSSP import *
import numpy as np
from mpl_toolkits import mplot3d
import matplotlib.pyplot as plt
import json
parser = PDBParser(PERMISSIVE=1)
def dsspList(dsspFile,id):
    def tupleToList(dict): #turns the tuple definitions of the key into lists for better manipulation
        returnDict = {}
        for key in dict:
            returnDict[key] = list(dict[key])
        return returnDict
    def dsspDictionarty(dsspFile): #creates the dictionary for each resiude containing all of the values
        returnDict = {}
        returnDict = make_dssp_dict(dsspFile)
        return returnDict
    def listRemove(dsspFile): #removes values from the dictionary so that we are just left with desired values
        dict = dsspFile
        for key in dict:
            del dict[key][5:]
            dict[key].pop(0)
        return dict
    dssp_tuple = dsspDictionarty(dsspFile)
    temp_dict = dssp_tuple[0]
    dssp_dict = tupleToList(temp_dict)
    dssp_dict = listRemove(dssp_dict)
    return dssp_dict[id] # id required should be in this format ('A', (' ', 14, ' ')), 'A' is the name of the chain
def prolineSuperimpose(targetResidue, movingResidue):
    sup = Superimposer()
    proline1 = targetResidue
    proline2 = movingResidue
    fixedAtoms = [atom for atom in proline1 if atom.name in ["C", "N", "O", "CA"]]
    movingAtoms = [atom for atom in proline2 if atom.name in ["C", "N", "O", "CA"]]
    sup.set_atoms(fixedAtoms,movingAtoms)
    sup.apply(proline2)
    return sup.rms

def conformationPlot(fileName, fileDirectory,dsspFile):
    def chiCalc(residue):
        N = residue["N"]
        CA = residue["CA"]
        CB = residue["CB"]
        CG = residue["CG"]
        v1 = N.get_vector()
        v2 = CA.get_vector()
        v3 = CB.get_vector()
        v4 = CG.get_vector()
        angle = calc_dihedral(v1,v2,v3,v4)
        return angle

    def prolineList(fileName,fileDirectory):
        returnDict = {}
        chiAngle = []
        structure = parser.get_structure(fileName, fileDirectory)
        model = structure[0]
        for chain in model:
            resID = []
            for residue in chain:
                if residue.get_resname() == "PRO":
                    resID.append(residue.id)
                    chiAngle.append(chiCalc(residue))
                else:
                    pass
            returnDict[chain.id] = resID
        return (returnDict,chiAngle)

    def idCreation(dict):
        keyList = []
        for key in dict:
            for i in dict[key]:
                addition = (key,i)
                keyList.append(addition)
        return keyList

    def psiPhiChi (fileName, fileDirectory, dsspFile):
        listVal = prolineList(fileName, fileDirectory)
        idList = idCreation(listVal[0])
        allChi = listVal[1]
        allPhi = []
        allPsi = []
        for id in idList:
            allPsi.append((dsspList(dsspFile,id))[3])
            allPhi.append((dsspList(dsspFile,id))[2])
        return(allPhi,allPsi,allChi)

    fig = plt.figure()
    ax = plt.axes(projection='3d')
    dataTuple = psiPhiChi(fileName, fileDirectory, dsspFile)
    zdata = np.array(dataTuple[2])
    ydata = np.array(dataTuple[1])
    xdata = np.array(dataTuple[0])
    ax.scatter3D(xdata,ydata,zdata,c=zdata,cmap='Greens')
    plt.show()

conformationPlot("7dwy","test/7dwy.pdb","test/DSSP7dwy.dssp")

