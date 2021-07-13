from PDBtoDSSP import pdb_to_dssp
from Bio.PDB import PDBParser, Superimposer, calc_dihedral
from Bio.PDB.DSSP import *
import numpy as np
import math
from mpl_toolkits import mplot3d
import matplotlib.pyplot as plt
import json
parser = PDBParser(PERMISSIVE=1)
def dsspList(dsspFile,id):

    def tupleToList(dsspDict): #turns the tuple definitions of the key into lists for better manipulation
        returnDict = {}
        for key in dsspDict:
            returnDict[key] = list(dsspDict[key])
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

def backboneSuperimpose(residue1, residue2): #Fucntion will superimpose two given residues
    sup = Superimposer()
    targetResidue = residue1 #the residue that will not move
    movingResidue = residue2 #the residue that will be moved. Will be superimposed over the target residue
    fixedAtoms = [atom for atom in targetResidue if atom.name in ["C", "N", "O", "CA"]]
    movingAtoms = [atom for atom in movingResidue if atom.name in ["C", "N", "O", "CA"]]
    sup.set_atoms(fixedAtoms,movingAtoms)
    sup.apply(movingResidue)
    return sup.rms

def chiCalc(residue): # calculates the chi angle of a giver residue
    N = residue["N"]
    CA = residue["CA"]
    CB = residue["CB"]
    CG = residue["CG"]
    v1 = N.get_vector()
    v2 = CA.get_vector()
    v3 = CB.get_vector()
    v4 = CG.get_vector()
    angle = calc_dihedral(v1,v2,v3,v4)
    angle = angle*(180/math.pi)
    return angle

def prolineDict(fileName,fileDirectory): #takes a pdb file and returns a dictionary of prolines with the key being the chain. This function is also repurposed to calculate the chi angles for all of the prolines in the file.
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
        returnDict[chain.id] = resID #eg. {'A': [(id1),(id2),(id3)]]
    return (returnDict,chiAngle)

def idCreation(residueDictionary): #takes the return dictioary from the function prolineDict and creates a new list of tuples which can be accepted by the dsspList function
    keyList = []
    for key in residueDictionary:
        for i in residueDictionary[key]:
            addition = (key,i) #this would look like ('A",("",80,"")) which is the id that is accepted by the dsspList function
            keyList.append(addition)
    return keyList

def psiPhiChi (fileName, fileDirectory, dsspFile): #creates a list of psi and phi and chi angles of all the prolines in the protein
    listVal = prolineDict(fileName, fileDirectory)
    idList = idCreation(listVal[0])
    allChi = listVal[1]
    allPhi = []
    allPsi = []
    for id in idList:
        allPsi.append((dsspList(dsspFile,id))[3])
        allPhi.append((dsspList(dsspFile,id))[2])
    return(allPhi,allPsi,allChi)

def conformationPlot(fileName, fileDirectory,dsspFile): #This takes the results from psiPhiChi and plots them on a 3d axis
    fig = plt.figure()
    ax = plt.axes(projection='3d')
    dataTuple = psiPhiChi(fileName, fileDirectory, dsspFile)
    zdata = np.array(dataTuple[2])
    ydata = np.array(dataTuple[1])
    xdata = np.array(dataTuple[0])
    ax.scatter3D(xdata,ydata,zdata,c=zdata,cmap='Greens')
    plt.xlabel("Phi Angle")
    plt.ylabel("Psi Angle")
    ax.set_zlabel("Chi Angle")
    ax.set_zlim(-180,180)
    ax.set_xlim(-180, 180)
    ax.set_ylim(-200, 200)
    plt.show()

def compatiblePositions(fileName,fileDirectory,dsspFile,dataFile):  #scans all the prolines in the protein and compares it to a list of comformations from a database and returns the ids of the prolines that match the acceptable conformations

    def toNum (value): #turns an exponential into a decimal
        valWO = value.split('e')
        ret_val = format(((float(valWO[0]))*(10**int(valWO[1]))), '.8f')
        return ret_val

    def angleCutoff(conformationDataFile): #scans the database for acceptable conformations. The probability of a conformation has to be greater than 0.01
        anglesCutoff = []
        lines = []
        with open(conformationDataFile) as file:
            lines = file.readlines()
            for line in lines:
                line = line.split()
                prob = float(line[2])
                try:
                    prob = float(toNum(str(prob)))
                except IndexError:
                    pass
                if prob>0.01:
                    anglesCutoff.append([float(line[0]),float(line[1])])
        return anglesCutoff

    angleList = angleCutoff(dataFile)
    listCompatible = []
    listVal = prolineDict(fileName, fileDirectory)
    idList = idCreation(listVal[0])

    for id in idList: #compares the phi and psi angles for all the prolines and compares it with the list of acceptable phi and psi conformations.
        phiPsi = []
        tempList = dsspList(dsspFile,id)
        phiPsi.extend([tempList[2],tempList[3]])
        if phiPsi in angleList:
            listCompatible.append(id)
    return listCompatible

print(compatiblePositions("7dwy","test/7dwy.pdb","test/DSSP7dwy.dssp","test/rama8000-transpro.data"))




