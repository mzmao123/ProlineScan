from Bio.PDB import PDBParser, Superimposer, calc_dihedral
from Bio.PDB.DSSP import make_dssp_dict
from Bio.PDB.SASA import *

import numpy as np
import math
from mpl_toolkits import mplot3d
import matplotlib.pyplot as plt
import json
import os

from PDBtoDSSP import pdb_to_dssp
from Bio.PDB.PDBIO import PDBIO
from collections import namedtuple
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

def dsspListWORemove(dsspFile,id):
    def tupleToList(dsspDict): #turns the tuple definitions of the key into lists for better manipulation
        returnDict = {}
        for key in dsspDict:
            returnDict[key] = list(dsspDict[key])
        return returnDict
    def dsspDictionarty(dsspFile): #creates the dictionary for each resiude containing all of the values
        returnDict = {}
        returnDict = make_dssp_dict(dsspFile)
        return returnDict

    dssp_tuple = dsspDictionarty(dsspFile)
    temp_dict = dssp_tuple[0]
    dssp_dict = tupleToList(temp_dict)

    return dssp_dict[id]

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

def toNum (value): #turns an exponential into a decimal
    valWO = value.split('e')
    ret_val = format(((float(valWO[0]))*(10**int(valWO[1]))), '.8f')
    return ret_val

class backboneCompatibility():
    def __init__(self,contourPlot,cutoff):
        self.angleDict = {}
        self.anglesCutoff = []
        with open(contourPlot) as file:
            self.lines = file.readlines()
            for line in self.lines:
                line = line.split()
                prob = float(line[2])
                try:
                    prob = float(toNum(str(prob)))
                except IndexError:
                    pass
                angVal = (float(line[0]),float(line[1]))
                self.angleDict[angVal] = prob
    # scans the database for acceptable conformations. The probability of a conformation has to be greater than 0.01
            for line in self.lines:
                line = line.split()
                prob = float(line[2])
                try:
                    prob = float(toNum(str(prob)))
                except IndexError:
                    pass
                if prob > cutoff:
                    self.anglesCutoff.append([float(line[0]), float(line[1])])

    def returnProb(phiPsi,self):
        return self.angleDict[phiPsi]


    def compatiblePositions(self,fileName,pdbFile,dsspFile):  # scans all the prolines in the protein and compares it to a list of comformations from a database and returns the ids of the prolines that match the acceptable conformations
        aminoList = ["ALA","ARG","ASN","ASP","CYS","GLN","GLU","GLY","HIS","LLE","ILE","LEU","LYS","MET","PHE","PRO","PYL","SER","SEC","THR","TRP","TYR","VAL"]
        struct = parser.get_structure(fileName,pdbFile)

        angleList = self.anglesCutoff
        phiList = []
        psiList = []

        for lis in angleList:
            phiList.append(lis[0])
            psiList.append(lis[1])

        idList = []
        idListWithModelNum = []
        listCompatible = []
        #model = struct[0]
        for model in struct:
            for chain in model:
                for residue in chain:
                    if residue.get_resname() in aminoList:
                        addition = (chain.id,residue.id)
                        secondAddition = (model.id,chain.id,residue.id)
                        idList.append(addition)
                        idListWithModelNum.append(secondAddition)

        for i in range(len(idList)):  # compares the phi and psi angles for all the prolines and compares it with the list of acceptable phi and psi conformations.
            phiPsi = []
            id = idList[i]
            tempList = dsspList(dsspFile, id)
            secondID = idListWithModelNum[i]
            phiPsi.extend([round(tempList[2]), round(tempList[3])])
            phi = round(tempList[2])
            psi = round(tempList[3])
            phiRange = set([phi-1,phi,phi+1])
            psiRange = set([psi-1,psi,psi+1])
            psiSet = set(psiList)
            phiSet = set(phiList)
            phiIntersectVal = phiSet.intersection(phiRange)
            psiIntersectVal = psiSet.intersection(psiRange)

            for val in phiIntersectVal:
                for val2 in psiIntersectVal:
                    if [val,val2] in angleList:
                        if secondID not in listCompatible:
                            listCompatible.append(secondID)
                        else:
                            pass
                        break
        return listCompatible

    def conformationPlot(fileName, fileDirectory,dsspFile):  # This takes the results from psiPhiChi and plots them on a 3d axis
        dataTuple = psiPhiChi(fileName, fileDirectory, dsspFile)
        return dataTuple

class prolineConformation():
    def __init__(self, fileName, dsspFile, fileDirectory):

        def prolineDict(fileName,fileDirectory):  # takes a pdb file and returns a dictionary of prolines with the key being the chain. This function is also repurposed to calculate the chi angles for all of the prolines in the file.
            returnDict = {}
            structure = parser.get_structure(fileName, fileDirectory)
            model = structure[0]
            for chain in model:
                resID = []
                for residue in chain:
                    if residue.get_resname() == "PRO":
                        resID.append(residue.id)
                    else:
                        pass
                returnDict[chain.id] = resID  # eg. {'A': [(id1),(id2),(id3)]]
            return (returnDict)
        def idCreation(residueDictionary):  # takes the return dictioary from the function prolineDict and creates a new list of tuples which can be accepted by the dsspList function
            keyList = []
            for key in residueDictionary:
                for i in residueDictionary[key]:
                    addition = (key,i)  # this would look like ('A",("",80,"")) which is the id that is accepted by the dsspList function
                    keyList.append(addition)
            return keyList
        idList = idCreation(prolineDict(fileName,fileDirectory))
        self.conformDict = {}
        self.psiList = []
        self.phiList = []
        for id in idList:
            angleList = dsspList(dsspFile, id)
            key = (angleList[2],angleList[3])
            self.conformDict[key] = id
            self.phiList.append(angleList[2])
            self.psiList.append(angleList[3])

    def similarPair(phiPsi,self):
        phi = phiPsi[0]
        psi = phiPsi[1]
        diff = 10000
        returnPhi = 0
        returnPsi = 0

        for i in range(len(self.phiList)):
            distance = abs((self.phiList[i]-phi)+(self.psiList[i]-psi))
            if (distance)<diff:
                diff = distance
                returnPhi = self.phiList[i]
                returnPsi = self.psiList[i]
        return [returnPhi,returnPsi]

def mutateSite(pdbFile,targetResidueFullID, targName, referenceStructure, refName, modelNum, chainName, resNum, colNum, conNum): #reference structure file is the pdb file for the structure that you want to extract the proline that you are going to mutate onto the target structure.

    structure = parser.get_structure(targName,pdbFile)
    newStructure = structure.copy()
    orgStruct = structure
    chain = newStructure[targetResidueFullID[1]][targetResidueFullID[2]]
    res = chain[targetResidueFullID[3]]
    count = 0
    for residue in chain:
        if residue != res:
            count+=1
        else:
            break
    resId = res.id
    chain.detach_child(resId)
    referenceStructure = parser.get_structure(refName,referenceStructure)
    replacePro = referenceStructure[modelNum][chainName][resNum]
    replacePro._id=resId
    chain.insert(count,replacePro)
    io = PDBIO()
    io.set_structure(newStructure)
    io.save("structureFile.pdb")
    costTuple = replacementCost(orgStruct,newStructure,chain.id,resId,chain[resId],targetResidueFullID[1],colNum,conNum)
    return costTuple

def distanceBetweenResidues(residue,discol,discon):
    collisions = 0
    mainChainAtoms = ["C", "N", "O", "CA"]
    mutatedResidueSideChain = []
    for atom in residue:
        if atom.get_name() not in mainChainAtoms:
            mutatedResidueSideChain.append(atom)
    contacts = 0
    chain = residue.get_parent()
    id = residue.get_id()[1]
    for res in chain:
        if res != residue:
            distance = 0
            for atom in res:
                for mutateAtom in mutatedResidueSideChain:
                    distance = abs(mutateAtom-atom)
                    if abs(distance)<=discon:
                        contacts+=1
                    if abs(distance)<=discol:
                        collisions+=1
    return [collisions,contacts]

def replacementCost(orgStrcut,newStructure, chainId, resId, mutatedResidue, modelID, distCollision, distContact):
    io = PDBIO()
    io.set_structure(newStructure)
    io.save("modifiedStructure.pdb")
    modifiedStructureDSSP = pdb_to_dssp("modifiedStructure.pdb","https://www3.cmbi.umcn.nl/xssp/")
    file = open("modifiedStructure.dssp", "w")
    file.write(modifiedStructureDSSP)
    file.close()
    id = (chainId, resId)
    listVal = dsspListWORemove("modifiedStructure.dssp",id)

    resNum = resId[1]
    orgResidue = orgStrcut[modelID][chainId][resNum]
    sr = ShrakeRupley()
    sr.compute(orgResidue, level='R')
    asa = orgResidue.sasa
    n_collisions = distanceBetweenResidues(mutatedResidue, distCollision, distContact)[0]

    n_contacts_wt = distanceBetweenResidues(orgResidue, distCollision, distContact)[1]
    n_contacts_pro = distanceBetweenResidues(mutatedResidue, distCollision, distContact)[1]

    ACC = listVal[2]
    mutate_Cost = namedtuple("Mutate_Cost",['phi', 'psi', 'SASA', 'SSE', 'ACC', 'n_collisions', 'n_contacts_wt', 'n_contacts_pro'])
    costTuple = mutate_Cost(listVal[3], listVal[4], asa, listVal[1], ACC, n_collisions, n_contacts_wt, n_contacts_pro) #phi,psi,asa,sse...
    return costTuple

    #os.remove("structureFile.dssp")
    #os.remove("structureFile.pdb")


'''
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
'''
