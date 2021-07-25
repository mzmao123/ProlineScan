from PDBtoDSSP import pdb_to_dssp
from Bio.PDB import PDBParser, Superimposer, calc_dihedral
from Bio.PDB.DSSP import make_dssp_dict
from Bio.PDB.SASA import *
import numpy as np
import math
from mpl_toolkits import mplot3d
import matplotlib.pyplot as plt
import json
from PDBtoDSSP import pdb_to_dssp
from Bio.PDB.PDBIO import PDBIO
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
def toNum (value): #turns an exponential into a decimal
    valWO = value.split('e')
    ret_val = format(((float(valWO[0]))*(10**int(valWO[1]))), '.8f')
    return ret_val


'''
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
'''

class backboneCompatibility():
    def __init__(self,contourPlot):
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
                if prob > 0.01:
                    self.anglesCutoff.append([float(line[0]), float(line[1])])

    def returnProb(phiPsi,self):
        return self.angleDict[phiPsi]


    def compatiblePositions(self,struct,dsspFile):  # scans all the prolines in the protein and compares it to a list of comformations from a database and returns the ids of the prolines that match the acceptable conformations

        angleList = self.anglesCutoff
        phiList = []
        psiList = []

        for lis in angleList:
            phiList.append(lis[0])
            psiList.append(lis[1])

        idList = []
        listCompatible = []
        model = struct[0]

        for chain in model:
            for residue in chain:
                addition = (chain.id,residue.id)
                idList.append(addition)

        for id in idList:  # compares the phi and psi angles for all the prolines and compares it with the list of acceptable phi and psi conformations.
            phiPsi = []
            tempList = dsspList(dsspFile, id)
            phiPsi.extend([round(tempList[2]), round(tempList[3])])
            phi = round(tempList[2])
            psi = round(tempList[3])

            #if phiPsi in angleList or [phiPsi[0]+1,phiPsi[1]] in angleList or[phiPsi[0],phiPsi[1]+1] in angleList or [phiPsi[0]+1,phiPsi[1]+1] in angleList or [phiPsi[0]-1,phiPsi[1]] in angleList or [phiPsi[0],phiPsi[1]-1] in angleList or [phiPsi[0]-1,phiPsi[1]-1] in angleList or [phiPsi[0]-1,phiPsi[1]+1] in angleList or [phiPsi[0]+1,phiPsi[1]-1] in angleList:
            phiRange = set([phi-1,phi,phi+1])
            psiRange = set([psi-1,psi,psi+1])
            psiSet = set(psiList)
            phiSet = set(phiList)
            phiIntersectVal = phiSet.intersection(phiRange)
            psiIntersectVal = psiSet.intersection(psiRange)

            for val in phiIntersectVal:
                for val2 in psiIntersectVal:
                    if [val,val2] in angleList:
                        if id not in listCompatible:
                            listCompatible.append(id)
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

def mutateSite(structure,residueFullID):
    newStructure = structure.copy()
    orgStruct = structure
    chain = newStructure[residueFullID[1]][residueFullID[2]]
    residue  = chain[residueFullID[3]]
    resId = residue.id
    replaceID = resId[1]
    chain.detach_child(resId)
    libraryStructure = parser.get_structure("7dwy","test/7dwy.pdb")
    replacePro = libraryStructure[0]["A"][809]
    replacePro.id = resId
    chain.insert(replaceID,replacePro)
    replacementCost(orgStruct,newStructure,chain.id,resId,chain[resId],residueFullID[1])
    return chain[resId]

def distanceBetweenResidues(residue):
    collisions = 0
    contacts = 0
    chain = residue.get_parent()
    id = residue.get_id()[1]
    for res in chain:
        if res != residue:
            try:
               distance = res["CA"]-residue["CA"]
            except KeyError:
                continue
            if abs(distance)<=4.5:
                contacts+=1
            if abs(distance)<=2.8:
                collisions+=1

    return [collisions,contacts]

def replacementCost(orgStrcut,newStructure, chainId, resId, mutatedResidue, modelID):
    io = PDBIO()
    io.set_structure(newStructure)
    io.save("structureFile.pdb")
    import os
    modifiedStructureDSSP = pdb_to_dssp("structureFile.pdb","https://www3.cmbi.umcn.nl/xssp/")
    file = open("structureFile.dssp", "w")
    file.write(modifiedStructureDSSP)
    file.close()
    structureDict = make_dssp_dict("structureFile.dssp")
    id = (chainId, resId)
    listVal = dsspListWORemove("structureFile.dssp",id)
    print(listVal)
    sr = ShrakeRupley()
    sr.compute(mutatedResidue, level='R')
    asa = mutatedResidue.sasa
    n_collisions = distanceBetweenResidues(mutatedResidue)[0]
    resNum = resId[1]
    n_contacts_wt = distanceBetweenResidues(orgStrcut[modelID][chainId][resNum])[1]
    n_contacts_pro = distanceBetweenResidues(mutatedResidue)[1]
    costTuple = (listVal[3],listVal[4],asa,listVal[1],n_collisions,n_contacts_wt,n_contacts_pro) #phi,psi,asa,sse...
    print(costTuple)
    os.remove("structureFile.dssp")
    os.remove("structureFile.pdb")