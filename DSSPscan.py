from Bio.PDB import PDBParser, Superimposer, calc_dihedral
from Bio.PDB.DSSP import make_dssp_dict

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
                angVal = (int(line[0]),int(line[1]))
                self.angleDict[angVal] = prob

    def returnProb(self, phiPsi):
        phi = round(phiPsi[0])
        psi = round(phiPsi[1])
        nb_pairs = []
        for i in (0, -1, 1):
            for j in (0, -1, 1):
                phi_nb = phi + i
                psi_nb = psi + j
                nb_pairs.append((phi_nb, psi_nb))
        probability = 0.0
        for phi_psi_pair in nb_pairs:
            if phi_psi_pair in self.angleDict:
                probability = self.angleDict[phi_psi_pair]
                break
        return probability

    def compatiblePositions(self, chain, dssp_dict, prob_cutoff=0.01):  # scans all the prolines in the protein and compares it to a list of comformations from a database and returns the ids of the prolines that match the acceptable conformations
        amino_acids = set("ALA ARG ASN ASP CYS GLN GLU GLY HIS LLE ILE LEU LYS MET PHE PRO PYL SER SEC THR TRP TYR VAL".split())
        listCompatible  = []
        for res in chain.get_residues():
            if res.name not in amino_acids: continue
            aa, sse, asa, phi, psi = dssp_dict[(chain.id, res.id)]
            probability = self.returnProb(phi, psi)
            if probability > prob_cutoff:
                listCompatible.append((resid, (sse, asa, phi, psi, probability)))
        return listCompatible

    def conformationPlot(fileName, fileDirectory,dsspFile):  # This takes the results from psiPhiChi and plots them on a 3d axis
        dataTuple = psiPhiChi(fileName, fileDirectory, dsspFile)
        return dataTuple

class RepresentiveProlines:
    def __init__(self, filelist):
        self.prolines = {}
        for fn in filelist:
            p, phi, psi = fn.split('.')[0].split('_')
            phi = int(phi)
            psi = int(psi)
            proline = parser.get_structure('', fn).get_residues()[0]
            if (phi, psi) in self.prolines:
                self.prolines[(phi, psi)].append(proline)
            else:
                self.prolines[(phi, psi)] = [proline]

    def get(self, phi, psi):
        dist = 1.0E9
        best_match_phi_psi = (0.0, 0.0)
        for rep_phi, rep_psi in self.prolines:
            rep_dist = abs(phi-rep_phi)+ abs(psi-rep_psi)
            if rep_dist < dist:
                best_match_phi_psi = (rep_phi, rep_psi)
                dist = rep_dist
        return self.prolines[best_match_phi_psi]

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
        diff = 1
        returnPhi = 0
        returnPsi = 0

        for i in range(len(self.phiList)):
            distance = abs((self.phiList[i]-phi)+(self.psiList[i]-psi))
            if (distance)<diff:
                diff = distance
                returnPhi = self.phiList[i]
                returnPsi = self.psiList[i]
        return [returnPhi,returnPsi]

def sidechain_compatibility(model, chain_id, res_id, dssp_dict, rep_prolines, collision_th, contact_dmin, contact_dmax):
    aa, sse, asa, phi, psi = dssp_dict[(chain_id, res_id)]
    # number of contacts in WT
    n_contacts_wt = cnt_sidechain_contacts(model, chain_id, res_id, dist_range=(contact_dmin, contact_dmax))
    #
    closest_prolines = rep_prolines.get((phi, psi))
    min_n_collision = 1000
    for pro in closest_prolines:
        working_model = model.copy()
        # backbone superimposition
        fixed_res = working_model[chain_id][res_id]
        moving_res = pro.copy()
        rms = backboneSuperimpose(fixed_res, moving_res)
        # replace the residue
        working_model[chain_id].detach_child(res_id)
        moving_res.id = res_id
        working_model[chain_id].add(moving_res)
        # collision, contacts afer the replacement
        n_collision = cnt_sidechain_contacts(working_model, chain_id, res_id, dist_range=(0,collision_th))
        if n_collision < min_n_collision:
            min_n_collision = n_collision
            n_contacts_pro = cnt_sidechain_contacts(working_model, chain_id, res_id, dist_range=(contact_dmin,contact_dmax))
    return (n_collision, n_contacts_wt, n_contacts_pro)

def cnt_sidechain_contacts(model, chain_id, res_id, dist_range):
    # get the atoms of interest
    main_chain_atoms = set('C N O CA'.split())
    res = model[chain_id][res_id]
    sc_atoms = [atom for atom in res if atom.name not in main_chain_atoms]
    other_atoms = [atom for atom in model.get_atoms() if atom.parent != res]
    # count number of atoms in the distance range
    squared_dmin = dist_range[0]*dist_range[0]
    squared_dmax = dist_range[1]*dist_range[1]
    sqdist_mat = get_sqdist_mat(other_atoms, sc_atoms)
    sqdist_min = numpy.min(sqdist_mat, axis=1)
    result = [(atom, min_dist) for atom, min_dist in zip(other_atoms, sqdist_min) if min_dist > squared_dmin and min_dist < squared_dmax]
    return len(result)

def get_sqdist_mat(atoms_a, atoms_b):
    coordsa = np.array([atom.coord for atom in atoms_a])
    coordsb = np.array([atom.coord for atom in atoms_b])
    return sqdist(coordsa, coordsb)

def sqdist(xyza, xyzb):
    ''' Get the distance matrix between coords array xyza and xyzb.

    Input: 
        xyza: [[xa1, ya1, za1], [xa2, ya2, za2], ...]
        xyzb: [[xb1, yb1, zb1], [xb2, yb2, zb2], ...]

    Output:
        distmatrix: (an x bn)
        [[D_a1_b1, D_a1_b2, D_a1_b3, ..., D_a1_bn], 
         [D_a2_b1, D_a2_b2, D_a2_b3, ..., D_a2_bn], 
         .
         .
         .
         [D_an_b1, D_an_b2, D_an_b3, ..., D_an_bn], 
    '''
    sizea = xyza.shape[0]
    sizeb = xyzb.shape[0]
    mat_a = xyza.reshape(sizea, 1, 3)
    mat_a = mat_a.repeat(sizeb, axis=1)
    # mat_a:
    # [[[xa1, ya1, za1], [[xa1, ya1, za1], ...],
    #  [[xa2, ya2, za2], [[xa2, ya2, za2], ...], 
    #  .
    #  .
    #  .
    #  [[xan, yan, zan], [[xan, yan, zan], ...]]
    mat_b = xyzb.reshape(1, sizeb, 3)
    mat_b = mat_b.repeat(sizea, axis=0)
    # mat_b:
    # [[[xb1, yb1, zb1], [xb2, yb2, zb2], ...],
    #  [[xb1, yb1, zb1], [xb2, yb2, zb2], ...],
    #  .
    #  .
    #  .
    #  [[xb1, yb1, zb1], [xb2, yb2, zb2], ...]]
    dist = mat_a - mat_b
    dist = numpy.sum(dist * dist, axis=2)
    return dist

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
    res2ID = (resId[0],resId[1]+1,resId[2])
    res3ID = (resId[0],resId[1]-1,resId[2])

    listVal = dsspListWORemove("modifiedStructure.dssp",id)
    try:
        listVal2 = dsspListWORemove("modifiedStructure.dssp",res2ID)
    except:
        listVal2 = listVal
    try:
        listVal3 = dsspListWORemove("modieifiedStructure.dssp",res3ID)
    except:
        listVal3 = listVal
    SSE1 = listVal2[1]
    SSE2 = listVal3[1]
    SSE = ""
    if SSE1 == "-":
        SSE = SSE2
    else:
        SSE = SSE1

    resNum = resId[1]
    orgResidue = orgStrcut[modelID][chainId][resNum]

    n_collisions = distanceBetweenResidues(mutatedResidue, distCollision, distContact)[0]

    n_contacts_wt = distanceBetweenResidues(orgResidue, distCollision, distContact)[1]
    n_contacts_pro = distanceBetweenResidues(mutatedResidue, distCollision, distContact)[1]

    ACC = listVal[2]
    asa = ACC


    mutate_Cost = namedtuple("Mutate_Cost",['phi', 'psi', 'SASA', 'SSE', 'ACC', 'n_collisions', 'n_contacts_wt', 'n_contacts_pro'])
    costTuple = mutate_Cost(listVal[3], listVal[4], asa, SSE, ACC, n_collisions, n_contacts_wt, n_contacts_pro) #phi,psi,asa,sse...
    return costTuple

    #os.remove("structureFile.dssp")
    #os.remove("structureFile.pdb")

def mutateSiteWhole(pdbFile, structName, referenceStructure, refName, modelNum, chainName, resNum, colNum, conNum, compatible): #reference structure file is the pdb file for the structure that you want to extract the proline that you are going to mutate onto the target structure.

    structure = parser.get_structure(structName,pdbFile)
    referenceStructure = parser.get_structure(refName, referenceStructure)
    replaceProMaster = referenceStructure[modelNum][chainName][resNum]

    newStructure = structure.copy()
    orgStruct = structure


    for id in compatible:
        replacepro = replaceProMaster.copy()
        chain = newStructure[id[0]][id[1]]
        res = chain[id[2]]
        count = 0
        for residue in chain:
            if residue != res:
                count+=1
            else:
                break
        resId = res.id
        chain.detach_child(resId)

        replacepro._id=resId

        chain.insert(count,replacepro)

        residue = chain[resId]
    io = PDBIO()
    io.set_structure(newStructure)
    io.save("structureFile.pdb")
    costDict = replacementCostDict(orgStruct,newStructure,colNum,conNum, compatible)
    return costDict

def replacementCostDict(orgStrcut,newStructure, distCollision, distContact, compatible):
    costDict = {}
    io = PDBIO()
    io.set_structure(newStructure)
    io.save("modifiedStructure.pdb")
    modifiedStructureDSSP = pdb_to_dssp("modifiedStructure.pdb", "https://www3.cmbi.umcn.nl/xssp/")
    file = open("modifiedStructure.dssp", "w")
    file.write(modifiedStructureDSSP)
    file.close()
    for pos in compatible:
        modelID = pos[0]
        chainId = pos[1]
        resId = pos[2]
        mutatedResidue = newStructure[modelID][chainId][resId]
        id = (pos[1],resId)
        res2ID = (resId[0], resId[1] + 1, resId[2])
        res3ID = (resId[0], resId[1] - 1, resId[2])

        listVal = dsspListWORemove("modifiedStructure.dssp", id)
        try:
            listVal2 = dsspListWORemove("modifiedStructure.dssp", res2ID)
        except:
            listVal2 = listVal
        try:
            listVal3 = dsspListWORemove("modieifiedStructure.dssp", res3ID)
        except:
            listVal3 = listVal
        SSE1 = listVal2[1]
        SSE2 = listVal3[1]
        SSE = ""
        if SSE1 == "-":
            SSE = SSE2
        else:
            SSE = SSE1

        resNum = resId[1]
        orgResidue = orgStrcut[modelID][chainId][resNum]
        sr = ShrakeRupley()
        sr.compute(orgResidue, level='R')
        asa = orgResidue.sasa
        n_collisions = distanceBetweenResidues(mutatedResidue, distCollision, distContact)[0]

        n_contacts_wt = distanceBetweenResidues(orgResidue, distCollision, distContact)[1]
        n_contacts_pro = distanceBetweenResidues(mutatedResidue, distCollision, distContact)[1]

        ACC = listVal[2]

        mutate_Cost = namedtuple("Mutate_Cost", ['phi', 'psi', 'SASA', 'SSE', 'ACC', 'n_collisions', 'n_contacts_wt','n_contacts_pro'])
        costTuple = mutate_Cost(listVal[3], listVal[4], asa, SSE, ACC, n_collisions, n_contacts_wt,n_contacts_pro)  # phi,psi,asa,sse...
        costDict[pos] = costTuple
    return costDict
