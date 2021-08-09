#!/usr/bin/env python3

import sys
import argparse
from PDBtoDSSP import pdb_to_dssp
from DSSPscan import backboneCompatibility as BC
from DSSPscan import mutateSite, mutateSiteWhole
import json

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description = "Check proline compatible positions in a given strucuture and return the potential cost of mutating a site to a proline"
    )
    parser.add_argument('pdbfile',
                        help = 'Input the structure in a pdb file format')

    parser.add_argument('-t','--targetLocation', nargs = "+",
                        help = "list the values with spaces in between containing name, modelnum, chain name and residuenum eg. '7dwy' 0 'A' 809", default = "N"
                        )
    parser.add_argument('-n','--structureName',
                        help = "if you are inputing the whole structure and not just a target location, enter the PDB ID of the structure eg. '7dwy' ")
    parser.add_argument('-r','--referenceStruct',
                        help = 'Input the pdb file for the structure from which the proline will come from, default will be the 7dwy.pdb structure',
                        default = 'test/7dwy.pdb')
    parser.add_argument('-l','--refPosTup', nargs="+",
                        help = "list the values with spaces in between containing name, modelnum, chain name and residuenum eg. '7dwy' 0 'A '809. default is ['7dwy',0,'A',809]",
                        default = ['7dwy',0,'A',809], type = list )
    parser.add_argument('-p','--prob_cutoff',
                        help='Probability cutoff for the compatible mainchain phi-psi angles (default: 0.01)',
                        default=0.01, type=float)
    parser.add_argument('-c','--dist_collision', help='collision distance cutoff, the default value is 2.8 angstroms',
                        default=2.8, type=float)
    parser.add_argument('-i','--dist_contact',
                        help='collision distance cutoff, the default value is 4.5 angstroms',
                        default=4.5, type=float)
    parser.add_argument('-d','--countourdata',
                        help = 'the data for the probabilities of given phi psi angles, the defaul data file is rama8000-transpro.data',
                        default = 'test/rama8000-transpro.data')
    args = parser.parse_args()

    pdbFile = args.pdbfile

    fullID = args.targetLocation

    referenceStructure = args.referenceStruct
    refModelNum = args.refPosTup[1]
    refChainName = args.refPosTup[2]
    refResID = int(args.refPosTup[3])
    name = args.refPosTup[0]

    countourPlot = args.countourdata
    cutoff = args.prob_cutoff
    contactDiscance = args.dist_contact
    collisionDistance = args.dist_collision

    dsspFile = pdb_to_dssp(pdbFile, "https://www3.cmbi.umcn.nl/xssp/")
    file = open("targetStructure.dssp", "w")
    file.write(dsspFile)
    file.close()

    if fullID != "N":
        pFileName = fullID[0]
        targModelNum = int(fullID[1])
        targetChainName = fullID[2]
        targetResID = int(fullID[3])
        targetLocation = (fullID[0], int(fullID[1]), fullID[2], (' ', int(fullID[3]), ' '))

        resultBC = BC(countourPlot,cutoff)
        compatiblePos = BC.compatiblePositions(resultBC,pFileName,pdbFile,"targetStructure.dssp")
        ID = (targModelNum,targetChainName, targetResID)
        if ID not in compatiblePos:
            print("WARNING: This residue is likely not compatible with a proline")
        costTuple = mutateSite(pdbFile, targetLocation, pFileName, referenceStructure, name, refModelNum, refChainName,refResID, collisionDistance, contactDiscance)
        print(costTuple)

    if fullID == "N":
        structName = args.structureName
        resultBC = BC(countourPlot, cutoff)
        compatiblePos = BC.compatiblePositions(resultBC, structName, pdbFile, "targetStructure.dssp")
        costDict = mutateSiteWhole(pdbFile, structName, referenceStructure, name, refModelNum, refChainName,refResID, collisionDistance, contactDiscance, compatiblePos)
        print(costDict)
        #with open("totalCostDictionary.txt","w") as file:
        #    file.write(json.dumps(costDict))

#python proline_scan.py test/7dwyTestMutateSite.pdb -t "7dwy" 0 "A" 809
#python proline_scan.py test/7dwyTestMutateSite.pdb -n "7dwyTestMutateSite"
