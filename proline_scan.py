#!/usr/bin/env python3

import sys
import argparse
from PDBtoDSSP import pdb_to_dssp
from DSSPscan import backboneCompatibility as BC
from DSSPscan import mutateSite

'''
def parse():
    parser = argparse.ArgumentParser(description="Check proline compatbible positions in a given structure and the potential cost of mutate them to proline") 
    parser.add_argument('pdbfile', help='Input structure (in pdb format)')
    parser.add_argument('--prob_cutoff', help='Probability cutoff for the compatible mainchain phi-psi angles (default: 0.01)', default=0.01, type=float)
    parser.add_argument('--dist_collision', help='XXX', default=2.8, type=float)
    parser.add_argument('--dist_contact_0', help='XXX', default=3.0, type=float)
    parser.add_argument('--dist_contact_1', help='XXX', default=4.5, type=float)
    parser.add_argument('--sse_window', help='XXX', default=3, type=int)
    parser.add_argument('output', help='Output file name')
    return parser.parse_args()

def main():
    para = parse()
    #
    print('#%s' %(' '.join(sys.argv)))
    #
    print(para.pdbfile)
    print(para.prob_cutoff)
    print(para.dist_collision)
    print(para.dist_contact_1)
    print(para.dist_contact_0)
    print(para.sse_window)
    print(para.output)
'''
if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description = "Check proline compatible positions in a given strucuture and return the potential cost of mutating a site to a proline"
    )
    parser.add_argument('pdbfile',
                        help = 'Input the structure in a pdb file format')

    parser.add_argument('-t','--targetLocation', nargs = "+", required = True,
                        help = "list the values with spaces in between containing name, modelnum, chain name and residuenum eg. '7dwy' 0 'A' 809",
                        )
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
    pFileName = fullID[0]
    targModelNum = int(fullID[1])
    targetChainName = fullID[2]
    targetResID = int(fullID[3])

    referenceStructure = args.referenceStruct
    refModelNum = args.refPosTup[1]
    refChainName = args.refPosTup[2]
    refResID = int(args.refPosTup[3])
    name = args.refPosTup[0]

    countourPlot = args.countourdata
    cutoff = args.prob_cutoff
    contactDiscance = args.dist_contact
    collisionDistance = args.dist_collision

    dsspFile = pdb_to_dssp(pdbFile,"https://www3.cmbi.umcn.nl/xssp/")
    file = open("targetStructure.dssp", "w")
    file.write(dsspFile)
    file.close()


    resultBC = BC(countourPlot,cutoff)
    compatiblePos = BC.compatiblePositions(resultBC,pFileName,pdbFile,"targetStructure.dssp")

    targetList = args.targetLocation
    targetLocation = (targetList[0],int(targetList[1]),targetList[2],(' ',int(targetList[3]),' '))
    print(compatiblePos)
    ''' this is what will create a cost tuple for each compatible position
    costDict = {}
    for id in compatiblePos:
        modelNum = id[0]
        chainName = id[1]
        resID = id[2]
        targLoc = (pFileName,modelNum,chainName,resID)
        print(id)
        print(targLoc)
        costDict[id] = mutateSite(pdbFile,targLoc,pFileName,referenceStructure,name,refModelNum,refChainName,refResID,collisionDistance,contactDiscance)
        print(costDict[id])
    f = open("totalCostDictionary.txt","w")
    f.write(costDict)
    f.close()
    
    '''
    costTuple = mutateSite(pdbFile,targetLocation,pFileName,referenceStructure, name, refModelNum, refChainName, refResID,collisionDistance,contactDiscance)
    print(costTuple)

#python proline_scan.py test/7dwyTestMutateSite.pdb -t "7dwy" 0 "A" 809
