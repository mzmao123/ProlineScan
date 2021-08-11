import os
import math

import numpy as np
from Bio.PDB import PDBParser, Superimposer

from PDBtoDSSP import pdb_to_dssp
from Bio.PDB.PDBIO import PDBIO
from collections import namedtuple
parser = PDBParser(PERMISSIVE=1, QUIET=True)

AA3to1 = {'ALA': 'A', 'CYS': 'C', 'ASP': 'D', 'GLU': 'E', 'PHE': 'F',
          'GLY': 'G', 'HIS': 'H', 'ILE': 'I', 'LYS': 'K', 'LEU': 'L',
          'MET': 'M', 'ASN': 'N', 'PRO': 'P', 'GLN': 'Q', 'ARG': 'R',
          'SER': 'S', 'THR': 'T', 'VAL': 'V', 'TRP': 'W', 'TYR': 'Y',
          'MSE': 'M', 'UNK': 'X'}

def backboneSuperimpose(residue1, residue2): #Fucntion will superimpose two given residues
    sup = Superimposer()
    targetResidue = residue1 #the residue that will not move
    movingResidue = residue2 #the residue that will be moved. Will be superimposed over the target residue
    fixedAtoms = [atom for atom in targetResidue if atom.name in ["C", "N", "O", "CA"]]
    movingAtoms = [atom for atom in movingResidue if atom.name in ["C", "N", "O", "CA"]]
    sup.set_atoms(fixedAtoms,movingAtoms)
    sup.apply(movingResidue)
    return sup.rms

BCInfo = namedtuple('BCInfo','resid sse asa phi psi prob'.split())
amino_acids = set("ALA ARG ASN ASP CYS GLN GLU GLY HIS LLE ILE LEU LYS MET PHE PRO PYL SER SEC THR TRP TYR VAL".split())

class BackboneCompatibility:
    def __init__(self, rama8000_fn, comment='#'):
        self.rama8000 = {}
        with open(rama8000_fn) as infile:
            for line in infile.readlines():
                cleaned_line = line.strip()
                if not cleaned_line.startswith('#'):
                    phi, psi, prob = cleaned_line.split()
                    self.rama8000[(int(float(phi)), int(float(psi)))] = float(prob)

    def prob(self, phi, psi):
        nb_pairs = []
        for i in (0, -1, 1):
            for j in (0, -1, 1):
                phi_nb = round(phi) + i
                psi_nb = round(psi) + j
                nb_pairs.append((phi_nb, psi_nb))
        probability = 0.0
        for phi_psi_pair in nb_pairs:
            if phi_psi_pair in self.rama8000:
                probability = self.rama8000[phi_psi_pair]
                break
        return probability

    def compatible_sites(self, chain, dssp_dict, prob_cutoff=0.01):  # scans all sides in the given chain 
        sites = []
        for res in chain.get_residues():
            if res.resname not in amino_acids: continue
            aa, sse, asa, phi, psi = dssp_dict[(chain.id, res.id)][:5]
            probability = self.prob(phi, psi)
            if probability > prob_cutoff:
                sites.append(BCInfo(res.id, sse, asa, phi, psi, probability))
        return sites

class RepresentiveProlines:
    def __init__(self, filelist):
        self.prolines = {}
        for fn in filelist:
            p, phi, psi = os.path.basename(fn).split('.')[0].split('_')  # assumed the file named as 'XX_{phi}_{psi}.XX' - should be improved later
            phi = int(phi)
            psi = int(psi)
            proline = next(parser.get_structure('', fn).get_residues())
            # print(proline)
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

def mutate_to_proline(model, chain_id, res_id, dssp_dict, rep_prolines, collision_th):
    aa, sse, asa, phi, psi = dssp_dict[(chain_id, res_id)][:5]
    closest_prolines = rep_prolines.get(phi, psi)  # Two of proline of different side chain conformations returned,  
    min_n_collision = 100
    min_collisions = None
    mutated_model = None
    for pro in closest_prolines:
        # Try both conformation, keep the one introduce minimum number of collisions
        working_model = model.copy()
        # backbone superimposition
        fixed_res = working_model[chain_id][res_id]
        moving_res = pro.copy()
        rms = backboneSuperimpose(fixed_res, moving_res)
        # replace the residue
        i = count_to_res(fixed_res)
        working_model[chain_id].detach_child(res_id)
        moving_res.id = res_id
        working_model[chain_id].insert(i, moving_res)
        # choose the proline conformation leads to fewer collisions
        collisions = sidechain_contacted_atoms(working_model, chain_id, res_id, dist_range=(0,collision_th))
        if len(collisions) < min_n_collision:
            min_collisions = collisions
            min_n_collision = len(collisions)
            mutated_model = working_model
    return (mutated_model, min_collisions)

def sidechain_contacted_atoms(model, chain_id, res_id, dist_range):
    dmin, dmax = dist_range
    main_chain_atoms = set('C N O CA'.split())
    res = model[chain_id][res_id]
    if res.resname != 'GLY':
        # get the atoms of interest 1: side chain atoms of selected residue
        sc_atoms = [atom for atom in res if atom.name not in main_chain_atoms]
        i = count_to_res(res)
        res_minus = res.parent.child_list[i-1]
        # get the atoms of interest 2: atoms on other resiudes
        other_atoms = []
        for atom in model.get_atoms():
            if (atom.parent == res_minus and atom.name in ['C', 'CA']) or atom.parent == res:
                # skip C/CA at i-1 (distance to proline_i is legally < 2.8A) and self
                continue
            else:
                other_atoms.append(atom)
        # select other_atoms within the distance range of sc_atoms
        squared_dmin = dmin * dmin
        squared_dmax = dmax * dmax
        sqdist_mat = get_sqdist_mat(other_atoms, sc_atoms)
        sqdist_min = np.min(sqdist_mat, axis=1)
        # select by the minimum distance for each other_atoms to sc_atoms
        result = [(atom, min_dist) for atom, min_dist in zip(other_atoms, sqdist_min) if min_dist > squared_dmin and min_dist < squared_dmax]
    else:  # GLY has no side chain atoms
        result = []
    return result

# auxiliary functions
def count_to_res(res):
    for i, r in enumerate(res.parent):
        if r == res:
            break
    return i

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
    dist = np.sum(dist * dist, axis=2)
    return dist

