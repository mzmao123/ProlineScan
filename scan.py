#!/usr/bin/env python3

import sys
import os.path
import argparse

from Bio.PDB.DSSP import make_dssp_dict
from PDBtoDSSP import pdb_to_dssp
import util

from Bio.PDB import PDBParser

AA3to1 = {'ALA': 'A', 'CYS': 'C', 'ASP': 'D', 'GLU': 'E', 'PHE': 'F',
          'GLY': 'G', 'HIS': 'H', 'ILE': 'I', 'LYS': 'K', 'LEU': 'L',
          'MET': 'M', 'ASN': 'N', 'PRO': 'P', 'GLN': 'Q', 'ARG': 'R',
          'SER': 'S', 'THR': 'T', 'VAL': 'V', 'TRP': 'W', 'TYR': 'Y',
          'MSE': 'M', 'UNK': 'X'}

def parse():
    parser = argparse.ArgumentParser(description="Check proline compatbible positions in a given structure and the potential cost of mutate them to proline") 
    parser.add_argument('pdbfile', help='Input structure (in pdb format)')
    parser.add_argument('chains', help='Chain ID(s) for the chain to be checked for proline compatibility (example: "AB")')
    parser.add_argument('--prob_cutoff', help='Probability cutoff for the compatible mainchain phi-psi angles (default: 0.01)', default=0.01, type=float)
    parser.add_argument('--dist_collision', help='Atom distance threshold below which will be considered as a collision (default: 2.8)', default=2.8, type=float)
    parser.add_argument('--dmin_contact', help='Minimum distance to be considered as a contact (default: 3.0)', default=3.0, type=float)
    parser.add_argument('--dmax_contact', help='Maximum distance to be considered as a contact (default: 4.5)', default=4.5, type=float)
    parser.add_argument('--max_collisions', help='Maximum allowed collision between side chain of introduced proline and the rest part of the protein (default: 1)', default=1, type=int)
    parser.add_argument('-o', '--output', help='output selected fasta sequences (default: stdout).', default=sys.stdout, type=argparse.FileType('wt'))
    return parser.parse_args()

def main():
    para = parse()
    #
    print('#%s' %(' '.join(sys.argv)), file=para.output)
    print("ChainID\tResID\tAA\tASA\tSSE\tBCProb\tN_collision\tdN_contacts\tN_contacs_wt\tN_contacts_pro\tPHI\tPSI", file=para.output)

    # A few parameters
    prob_cutoff = para.prob_cutoff
    collision_th = para.dist_collision
    contact_dmin = para.dmin_contact
    contact_dmax = para.dmax_contact

    # Read structure
    pdb_parser = PDBParser(PERMISSIVE=1, QUIET=True)
    model = pdb_parser.get_structure('', para.pdbfile)[0]

    # SSE assignment by DSSP
    dssp_file = para.pdbfile + '.dssp'
    if not os.path.isfile(dssp_file):
        dssp_out = pdb_to_dssp(para.pdbfile, "https://www3.cmbi.umcn.nl/xssp/")
        with open(dssp_file, 'w') as tempfile:
            tempfile.write(dssp_out)
    dssp_dict = make_dssp_dict(dssp_file)[0]

    # Scanning
    for chain_id in para.chains:
        chain = model[chain_id]
        # Backbone compatibility
        proline_bc = util.BackboneCompatibility('data/rama8000-transpro.data')
        bbcompatible_sites = proline_bc.compatible_sites(chain, dssp_dict, prob_cutoff)
        # Sidechain compatibility
        proline_conformations = './data/P_-58_37.1.pdb ./data/P_-58_37.2.pdb ./data/P_-60_140.1.pdb ./data/P_-60_140.2.pdb'.split()
        rep_prolines = util.RepresentiveProlines(proline_conformations)
        for res_id, sse, asa, phi, psi, prob in bbcompatible_sites:
            res = model[chain_id][res_id]
            #if res.resname == 'PRO':  # skip side chain compatibility check for prolines
            #    continue
            contacts_wt = util.sidechain_contacted_atoms(model, chain_id, res_id, dist_range=(contact_dmin,contact_dmax))
            mutated_model, collisions = util.mutate_to_proline(model, chain_id, res_id, dssp_dict, rep_prolines, collision_th)
            contacts_pro = util.sidechain_contacted_atoms(mutated_model, chain_id, res_id, dist_range=(contact_dmin,contact_dmax))
            # output
            if len(collisions) <= para.max_collisions:
                aa = AA3to1[res.resname]
                n_collisions = len(collisions)
                n_contacts_wt = len(contacts_wt)
                n_contacts_pro = len(contacts_pro)
                if para.output != sys.stdout:
                    print(f'{chain_id}\t{res_id[1]}\t{aa}\t{asa}\t{sse}\t{prob:.3f}\t{n_collisions}\t{n_contacts_pro-n_contacts_wt}')
                print(f'{chain_id}\t{res_id[1]}\t{aa}\t{asa}\t{sse}\t{prob:.3f}\t{n_collisions}\t{n_contacts_pro-n_contacts_wt}\t{n_contacts_wt}\t{n_contacts_pro}\t{phi}\t{psi}', file=para.output)


if __name__ == "__main__":
    main()
