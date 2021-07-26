#!/usr/bin/env python3

import sys
import argparse

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
    
if __name__ == "__main__":
    main()
