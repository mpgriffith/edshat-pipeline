#!/usr/bin/env python3

import ete3
import argparse
import sys
import os

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('tree')
    parser.add_argument('-o', '--output', default=None)
    return parser.parse_args()

def main():
    args = parse_args()
    ete3.parser.newick.set_float_format('%0.20f')
    tree = ete3.Tree(args.tree, format=1)
    if len(tree) > 2:
        tree.set_outgroup(tree.get_midpoint_outgroup())
    if args.output:
        tree.write(outfile=args.output, format=1)
    else:
        print(tree.write(format=1))

if __name__ == "__main__":
    main()