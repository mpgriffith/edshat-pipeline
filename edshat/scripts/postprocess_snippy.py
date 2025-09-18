#!/usr/bin/env python3

import os
from pathlib import Path
import subprocess
import argparse
import glob
import itertools
from Bio import SeqIO, AlignIO
import pandas as pd
import numpy as np
missing_chars = set(['?', '-', 'N', 'n', 'X'])


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('snippy_folder')
    parser.add_argument('-c', '--redo_core_pct', default=95)
    parser.add_argument('-p', '--prefix', default=None)
    parser.add_argument('-r', '--reference')
    return parser.parse_args()

def run_snp_dists(aln, output, linear=True):
    snp_dists_cmd = ['snp-dists', '-b', '-c', aln]
    with open(output, 'w') as f:
        subprocess.check_call(snp_dists_cmd, stdout=f)
    if linear:
        df = pd.read_csv(output, index_col=0)
        keep = np.triu(np.ones(df.shape)).astype('bool').reshape(df.size)
        df = df.stack()[keep]
        df.to_csv(output.replace('.csv', '_linear.csv'))
    

def run_snp_sites(aln, output):
    snp_sites_cmd = ['snp-sites', aln]
    with open(output, 'w') as f:
        subprocess.check_call(snp_sites_cmd, stdout=f)

def annotate_vcf(snippy_folder, reference, prefix=None):
    pass

def redo_core(snippy_folder, prefix, reference=None, min_proportion=.95):
    full_aln = snippy_folder / Path(f"{prefix}.full.aln")
    aln = AlignIO.read(full_aln, 'fasta')
    present_in = []
    for i in range(aln.get_alignment_length()):
        res = [aln[g][i] for g in range(len(aln))]
        present = len([1 for v in res if v not in missing_chars])
        present_in.append(present)
    l = [idx for idx, value in enumerate(present_in) if float(value) >= min_proportion ]
    ints = [list(g) for _,g in itertools.groupby(l, key=lambda n,c=itertools.count(): n-next(c))]
    ints = [[i[0], i[-1]] for i in ints]
    core_aln = []
    for i in ints:
        if not core_aln:
            core_aln = aln[:, i[0]:i[1]+1]
        else:
            core_aln += aln[:, i[0]: i[1]+1]
    AlignIO.write(core_aln, snippy_folder / Path(f"{prefix}.core.aln"), 'fasta')
    run_snp_dists(snippy_folder  / Path(f"{prefix}.core.aln"), snippy_folder / Path(f"{prefix}_SNP_matrix_core.csv"))
    run_snp_sites(snippy_folder / Path(f"{prefix}.core.aln"), snippy_folder / Path(f"{prefix}.snps.aln"))

def postprocess_snippy(snippy_folder, prefix=None):
    try:
        if not prefix:
            alns = glob.glob(os.path.join(snippy_folder, '*.aln'))
            alns = [l for l in alns if '.full.aln' not in l]
            if len(alns) != 1:
                raise Exception("snippy-folder doesn't contain a alignment")
            prefix = os.path.basename(alns[0]).rsplit('.', 1)[0]
            aln = alns[0]
        else:
            alns = glob.glob(os.path.join(snippy_folder, f'{prefix}.aln'))
            alns = [l for l in alns if '.full.aln' not in l]
            
            if len(alns) != 1:
                raise Exception
            aln = alns[0]
    except:
        raise Exception("snippy-folder doesn't contain a alignment")
    
    full_aln = Path(f"{snippy_folder}/{prefix}.full.aln")
    full_seqs = list(SeqIO.parse(full_aln, 'fasta'))
    full_seqs = [s for s in full_seqs if s.id != 'Reference']
    snp_seqs = list(SeqIO.parse(aln, 'fasta'))
    snp_seqs = [s for s in snp_seqs if s.id != 'Reference']
    SeqIO.write(full_seqs, full_aln, 'fasta')
    SeqIO.write(snp_seqs, aln, 'fasta')

    if len(full_seqs) >= 20:
        redo_core(snippy_folder, prefix, min_proportion=.95)
    
    else:
        run_snp_dists(aln, f"{snippy_folder}/{prefix}_SNP_matrix_core.csv")
        run_snp_sites(aln, f"{snippy_folder}/{prefix}.snps.aln")

    run_snp_dists(full_aln, f"{snippy_folder}/{prefix}_SNP_matrix_full.csv")


def main():
    args = parse_args()
    postprocess_snippy(Path(args.snippy_folder), prefix=args.prefix)
    if args.reference:
        annotate_vcf(Path(args.snippy_folder), Path(args.reference), prefix=args.prefix)


if __name__ == '__main__':
    main()