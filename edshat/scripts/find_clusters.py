#!/usr/bin/env python3

import scipy.cluster.hierarchy as sch
from skbio import DistanceMatrix
import argparse
import sys
import itertools
from string import ascii_uppercase
import collections
import csv

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('matrix')
    parser.add_argument('-m', '--method', default='single')
    parser.add_argument('-c', '--criterion', default='distance')
    parser.add_argument('-t', '--cutoff', default=15, type=int)
    parser.add_argument('-i', '--include_singletons', action='store_true', default=False)
    parser.add_argument('-o', '--output', default=None)
    parser.add_argument('-f', '--fastani_matrix', action='store_true', default=False)
    parser.add_argument('-s', '--similarity_matrix', action='store_true', default=False)
    return parser.parse_args()


def make_matrix(matrix_file, fastani=False, similarity=False, make_symmetrical=False):
    delimiter = ',' if not fastani else '\t'
    with open(matrix_file) as f:
        reader = csv.reader(f, delimiter=delimiter)
        rows = [r for r in reader if r]
    if fastani:
        data = []
        m = {}
        for r in rows:
            # print(r)
            s1 = r[0].split('/')[-1].split('.')[0]
            s2 = r[1].split('/')[-1].split('.')[0]
            v = float(r[2])
            if s1 not in m: m[s1] = {}
            if s2 not in m: m[s2] = {}
            m[s1][s2] = 100.0 - v
            m[s2][s1] = 100.0 - v
            m[s1][s1] = 0.0
            m[s2][s2] = 0.0
        ids = sorted(m.keys())
        for s1 in ids:
            r = []
            for s2 in ids:
                v = m[s1][s2] if s2 in m[s1] else 100.0
                r.append(v)
            data.append(r)
    elif similarity:
        data = []
        ids = [str(v) for v in rows[0][1:]]
        m = {}
        for r in rows[1:]:
            m[r[0]] = {ids[i]: v for i, v in enumerate(r[1:])}
        for i1, i2 in itertools.combinations_with_replacement(ids, 2):
            v1 = float(m[i1][i2]) if m[i1][i2] else 0
            v2 = float(m[i2][i1]) if m[i2][i1] else 0
            v = 100 - min([v1, v2])
            m[i1][i2] = v
            m[i2][i1] = v
            m[i1][i1] = 0
            m[i2][i2] = 0
        for i1 in ids:
            r = []
            for i2 in ids:
                v = m[i1][i2]
                r.append(v)
            data.append(r)  
    else:        
        data = [r[1:] for r in rows[1:]]
        ids = rows[0][1:]
    #print(data)
    dm = DistanceMatrix(data, ids)
    return dm

def make_clusters(dm, method='single', criterion='distance', cutoff=15, include_singletons=False):
    #fix matrix file
    ids = dm.ids
    Z = sch.linkage(dm.condensed_form(), method=method)
    cluster_list = sch.fcluster(Z, criterion=criterion, t=cutoff)
    clusters = []
    all_clusters_count = 0 if include_singletons else 1
    cluster_counts = {c: count for c, count in collections.Counter(cluster_list).items() if count > all_clusters_count}
    for c in sorted(cluster_counts, key=lambda i: cluster_counts[i], reverse=True):
        inds = [i for i, clus in enumerate(cluster_list) if clus == c]
        samples = [ids[i] for i in inds]
        clusters.append(samples)
    return clusters


def iter_all_strings():
    size = 1
    while True:
        for s in itertools.product(ascii_uppercase, repeat=size):
            yield ''.join(s)
        size += 1
    
def main():
    args = parse_args()
    matrix_file = args.matrix
    method = args.method
    crit = args.criterion
    t = args.cutoff
    if args.similarity_matrix:
        t = 100 - t
    matrix = make_matrix(matrix_file, fastani=args.fastani_matrix, similarity=args.similarity_matrix)
    clusters = make_clusters(matrix, method=method, criterion=crit, cutoff=t, include_singletons=args.include_singletons)

    output = open(args.output, 'w') if args.output else sys.stdout
    writer = csv.writer(output)
    clusters = sorted(clusters, key=lambda c: sorted(c)[0])
    for s, c in zip(iter_all_strings(), clusters):
        for i in c:
            writer.writerow([i, s])



if __name__ == "__main__":
    main()
