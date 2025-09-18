#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import subprocess
from pathlib import Path
import argparse
import pandas as pd
import sys
import multiprocessing as mp
import itertools
import math

def parse_args():
    parser = argparse.ArgumentParser(description="Run ska compare for multiple SKF files.")
    # parser.add_argument('reference_skf')
    parser.add_argument('query_skfs', nargs='+', help='List of query SKF files to compare against the reference.')
    parser.add_argument('-r', '--reference', default=None)
    parser.add_argument('-p', '--prev_results', default=None, type=str, help='Path to previous results CSV file to skip already processed queries.')
    parser.add_argument('-o', '--output', default=None)
    parser.add_argument('-t', '--threads', default=None)
    return parser.parse_args()

def ska_comp(reference_skf, query_skfs, suffix="_ska_compare.tsv"):
    ska_cmd = ["ska", "compare", '-q', reference_skf] + query_skfs
    print(reference_skf.stem, len(query_skfs))
    try:
        result = subprocess.run(ska_cmd, capture_output=True, text=True, check=True)
        lines = result.stdout.split('\n')
        df = pd.DataFrame([l.split('\t') for l in lines[1:]], columns = lines[0].split('\t'))
        df['Reference'] = reference_skf.stem
        prev_results = reference_skf.parent / f"{reference_skf.stem}{suffix}"
        if prev_results.is_file():
            prev_df = pd.read_csv(prev_results)
            df = pd.concat([prev_df, df], ignore_index=True)
        # df['Jaccard Index'] = (df['Matches'].astype(int)) / (df['Matches'].astype(int) + df['Kmers unique to Subject'].astype(int) + df['Kmers unique to Query'].astype(int))
        # df['Mash-like distance'] =  df['Jaccard Index'].map(lambda j: (-1/(2 * 15 +1))* math.log(2 * j /(1+j)) if j != 0 else  1)
        # df['SNP Distance'] = df['SNPs'] / df['Mathces']
        df.to_csv(prev_results, index=False)
        return df
    except subprocess.CalledProcessError as e:
        print(f"Error comparing {reference_skf} with : {e.stderr}")
        return None

def greedy_batch_queries(skf_files, missing_pairs):
    batches = []
    missing_by_query = {q: set() for q in skf_files}

    for a, b in missing_pairs:
        missing_by_query[a].add(b)
        missing_by_query[b].add(a)

    used_pairs = set()
    missing_pairs = set(missing_pairs)
    print(missing_pairs)
    print(missing_by_query)
    while missing_pairs - used_pairs:
#        for q in missing_by_query:
#            if len(missing_by_query[q]) == 0:
        print(missing_by_query)
        print(missing_pairs - used_pairs)                
        best_query = max(missing_by_query, key=lambda q: len(missing_by_query[q]))
        targets = list(missing_by_query[best_query])
        
        if not targets:
            del missing_by_query[best_query]
            continue

        batch = (best_query, targets)
        batches.append(batch)

        for t in targets:
            used_pairs.add(tuple(sorted((best_query, t))))
            if t in missing_by_query:
                missing_by_query[t].discard(best_query)
        missing_by_query[best_query].clear()

    return batches

def run_ska_compare(query_skfs, batches, threads=8):
    """
    Run ska compare for a reference SKF file against a list of query SKF files.
    Skip comparisons for query files that already have data.

    :param reference_skf: Path to the reference SKF file.
    :param query_skfs: List of paths to query SKF files.
    :param results_dir: Directory to store comparison results.
    """
   # combs = list(itertools.combinations(query_skfs.keys(), 2))
    skf_batches = []
    for s, targets in batches:
        skf_batches += [(query_skfs[s], [query_skfs[t] for t in targets[i: min(len(targets), i+ 50)]]) for i in range(0, len(targets), 50)]
    print('ska compares to do', len(skf_batches))
    with mp.Pool(threads) as pool:
        results = pool.starmap(ska_comp, skf_batches)
    # Combine results into a single DataFrame
    results = [res for res in results if res is not None]
    if not results:
        return pd.DataFrame()
    combined_results = pd.concat([res for res in results if res is not None], ignore_index=True)

    # Save combined results to the output file
    return combined_results

def get_query_skfs(query_skfs, prev_results_suffix='_ska_compare.tsv'):
    if len(query_skfs) == 1 and query_skfs[0] == '*':
        query_skfs = list(Path('.').glob('*.skf'))
    elif len(query_skfs) == 1 and not query_skfs[0].endswith('.skf'):
        query_skfs = [Path(line.strip()) for line in open(query_skfs[0])]
    else:
        query_skfs = [Path(skf) for skf in query_skfs]
    query_skfs = {skf.stem: skf for skf in query_skfs}

    prev_results = {}
    for q, f in query_skfs.items():
        s = f.stem
        prf = f.parent / f'{s}{prev_results_suffix}'
        if prf.is_file():
            prev_results[s] = prf
    return query_skfs, prev_results

def get_existing_data(query_skfs, prev_files):
    seen = set()
    prev_dfs ={}
    for s, f in prev_files.items():
        df = pd.read_csv(f)
        print(df)
        prev_results_time = f.stat().st_mtime
        ref_time = query_skfs[s].stat().st_mtime
        if ref_time > prev_results_time:
            continue
        keep_rows = []
        for idx, row in df.iterrows():
            query = row['Subject']
            if query not in query_skfs:
                continue
            query_time = query_skfs[query].stat().st_mtime
            if query_time > prev_results_time:
                continue
            pair = tuple(sorted((row['Reference'], row['Subject'])))
            seen.add(pair)
            keep_rows.append(idx)
        df = df.loc[keep_rows]
        df.to_csv(f, index=False)
        prev_dfs[s] = df
    return seen, prev_dfs

def df_column_switch(df, col1, col2):
    #i = list(df.columns)
   # a, b = i.index(column1), i.index(column2)i
    cols = df.columns
    col1_idx = cols.get_loc(col1)
    col2_idx = cols.get_loc(col2)
    df[[cols[col1_idx], cols[col2_idx]]] = df[[cols[col2_idx], cols[col1_idx]]]
    #i[b], i[a] = i[a], i[b]
    #df = df[i]
    return df

def add_switched_rows(result_df):
    switch_df = df_column_switch(result_df, 'Reference', 'Subject')
    switch_df = df_column_switch(switch_df, 'Kmers unique to Subject', 'Kmers unique to Query')
    switch_df = df_column_switch(switch_df, '% kmers in Subject matching', '% kmers in Query matching')
    switch_df = df_column_switch(switch_df, '%ID of Subject kmers', '%ID of Query kmers')
    switch_df = df_column_switch(switch_df, 'Ns in Subject', 'Ns in Query')
    df = pd.concat([result_df, switch_df])
    return df

def switch_reference_query(df):
    switch_cols = {'Reference': 'Subject', 
                    'Kmers unique to Subject': 'Kmers unique to Query', 
                    '% kmers in Subject matching': '% kmers in Query Matching', 
                    '%ID of Subject kmers': '%ID of Query kmers',
                    'Ns in Subject': 'Ns in Query'}
    col_order = df.columns.tolist()
    df = df.rename(columns={**switch_cols, **{v:k for k,v in switch_cols.items()}})
    df = df.reindex(col_order, axis=1)
    return df


def write_results_files(all_results, query_skfs, results_file):
    isos = [Path(f).stem for f in query_skfs]
    combs = list(itertools.product(isos, isos))
    all_results = pd.concat([all_results, switch_reference_query(all_results)])
    just_comb_results = all_results[pd.Series(list(zip(all_results['Reference'], all_results['Subject']))).isin(combs)]
    just_comb_results.to_csv(results_file, index=False)


    for s in query_skfs:
        df = all_results[(all_results['Reference'] == s)]
       # df_query = all_results[(all_results['Subject'] == s)]
       # df_column_switch(df_query, 'Reference', 'Subject')
       # df_column_switch(df_query, 'Kmers unique to Subject', 'Kmers unique to Query')
       # df_column_switch(df_query, '% kmers in Subject matching', '% kmers in Query matching')
       # df_column_switch(df_query, '%ID of Subject kmers', '%ID of Query kmers')
       # df_column_switch(df_query, 'Ns in Subject', 'Ns in Query')
       # df = pd.concat([df, df_query])
        df = df.drop_duplicates(subset=['Reference', 'Subject'])
        df.to_csv(query_skfs[s].parent / f'{s}_ska_compare.tsv', index=False)
    
def get_previous_results(prev_results):
    if not prev_results:
        return pd.DataFrame()
    dfs = []
    for f in prev_results.values():
        df = pd.read_csv(f)
        dfs.append(df)
    if len(dfs) > 0:
        return pd.concat(dfs, ignore_index=True)
    else:
        return pd.DataFrame()


if __name__ == "__main__":
    # Example usage
    args = parse_args()
    threads = int(args.threads) if args.threads else 8
    query_skfs, prev_results = get_query_skfs(args.query_skfs)
    existing_pairs, prev_dfs = get_existing_data(query_skfs, prev_results)
    all_pairs_needed = list(itertools.combinations_with_replacement(sorted(query_skfs.keys()), 2))
    pairs_to_do = [pair for pair in all_pairs_needed if pair not in existing_pairs]
    all_combs = list(itertools.product(query_skfs.keys(), repeat=2))
    batches = greedy_batch_queries(query_skfs, pairs_to_do)
    #reference_skf = Path(args.reference)
    #reference_skf = Path(args.reference_skf)
    # prev_results = args.prev_results
    # if prev_results:
    #     prev_results_time = os.stat(prev_results).st_mtime
    #     if prev_results_time <= reference_skf.stat().st_mtime:
    #         prev_results = None
    #     else:
    #         prev_results = pd.read_csv(prev_results)
    # else:
    #     prev_results = None
    #     prev_results_time = None
    results = run_ska_compare(query_skfs, batches, threads=threads)
    if len(prev_dfs) > 0:
        prev_results_df = pd.concat(prev_dfs.values(), ignore_index=True)
    else:
        prev_results_df = pd.DataFrame()
    if len(prev_results_df) > 0:
        results = pd.concat([results, prev_results_df])

    output_file = args.output if args.output else "ska_combined_results.csv"

    write_results_files(results, query_skfs, output_file)
    #results.to_csv(output_file, index=False)
    print(f"Results saved to {output_file}")
