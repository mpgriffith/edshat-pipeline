#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import subprocess
from pathlib import Path
import argparse
import pandas as pd
import sys
import multiprocessing as mp

def parse_args():
    parser = argparse.ArgumentParser(description="Run ska compare for multiple SKF files.")
    parser.add_argument('reference_skf')
    parser.add_argument('query_skfs', nargs='+', help='List of query SKF files to compare against the reference.')
    parser.add_argument('-p', '--prev_results', default=None, type=str, help='Path to previous results CSV file to skip already processed queries.')
    parser.add_argument('-o', '--output', default=None)
    return parser.parse_args()

def ska_comp(reference_skf, query_skf):
    ska_cmd = ["ska", "compare", reference_skf, query_skf]
    try:
        result = subprocess.run(ska_cmd, capture_output=True, text=True, check=True)
        df = pd.read_csv(result.stdout, sep='\t')
        df['Reference'] = os.path.splitext(os.path.basename(reference_skf))[0]
        return df
    except subprocess.CalledProcessError as e:
        print(f"Error comparing {reference_skf} with {query_skf}: {e.stderr}")
        return None


def run_ska_compare(reference_skf, query_skfs, threads=8, prev_results=None, prev_results_time=None):
    """
    Run ska compare for a reference SKF file against a list of query SKF files.
    Skip comparisons for query files that already have data.

    :param reference_skf: Path to the reference SKF file.
    :param query_skfs: List of paths to query SKF files.
    :param results_dir: Directory to store comparison results.
    """
    queries_to_do = []
    for query_skf in query_skfs:
        query_skf_time = query_skf.stat().st_mtime
        if prev_results is not None:
            # Check if the query SKF has already been processed
            query_name = query_skf.stem
            if query_name in prev_results['Query'].values and  query_skf_time <= prev_results_time:
                continue
        queries_to_do.append(query_skf)

    with mp.Pool(threads) as pool:
        results = pool.starmap(ska_comp, [(reference_skf, query_skf) for query_skf in queries_to_do])

    # Combine results into a single DataFrame
    combined_results = pd.concat([res for res in results if res is not None], ignore_index=True)

    # Save combined results to the output file
    if prev_results is not None:
        combined_results = pd.concat([prev_results, combined_results], ignore_index=True)
    return combined_results

def get_query_skfs(query_skfs):
    if len(query_skfs) == 1 and query_skfs[0] == '*':
        query_skfs = list(Path('.').glob('*.skf'))
    elif len(query_skfs) == 1 and not query_skfs[0].endswith('.skf'):
        query_skfs = [Path(line.strip()) for line in open(query_skfs[0])]
    else:
        query_skfs = [Path(skf) for skf in query_skfs]
    return query_skfs

if __name__ == "__main__":
    # Example usage
    args = parse_args()
    query_skfs = get_query_skfs(args.query_skfs)
    reference_skf = Path(args.reference_skf)
    prev_results = args.prev_results
    if prev_results:
        prev_results_time = os.stat(prev_results).st_mtime
        if prev_results_time <= reference_skf.stat().st_mtime:
            prev_results = None
        else:
            prev_results = pd.read_csv(prev_results)
    else:
        prev_results = None
        prev_results_time = None
    results = run_ska_compare(args.reference_skf, query_skfs, prev_results=prev_results, prev_results_time=prev_results_time)
    output_file = args.output if args.output else f"{reference_skf.stem}_combined_results.csv"
    results.to_csv(output_file, index=False)
    print(f"Results saved to {output_file}")