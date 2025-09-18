import argparse
import yaml
import pandas as pd
import sys
import tempfile
from pathlib import Path
import re


from . import config_manager
from . import runner

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_arugment('--config', default=None, nargs='+')
    parser.add_argument('--config_file', default=None, nargs='+')


    run_parser = argparse.ArgumentParser(add_help=False)
    run_parser.add_argument('-o', '--output_dir', default=None)
    run_parser.add_argument('-j', '--threads', default=12)
    run_parser.add_argument('-i', '--samples', nargs='+')
    run_parser.add_argument('--reads_dir', default=None)
    run_parser.add_argument('-s', '--sample_dir', default=None)
    run_parser.add_argument('-n', '--name', default=None)
    run_parser.add_argument('--snakemake_options', default=None, nargs='+')
    run_parser.add_argument('--dry-run', action='store_true', default=False)

    subparsers = parser.add_subparsers(dest='cmd')
    
    config = subparsers.add_parser('config', help='View or edit the user configuration.')
    config.add_argument('variable', nargs='?', default=None, 
                        help="The configuration variable to set (e.g., 'database_paths;kraken'). Use ';' for nested keys. Leave blank to show current config.")
    config.add_argument('value', nargs='?', default=None, 
                        help="The value to set for the variable.")


    run = subparsers.add_parser('run', parents=[run_parser])
    run.add_argument('--only_per_sample', action='store_true', default=False)
    run.add_argument('-r', '--reference', default=None)
    run.add_argument('--reference_metric', default='N50', help='Metric to choose reference')
    run.add_argument('--species', default=None, help='Expected Species (or function/mapping)')
    run.add_argument('--targets', default=None, nargs='+')
    run.add_argument('--split_by', default=None)
    run.add_argument('--use_all', default=False)
    run.add_argument('-d', '--prev_data', default=None)
    run.add_argument('--filter', nargs='+')
    run.add_argument('--logic', default='AND', choices=['AND', 'OR'])
    
    run.add_argument('--output-columns', nargs='+')
    run.add_argument('--previous_reads', default=None)
    run.add_argument('--duplicate-columns', nargs='+')

    download = subparsers.add_parser('download')
    download.add_argument('--db_dir', default=None)
    download.add_argument('--kraken_db', default=None)
    download.add_argument('--kraken_size', default=None)
    download.add_argument('--dry-run', action='store_true', default=False)
    download.add_argument('-t', '--threads', default=4)

    return parser.parse_args()
    
def get_species(samples, species, config):
    if type(species) == str and Path(species).isfile():
        species = yaml.load(open(species), Loader=yaml.FullLoader)
    else:
        species = config['species_mapping']
    
    print(species)
    def species_map(name):
        
        match = re.match(species['mapping'], name)
        print(name, match, species['mapping'].strip("'").strip('"'))
        if match:
            species_name = species['species'].get(match.group(1), None)
        else:
            species_name = None
        return species_name
    species_samples = samples['Sample'].map(species_map)
    return species_samples

def find_reads_for_sample(sample_name: str, base_reads_dir: Path) -> dict:
    """
    Finds paired-end FASTQ files for a given sample name in a base directory.

    Args:
        sample_name: The name of the sample.
        base_reads_dir: The top-level directory containing sample subfolders.

    Returns:
        A dictionary with sample info if reads are found, otherwise None.
    """
    sample_dir = base_reads_dir / sample_name
    if not sample_dir.is_dir():
        print(f"⚠️ Warning: Could not find directory for sample '{sample_name}' in {base_reads_dir}")
        return None

    # Use a more specific glob pattern first, then a general one as a fallback
    reads = sorted(list(sample_dir.glob('*_R[12]_*.f*q.gz')))
    if len(reads) < 2:
        reads = sorted(list(sample_dir.glob('*_[12].f*q.gz')))
    if len(reads) < 2:
        reads = sorted(list(sample_dir.glob('*.f*q.gz')) + list(sample_dir.glob('*.fastq')))
    
    if len(reads) >= 2:
        return {'Sample': sample_name, 'Read1': str(reads[0]), 'Read2': str(reads[1])}
    else:
        print(f"⚠️ Warning: Found < 2 read files for sample '{sample_name}'. Skipping.")
        return None
    
def prepare_sample_sheet(samples_input: list, reads_dir: str, species: str = None, config: dict = None) -> pd.DataFrame:
    """
    Prepares a standardized sample sheet DataFrame from various input types.

    Supports four methods for specifying samples:
    1.  A CSV/TSV sample sheet file with 'Sample', 'Read1', 'Read2' columns.
    2.  A directory containing one subdirectory per sample.
    3.  A plain text file listing one sample name per line.
    4.  A list of sample names provided directly on the command line.
    """
    base_reads_dir = Path(reads_dir) if reads_dir else Path('Reads')
    sample_data_rows = []

    # Case: No input provided, so scan the default reads directory
    if not samples_input:
        print(f"ℹ️ No samples provided via -i/--samples. Scanning default directory: '{base_reads_dir}'")
        if not base_reads_dir.is_dir():
            print(f"❌ Error: Default reads directory '{base_reads_dir}' not found. Please specify samples or a reads directory.")
            sys.exit(1)
        # Treat this case as providing a directory
        samples_input = [str(base_reads_dir)]

    input_path = Path(samples_input[0])

    # --- Determine Input Type ---
    
    # 1. Input is a CSV/TSV sample sheet file
    if len(samples_input) == 1 and input_path.is_file() and input_path.suffix.lower() in ['.csv', '.tsv']:
        print(f"✅ Detected sample sheet file: {input_path}")
        separator = '\t' if input_path.suffix.lower() == '.tsv' else ','
        samples_df = pd.read_csv(input_path, sep=separator)
        # Basic validation
        required_cols = ['Sample', 'Read1']
        if not all(col in samples_df.columns for col in required_cols):
            print(f"❌ Error: Sample sheet must contain at least the columns: {required_cols}")
            sys.exit(1)

    # 2. Input is a directory containing sample subdirectories
    elif len(samples_input) == 1 and input_path.is_dir():
        print(f"✅ Detected directory of samples: {input_path}")
        for sample_dir in input_path.iterdir():
            if sample_dir.is_dir():
                sample_info = find_reads_for_sample(sample_dir.name, input_path)
                if sample_info:
                    sample_data_rows.append(sample_info)
        samples_df = pd.DataFrame(sample_data_rows)

    # 3. Input is a plain text file with sample names
    elif len(samples_input) == 1 and input_path.is_file():
        print(f"✅ Detected file with sample list: {input_path}")
        with open(input_path, 'r') as f:
            sample_names = [line.strip() for line in f if line.strip()]
        for name in sample_names:
            sample_info = find_reads_for_sample(name, base_reads_dir)
            if sample_info:
                sample_data_rows.append(sample_info)
        samples_df = pd.DataFrame(sample_data_rows)

    # 4. Input is a list of sample names from the command line
    else:
        for name in samples_input:
            sample_info = find_reads_for_sample(name, base_reads_dir)
            if sample_info:
                sample_data_rows.append(sample_info)
        samples_df = pd.DataFrame(sample_data_rows)

    # --- Final Processing ---

    if samples_df.empty:
        print("❌ Error: No valid samples were found from the provided input. Exiting.")
        sys.exit(1)

    # Add species information and set the index
    samples_df['species'] = species if isinstance(species, str) else get_species(samples_df, species, config)
    samples_df = samples_df.set_index('Sample', drop=False)
    
    print("\n--- Successfully loaded samples ---")
    print(samples_df[['Sample', 'Read1', 'Read2', 'species']].head())
    print("---------------------------------\n")
    return samples_df

def filter_prev_samples(prev_samples, filters, logic='AND'):
    if not filters:
        return prev_samples
    filter_funcs = []
    for f in filters:
        if '>=' in f:
            col, val = f.split('>=')
            filter_funcs.append(lambda df: df[col] >= val)
        elif '>' in f:
            col, val = f.split('>')
            filter_funcs.append(lambda df: df[col] > val)
        elif '<=' in f:
            col, val = f.split('>')
            filter_funcs.append(lambda df: df[col] <= val)
        elif '<' in f:
            col, val = f.split('<')
            filter_funcs.append(lambda df: df[col] < val)
        elif '=' in f:
            col, val = f.split('=')
            filter_funcs.append(lambda df: df[col] == val)
        else:
            raise ValueError(f"Invalid filter format: {f}")
    
    if logic == 'AND':
        combined_filter = lambda df: pd.concat([func(df) for func in filter_funcs], axis=1).all(axis=1)
    elif logic == 'OR':
        combined_filter = lambda df: pd.concat([func(df) for func in filter_funcs], axis=1).any(axis=1)
    else:
        raise ValueError(f"Invalid logic: {logic}")
    
    filtered_samples = prev_samples[combined_filter(prev_samples)]
    return filtered_samples

def handle_run_command(args):
    """
    Handles the integrated 'run' command, supporting both simple and 
    surveillance-style workflows.
    """
    # --- 1. Initial Setup and Config Loading ---
    config = config_manager.get_config()
    
    if args.name:
        config['run_name'] = args.name
        run_output_dir = Path(config.get('output_dir', '.')) / args.name
        config['isolate_dir'] = run_output_dir / 'isolates'
        config['set_dir'] = run_output_dir / 'comparative'
    
    # --- 2. Per-Sample Analysis (Always Runs) ---
    print("--- Stage 1: Running Per-Sample Analysis ---")
    
    new_samples_df = prepare_sample_sheet(args.samples, args.reads_dir, args.species, config=config)
    with tempfile.NamedTemporaryFile(mode='w', suffix=".csv", delete=False) as tmp:
        new_samples_df.to_csv(tmp.name)
        per_sample_config = config.copy()
        per_sample_config['samples'] = tmp.name

    runner.run_snakemake(
        snakefile=config_manager.WORKFLOW_DIR / 'Snakefile',
        config=per_sample_config,
        threads=args.threads,
        workdir=Path(config.get('workdir', '.')),
        dry_run=args.dry_run
    )

    # --- 3. Conditional Comparative Analysis ---
    if args.only_per_sample:
        print("\n✅ Per-sample analysis complete. Skipping comparative step as requested.")
        return

    print("\n--- Stage 2: Running Comparative Analysis ---")
    
    # Get the list of explicitly NEW samples
    new_sample_ids = list(new_samples_df.index)
    
    # The scope for comparison starts with the new samples.
    final_sample_list = new_sample_ids[:] # Start with a copy
    
    if args.prev_data:
        print(f"Including previous data from {args.prev_data}")
        prev_metrics_df = pd.read_csv(args.prev_data)
        final_sample_list.extend(prev_metrics_df['SpecID'].tolist())
        if Path(args.prev_data).is_file():
            if args.prev_data.endswith('xlsx') or args.prev_data.endswith('xls'):
                other_samples = pd.read_excel(args.prev_data, index_col=0)
            else:
                other_samples = pd.read_csv(args.prev_data, index_col=0)
            filtered_samples = filter_prev_samples(other_samples, args.filter, args.logic)
            print(filtered_samples)
            print(samples)
            all_samples = pd.concat([samples, filtered_samples])
            all_samples = all_samples[~all_samples.index.duplicated(keep='last')]
            final_sample_list = list(all_samples.index)

    comparative_config = config.copy()
    comparative_config['new_samples'] = new_sample_ids

    final_sample_list = sorted(list(set(final_sample_list)))
    final_sample_df = prepare_sample_sheet(final_sample_list, args.reads_dir, args.species, config=config)
    with tempfile.NamedTemporaryFile(mode='w', suffix=".csv", delete=False) as tmp:
        final_sample_df.to_csv(tmp.name)
        comparative_config['analysis_samples'] = tmp.name

    comparative_config['isolate_dir'] = Path(config.get('workdir', '.')) / 'isolates'
    comparative_config['previous_data'] = args.prev_data
    print(f"Running comparative analysis on {len(final_sample_list)} total samples.")
    runner.run_snakemake(
        snakefile=config_manager.WORKFLOW_DIR / 'comparative.smk',
        config=comparative_config,
        threads=args.threads,
        workdir=Path(config.get('workdir', '.')),
        dry_run=args.dry_run
    )

def handle_config_command(args):
    # If both variable and value are provided, set the config
    if args.variable and args.value:
        config_manager.set_config(args)
    # If only a variable is provided, it's an error
    elif args.variable and not args.value:
        print("❌ Error: A value must be provided when specifying a variable.")
    # If no arguments are provided, show the current config
    
    else:
        config = config_manager.get_config()
        print("# --- Current EDSHAT Configuration ---")
        print("# Default values are merged with your user-specific settings.")
        print(f"# User config file: {USER_CONFIG_FILE}\n")
        print(yaml.dump(config, default_flow_style=False, sort_keys=False))


def handle_download_command(args):
    config = config_manager.get_config(args)
    db_dir = args.db_dir if args.db_dir else Path(config_manager.SCRIPT_DIR) /config['database_dir']


    config['database_dir'] = str(db_dir)
    kraken_db = args.kraken_db if args.kraken_db else 'PlusPFP'
    kraken_size = '-' +  str(args.kraken_size) if args.kraken_size else ''
    config['kraken_db'] = kraken_db + kraken_size
    runner.run_snakemake(snakefile=config_manager.WORKFLOW_DIR / 'rules' / 'download_db.smk',
                         config=config,
                         threads=args.threads,
                         workdir=Path(config.get('workdir', '.')),
                         dry_run=args.dry_run,
                         targets=['download']
        )
    
    if args.db_dir:
        args.variable = 'database_dir'
        args.value = Path(args.db_dir).resolve()
        config_manager.set_config(args)

def main():
    args = parse_args()
    if args.cmd == 'run':
        handle_run_command(args)
    elif args.cmd == 'config':
        handle_config_command(args)
    elif args.cmd == 'combine':
        handle_combine_command(args)
    elif args.cmd == 'download':
        handle_download_command(args)


if __name__ == '__main__':
    main()