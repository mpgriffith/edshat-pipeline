import argparse
import yaml
import pandas as pd
import sys
import tempfile
import datetime
from pathlib import Path
import re
import logging

from . import config_manager
from . import runner

logging.basicConfig(level=logging.INFO, format='%(levelname)s: %(message)s', stream=sys.stdout)

def parse_args():
    parser = argparse.ArgumentParser(
        description="EDS-HAT: A pipeline for hospital-acquired transmission analysis."
    )
    parser.add_argument('--config', default=None, nargs='+')
    parser.add_argument('--config_file', default=None, nargs='+')


    run_parser = argparse.ArgumentParser(add_help=False)
    run_parser.add_argument('-o', '--output_dir', default=None)
    run_parser.add_argument('-j', '--threads', default=12)
    run_parser.add_argument('-i', '--samples', nargs='+')
    run_parser.add_argument('--reads_dir', default=None)
    run_parser.add_argument('-s', '--sample_dir', default=None)
    run_parser.add_argument('-n', '--name', default=None)
    run_parser.add_argument('--snakemake_options', default=None, nargs='+')
    run_parser.add_argument('--dry-run', action='store_true', default=False, help="Perform a dry run of the pipeline.")

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
    """Determines the species for each sample based on a mapping configuration."""
    species_config = config['species_mapping']
    if isinstance(species, str) and Path(species).is_file():
        with open(species, 'r') as f:
            species_config = yaml.safe_load(f)

    def species_map(name):
        
        match = re.match(species_config.get('mapping', ''), name)
        if match:
            species_name = species_config['species'].get(match.group(1), None)
            # Assumes the first captured group is the species slug
            species_name = species_config.get('species', {}).get(match.group(1), None)
        else:
            species_name = None
        return species_name or "Unknown"
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
        logging.warning(f"Could not find directory for sample '{sample_name}' in {base_reads_dir}")
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
        logging.warning(f"Found < 2 read files for sample '{sample_name}'. Skipping.")
        return None

def _prepare_from_sheet(file_path: Path) -> pd.DataFrame:
    """Load samples from a CSV or TSV file."""
    logging.info(f"Detected sample sheet file: {file_path}")
    separator = '\t' if file_path.suffix.lower() == '.tsv' else ','
    samples_df = pd.read_csv(file_path, sep=separator)
    required_cols = ['Sample', 'Read1']
    if not all(col in samples_df.columns for col in required_cols):
        logging.error(f"Sample sheet must contain at least the columns: {required_cols}")
        sys.exit(1)
    return samples_df

def _prepare_from_dir(directory: Path) -> pd.DataFrame:
    """Scan a directory for sample subdirectories and find their reads."""
    logging.info(f"Detected directory of samples: {directory}")
    sample_data_rows = []
    for sample_dir in directory.iterdir():
        if sample_dir.is_dir():
            sample_info = find_reads_for_sample(sample_dir.name, directory)
            if sample_info:
                sample_data_rows.append(sample_info)
    return pd.DataFrame(sample_data_rows)

def _prepare_from_list_file(file_path: Path, base_reads_dir: Path) -> pd.DataFrame:
    """Load sample names from a text file and find their reads."""
    logging.info(f"Detected file with sample list: {file_path}")
    sample_data_rows = []
    with open(file_path, 'r') as f:
        sample_names = [line.strip() for line in f if line.strip()]
    for name in sample_names:
        sample_info = find_reads_for_sample(name, base_reads_dir)
        if sample_info:
            sample_data_rows.append(sample_info)
    return pd.DataFrame(sample_data_rows)

def _prepare_from_cli_list(sample_names: list, base_reads_dir: Path) -> pd.DataFrame:
    """Find reads for a list of sample names given on the command line."""
    logging.info("Processing sample names from command line.")
    sample_data_rows = []
    for name in sample_names:
        sample_info = find_reads_for_sample(name, base_reads_dir)
        if sample_info:
            sample_data_rows.append(sample_info)
    return pd.DataFrame(sample_data_rows)

def prepare_sample_sheet(
    samples_input: list,
    reads_dir: str,
    species: str,
    config: dict
) -> pd.DataFrame:
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
    samples_df = pd.DataFrame()

    # Case: No input provided, so scan the default reads directory
    if not samples_input:
        logging.info(f"No samples provided via -i/--samples. Scanning default directory: '{base_reads_dir}'")
        if not base_reads_dir.is_dir():
            logging.error(f"Default reads directory '{base_reads_dir}' not found. Please specify samples or a reads directory.")
            sys.exit(1)
        samples_df = _prepare_from_dir(base_reads_dir)
    else:
        input_path = Path(samples_input[0])

        # --- Determine Input Type ---
        if len(samples_input) == 1 and input_path.is_file() and input_path.suffix.lower() in ['.csv', '.tsv']:
            samples_df = _prepare_from_sheet(input_path)
        elif len(samples_input) == 1 and input_path.is_dir():
            samples_df = _prepare_from_dir(input_path)
        elif len(samples_input) == 1 and input_path.is_file():
            samples_df = _prepare_from_list_file(input_path, base_reads_dir)
        else:
            samples_df = _prepare_from_cli_list(samples_input, base_reads_dir)

    # --- Final Processing ---

    if samples_df.empty:
        logging.error("No valid samples were found from the provided input. Exiting.")
        sys.exit(1)

    # Add species information and set the index
    if species and not Path(species).is_file():
        samples_df['species'] = species
    else:
        samples_df['species'] = get_species(samples_df, species, config)

    samples_df = samples_df.set_index('Sample', drop=False)
    
    logging.info("Successfully loaded samples:")
    logging.info("\n" + samples_df[['Sample', 'Read1', 'Read2', 'species']].head().to_string())
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

def _run_per_sample_analysis(args, config, new_samples_df):
    """Run the per-sample part of the Snakemake workflow."""
    logging.info("--- Stage 1: Running Per-Sample Analysis ---")
    with tempfile.NamedTemporaryFile(mode='w', suffix=".csv", delete=False, newline='') as tmp:
        new_samples_df.to_csv(tmp.name, index=False)
        per_sample_config = config.copy()
        per_sample_config['samples'] = tmp.name

    runner.run_snakemake(
        snakefile=config_manager.WORKFLOW_DIR / 'Snakefile',
        config=per_sample_config,
        threads=args.threads,
        workdir=Path(config.get('workdir', '.')),
        dry_run=args.dry_run
    )
    return Path(tmp.name)

def _run_comparative_analysis(args, config, new_samples_df, per_sample_sheet_path):
    """Run the comparative part of the Snakemake workflow."""
    logging.info("\n--- Stage 2: Running Comparative Analysis ---")
    config['run_name'] = args.name if args.name else datetime.datetime.now().strftime("%Y%m%d_%H%M%S")

    # Get the list of explicitly NEW samples
    new_sample_ids = list(new_samples_df.index)
    
    # The scope for comparison starts with the new samples.
    final_sample_list = new_sample_ids[:] # Start with a copy
    
    if args.prev_data:
        logging.info(f"Including previous data from {args.prev_data}")
        if Path(args.prev_data).is_file():
            if args.prev_data.endswith('xlsx') or args.prev_data.endswith('xls'):
                other_samples = pd.read_excel(args.prev_data, index_col=0)
            else:
                other_samples = pd.read_csv(args.prev_data, index_col=0)
            filtered_samples = filter_prev_samples(other_samples, args.filter, args.logic)
            logging.info(f"Found {len(filtered_samples)} samples in previous data after filtering.")
            all_samples = pd.concat([new_samples_df, filtered_samples])
            all_samples = all_samples[~all_samples.index.duplicated(keep='last')]
            final_sample_list = list(all_samples.index)

    comparative_config = config.copy()
    comparative_config['isolate_dir'] = str(Path(config.get('workdir', '.')) / 'isolates')
    comparative_config['set_dir'] = str(Path(config.get('workdir', '.')) / 'sets' / config['run_name'])
    comparative_config['workdir'] = str(Path(config.get('workdir', '.')))

    comparative_config['new_samples'] = new_sample_ids
    comparative_config['analysis_samples'] = final_sample_list

    final_sample_list = sorted(list(set(final_sample_list)))
    final_sample_df = prepare_sample_sheet(final_sample_list, args.reads_dir, args.species, config=config)
    with tempfile.NamedTemporaryFile(mode='w', suffix=".csv", delete=False, newline='') as tmp:
        final_sample_df.to_csv(tmp.name, index=False)
        comparative_config['samples'] = tmp.name
    

    comparative_config['previous_data'] = args.prev_data
    logging.info(f"Running comparative analysis on {len(final_sample_list)} total samples.")
    runner.run_snakemake(
        snakefile=config_manager.WORKFLOW_DIR / 'comparative.smk',
        config=comparative_config,
        threads=args.threads,
        workdir=Path(config.get('workdir', '.')),
        dry_run=args.dry_run
    )

def handle_run_command(args, config):
    """
    Handles the integrated 'run' command, supporting both simple and
    surveillance-style workflows.
    """
    if args.name:
        config['run_name'] = args.name
        run_output_dir = Path(config.get('output_dir', '.')) / args.name
        config['isolate_dir'] = run_output_dir / 'isolates'
        config['set_dir'] = run_output_dir / 'comparative'

    new_samples_df = prepare_sample_sheet(args.samples, args.reads_dir, args.species, config=config)
    per_sample_sheet_path = _run_per_sample_analysis(args, config, new_samples_df)

    if not args.only_per_sample:
        _run_comparative_analysis(args, config, new_samples_df, per_sample_sheet_path)
    else:
        logging.info("\nPer-sample analysis complete. Skipping comparative step as requested.")

    # Clean up the temporary file
    per_sample_sheet_path.unlink()

def handle_config_command(args, config):
    # If both variable and value are provided, set the config
    if args.variable and args.value:
        config_manager.set_config(args)
    # If only a variable is provided, it's an error
    elif args.variable and not args.value:
        logging.error("A value must be provided when specifying a variable.")
    # If no arguments are provided, show the current config
    else:
        logging.info("--- Current EDSHAT Configuration ---")
        logging.info("Default values are merged with your user-specific settings.")
        logging.info(f"User config file: {config_manager.USER_CONFIG_FILE}\n")
        sys.stdout.write(yaml.dump(config, default_flow_style=False, sort_keys=False))

def handle_download_command(args, config):
    # If --db_dir is specified, it overrides the config value for this run.
    if args.db_dir:
        config['database_dir'] = str(Path(args.db_dir).resolve())

    kraken_db = args.kraken_db if args.kraken_db else 'PlusPFP'
    kraken_size = '-' +  str(args.kraken_size) if args.kraken_size else ''
    config['kraken2_db'] = kraken_db + kraken_size
    config['workdir'] = str(Path(config.get('workdir', '.')))
    runner.run_snakemake(snakefile=config_manager.WORKFLOW_DIR / 'rules' / 'download_db.smk',
                         config=config,
                         threads=args.threads,
                         workdir=Path(config.get('workdir', '.')),
                         dry_run=args.dry_run
        )
    
    if args.db_dir:
        args.variable = 'database_dir'
        args.value = Path(args.db_dir).resolve()
        config_manager.set_config(args)

def main():
    args = parse_args()
    config = config_manager.get_config(args)
    if args.cmd == 'run':
        handle_run_command(args, config)
    elif args.cmd == 'config':
        handle_config_command(args, config)
    elif args.cmd == 'download':
        handle_download_command(args, config)


if __name__ == '__main__':
    main()