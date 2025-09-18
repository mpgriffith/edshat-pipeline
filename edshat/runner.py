# edshat/runner.py

import snakemake
import sys
from pathlib import Path
from typing import Dict, List
import tempfile
import yaml

def run_snakemake(
    snakefile: Path,
    config: Dict,
    workdir: Path,
    threads: int,
    dry_run: bool = False,
    unlock: bool = False,
    targets: List[str] = None
):
    """
    Constructs and executes a Snakemake command using its Python API.

    Args:
        snakefile (Path): Path to the Snakefile to execute.
        config (Dict): Configuration dictionary for the workflow.
        workdir (Path): The working directory for Snakemake.
        threads (int): Number of cores to provide to Snakemake.
        dry_run (bool): If True, run Snakemake in dry-run mode.
        unlock (bool): If True, unlock the working directory.
    """
    # Base arguments for the Snakemake API
    snakemake_args = [
        "--snakefile", str(snakefile),
        "--directory", str(workdir),
        "--cores", str(threads),
       # "--use-conda",      # Assuming you want to use conda environments
        "--reason"          # Print the reason for each rule execution
    ]

    # Add conditional flags
    if dry_run:
        snakemake_args.append("--dry-run")
    
    if unlock:
        snakemake_args.append("--unlock")
        # For unlock, we don't need other configs, so we can run and exit early
        print("🔑 Unlocking working directory...")
        snakemake.main(snakemake_args)
        return

    # Write config dictionary to a temporary YAML file
    with tempfile.NamedTemporaryFile(mode="w", suffix=".yaml", delete=False) as tmp_config:
        yaml.dump(config, tmp_config)
        configfile_path = tmp_config.name

    # Add the configfile argument to Snakemake
    snakemake_args.extend(["--configfile", configfile_path])

    if targets:
        snakemake_args.extend(targets)
    print("🚀 Executing Snakemake workflow with the following arguments:")
    print(f"   snakemake {' '.join(snakemake_args)}\n")

    try:
        # Execute Snakemake. This will block until the workflow is complete.
        status = snakemake.main(snakemake_args)
        if status:
            print("✅ Workflow finished successfully.")
        else:
            print("❌ Workflow finished with errors.", file=sys.stderr)
            sys.exit(1)
            
    except SystemExit as e:
        # Snakemake calls sys.exit() on certain errors, which raises SystemExit
        print(f"❌ Snakemake execution failed with exit code {e.code}.", file=sys.stderr)
        sys.exit(e.code)