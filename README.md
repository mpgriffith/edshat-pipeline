# EDS-HAT: Enhanced Detection System for Hospital Acquired Transmission

EDS-HAT is a modular, reproducible bioinformatics pipeline for the detection and analysis of hospital-acquired transmission events using whole genome sequencing data. It automates quality control, assembly, annotation, AMR gene detection, species identification, SNP analysis, and clustering for bacterial isolates.

## Features

- **Automated assembly and QC**: Uses Unicycler and QUAST for genome assembly and quality assessment.
- **Annotation**: Supports Prokka or Bakta for genome annotation.
- **AMR gene detection**: Integrates AMRFinder for antimicrobial resistance gene identification.
- **Species identification**: Uses Kraken2 and MLST for species and sequence type assignment.
- **SNP and cluster analysis**: Computes SNP distances (SKA, Snippy), builds phylogenies, and clusters isolates.
- **Flexible configuration**: YAML-based config files for easy customization.
- **Batch and comparative analyses**: Supports combining new and previous datasets for longitudinal surveillance.

## Installation

1. **Clone the repository**  
   ```bash
   git clone <your-edshat-repo-url>
   cd edshat
   ```

2. **Install dependencies**  
   EDS-HAT requires [Snakemake](https://snakemake.readthedocs.io), Python 3, and the following tools in your `$PATH`:  
   - unicycler
   - quast
   - prokka or bakta
   - amrfinder
   - mlst
   - kraken2
   - ska
   - snippy
   - raxml-ng

   You can install dependencies via conda:
   ```bash
   conda env create -f environment.yaml
   conda activate edshat
   ```

3. **Configure reference databases**  
   Update `workflow/config/config.yaml` with paths to your reference databases (e.g., Kraken2 DB).

## Usage

### 1. Prepare your data

- Place paired-end FASTQ files in a directory (e.g., `Reads/Sample1/`).
- Optionally, prepare a CSV sample sheet with columns: `Sample,Read1,Read2`.

### 2. Run the pipeline

Basic run:
```bash
python edshat-pipeline run -i Sample1 Sample2 -o results/
```

With a sample sheet:
```bash
python edshat-pipeline run -i samples.csv -o results/
```

Key options:
- `-j/--threads`: Number of threads (default: 12)
- `--reads_dir`: Directory containing reads (default: Reads/)
- `--config`: Additional config YAML(s)
- `--species`: Expected species or mapping file
- `--targets`: Pipeline targets (assembly, annotation, amr, species, metrics, ska, snps, tree, etc.)

### 3. Combine with previous data

To add new isolates to an existing dataset and re-cluster:
```bash
python edshat-pipeline combine -i new_samples.csv -d previous_metrics_or_data.csv --output_dir results/
```

### 4. Output

- Per-sample results in `results/isolates/<sample>/`
- Set-level comparative results in `results/sets/<set_name>/`
- Key outputs:
  - Assembly FASTA, annotation files, AMR results, MLST, Kraken2 reports
  - SNP distance matrices, cluster assignments, phylogenetic trees
  - Metrics and summary CSVs for downstream analysis

## Configuration

- Main config: [`workflow/config/config.yaml`]
- Species mapping: [`workflow/config/species_map.yaml`]
- Genome stats: [`workflow/config/genome_stats.yaml`]


## Citation

If you use EDS-HAT in your work, please cite the relevant tools and this pipeline.

## License

See LICENSE file for details.

## Contact

For questions or issues, please open an issue on the repository or contact the authors.
