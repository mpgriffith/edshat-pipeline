include: "common.smk"

import yaml
from pathlib import Path

rule annotate:
    input: isolate_dir + "/{sample}/{sample}.fasta"
    output: isolate_dir + "/{sample}/annotation/{sample}.faa",
            isolate_dir + "/{sample}/annotation/{sample}.gbk"
    threads: 8
    log: isolate_dir + "/{sample}/logs/annotation.log"
    benchmark: isolate_dir + "/{sample}/benchmarks/annotation.tsv"
    params:
        out_dir = isolate_dir + "/{sample}/annotation"
    
    run:
        if config.get('annotation', 'prokka') == 'bakta':
            db_arg = f"--db {config['database_dir']}/bakta" if 'bakta' in config.get('database_dir', '') else ''
            shell("bakta {db_arg} --outdir {params.out_dir} --threads {threads} --force {input} --prefix {wildcards.sample} > {log} 2>&1")
        else:
            shell("prokka --prefix {wildcards.sample} --outdir {params.out_dir} --cpus {threads} --force --locustag {wildcards.sample} {input} > {log} 2>&1")

def get_mlst_schema(wildcards):
    species = get_species(wildcards)
    if species == 'Mycobacterium abscessus':
        return '--scheme mabscessus'
    if species == 'Acinetobacter baumannii':
        return '--scheme abaumannii2'
    else:
        return ''

rule mlst:
    input: fasta=isolate_dir + "/{sample}/{sample}.fasta",
              species=isolate_dir + "/{sample}/{sample}_species.txt"
    output: isolate_dir + "/{sample}/{sample}_mlst.tsv"
    params: schema = get_mlst_schema
    benchmark: isolate_dir + "/{sample}/benchmarks/mlst.tsv"
    shell:
        "mlst {params.schema} {input} > {output}"

amrfinder_species = {}


def get_amrfinder_species_arg(wildcards):
    species = get_species(wildcards)
    if species in amrfinder_species:
        return f"--species {amrfinder_species[species]}"
    else:
        return ''

rule amrfinder:
    input: fasta=isolate_dir + "/{sample}/{sample}.fasta", 
              annotation=isolate_dir + "/{sample}/annotation/{sample}.gbk",
              species = isolate_dir + "/{sample}/{sample}_species.txt"
    output: isolate_dir + "/{sample}/{sample}.amrfinder.tsv"
    log: isolate_dir + "/{sample}/logs/amrfinder.log"
    benchmark: isolate_dir + "/{sample}/benchmarks/amrfinder.tsv"
    params:
        species = get_amrfinder_species_arg,
        annot_dir = isolate_dir + "/{sample}/annotation",
        annot_arg = '-a bakta' if config.get('annotation', 'prokka') == 'bakta' else '-a prokka',
    shell:
        "amrfinder --plus -n {input.fasta} -p {params.annot_dir}/{wildcards.sample}.faa -g {params.annot_dir}/{wildcards.sample}.gff --output {output} {params.species} {params.annot_arg} --threads {threads}"


def pass_species_check(kraken_results, expected):
    exp_genus = expected.split()[0]
    exp_is_genus = len(expected.split()) == 1
    if kraken_results['Species #2'] != 'unclassified' and kraken_results['Species #1'].split()[0] != exp_genus and exp_is_genus:
        if kraken_results['Species #2'].split()[0] != exp_genus and kraken_results['% Species #2'] >= 10:
            return False
        else:
            return True
    elif kraken_results['Species #1'] == expected:
        if kraken_results['Species #2'] != 'unclassified' and kraken_results['Species #2'].split()[0] != exp_genus and kraken_results['% Species #2'] >= 10:
            return False
        else:
            return True
    elif kraken_results['Species #1'] == 'unclassified' and kraken_results['Species #2'].split()[0] == exp_genus and exp_is_genus:
        return True
    elif kraken_results['Species #1'] == 'unclassified' and kraken_results['Species #2'] == expected:
        return True
    else:
        return False

def pass_assembly_check(quast_metrics, expected, found_species=None):
    map_file = config.get('genome_stats_file', None)
    map_file = workflow.source_path(map_file) if map_file else None
    print(map_file)
    if map_file is None:
        genome_size = config.get('expected_genome_size', 5e6)  # Default to 5 Mb if not specified
        genome_sizes = {}
    else:
        genome_size = config.get('expected_genome_size', 5e6)
        with open(map_file, 'r') as f:
            genome_sizes =  yaml.safe_load(f)
        exp_genus = expected.split()[0]
    exp_is_genus = len(expected.split()) == 1
    print(genome_sizes)
    genome_size = genome_sizes.get(expected, genome_size)
    if exp_is_genus and found_species.split()[0] == exp_genus:
        if found_species in genome_sizes:
            genome_size = genome_sizes[found_species]
        elif exp_genus in genome_sizes:
            genome_size = genome_sizes[exp_genus]
    elif found_species in genome_sizes:
        genome_size = genome_sizes[found_species]
    elif exp_genus in genome_sizes:
        genome_size = genome_sizes[exp_genus]
    
    
    length = quast_metrics['Total Length']
    num_contigs = quast_metrics['NumContigs']
    depth = quast_metrics.get('Coverage', 36)

    contig_qc = config.get('min_contigs', 350)
    depth_qc = config.get('min_depth', 35)
    if length < 0.8 * genome_size or length > 1.2 * genome_size:
        return False
    else:
        return True if num_contigs <= contig_qc and depth >= depth_qc else False

rule checkm:
    input: isolate_dir + "/{sample}/{sample}.fasta"
    output: isolate_dir +"}/{sample}/checkm/checkm.tsv"
    params: out_dir=isolate_dir +"/{sample}/checkm",
           sample_dir=isolate_dir +"/{sample}"
    log: isolate_dir +"/{sample}/logs/checkm.log"
    benchmark: isolate_dir + "/{sample}/benchmarks/checkm.tsv"
    shell:
        "checkm lineage_wf -x fasta {params.sample_dir} {params.out_dir} -t {threads} --tab_table -f {output} 2> {log}"

rule metrics:
    #TODO: fix metrics to our columns
    #      kraken results to top 2 species/genus?
    input:  quast=isolate_dir + "/{sample}/quast/report.tsv",
            kraken=isolate_dir + "/{sample}/{sample}_kraken.tab",
            mlst=isolate_dir + "/{sample}/{sample}_mlst.tsv",
            gtdbtk=isolate_dir + "/{sample}/gtdbtk/gtdbtk.bac120.summary.tsv",
            species=isolate_dir + "/{sample}/{sample}_species.txt"
    output: isolate_dir + "/{sample}/{sample}.metrics.csv"
    run:
        import pandas as pd
        quast_df = pd.read_csv(input.quast, sep='\t', index_col=0)
        kraken_df = pd.read_csv(input.kraken, sep='\t', header=None, names=['pct_reads', 'num_reads', 'num_reads_parent', 'rank', 'taxid', 'tax_name'])
        kraken_top_2 = {''}
        mlst_df = pd.read_csv(input.mlst, sep='\t', header=None)
        gtdbtk_df = pd.read_csv(input.gtdbtk, sep='\t')
        cols = ['SpecID', 'Species #1', '% Species #1', 'Species #2', '% Species #2', 'Isolation Flag', 'NumContigs', 'Total Length', 'N50', 'Coverage', 'Resequencing Flag', 'ST', 'rMLST', 'species']
        species = next(open(input.species, 'r')).strip()
        # Extract relevant metrics
        metrics = {
            'SpecID': wildcards.sample,
            'NumContigs': quast_df.loc['# contigs (>= 0 bp)', wildcards.sample],
            'Total Length': quast_df.loc['Total length (>= 0 bp)', wildcards.sample],
            'N50': quast_df.loc['N50', wildcards.sample],
            'Coverage': quast_df.loc['Avg. coverage depth', wildcards.sample] if 'Avg. coverage depth' in quast_df.index else 36,
            'ST': mlst_df[2].iloc[0] if not mlst_df.empty else None,
            'rMLST': None,
            'species': species,
        
        }
        
        kraken_sp_df = kraken_df[kraken_df['rank'].str.contains('S|U')].sort_values('pct_reads', ascending=False)
        for i in range(0,2):
            if i < len(kraken_sp_df):
                metrics['Species #' + str(i+1)] = kraken_sp_df.iloc[i]['tax_name'].strip()
                metrics['% Species #' + str(i+1)] = kraken_sp_df.iloc[i]['pct_reads']
            else:
                metrics['Species #' + str(i+1)] = 'unclassified'
                metrics['% Species #' + str(i+1)] = 0.0
        exp_species = samples.loc[wildcards.sample]['species']
        print(metrics)
        sp = get_species(wildcards)
        metrics['Isolation Flag'] =  pass_species_check(metrics, exp_species)
        metrics['Resequencing Flag'] = pass_assembly_check(metrics, exp_species, found_species=sp)
        # Create DataFrame and save to CSV
        metrics_df = pd.DataFrame([metrics])
        metrics_df.to_csv(output[0], columns=cols, index=False)
