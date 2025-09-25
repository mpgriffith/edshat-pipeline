
isolate_dir = config.get('isolate_dir', os.path.join(config.get('workdir', os.path.abspath(os.curdir)), 'isolates'))
db_dir = config.get('database_dir', )

rule gtdbtk:
    input: isolate_dir  + "/{sample}/{sample}.fasta"
    output: isolate_dir  + "/{sample}/gtdbtk/gtdbtk.bac120.summary.tsv"
    log: isolate_dir  + "/{sample}/logs/gtdbtk.log"
    params: sample_dir = isolate_dir  + "/{sample}",
            gtdbtk_data_dir = config['database_dir'] + '/gtdbtk'
    shell:
        """
        export GTDBTK_DATA_PATH={params.gtdbtk_data_dir}
        gtdbtk classify_wf --genome_dir {params.sample_dir} --out_dir {params.sample_dir}/gtdbtk -x fasta --skip_ani_screen
        """


checkpoint specieate:
    input: gtdbtk =isolate_dir + "/{sample}/gtdbtk/gtdbtk.bac120.summary.tsv"
    output: species=isolate_dir + "/{sample}/{sample}_species.txt"
    run:
        import pandas as pd
        import re
        df = pd.read_csv(input.gtdbtk, sep='\t')
        species = df['classification'].values[0].split('__')[-1]
        # remove species group designations
        species = re.sub(r'_[A-Z]$', '', species)
        shell("echo {species} > {output}")

def get_species(wildcards):
    species_file = checkpoints.specieate.get(sample=wildcards.sample).output.species
    with open(species_file) as f:
        species = f.read().strip()
    return species

def get_species_files(wildcards):
    species = get_species(wildcards)
    files = []
    files.append(isolate_dir + "/{sample}/{sample}_species.txt")
    genus = species.split()[0]
    if genus == 'Klebsiella':
        files.append(isolate_dir + "/{sample}/kleborate/{sample}_kleborate.tsv")
    # if genus in ['Mycobacterium', 'Mycolicibacterium']:
    #     files.append(isolate_dir  + "/{sample}/{sample}_subspecies.txt")
    return files

rule kleborate:
    input: fasta = isolate_dir  + "/{sample}/{sample}.fasta",
           species_file = isolate_dir + "/{sample}/{sample}_species.txt"
    output: kleborate=isolate_dir  + "/{sample}/kleborate/{sample}_kleborate.tsv"
    params: out_dir = isolate_dir + "/{sample}/kleborate"
    log: isolate_dir  + "/{sample}/logs/{sample}_kleborate.log"
    run:
        sp = get_species(wildcards)
        genus = sp.split()[0]
        preset = '--preset kpsc'
        complex = 'klebsiella_pneumo_complex'
        if sp == 'Klebsiella pneumoniae':
            preset = '--preset kpsc'
            complex = 'klebsiella_pneumo_complex'
        elif sp == 'Klebsiella oxytoca':
            preset = '--preset kosc'
            complex = 'klebsiella_oxytoca_complex'
        elif sp == 'Escherichia coli':
            preset = '--preset escherichia'
            comples = 'escherichia'

        shell(f"kleborate -a {{input.fasta}} --outdir {{params.out_dir}} {preset} --trim_headers > {{log}} 2>&1")
        shell("cp {{params.out_dir}}/{}_output.txt {{output.kleborate}}".format(complex))       

rule species_summary:
    input: get_species_files
    output: isolate_dir  + "/{sample}/{sample}_species_summary.tsv"
    run:
        for f in input:
            shell("echo {f} >> {output}")
