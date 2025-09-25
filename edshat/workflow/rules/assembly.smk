include: "common.smk"


#todo: update for long/hybrid assemblies
rule unicycler:
    input: reads=get_input_reads
    output:
        assembly = isolate_dir  + "/{sample}/{sample}.fasta"
    params:
        out_dir= isolate_dir + "/{sample}/unicycler"
    resources:
    threads: 8
    benchmark: isolate_dir + "/{sample}/benchmarks/unicycler.tsv"
    log: isolate_dir  + "/{sample}/unicycler/unicycler.log"
    shell:
        """
        unicycler -1 {input.reads[0]} -2 {input.reads[1]} -o {params.out_dir} -t {threads}
        cp {params.out_dir}/assembly.fasta {output.assembly}
        """

rule quast:
    input: fasta=isolate_dir  + "/{sample}/{sample}.fasta",
            reads=get_input_reads
    output: isolate_dir  + "/{sample}/quast/report.tsv"
    threads: 8
    log: isolate_dir  + "/{sample}/quast/quast.log"
    benchmark: isolate_dir + "/{sample}/benchmarks/quast.tsv"
    params:
        out_dir=isolate_dir  + "/{sample}/quast",
        read_args = lambda w, input:  f"-1 {input.reads[0]} -2 {input.reads[1]}" if len(input.reads) > 1 else f"-1 {input.reads[0]}"
    shell:
        "quast -o {params.out_dir} {params.read_args} {input.fasta} --threads {threads} > {log} 2>&1"

rule ska_fastq:
    input: get_input_reads
    output: isolate_dir  + "/{sample}/{sample}.skf"
    log: isolate_dir  + "/{sample}/logs/ska_fastq.log"
    params:
        prefix=isolate_dir  + "/{sample}/{sample}"
    benchmark: isolate_dir + "/{sample}/benchmarks/ska_fastq.tsv"
    shell:
        "ska fastq -o {params.prefix} {input} > {log} 2>&1"

rule kraken2:
    input: get_input_reads
    output: report= isolate_dir  + "/{sample}/{sample}_kraken.tab",
            out = isolate_dir  + "/{sample}/{sample}_kraken.out"
    params:
        paired = '--paired',
        db = config['database_dir'] + '/kraken2/' 
    benchmark: isolate_dir + "/{sample}/benchmarks/kraken2.tsv"
    shell:
        "kraken2 --db {params.db} {params.paired} --report {output.report} --output {output.out} {input}"
