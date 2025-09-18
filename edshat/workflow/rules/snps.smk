


include: "common.smk"

wildcard_constraints:
    sample="^/"

isolate_dir = config.get('isolate_dir', os.path.join(config['workdir'], 'isolates'))
set_dir = config['set_dir']
script_dir =  os.path.abspath(srcdir('../scripts'))

rule ska_compare:
    input: 
           query_skas =lambda wildcards: get_ska_skfs(wildcards)
    output: set_dir + "/{set}/{set}_ska.distances.csv"
    threads: 8
    shell:
        """
        {script_dir}/ska_compare_to_dist.py -t {threads} -o {output} {input.query_skas}
        """


rule ska_matrix:
    input: set_dir + "/{set}/{set}_ska.distances.csv"
    output: set_dir + "/{set}/{set}_ska.matrix.csv"
    benchmark: set_dir + "/benchmarks/{set}_combine_snps.tsv"
    run:
        import pandas as pd
        import numpy as np

        ska_df = pd.read_csv(input[0])
        ska_df = ska_df.drop_duplicates(subset=['Subject', 'Reference'])
        matrix = ska_df.pivot(index='Subject', columns='Reference', values='SNPs')
        matrix = matrix.fillna(0)
        matrix = matrix.add(matrix.T, fill_value=0)
        matrix.to_csv(output[0], index=True, header=True)
    

rule snippy_sample:
    input: reference=lambda wildcards: get_snippy_input(wildcards, itype="reference"),
            reads=lambda wildcards: get_input_reads(wildcards)
    output: isolate_dir + "/{sample}/snippy_{ref}/{sample}/snps.tab"
    benchmark: isolate_dir + "/{sample}/benchmarks/{sample}_snippy_{ref}.tsv"
    params:
        output_directory=isolate_dir + "/{sample}/snippy_{ref}/{sample}",
        sample="{sample}",
        read_param = lambda wildcards, input: "-R1 " + input.reads[0] + " -R2 " + input.reads[1] if len(input.reads) == 2 else "--contigs " + input.reads[0]
    shell:
        "snippy --force --cpus {threads} --outdir {params.output_directory} --reference {input.reference} {params.read_param}"


rule snippy_core:
    input: ref = lambda wildcards: get_snippy_input(wildcards, itype="reference"),
            snippy_res = lambda wildcards: get_snippy_input(wildcards, itype="snippy_results")
    
    output: set_dir + "/{set}/{set}.aln"
    log: set_dir + "/{set}/logs/snippy_core.log"
    benchmark: set_dir + "/benchmarks/{set}_snippy_core.tsv"
    params: out_dir=set_dir + "/{set}",
            snippy_folds = lambda wildcards: get_snippy_input(wildcards, itype="snippy_folds")
    shell:
        """
        mkdir -p {params.out_dir}
        cd {params.out_dir}
        snippy-core --prefix {wildcards.set} --ref {input.ref} {params.snippy_folds}   2> {log}
        """

rule snippy_postprocess:
    input: aln=set_dir + "/{set}/{set}.aln",
            ref = lambda wildcards: get_snippy_input(wildcards, itype='reference', suffix='gbk')

    output: matrix=set_dir + "/{set}/{set}_SNP_matrix_core.csv",
            aln=set_dir + "/{set}/{set}.snps.aln"
    params: snippy_dir = set_dir + "/{set}"
    log: set_dir + "/{set}/logs/snippy_postprocess.log"
    benchmark: set_dir + "/benchmarks/{set}_snippy_post.tsv"
    shell:
        """
        {script_dir}/postprocess_snippy.py -r {input.ref} {params.snippy_dir}  2> {log}
        """

rule raxml:
    input: set_dir + "/{set}/{set}.snps.aln"
    output: set_dir + "/{set}/{set}.raxml.tree"
    params: out_dir=set_dir + "/{set}"
    benchmark: set_dir + "/benchmarks/{set}_tree.tsv"
    threads: 8
    params: model = "GTR+ASC_LEWIS",
            bootstraps = 100
    shell:
        """
        cd {set_dir}/{wildcards.set}
        raxml-ng --all --msa {input} --model GTR+ASC_LEWIS --bs-trees 100 --threads {threads} --redo --prefix {wildcards.set}
        cp {wildcards.set}.raxml.besttree  {output}
        """

