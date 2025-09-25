


include: "common.smk"

wildcard_constraints:
    sample="^/"

isolate_dir = config.get('isolate_dir', 'isolates')
set_dir = config['set_dir']
script_dir =  os.path.abspath(srcdir('../scripts'))

rule ska_compare:
    input: 
           query_skas =lambda wildcards: get_ska_skfs(wildcards)
    output: set_dir + "/{set}/{set}_ska.distances.csv"
    threads: 8
    benchmark: set_dir + "/benchmarks/{set}_ska.tsv"
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
    

# checkpoint snippy_sample:
#     input: reference=isolate_dir + "/{ref}/{ref}.fasta",
#             reads=lambda wildcards: get_input_reads(wildcards)
#     output: isolate_dir + "/{sample}/snippy_{ref}/{sample}/snps.tab"
#     benchmark: isolate_dir + "/{sample}/benchmarks/{sample}_snippy_{ref}.tsv"
#     params:
#         output_directory=isolate_dir + "/{sample}/snippy_{ref}/{sample}",
#         sample="{sample}",
#         read_param = lambda wildcards, input: "-R1 " + input.reads[0] + " -R2 " + input.reads[1] if len(input.reads) == 2 else "--contigs " + input.reads[0]
#     shell:
#         "snippy --force --cpus {threads} --outdir {params.output_directory} --reference {input.reference} {params.read_param}"


def get_all_snippy_results(wildcards):
    """
    This function is evaluated *after* the snippy_sample checkpoint has been completed.
    It gathers all the required snps.tab files for a given set.
    """
    ska_sets, snippy_sets = get_sets()
    st_sets_map = get_st_sets()
    
    # This function is now just for getting the folder paths for the shell command
    # The actual dependency is handled by the `input` function of the rule.
    if wildcards.set in st_sets_map:
        ref = st_sets_map[wildcards.set][0]
        samples_in_set = st_sets_map[wildcards.set]
        # Return the parent directories for the snippy-core command
        return [str(Path(f).parent) for f in expand(isolate_dir + "/{sample}/snippy_{ref}/{sample}/snps.tab", sample=samples_in_set, ref=ref)]
    return []

def get_snippy_core_input(wildcards):
    """
    Input function for snippy_core.
    Determines the samples and reference for a given set, then gets the corresponding
    output from the snippy_sample checkpoint.
    """
    st_sets_map = get_st_sets()
    samples_in_set = st_sets_map.get(wildcards.set, [])
    if samples_in_set:
        ref = samples_in_set[0] # The reference is the first sample in the set
        # Use expand to generate the list of expected output files.
        # Snakemake will see this list and know which checkpoint jobs to run.
        return {'reference': isolate_dir + '/' + ref + '/' + ref + '.fasta', 
                'samples':expand(isolate_dir + "/{sample}/snippy_{ref}/{sample}/snps.tab",
                      sample=samples_in_set,
                      ref=ref)}
    return {'reference': '', 'samples': []}

# rule snippy_core:
#     input: get_snippy_core_input
#     output: set_dir + "/{set}/{set}.aln"
#     log: set_dir + "/{set}/logs/snippy_core.log"
#     benchmark: set_dir + "/benchmarks/{set}_snippy_core.tsv"
#     params: out_dir=set_dir + "/{set}", snippy_folds = get_all_snippy_results
#     shell:
#         """
#         snippy-core --prefix {params.out_dir}/{wildcards.set} {params.snippy_folds} > {log} 2>&1
#         """

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
    params: out_dir=set_dir + "/{set}",
            model = "GTR+ASC_LEWIS",
            bootstraps = 100,
            aln = "{set}.snps.aln",
            output_stem = "{set}.raxml.tree"
    benchmark: set_dir + "/benchmarks/{set}_tree.tsv"
    threads: 8
    shell:
        """
        cd {set_dir}/{wildcards.set}
        raxml-ng --all --msa {params.aln} --model GTR+ASC_LEWIS --bs-trees 100 --threads {threads} --redo --prefix {wildcards.set}
        cp {wildcards.set}.raxml.bestTree  {wildcards.set}.raxml.tree
        """
