import datetime


include: "rules/common.smk"
include: "rules/snps.smk"
include: "rules/clusters.smk"
include: "rules/annotate.smk"

workdir: config['workdir']
set_dir = config['set_dir']


wildcard_constraints:
    sample= "[^/]+",
    ref = "[^/]+"

ska_sets, sets = get_sets()

script_dir =  os.path.abspath(srcdir('../scripts'))
isolate_dir = config.get('isolate_dir', os.path.join(config['workdir'], 'isolates'))
set_name = config.get('run_name', Path(set_dir).stem)

def get_comp_outputs(wildcards, targets=None):
    ska_sets, sets = get_sets()

    if targets == None:
        targets = ['ska', 'snippy', 'snps','clusters', 'tree', 'pangenome', 'combine']
    outputs = []
    #set_dir = config['output_directory']
    if 'ska' in targets:
        ska_dfs = expand(set_dir +'/{name}/{name}_ska.matrix.csv', name=ska_sets.keys())
        outputs += ska_dfs
    if 'snippy' in targets:
        snippy_results = expand(set_dir + "/{name}/{name}_SNP_matrix_core.csv", name=get_st_sets_names(st_sets))
        outputs += snippy_results
    if 'clusters' in targets:
        cluster_files = expand(set_dir + '/clusters/{name}_clusters.csv', name=ska_sets.keys())
        outputs += cluster_files
    # if 'snps' in targets:
    #     for sp, sets in sets.items():
    #         for st, samples in sets.items():
    #             ref = samples[0]
    #             for sample in samples:
    #                 outputs += [isolate_dir + '/' + sample + '/' + 'snippy_' + ref + '/' + sample + '/snps.tab']
    # if 'tree' in targets:
        trees = expand(set_dir + "/{name}/{name}.raxml.tree", name=get_st_sets_names(st_sets, min_len=4))
        outputs += trees
    if 'pangenome' in targets:
        pass
    if 'recomb' in targets:
        recomb = expand(set_dir + "/{name}/{name}_recomb.gff", name=get_st_sets_names(st_sets))
    if 'combine' in targets:
        cluster_data = set_dir + "/EDSHAT_" + set_name + "_New_cluster_data.csv"
        cluster_snps = set_dir + "/EDSHAT_" + set_name + "_New_cluster_SNP_distances.csv"
        outputs += [cluster_data, cluster_snps]
    print(outputs)
    return outputs


rule all:
    input: lambda wildcards: get_comp_outputs(wildcards)

rule snippy:
    input: lambda wildcards: get_comp_outputs(wildcards, targets=['snippy'])

rule snps:
    input: lambda wildcards: get_comp_outputs(wilcards, targets=['snps'])

rule tree:
    input: lambda wildcards: get_comp_outputs(wildcards, targets=['tree'])

rule ska_snps:
    input: lambda wildcards: get_comp_outputs(wildcards, targets=['ska'])

rule pangenome:
    input: lambda wildcards: get_comp_outputs(wildcards, targets=['pangenome'])


rule recomb:
    input: lambda wildcards: get_comp_outputs(wildcards, targets=['recomb'])

rule combine_data:
    input:
    

checkpoint snippy_sample:
    input: reference=isolate_dir + "/{ref}/{ref}.fasta",
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
    input: unpack(get_snippy_core_input)
    output: set_dir + "/{set}/{set}.aln"
    log: set_dir + "/{set}/logs/snippy_core.log"
    benchmark: set_dir + "/benchmarks/{set}_snippy_core.tsv"
    params: out_dir=set_dir + "/{set}", snippy_folds = get_all_snippy_results
    shell:
        """
        snippy-core --ref {input.reference} --prefix {params.out_dir}/{wildcards.set} {params.snippy_folds} > {log} 2>&1
        """