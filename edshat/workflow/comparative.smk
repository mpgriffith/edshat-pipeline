import datetime


include: "rules/common.smk"
include: "rules/snps.smk"
include: "rules/clusters.smk"
include: "rules/annotation.smk"

workdir: config['workdir']
set_dir = config['set_dir']


wildcard_constraints:
    sample= "[^/]+",
    ref = "[^/]+"

ska_sets, sets = get_sets()

script_dir =  os.path.abspath(srcdir('../scripts'))
isolate_dir = config.get('isolate_dir', os.path.join(config['workdir'], 'isolates'))
set_name = config.get('set_name', Path(set_dir).stem)

def get_comp_outputs(wildcards, targets=None):
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
    if 'snps' in targets:
        pass
    if 'tree' in targets:
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
    
