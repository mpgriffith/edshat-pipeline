include: "common.smk"

set_dir = config['set_dir']
ska_sets, st_sets = get_sets()
set_name = config.get('set_name', Path(set_dir).stem)

rule combine_snps:
    input: ska_snps=set_dir + "/{set}/{set}_ska.distances.csv",
            snippy_snps=lambda wildcards: expand(set_dir + "/{st_name}/{st_name}_SNP_matrix_core.csv", st_name=get_st_species(wildcards.set))
    output: set_dir + "/{set}_combined_SNP_dist.csv",
            set_dir + "/{set}_combined_SNP_matrix.csv"
    benchmark: set_dir + "/benchmarks/{set}_combine_snps.tsv"
    run:
        import pandas as pd
        import itertools
        import numpy as np
        ska_df = pd.read_csv(input.ska_snps)
        snippy_ms = [pd.read_csv(f, index_col=0) for f in input.snippy_snps]
        ska_df = ska_df.drop_duplicates(subset=['Subject', 'Reference'])
        avg_matches = ska_df['Matches'].mean()
        for idx, row in ska_df.iterrows():
            i1 = row['Subject']
            i2 = row['Reference']
            matches = row['Matches']
            if matches <= avg_matches / 2:
                row['SNPs'] = 100000
            for m in snippy_ms:
                if i1 in m.index and i2 in m.index:
                    snippy_snps = m.loc[i1, i2]
                    if snippy_snps < row['SNPs']:
                        ska_df.loc[idx, 'SNPs'] = snippy_snps
        ska_df.to_csv(output[0])
        #todo: make a distance matrix from this
        matrix = ska_df.pivot(index='Subject', columns='Reference', values='SNPs')
        matrix = matrix.fillna(0)
        matrix = matrix.add(matrix.T, fill_value=0)
        matrix.to_csv(output[0], index=True, header=True)

        matrix.to_csv(output[1])


def get_species_snp_params(wildcards, param='linkage'):
    species_d = config.get('species', {})
    species_spec = species_d.get('species_specific', None)
    species_name = wildcards.set.replace('_', ' ')
    if not species_spec:
       return config[param]
    if not species_name in species_spec:
       return config[param]
    return species_spec[species_name].get(param, config[param])


rule get_clusters:
    input: set_dir + "/{set}_combined_SNP_matrix.csv"
    output: set_dir + "/clusters/{set}_clusters.csv"
    params: method = lambda wildcards: get_species_snp_params(wildcards, param='linkage'),
            threshold = lambda wildcards: get_species_snp_params(wildcards, param='SNP_threshold')
    benchmark: set_dir + "/benchmarks/{set}_clusters.tsv"
    shell:
        """
        {script_dir}/find_clusters.py -m {params.method} -t {params.threshold} -o {output} {input}
        """

rule combine_clusters:
    input:
        cluster_data = expand(set_dir + "/clusters/{set}_clusters.csv", set=get_sets()[0].keys()),
        snp_matrices = expand(set_dir + "/{set}_combined_SNP_matrix.csv", set=get_sets()[0].keys()) # Pass the list of SNP matrices via config
    output:
        report_csv = set_dir + "/EDSHAT_" + set_name + "_New_cluster_data.csv",
        snp_report_csv = set_dir + "/EDSHAT_" + set_name + "New_cluster_SNP_distances.csv"
    params:
        output_columns = config.get('output_cols', [])
    run:
        cluster_data = pd.concat([pd.read_csv(f, header=None, index_col=0, names=['Cluster']).assign(Set=f.stem.replace('_clusters', '')) for f in cluster_files]) if len(cluster_files) > 0 else pd.DataFrame()
        if cluster_data.empy:
            pass
        cluster_data['ClusterName'] = cluster_df['Set'] + '_' + cluster_df['Cluster'].astype(str)

        snp_ms = {
            Path(f).name.split('_')[0].replace('-', ' '): pd.read_csv(f, index_col=0) 
            for f in input.snp_matrices
        }

        # Filtering logic
        # filter to clusters with new isolates
        new_samples = set(config.get('new_samples', []))
        new_sample_clusters = cluster_df[cluster_data.index.isin(new_samples)].ClusterName
        new_cluster_data = cluster_data[cluster_data['ClusterName'.isin(new_sample_clusters)]]
        # SNP calculation logic
        snp_rows = []
        for cluster, cdf in new_cluster_data.groupby('Cluster'):
            sp = cdf['Species'].mode()[0]
            if sp in snp_ms:
                snp_m = snp_ms[sp]
                cluster_samples = cdf.index.tolist()
                snps = []
                for s1, s2 in itertools.combinations(cluster_samples, 2):
                    if s1 in snp_m.index and s2 in snp_m.columns:
                        snps.append(snp_m.loc[s1, s2])
                        snp_rows.append([s1, s2, sp, cluster, snp_m.loc[s1, s2]])
                
                if snps:
                    cluster_data.loc[cdf.index, 'SNP Range'] = f'{min(snps)}-{max(snps)}'

        snp_df = pd.DataFrame(snp_rows, columns=['Sample 1', 'Sample 2', 'Species', 'Cluster', 'SNPs'])
        
        # Select output columns and save
        # out_cols = [c for c in params.output_columns if c in cluster_data.columns]
        if config.get('previous_data', None):
            prev_data  = config['previous_data']
            if prev_data.endswith('xlsx') or prev_data.endswith('xls'):
                prev_data = pd.read_excel(prev_data, index_col=0)
            else:
                prev_data = pd.read_csv(prev_data, index_col=0)
            db_data = prev_data[prev_data.index.isin(new_clusters.index)]
            db_data['Species'] = cluster_data.index.map(lambda x: new_clusters.loc[x, 'Set'].replace('_', ' ') if x in new_clusters.index else None)
            db_data['Cluster'] = cluster_data.index.map(lambda x: new_clusters.loc[x, 'ClusterName'] if x in new_clusters.index else None)
            db_data['NewIsolate'] = cluster_data.index.map(lambda x: 1 if x in samples.index else 0)
            db_data['DateAnalyzed'] = datetime.now().strftime("%Y-%m-%d")
            db_data['SNP Range'] = cluster_data.index.map(lambda x: )
            dup_cols = config.get('duplicate_columns', [])
            remove_clusters = []
            if dup_cols:
                for cluster, cdf in db_data.groupby('Cluster'):
                    cdf_dd = cdf.drop_duplicates(subset=dup_cols)
                    if len(cdf_dd) == 1:
                        remove_clusters.append(cluster)
                db_data = db_data[~db_data['Cluster'].isin(remove_clusters)]
            output_cols = config.get('output_cols', db_data.columns.tolist())
            db_data[output_cols].to_csv(output.report_csv, index=True)

        else:    
            cluster_data.to_csv(output.report_csv, index=True)
        snp_df.to_csv(output.snp_report_csv, index=False)