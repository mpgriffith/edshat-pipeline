import os
import pandas as pd
import collections
from pathlib import Path

isolate_dir = config.get('isolate_dir', os.path.join(config.get('workdir', os.path.abspath(os.curdir)), 'isolates'))

samples = pd.read_csv(config['samples'], index_col=0)

#script_dir =  os.path.abspath(srcdir('scripts'))

def get_input_reads(wildcards, type='illumina'):
    reads = samples.loc[wildcards.sample]
    if type == 'illumina':
        return [reads['Read1'], reads['Read2']]
    elif type == 'nanopore':
        return reads['nanopore']
    else:
        raise ValueError("Invalid read type specified.")


def get_st_sets_names(st_sets, min_len=0):
    st_set_names = []
    for species in st_sets.keys():
        for st_set in st_sets[species]:
            if len(st_sets[species][st_set]) > min_len:
                st_set_names.append(st_set)
    return st_set_names

def get_st_sets(st_sets=None):
    if not st_sets:
        ska_sets, st_sets = get_sets()
    st_sets_collapse = {}
    for species in st_sets:
        st_sets_collapse.update(st_sets[species])

    return st_sets_collapse

def get_st_species(species, st_sets=None):
    if not st_sets:
        ska_sets, st_sets = get_sets()
    return list(st_sets[species].keys())

def get_ska_skfs(wildcards, ska_sets=None):
    if not ska_sets:
        ska_sets, st_sets = get_sets()
    
    samples = ska_sets[wildcards.set]
    return expand(config['isolate_dir'] + "/{sample}/{sample}.skf", sample=samples)

def get_snippy_input(wildcards, itype="reference", suffix='fasta'):
    ska_sets, snippy_sets = get_sets()
    st_sets = get_st_sets()
    # reference = st_sets[wildcards.set][0]
    reference = wildcards.ref if hasattr(wildcards, 'ref') else st_sets[wildcards.set][0]
    if itype == "reference":
        return config['isolate_dir'] + "/" + reference + "/annotation/" + reference + "." + suffix
    elif itype == "snippy_folds":
        return expand(config['isolate_dir'] + "/{sample}/snippy_" + reference + "/{sample}", sample=st_sets[wildcards.set])
    elif itype == "snippy_results":
        return expand(config['isolate_dir'] + "/{sample}/snippy_" + reference + "/{sample}/snps.tab", sample=st_sets[wildcards.set])
    else:
        raise ValueError("Invalid input type specified for snippy input.")


def get_sets():
    """
    Reads metrics files for an EXPLICIT list of samples provided in the config
    to create Ska and Snippy comparison sets.
    """
    isolate_dir = Path(config.get('isolate_dir'))
    
    samples_to_analyze = config.get('analysis_samples', [])
    if not samples_to_analyze:
        print("Warning: No samples provided in 'analysis_samples' config. No sets will be made.")
        return {}, {}

    # Construct paths to only the required metrics files
    metrics_files = [isolate_dir / s / f"{s}.metrics.csv" for s in samples_to_analyze]
    
    # Read only the files that actually exist, gracefully skipping missing ones
    df_list = []
    for f in metrics_files:
        if f.exists():
            df_list.append(pd.read_csv(f))
        else:
            print(f"Warning: Could not find metrics file for sample: {f.parent.name}")

    if not df_list:
        return {}, {}
    

    all_metrics_df = pd.concat(df_list, ignore_index=True)
    pass_metrics_df = all_metrics_df[all_metrics_df['PassQC'] == True]

    # --- 1. Create Ska Sets (grouped by species) ---
    ska_sets = collections.defaultdict(list)
    for _, row in pass_metrics_df.iterrows():
        species_name = str(row['species']).replace(' ', '_')
        ska_sets[species_name].append(row['SpecID'])
    
    # --- 2. Create Snippy Sets (configurable grouping) ---
    snippy_grouping = config.get('snippy_grouping_method', 'ST').upper()
    snippy_sets = collections.defaultdict(list)

    if snippy_grouping == 'ST':
        # Group by both species and ST for high-resolution comparison
        # Clean up ST column: fill missing values and remove '.0' from numbers
        pass_metrics_df['ST'] = pass_metrics_df['ST'].fillna('NoST').astype(str).str.replace(r'\.0$', '', regex=True)
        
        for group_keys, group_df in all_metrics_df.groupby(['species', 'ST']):
            species, st = group_keys
            species_name = str(species).replace(' ', '_')
            set_name = f"{species_name}_ST{st}" # e.g., Klebsiella_pneumoniae_ST11
            if not species_name in snippy_sets: snippy_sets[species_name] = {}
            snippy_sets[species_name][set_name] = group_df['SpecID'].tolist()
    else: # If method is 'SPECIES' or anything else, group by species
        snippy_sets = ska_sets


    new_samples_set = set(config.get('new_samples', []))
    
    if not new_samples_set:
        # no new samples pass, assume do all samples
        new_sample_set = samples

    active_ska_sets = {name: samples for name, samples in ska_sets.items()
        if not new_samples_set.isdisjoint(samples)}
    
    active_snippy_sets = {}
    for species, st_sets in snippy_sets.items():
        for st_name, samples in st_sets.items():
            if not new_samples_set.isdisjoint(samples):
                if not species in active_snippy_sets:
                    active_snippy_sets[species] = {}
                active_snippy_sets[species][st_name] = samples
    

    # A comparison requires at least two samples.
    ska_sets_filtered = {k: v for k, v in active_ska_sets.items() if len(v) > 1}
    snippy_sets_filtered = {}
    for species, st_sets in active_snippy_sets.items():
        for st_name, samples in st_sets.items():
            if len(samples) > 1:
                if not species in snippy_sets_filtered:
                    snippy_sets_filtered[species] = {}
                snippy_sets_filtered[species][st_name] = samples
    
    print("🔬 Active Species sets found (by species with new isolates):")
    for name, samples in ska_sets_filtered.items():
        print(f"  - {name}: {len(samples)} samples")

    print(f"🔬 Active Subgroup sets found (by {config.get('snippy_grouping_method', 'ST').upper()} with new isolates):")
    for name, samples in snippy_sets_filtered.items():
        print(f"  - {name}: {len(samples)} samples")
        
    # Return the ACTIVATED and filtered dictionaries
    return dict(ska_sets_filtered), dict(snippy_sets_filtered)
