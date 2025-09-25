
from pathlib import Path


db_dir = config['database_dir']

kraken_db_urls = {"Standard": "https://genome-idx.s3.amazonaws.com/kraken/k2_standard_20250714.tar.gz",
            "PlusPF": "https://genome-idx.s3.amazonaws.com/kraken/k2_pluspf_20250714.tar.gz",
            "PlusPFP": "https://genome-idx.s3.amazonaws.com/kraken/k2_pluspfp_20250714.tar.gz",
            "Standard-8": "https://genome-idx.s3.amazonaws.com/kraken/k2_standard_08_GB_20250714.tar.gz",
            "PlusPF-8": "https://genome-idx.s3.amazonaws.com/kraken/k2_pluspf_08_GB_20250714.tar.gz",
            "PlusPFP-8": "https://genome-idx.s3.amazonaws.com/kraken/k2_pluspfp_08_GB_20250714.tar.gz",
            "Standard-16": "https://genome-idx.s3.amazonaws.com/kraken/k2_standard_16_GB_20250714.tar.gz",
            "PlusPF-16": "https://genome-idx.s3.amazonaws.com/kraken/k2_pluspf_16_GB_20250714.tar.gz",
            "PlusPFP-16": "https://genome-idx.s3.amazonaws.com/kraken/k2_pluspfp_16_GB_20250714.tar.gz"
            }

dbs = [db_dir + '/kraken2/taxonomy.dmp', db_dir + '/gtdbtk/metadata/metadata.txt']

rule download:
    input: dbs

rule download_kraken_db:
    output: db_dir + "/kraken2/taxonomy.dmp"
    params: db_name = lambda wildcards: config.get('kraken2_db', 'PlusPFP'),
            url = kraken_db_urls[config.get('kraken2_db', 'PlusPFP')],
            tar_output = db_dir + '/kraken2.tar.gz',
            kraken_dir = db_dir + '/kraken2'
    shell:
        """
        wget {params.url} -O {params.tar_output}
        cd {params.kraken_dir}
        tar -xvf {params.tar_output}
        rm {params.tar_output}
        """

gtdbtk_db_url = "https://data.ace.uq.edu.au/public/gtdb/data/releases/latest/auxillary_files/gtdbtk_package/full_package/gtdbtk_data.tar.gz"

rule download_gtdbtk:
    output: temp(db_dir + '/gtdbtk_data.tar.gz')
    shell:
        """
        wget {gtdbtk_db_url} -O {output}
        """
rule extract_gtdbtk:
    input: db_dir + "/gtdbtk_data.tar.gz"
    output: db_dir + '/gtdbtk/metadata/metadata.txt'
    params: db_dir = db_dir
    shell:
        """
        mkdir -p {params.db_dir}/gtdbtk
        cd {params.db_dir}/gtdbtk/
        tar -xvf {input} 
        """
