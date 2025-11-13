from os.path import join, abspath, basename
from os import makedirs

output_directory = 'db'
makedirs(output_directory, exist_ok=True)

# https://zenodo.org/records/15232972/files/S5.4.0.GTDB_r226.metapackage_20250331.smpkg.zb.tar.gz?download=1
singlem_metapackage = join(output_directory, 'S5.4.0.GTDB_r226.metapackage_20250331.smpkg')
singlem_metapackage_tgz = singlem_metapackage + '.zb.tar.gz'

sylph_package = join(output_directory, 'gtdb-r226-c200-dbv1.syldb')

#gtdb_bac120=join(output_directory, 'bac120_metadata_r226.tsv')
#gtdb_bac120_gz = gtdb_bac120 + '.gz'
gtdb_bac120_tax=join(output_directory, 'bac120_taxonomy_r226.tsv')
gtdb_bac120_tax_gz = gtdb_bac120_tax + '.gz'

#gtdb_ar53 = join(output_directory, 'ar53_metadata_r226.tsv')
#gtdb_ar53_gz = gtdb_ar53_tax + '.gz'
gtdb_ar53_tax = join(output_directory, 'ar53_taxonomy_r226.tsv')
gtdb_ar53_tax_gz = gtdb_ar53_tax + '.gz'

tools = ['singlem', 'sylph']

rule all:
    input:
        [join(output_directory, f'{tool}.done') for tool in tools],
        join(output_directory, 'gtdb.done'),
        join(output_directory, 'gtdb-ar.done')

rule singlem_download:
    output:
        done=touch(join(output_directory, 'singlem-download.done')),
        singlem_metapackage_tgz=singlem_metapackage_tgz
    log:
        join(output_directory, 'singlem-download.log')
    shell:
        "wget 'https://zenodo.org/records/15232972/files/S5.4.0.GTDB_r226.metapackage_20250331.smpkg.zb.tar.gz?download=1' -O {output.singlem_metapackage_tgz} &> {log}"

rule singlem_extract:
    input:
        done=join(output_directory, 'singlem-download.done'),
        singlem_metapackage_tgz=singlem_metapackage_tgz,
    output:
        done=touch(join(output_directory, 'singlem.done')),
        singlem_metapackage=directory(singlem_metapackage)
    log:
        abspath(join(output_directory, 'singlem-extract.log'))
    params:
        singlem_metapackage_basename = basename(singlem_metapackage)
    shell:
        "bash -c 'cd {output_directory} && tar -xzf {params.singlem_metapackage_basename}.zb.tar.gz && mv -v {params.singlem_metapackage_basename}.zb/payload_directory ../{output.singlem_metapackage}' &> {log}"

rule sylph_download:
    output:
        done=touch(join(output_directory, 'sylph.done')),
        sylph_package=sylph_package
    log:
        join(output_directory, 'sylph-download.log')
    shell:
        "wget 'http://faust.compbio.cs.cmu.edu/sylph-stuff/gtdb-r226-c200-dbv1.syldb' -O {output.sylph_package} &> {log}"

rule gtdb_download_bac120:
    output:
        done=touch(join(output_directory, 'gtdb_download.done')),
        tsv = gtdb_bac120_tax_gz
    log:
        join(output_directory, 'gtdb.log')
    shell:
        'wget https://data.gtdb.ecogenomic.org/releases/release226/226.0/bac120_taxonomy_r226.tsv.gz -O {output.tsv} &> {log}'
        #'wget https://data.gtdb.ecogenomic.org/releases/release226/226.0/bac120_metadata_r226.tsv.gz -O {output.tsv} &> {log}'

rule gtdb_extract_bac120:
    input:
        done=join(output_directory, 'gtdb_download.done'),
        tsv = gtdb_bac120_tax_gz
    output:
        done=touch(join(output_directory, 'gtdb.done')),
    log:
        join(output_directory, 'gtdb-extract.log')
    shell:
        'gzip -kd {input.tsv} &> {log}'

rule gtdb_download_ar53:
    output:
        done=touch(join(output_directory, 'gtdb_download_ar53.done')),
        tsv = gtdb_ar53_tax_gz
    log:
        join(output_directory, 'gtdb.log')
    shell:
        'wget https://data.gtdb.ecogenomic.org/releases/release226/226.0/ar53_taxonomy_r226.tsv.gz -O {output.tsv} &> {log}'
        #'wget https://data.gtdb.ecogenomic.org/releases/release226/226.0/ar53_metadata_r226.tsv.gz -O {output.tsv} &> {log}'

rule gtdb_extract_ar53:
    input:
        done=join(output_directory, 'gtdb_download_ar53.done'),
        tsv = gtdb_ar53_tax_gz
    output:
        done=touch(join(output_directory, 'gtdb-ar.done')),
    log:
        join(output_directory, 'gtdb-extract.log')
    shell:
        'gzip -kd {input.tsv} &> {log}'

