rule list_ancestries:
    ''' Make a tsv showing all unique ancestries in study table
    '''
    input:
        'output/ot_genetics_studies_table.{version}.tsv'
    output:
        tmpdir + '/ancestry_list.{version}.tsv'
    shell:
        'python scripts/list_ancestries.py '
        '--inf {output} '
        '--outf {output}'

rule get_postgap_data:
    ''' Download postgap data. This will be used as a temporary stop gap to
        extract LD table
    '''
    input:
        HTTPRemoteProvider().remote('https://storage.googleapis.com/postgap-data/postgap.20180615.txt.gz')
    output:
        tmpdir + '/postgap.20180615.txt.gz'
    shell:
        'cp {input} {output}'

rule get_1000g_variantID_lut:
    ''' Download rsid to variant id look up table for 1000G variants
    '''
    input:
        GSRemoteProvider().remote('genetics-portal-data/lut/1000g_rsid_to_variantid_lut.tsv.gz')
    output:
        tmpdir + '/1000g_rsid_to_variantid_lut.tsv.gz'
    shell:
        'cp {input} {output}'

rule parse_postgap:
    ''' Parses LD information from the postgap output. Temporary fix.
    '''
    input:
        in_postgap=tmpdir + '/postgap.20180615.txt.gz',
        in_lut=tmpdir + '/1000g_rsid_to_variantid_lut.tsv.gz'
    output:
        tmpdir + '/postgap_ld_table.temp.tsv.gz'
    shell:
        'python scripts/parse_ld_from_postgap_output.py '
        '--in_postgap {input.in_postgap} '
        '--in_lut {input.in_lut} '
        '--outf {output}'

rule merge_ld_to_loci:
    ''' Merges the postgap ld table to the top loci table and outputs into
        a standard format
    '''
    input:
        in_ld=tmpdir + '/postgap_ld_table.temp.tsv.gz',
        in_loci='output/ot_genetics_toploci_table.{version}.tsv'
    output:
        'output/ot_genetics_ld_table.{version}.tsv.gz'
    shell:
        'python scripts/merge_postgap_ld_to_top_loci.py '
        '--in_ld {input.in_ld} '
        '--in_loci {input.in_loci} '
        '--outf {output}'

rule ld_to_GCS:
    ''' Copy to GCS
    '''
    input:
        'output/ot_genetics_ld_table.{version}.tsv.gz'
    output:
        GSRemoteProvider().remote(
            '{gs_dir}/{{version}}/ot_genetics_ld_table.{{version}}.tsv.gz'.format(gs_dir=config['gs_dir'])
            )
    shell:
        'cp {input} {output}'
