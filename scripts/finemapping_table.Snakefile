rule convert_finemapping_to_standard:
    ''' Extract required fields from credible set file
    '''
    input:
        GSRemoteProvider().remote(config['credset'], keep_local=True) # DEBUG
    output:
        'output/ot_genetics_finemapping_table.{version}.tsv.gz'
    shell:
        'python scripts/format_finemapping_table.py '
        '--inf {input} '
        '--outf {output}'

rule finemap_to_GCS:
    ''' Copy to GCS
    '''
    input:
        'output/ot_genetics_finemapping_table.{version}.tsv.gz'
    output:
        GSRemoteProvider().remote(
            '{gs_dir}/{{version}}/ot_genetics_finemapping_table.{{version}}.tsv.gz'.format(gs_dir=config['gs_dir'])
            )
    shell:
        'cp {input} {output}'
