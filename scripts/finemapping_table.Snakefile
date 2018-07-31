rule convert_finemapping_to_standard:
    ''' Extract required fields from credible set file
    '''
    input:
        GSRemoteProvider().remote(config['credset'], keep_local=KEEP_LOCAL) # DEBUG
    output:
        'output/{version}/finemapping.tsv.gz'
    shell:
        'python scripts/format_finemapping_table.py '
        '--inf {input} '
        '--outf {output}'

rule finemap_to_GCS:
    ''' Copy to GCS
    '''
    input:
        'output/{version}/finemapping.tsv.gz'
    output:
        GSRemoteProvider().remote(
            '{gs_dir}/{{version}}/finemapping.tsv.gz'.format(gs_dir=config['gs_dir'])
            )
    shell:
        'cp {input} {output}'
